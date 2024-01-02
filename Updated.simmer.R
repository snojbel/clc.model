# Updated simmer with jobs

library(job)

# Resource initializations -------------------

Num.Res <- 25
res.Abund <-  50000

# Evenly distributed Resources

# SLC:
resource.freq.even.slc <- rep(1/Num.Res, times = Num.Res)                                      # res. freq. 
resource.prop.even.slc <- c(seq(from = -2.5, to = 2.5, length.out = Num.Res))            # res. property 
resource.freq.even.slc <- res.Abund*resource.freq.even.slc


# CLC:

resource.property.even.clc <- c(seq(from = -2.5, to = 2.5, length.out = Num.Res)) 

resource.frequency.even.clc <- rep(1/Num.Res, times = Num.Res)

resource.abundance.adults.even.clc     <- res.Abund                              # res. abundance of adults and juveniles
resource.abundance.juveniles.even.clc  <- res.Abund

resFreqMatrix.even.clc <- matrix(resource.frequency.even.clc, nrow=2, ncol=length(resource.frequency.even.clc ), byrow = TRUE)
resFreqMatrix.even.clc [1, ] <- resFreqMatrix.even.clc [1, ]*resource.abundance.adults.even.clc 
resFreqMatrix.even.clc [2, ] <- resFreqMatrix.even.clc [2, ]*resource.abundance.juveniles.even.clc 

resPropMatrix.even.clc <- matrix(resource.property.even.clc, nrow=2, ncol=length(resource.property.even.clc ), byrow = TRUE) 

rownames(resFreqMatrix.even.clc) <- c("Adult", "Juvenile")
colnames(resFreqMatrix.even.clc)  <- paste0("Resource ", 1:ncol(resFreqMatrix.even.clc))

rownames(resPropMatrix.even.clc)<-c("Adult", "Juvenile")
colnames(resPropMatrix.even.clc)  <- paste0("Resource ", 1:ncol(resPropMatrix.even.clc))

# Normal resources:

m <- 0 
s <- 1
N.resource.frequency <- c()
N.resource.property<- c(seq(from = -2.5, to = 2.5, length.out = Num.Res)) 

mid.add <- c()
midpoint <- c()

for(i in 1:(length(N.resource.property))){
  mid.add <- (N.resource.property[i+1]-N.resource.property[i])/2
  high.midpoint <- N.resource.property[(i)]+mid.add
  low.midpoint <- N.resource.property[(i)]-mid.add
  if(i == 1){
    N.resource.frequency[i] <- pnorm(high.midpoint, mean = m, sd = s) 
  }else if(i == length(N.resource.property)){
    low.midpoint <- N.resource.property[(i-1)] + (N.resource.property[i]-N.resource.property[i-1])/2
    N.resource.frequency[i] <- pnorm(low.midpoint, mean = m, sd = s, lower.tail = FALSE)
  }else{
    N.resource.frequency[i] <- pnorm(high.midpoint, mean = m, sd = s) - pnorm(low.midpoint, mean = m, sd = s) 
  }
}



resource.abundance.adults.norm.clc      <- res.Abund                              # res. abundance of adults and juveniles
resource.abundance.juveniles.norm.clc   <- res.Abund

# SLC:

resource.prop.norm.slc  <- N.resource.property             # res. property 
resource.freq.norm.slc  <- res.Abund*N.resource.frequency


# CLC:

resFreqMatrix.norm.clc  <- matrix(N.resource.frequency, nrow=2, ncol=length(N.resource.frequency), byrow = TRUE)

resFreqMatrix.norm.clc [1, ] <- resFreqMatrix.norm.clc [1, ]*resource.abundance.adults.norm.clc 
resFreqMatrix.norm.clc [2, ] <- resFreqMatrix.norm.clc [2, ]*resource.abundance.juveniles.norm.clc 

rownames(resFreqMatrix.norm.clc ) <- c("Adult", "Juvenile")
colnames(resFreqMatrix.norm.clc )  <- paste0("Resource ", 1:ncol(resFreqMatrix.norm.clc ))


resPropMatrix.norm.clc  <- matrix(N.resource.property, nrow=2, ncol=length(N.resource.property), byrow = TRUE) 


rownames(resPropMatrix.norm.clc )<-c("Adult", "Juvenile")
colnames(resPropMatrix.norm.clc )  <- paste0("Resource ", 1:ncol(resPropMatrix.norm.clc))


# Skewed resource distribution

# SLC

tot <- (Num.Res*(Num.Res+1))/2

x <- 1/tot
resource.freq.skew.slc <- c()

for (i in 1:Num.Res){ 
  resource.freq.skew.slc[i] <- i*x 
}                                     # res. freq. 

resource.prop.skew.slc <- c(seq(from = -2.5, to = 2.5, length.out = Num.Res))            # res. property 
resource.freq.skew.slc <- res.Abund*resource.freq.skew.slc


# CLC:

resource.property.skew.clc <- c(seq(from = -2.5, to = 2.5, length.out = Num.Res)) 


resource.frequency.skew.clc <- c()
for (i in 1:Num.Res){ 
  resource.frequency.skew.clc[i] <- i*x 
}    


resource.abundance.adults     <- res.Abund                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- res.Abund

resFreqMatrix.skew.clc <- matrix(resource.frequency.skew.clc, nrow=2, ncol=length(resource.frequency.skew.clc), byrow = TRUE)
resFreqMatrix.skew.clc[1, ] <- resFreqMatrix.skew.clc[1, ]*resource.abundance.adults
resFreqMatrix.skew.clc[2, ] <- resFreqMatrix.skew.clc[2, ]*resource.abundance.juveniles

resPropMatrix.skew.clc <- matrix(resource.property.skew.clc, nrow=2, ncol=length(resource.property.skew.clc), byrow = TRUE) 

rownames(resFreqMatrix.skew.clc) <- c("Adult", "Juvenile")
colnames(resFreqMatrix.skew.clc)  <- paste0("Resource ", 1:ncol(resFreqMatrix.skew.clc))

rownames(resPropMatrix.skew.clc)<-c("Adult", "Juvenile")
colnames(resPropMatrix.skew.clc)  <- paste0("Resource ", 1:ncol(resPropMatrix.skew.clc))


# Running simulations -----------------


# Even distribution -----------------------

# SLC

job::job(output.Even.SLC = {
  output.Even.SLC <- resourceCompetitionSLC(resProp=resource.prop.even.slc, iniP = 0, resFreq=resource.freq.even.slc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Even.SLC)
}, import = c(resource.prop.even.slc, resource.freq.even.slc, resourceCompetitionSLC))

# CLC

job::job(output.Even.CLC = {
  output.Even.CLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix.even.clc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Even.CLC)
}, import = c(resPropMatrix.even.clc, resFreqMatrix.even.clc, resourceCompetitionCLC))



# Normal Distribution --------------------------

# SLC

job::job(output.Norm.SLC = {
  output.Norm.SLC <- resourceCompetitionSLC(resProp=resource.prop.norm.slc, iniP = 0, resFreq=resource.freq.norm.slc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Norm.SLC)
}, import = c(resource.prop.norm.slc, resource.freq.norm.slc, resourceCompetitionSLC))

# CLC

job::job(output.Norm.CLC = {
  output.Norm.CLC <- resourceCompetitionCLC(resProp=resPropMatrix.norm.clc, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix.norm.clc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Norm.CLC)
}, import = c(resPropMatrix.norm.clc, resFreqMatrix.norm.clc, resourceCompetitionCLC))




# Skewed distribution ------------------------


# SLC

job::job(output.Skew.SLC = {
  output.Skew.SLC <- resourceCompetitionSLC(resProp=resource.prop.skew.slc, iniP = 0, resFreq=resource.freq.skew.slc, 
                                            resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Skew.SLC)
}, import = c(resource.prop.skew.slc, resource.freq.skew.slc, resourceCompetitionSLC))

# CLC

job::job(output.Skew.CLC = {
  output.Skew.CLC <- resourceCompetitionCLC(resProp=resPropMatrix.skew.clc, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix.skew.clc, 
                                            resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Skew.CLC)
}, import = c(resPropMatrix.skew.clc, resFreqMatrix.skew.clc, resourceCompetitionCLC))


# ------------------------------











