

# Simulation runner for the purpose of only running one sim and extracting all results. 

# Loading libaries

library(job)
library(gganimate)
library(gifski)
library(png)
library(ggmatplot)
library(tidyverse)
library(viridisLite)  # Color things
library(viridis)
library(ggplot2)      # Prettier plots
library(gridExtra)  
library(ggpubr)
library(cowplot)
library(extrafont) 
library(patchwork)
library(dplyr)
library(grid)




# Resource initializations -------------------


Num.Res <- 16
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



# Two resources 

resource.prop <- c(-1,1)
resource.frequency <- c(0.5, 0.5)
resource.frequency.as <- c(0.2, 0.8)

resFreqMatrix.2res <- matrix(resource.frequency, nrow=2, ncol=length(resource.frequency), byrow = TRUE)
resFreqMatrixAs.2res <- matrix(resource.frequency.as, nrow=2, ncol=length(resource.frequency.as), byrow = TRUE)


resFreqMatrix.2res[1, ] <- resFreqMatrix.2res[1, ]*res.Abund
resFreqMatrix.2res[2, ] <- resFreqMatrix.2res[2, ]*res.Abund

resFreqMatrixAs.2res[1, ] <- resFreqMatrixAs.2res[1, ]*res.Abund
resFreqMatrixAs.2res[2, ] <- resFreqMatrixAs.2res[2, ]*res.Abund

rownames(resFreqMatrix.2res) <- c("Adult", "Juvenile")
colnames(resFreqMatrix.2res)  <- paste0("Resource ", 1:ncol(resFreqMatrixAs.2res))

rownames(resFreqMatrixAs.2res) <- c("Adult", "Juvenile")
colnames(resFreqMatrixAs.2res)  <- paste0("Resource ", 1:ncol(resFreqMatrixAs.2res))

resPropMatrix.2res <- matrix(resource.prop, nrow=2, ncol=length(resource.prop), byrow = TRUE) 


rownames(resPropMatrix.2res)<-c("Adult", "Juvenile")
colnames(resFreqMatrix.2res)  <- paste0("Resource ", 1:ncol(resPropMatrix.2res))


# ------------------------

# Parameter initializations ----------------------------


popSize <- 10
sigma <- seq(from = 0.05, to = 1.25, length.out = 5)
im <-  0 
fmax <-  2
kA <-  0.5
kJ <-  0.5
mutProb <- 0.0005
mutVar <- 0.05
time.steps <- 10000
iniP <- 0
iniPJ <- 0
iniPA <- 0
nmorphs <-  1
threshold <-  0.005


# -------------------------

# Running simulations -----------------------------


# Even

job::job(even = {
  
  Total.species.SLC.single.even <- c()
  
  Total.species.CLC.even <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.even) <- sigma  #ADULTS
  colnames(Total.species.CLC.even) <- sigma #JUVENILES
  
  
  
  # SLC
  
  Total.species.SLC.even <- list()
  Total.endpoint.SLC.even <- list()
  
  for(r in 1:10) {
    
    print(paste0("loop ", r, " started"))
    
    Number.species.SLC.even <- c()
    endpoint.SLC.even <- list()
    
    for(i in 1:length(sigma)){
      
      
      outputSLC <- resourceCompetitionSLC(resProp=resource.prop.even.slc, iniP = iniP, resFreq=resource.freq.even.slc, resGen=matrix(c(sigma[i],sigma[i])),
                                          popSize = popSize, mutProb=mutProb, mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                          threshold = threshold)
      
      
      #Filter out similar "species" and collect number of species data
      
      final.data.SLC.even <- slc.groups(output = outputSLC)
      Number.species.SLC.even[i] <- nrow(final.data.SLC.even)
      
      #Collect endpoint data
      
      endpoint.SLC.even[[i]] <- final.data.SLC.even 
    }
    
    Total.species.SLC.even[[r]] <- Number.species.SLC.even
    Total.endpoint.SLC.even[[r]] <- endpoint.SLC.even
  }
  
  
  # Caluclating mean and SD of 10 runs
  
  
  
  Total.mean.SLC.even <- sapply(1:length(sigma), function(i) mean(sapply(Total.species.SLC.even, function(x) x[i])))
  
  array.data.SLC <- array(unlist(Total.species.SLC.even), dim = c(dim(Total.species.SLC.even[[1]]), length(Total.species.SLC.even)))
  
  Total.sd.SLC.even <- sapply(1:length(sigma), function(i) sd(sapply(Total.species.SLC.even, function(x) x[i])))
  
  
  # CLC
  
  print("clc start")
  
  Total.species.CLC.even <- list()
  
  Total.endpoint.CLC.even <- list()
  
  
  for(a in 1:10){
    print(paste0("loop ", a, " started"))
    
    species.CLC.even <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
    rownames(species.CLC.even) <- sigma  #ADULTS
    colnames(species.CLC.even) <- sigma #JUVENILES
    
    endpoint.CLC.even <- c()
    
    for(b in 1:length(sigma)){
      
      for(k in 1:length(sigma)){
        
        
        outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, resFreq=resFreqMatrix.even.clc, iniPA = iniPA, iniPJ = iniPJ, resGen=matrix(c(sigma[b],sigma[k])), 
                                            popSize = popSize, mutProb=mutProb, mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                            threshold = threshold)
        
        
        
        #Filter out similar "species"
        final.data.CLC.even <- clc.groups(output = outputCLC)
        
        # Collect Data
        
        species.CLC.even[b, k] <- nrow(final.data.CLC.even)
        final.data.CLC.even$Adult.gen <- c(rep(sigma[b], times = nrow(final.data.CLC.even)))
        final.data.CLC.even$Juv.gen <- c(rep(sigma[k], times = nrow(final.data.CLC.even)))
        endpoint.CLC.even <- rbind(endpoint.CLC.even, final.data.CLC.even) 
        
      }
      
    }
    Total.species.CLC.even[[a]] <- species.CLC.even
    Total.endpoint.CLC.even[[a]] <- endpoint.CLC.even
  }
  
  # Calculating mean of 10 runs
  
  # Combine matrices in the list into a 3D array
  array.data.CLC <- array(unlist(Total.species.CLC.even), dim = c(dim(Total.species.CLC.even[[1]]), length(Total.species.CLC.even)))
  
  
  # Calculate mean and standard deviation along the third dimension (across the list)
  Total.mean.CLC.even <- apply(array.data.CLC, c(1, 2), mean)
  Total.sd.CLC.even <- apply(array.data.CLC, c(1, 2), sd)
  
  
  
  job::export(list(Total.mean.CLC.even, Total.sd.CLC.even, Total.mean.SLC.even, Total.sd.SLC.even, Total.endpoint.SLC.even, Total.endpoint.CLC.even))
}, import = c(resPropMatrix.even.clc, resFreqMatrix.even.clc, resourceCompetitionCLC, resource.prop.even.slc, resource.freq.even.slc, resourceCompetitionSLC, clc.groups, slc.groups, sigma,
              popSize, im, fmax, kA, kJ, mutProb, mutVar, time.steps, iniP, iniPA, iniPJ, nmorphs, threshold))




# Normal 

job::job(norm = {
  
  Total.species.SLC.single.norm <- c()
  
  Total.species.CLC.norm <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.norm) <- sigma  #ADULTS
  colnames(Total.species.CLC.norm) <- sigma #JUVENILES
  
  
  
  # SLC
  
  Total.species.SLC.norm <- list()
  Total.endpoint.SLC.norm <- list()
  
  for(r in 1:10) {
    
    print(paste0("loop ", r, " started"))
    
    Number.species.SLC.norm <- c()
    endpoint.SLC.norm <- list()
    
    for(i in 1:length(sigma)){
      
      
      outputSLC <- resourceCompetitionSLC(resProp=resource.prop.norm.slc, iniP = iniP, resFreq=resource.freq.norm.slc, resGen=matrix(c(sigma[i],sigma[i])),
                                          popSize = popSize, mutProb=mutProb, mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                          threshold = threshold)
      
      
      #Filter out similar "species" and collect number of species data
      
      final.data.SLC.norm <- slc.groups(output = outputSLC)
      Number.species.SLC.norm[i] <- nrow(final.data.SLC.norm)
      
      #Collect endpoint data
      
      endpoint.SLC.norm[[i]] <- final.data.SLC.norm 
    }
    
    Total.species.SLC.norm[[r]] <- Number.species.SLC.norm
    Total.endpoint.SLC.norm[[r]] <- endpoint.SLC.norm
  }
  
  
  # Caluclating mean and SD of 10 runs
  
  
  
  Total.mean.SLC.norm <- sapply(1:length(sigma), function(i) mean(sapply(Total.species.SLC.norm, function(x) x[i])))
  
  array.data.SLC <- array(unlist(Total.species.SLC.norm), dim = c(dim(Total.species.SLC.norm[[1]]), length(Total.species.SLC.norm)))
  
  Total.sd.SLC.norm <- sapply(1:length(sigma), function(i) sd(sapply(Total.species.SLC.norm, function(x) x[i])))
  
  
  # CLC
  
  print("clc start")
  
  Total.species.CLC.norm <- list()
  
  Total.endpoint.CLC.norm <- list()
  
  
  for(a in 1:10){
    print(paste0("loop ", a, " started"))
    
    species.CLC.norm <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
    rownames(species.CLC.norm) <- sigma  #ADULTS
    colnames(species.CLC.norm) <- sigma #JUVENILES
    
    endpoint.CLC.norm <- c()
    
    for(b in 1:length(sigma)){
      
      for(k in 1:length(sigma)){
        
        
        outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.norm.clc, resFreq=resFreqMatrix.norm.clc, iniPA = iniPA, iniPJ = iniPJ, resGen=matrix(c(sigma[b],sigma[k])), 
                                            popSize = popSize, mutProb=mutProb, mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                            threshold = threshold)
        
        
        
        #Filter out similar "species"
        final.data.CLC.norm <- clc.groups(output = outputCLC)
        
        # Collect Data
        
        species.CLC.norm[b, k] <- nrow(final.data.CLC.norm)
        final.data.CLC.norm$Adult.gen <- c(rep(sigma[b], times = nrow(final.data.CLC.norm)))
        final.data.CLC.norm$Juv.gen <- c(rep(sigma[k], times = nrow(final.data.CLC.norm)))
        endpoint.CLC.norm <- rbind(endpoint.CLC.norm, final.data.CLC.norm) 
        
      }
      
    }
    Total.species.CLC.norm[[a]] <- species.CLC.norm
    Total.endpoint.CLC.norm[[a]] <- endpoint.CLC.norm
  }
  
  # Calculating mean of 10 runs
  
  # Combine matrices in the list into a 3D array
  array.data.CLC <- array(unlist(Total.species.CLC.norm), dim = c(dim(Total.species.CLC.norm[[1]]), length(Total.species.CLC.norm)))
  
  
  # Calculate mean and standard deviation along the third dimension (across the list)
  Total.mean.CLC.norm <- apply(array.data.CLC, c(1, 2), mean)
  Total.sd.CLC.norm <- apply(array.data.CLC, c(1, 2), sd)
  
  
  
  job::export(list(Total.mean.CLC.norm, Total.sd.CLC.norm, Total.mean.SLC.norm, Total.sd.SLC.norm, Total.endpoint.SLC.norm, Total.endpoint.CLC.norm))
}, import = c(resPropMatrix.norm.clc, resFreqMatrix.norm.clc, resourceCompetitionCLC, resource.prop.norm.slc, resource.freq.norm.slc, resourceCompetitionSLC, clc.groups, slc.groups, sigma,
              popSize, im, fmax, kA, kJ, mutProb, mutVar, time.steps, iniP, iniPA, iniPJ, nmorphs, threshold))



# Skewed


job::job(skew = {
  
  Total.species.SLC.single.skew <- c()
  
  Total.species.CLC.skew <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.skew) <- sigma  #ADULTS
  colnames(Total.species.CLC.skew) <- sigma #JUVENILES
  
  
  
  # SLC
  
  Total.species.SLC.skew <- list()
  Total.endpoint.SLC.skew <- list()
  
  for(r in 1:10) {
    
    print(paste0("loop ", r, " started"))
    
    Number.species.SLC.skew <- c()
    endpoint.SLC.skew <- list()
    
    for(i in 1:length(sigma)){
      
      
      outputSLC <- resourceCompetitionSLC(resProp=resource.prop.skew.slc, iniP = iniP, resFreq=resource.freq.skew.slc, resGen=matrix(c(sigma[i],sigma[i])),
                                          popSize = popSize, mutProb=mutProb, mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                          threshold = threshold)
      
      
      #Filter out similar "species" and collect number of species data
      
      final.data.SLC.skew <- slc.groups(output = outputSLC)
      Number.species.SLC.skew[i] <- nrow(final.data.SLC.skew)
      
      #Collect endpoint data
      
      endpoint.SLC.skew[[i]] <- final.data.SLC.skew 
    }
    
    Total.species.SLC.skew[[r]] <- Number.species.SLC.skew
    Total.endpoint.SLC.skew[[r]] <- endpoint.SLC.skew
  }
  
  
  # Caluclating mean and SD of 10 runs
  
  
  
  Total.mean.SLC.skew <- sapply(1:length(sigma), function(i) mean(sapply(Total.species.SLC.skew, function(x) x[i])))
  
  array.data.SLC <- array(unlist(Total.species.SLC.skew), dim = c(dim(Total.species.SLC.skew[[1]]), length(Total.species.SLC.skew)))
  
  Total.sd.SLC.skew <- sapply(1:length(sigma), function(i) sd(sapply(Total.species.SLC.skew, function(x) x[i])))
  
  
  # CLC
  
  print("clc start")
  
  Total.species.CLC.skew <- list()
  
  Total.endpoint.CLC.skew <- list()
  
  
  for(a in 1:10){
    print(paste0("loop ", a, " started"))
    
    species.CLC.skew <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
    rownames(species.CLC.skew) <- sigma  #ADULTS
    colnames(species.CLC.skew) <- sigma #JUVENILES
    
    endpoint.CLC.skew <- c()
    
    for(b in 1:length(sigma)){
      
      for(k in 1:length(sigma)){
        
        
        outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.skew.clc, resFreq=resFreqMatrix.skew.clc, iniPA = iniPA, iniPJ = iniPJ, resGen=matrix(c(sigma[b],sigma[k])), 
                                            popSize = popSize, mutProb=mutProb, mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                            threshold = threshold)
        
        
        
        #Filter out similar "species"
        final.data.CLC.skew <- clc.groups(output = outputCLC)
        
        # Collect Data
        
        species.CLC.skew[b, k] <- nrow(final.data.CLC.skew)
        final.data.CLC.skew$Adult.gen <- c(rep(sigma[b], times = nrow(final.data.CLC.skew)))
        final.data.CLC.skew$Juv.gen <- c(rep(sigma[k], times = nrow(final.data.CLC.skew)))
        endpoint.CLC.skew <- rbind(endpoint.CLC.skew, final.data.CLC.skew) 
        
      }
      
    }
    Total.species.CLC.skew[[a]] <- species.CLC.skew
    Total.endpoint.CLC.skew[[a]] <- endpoint.CLC.skew
  }
  
  # Calculating mean of 10 runs
  
  # Combine matrices in the list into a 3D array
  array.data.CLC <- array(unlist(Total.species.CLC.skew), dim = c(dim(Total.species.CLC.skew[[1]]), length(Total.species.CLC.skew)))
  
  
  # Calculate mean and standard deviation along the third dimension (across the list)
  Total.mean.CLC.skew <- apply(array.data.CLC, c(1, 2), mean)
  Total.sd.CLC.skew <- apply(array.data.CLC, c(1, 2), sd)
  
  
  
  job::export(list(Total.mean.CLC.skew, Total.sd.CLC.skew, Total.mean.SLC.skew, Total.sd.SLC.skew, Total.endpoint.SLC.skew, Total.endpoint.CLC.skew))
}, import = c(resPropMatrix.skew.clc, resFreqMatrix.skew.clc, resourceCompetitionCLC, resource.prop.skew.slc, resource.freq.skew.slc, resourceCompetitionSLC, clc.groups, slc.groups, sigma,
              popSize, im, fmax, kA, kJ, mutProb, mutVar, time.steps, iniP, iniPA, iniPJ, nmorphs, threshold))


# -------------------------------------------------

# Plotting Mean number of Species ----------------------------------------------

# Even

Total.mean.CLC.even <- even$Total.mean.CLC.even
Total.mean.SLC.even <- even$Total.mean.SLC.even

Total.sd.CLC.even <-  even$Total.sd.CLC.even
Total.sd.SLC.even <-  even$Total.sd.SLC.even

x <- as.factor(sigma)



df.CLC <- data.frame(
  Juvenile.trait = rep(x, each = length(x)),
  Adult.trait = rep(x, times = length(x)),
  Richness = as.vector(Total.mean.CLC.even),
  sd = as.vector(Total.sd.CLC.even),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df.SLC <- data.frame(
  Juvenile.trait = x,
  Adult.trait = x,
  Richness = as.vector(Total.mean.SLC.even),
  sd = as.vector(Total.sd.SLC.even),
  Cycle = rep("Simple", times = length(x))
)



df.combined <- rbind(df.CLC, df.SLC)


ggplot(df.combined, aes(x = Adult.trait, y = Richness, shape = Cycle, color = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 5) +
  #geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Number of species") +
  labs(color = "Juvenile Generalism", shape = "Life strategy") +
  ggtitle("Even Resource distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))#+
  #scale_color_manual(values = c("slateblue", "thistle"))


# Normal


Total.sd.CLC.norm <-  norm$Total.sd.CLC.norm
Total.sd.SLC.norm <-  norm$Total.sd.SLC.norm

x <- as.factor(sigma)



df.CLC <- data.frame(
  Juvenile.trait = rep(x, each = length(x)),
  Adult.trait = rep(x, times = length(x)),
  Richness = as.vector(Total.mean.CLC.norm),
  sd = as.vector(Total.sd.CLC.norm),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df.SLC <- data.frame(
  Juvenile.trait = x,
  Adult.trait = x,
  Richness = as.vector(Total.mean.SLC.norm),
  sd = as.vector(Total.sd.SLC.norm),
  Cycle = rep("Simple", times = length(x))
)



df.combined <- rbind(df.CLC, df.SLC)


ggplot(df.combined, aes(x = Adult.trait, y = Richness, shape = Cycle, color = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 5) +
  #geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Number of species") +
  labs(color = "Juvenile Generalism", shape = "Life strategy") +
  ggtitle("Normal Resource distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))#+
#scale_color_manual(values = c("slateblue", "thistle"))


# Skewed


Total.mean.CLC.skew <- skew$Total.mean.CLC.skew
Total.mean.SLC.skew <- skew$Total.mean.SLC.skew

Total.sd.CLC.skew <-  skew$Total.sd.CLC.skew
Total.sd.SLC.skew <-  skew$Total.sd.SLC.skew

x <- as.factor(sigma)



df.CLC <- data.frame(
  Juvenile.trait = rep(x, each = length(x)),
  Adult.trait = rep(x, times = length(x)),
  Richness = as.vector(Total.mean.CLC.skew),
  sd = as.vector(Total.sd.CLC.skew),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df.SLC <- data.frame(
  Juvenile.trait = x,
  Adult.trait = x,
  Richness = as.vector(Total.mean.SLC.skew),
  sd = as.vector(Total.sd.SLC.skew),
  Cycle = rep("Simple", times = length(x))
)



df.combined <- rbind(df.CLC, df.SLC)


ggplot(df.combined, aes(x = Adult.trait, y = Richness, shape = Cycle, color = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 5) +
  #geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Number of species") +
  labs(color = "Juvenile Generalism", shape = "Life strategy") +
  ggtitle("Skewed Resource distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))#+
#scale_color_manual(values = c("slateblue", "thistle"))


#--------------------------------

# Plotting Phenotype Endpoint --------------------------------------------------


#--------------------------------







