

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


a <- list(1,2,3)
b <- list("Hello", "Goodday")

matrix(c(a, b))


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
time.steps <- 50000
iniP <- 0
iniPJ <- 0
iniPA <- 0
nmorphs <-  1
threshold <-  0.005


# -------------------------

# Running simulation


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
    rownames(Total.species.CLC.even) <- sigma  #ADULTS
    colnames(Total.species.CLC.even) <- sigma #JUVENILES
    
    endpoint.CLC.even <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
    rownames(Total.species.CLC.even) <- sigma  #ADULTS
    colnames(Total.species.CLC.even) <- sigma #JUVENILES
    
    for(b in 1:length(sigma)){
      
      for(k in 1:length(sigma)){
        
        
        outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, resFreq=resFreqMatrix.even.clc, iniPA = iniPA, iniPJ = iniPJ, resGen=matrix(c(sigma[b],sigma[k])), 
                                            popSize = popSize, mutProb=mutProb, mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                            threshold = threshold)
        
        
        
        #Filter out similar "species"
        final.data.CLC.even <- clc.groups(output = outputCLC)
        
        # Collect Data
        
        species.CLC.even[b, k] <- nrow(final.data.CLC.even)
        endpoint.CLC.even[b, k] <- final.data.CLC.even 
        
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












