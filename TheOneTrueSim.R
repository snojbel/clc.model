

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
library(FamilyRank)



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



# Bimodal resources


m1 <- -1.25
m2 <-  1.25
s <- 0.5
Bi.resource.frequency <- c()
Bi.resource.property<- c(seq(from = -2.5, to = 2.5, length.out = Num.Res)) 



mid.add <- c()
midpoint <- c()

for(i in 1:(length(Bi.resource.property))){
  mid.add <- (Bi.resource.property[i+1]-Bi.resource.property[i])/2
  high.midpoint <- Bi.resource.property[(i)]+mid.add
  low.midpoint <- Bi.resource.property[(i)]-mid.add
  if(i == 1){
    Bi.resource.frequency[i] <- pnorm(high.midpoint, mean = m1, sd = s)/2 
  }else if(i == length(Bi.resource.property)){
    low.midpoint <- Bi.resource.property[(i-1)] + (Bi.resource.property[i]-Bi.resource.property[i-1])/2
    Bi.resource.frequency[i] <- pnorm(low.midpoint, mean = m2, sd = s, lower.tail = FALSE)
  }else if (Bi.resource.property[i]<0) {
    Bi.resource.frequency[i] <- (pnorm(high.midpoint, mean = m1, sd = s) - pnorm(low.midpoint, mean = m1, sd = s))/2
  }else{
    Bi.resource.frequency[i] <- (pnorm(high.midpoint, mean = m2, sd = s) - pnorm(low.midpoint, mean = m2, sd = s))/2
  }
}


resource.abundance.adults.binorm.clc      <- res.Abund                              # res. abundance of adults and juveniles
resource.abundance.juveniles.binorm.clc   <- res.Abund

# SLC:

resource.prop.binorm.slc  <- Bi.resource.property             # res. property 
resource.freq.binorm.slc  <- res.Abund*Bi.resource.frequency


# CLC:

resFreqMatrix.binorm.clc  <- matrix(Bi.resource.frequency, nrow=2, ncol=length(Bi.resource.frequency), byrow = TRUE)

resFreqMatrix.binorm.clc [1, ] <- resFreqMatrix.binorm.clc [1, ]*resource.abundance.adults.binorm.clc 
resFreqMatrix.binorm.clc [2, ] <- resFreqMatrix.binorm.clc [2, ]*resource.abundance.juveniles.binorm.clc 

rownames(resFreqMatrix.binorm.clc ) <- c("Adult", "Juvenile")
colnames(resFreqMatrix.binorm.clc )  <- paste0("Resource ", 1:ncol(resFreqMatrix.binorm.clc ))


resPropMatrix.binorm.clc  <- matrix(Bi.resource.property, nrow=2, ncol=length(Bi.resource.property), byrow = TRUE) 


rownames(resPropMatrix.binorm.clc )<-c("Adult", "Juvenile")
colnames(resPropMatrix.binorm.clc )  <- paste0("Resource ", 1:ncol(resPropMatrix.binorm.clc))


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
#sigma <- seq(from = 0.05, to = 0.8, length.out = 6)
sigma <- seq(from = 0.05, to = 0.2, length.out = 6) 
im <-  0.5 
fmax <-  2
kA <-  0.5
kJ <-  0.5
mutProb <- 0.005
mutVar <- 0.05
time.steps <- 50000
iniP <- 2
iniPJ <- 2
iniPA <- 2
nmorphs <-  1
threshold <-  0.005
maxTr = 3
minTr = -3

# For No mutation:
nmorphs <- 200

iniP <- runif(200, min = minTr, max = maxTr)


iniPA <- runif(200, min = minTr, max = maxTr)
iniPJ <- runif(200, min = minTr, max = maxTr)



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
  
  for(r in 1:3) {
    
    
    id <- 1
    
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
      final.data.SLC.even$ID <- c(rep(id, times = nrow(final.data.SLC.even)))
     
      
      endpoint.SLC.even[[i]] <- final.data.SLC.even 
      id <- id + 1 
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
  
  
  for(a in 1:3){
    print(paste0("loop ", a, " started"))
    
    
    id <- 1
    
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
        final.data.CLC.even$ID <- c(rep(id, times = nrow(final.data.CLC.even)))
        endpoint.CLC.even <- rbind(endpoint.CLC.even, final.data.CLC.even) 
        id <- id + 1
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
              popSize, im, fmax, kA, kJ, mutProb, mutVar, time.steps, iniP, iniPA, iniPJ, nmorphs, threshold, maxTr, minTr))




# Normal 

job::job(norm = {
  
  Total.species.SLC.single.norm <- c()
  
  Total.species.CLC.norm <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.norm) <- sigma  #ADULTS
  colnames(Total.species.CLC.norm) <- sigma #JUVENILES
  
  
  
  # SLC
  
  Total.species.SLC.norm <- list()
  Total.endpoint.SLC.norm <- list()
  
  for(r in 1:3) {
    
    
    id <- 1
    
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
      final.data.SLC.norm$ID <- c(rep(id, times = nrow(final.data.SLC.norm)))
      
      
      endpoint.SLC.norm[[i]] <- final.data.SLC.norm 
      id <- id + 1 
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
  
  
  for(a in 1:3){
    print(paste0("loop ", a, " started"))
    
    
    id <- 1
    
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
        final.data.CLC.norm$ID <- c(rep(id, times = nrow(final.data.CLC.norm)))
        endpoint.CLC.norm <- rbind(endpoint.CLC.norm, final.data.CLC.norm) 
        id <- id + 1
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
              popSize, im, fmax, kA, kJ, mutProb, mutVar, time.steps, iniP, iniPA, iniPJ, nmorphs, threshold, maxTr, minTr))





# Skewed


job::job(skew = {
  
  Total.species.SLC.single.skew <- c()
  
  Total.species.CLC.skew <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.skew) <- sigma  #ADULTS
  colnames(Total.species.CLC.skew) <- sigma #JUVENILES
  
  
  
  # SLC
  
  Total.species.SLC.skew <- list()
  Total.endpoint.SLC.skew <- list()
  
  for(r in 1:3) {
    
    
    id <- 1
    
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
      final.data.SLC.skew$ID <- c(rep(id, times = nrow(final.data.SLC.skew)))
      
      
      endpoint.SLC.skew[[i]] <- final.data.SLC.skew 
      id <- id + 1 
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
  
  
  for(a in 1:3){
    
    
    print(paste0("loop ", a, " started"))
    
    id <- 1
    
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
        final.data.CLC.skew$ID <- c(rep(id, times = nrow(final.data.CLC.skew)))
        endpoint.CLC.skew <- rbind(endpoint.CLC.skew, final.data.CLC.skew) 
        id <- id + 1
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
              popSize, im, fmax, kA, kJ, mutProb, mutVar, time.steps, iniP, iniPA, iniPJ, nmorphs, threshold, maxTr, minTr))


# Bimodal Normal

job::job(binorm= {
  
  Total.species.SLC.single.binorm <- c()
  
  Total.species.CLC.binorm <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.binorm) <- sigma  #ADULTS
  colnames(Total.species.CLC.binorm) <- sigma #JUVENILES
  
  
  
  # SLC
  
  Total.species.SLC.binorm <- list()
  Total.endpoint.SLC.binorm <- list()
  
  for(r in 1:3) {
    
    
    id <- 1
    
    print(paste0("loop ", r, " started"))
    
    Number.species.SLC.binorm <- c()
    endpoint.SLC.binorm <- list()
    
    for(i in 1:length(sigma)){
      
      
      outputSLC <- resourceCompetitionSLC(resProp=resource.prop.binorm.slc, iniP = iniP, resFreq=resource.freq.binorm.slc, resGen=matrix(c(sigma[i],sigma[i])),
                                          popSize = popSize, mutProb=mutProb, mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                          threshold = threshold)
      
      
      #Filter out similar "species" and collect number of species data
      
      final.data.SLC.binorm <- slc.groups(output = outputSLC)
      Number.species.SLC.binorm[i] <- nrow(final.data.SLC.binorm)
      
      #Collect endpoint data
      final.data.SLC.binorm$ID <- c(rep(id, times = nrow(final.data.SLC.binorm)))
      
      
      endpoint.SLC.binorm[[i]] <- final.data.SLC.binorm 
      id <- id + 1 
    }
    
    Total.species.SLC.binorm[[r]] <- Number.species.SLC.binorm
    Total.endpoint.SLC.binorm[[r]] <- endpoint.SLC.binorm
  }
  
  
  # Caluclating mean and SD of 10 runs
  
  
  
  Total.mean.SLC.binorm <- sapply(1:length(sigma), function(i) mean(sapply(Total.species.SLC.binorm, function(x) x[i])))
  
  array.data.SLC <- array(unlist(Total.species.SLC.binorm), dim = c(dim(Total.species.SLC.binorm[[1]]), length(Total.species.SLC.binorm)))
  
  Total.sd.SLC.binorm <- sapply(1:length(sigma), function(i) sd(sapply(Total.species.SLC.binorm, function(x) x[i])))
  
  
  # CLC
  
  print("clc start")
  
  Total.species.CLC.binorm <- list()
  
  Total.endpoint.CLC.binorm <- list()
  
  
  for(a in 1:3){
    print(paste0("loop ", a, " started"))
    
    id <- 1
    
    species.CLC.binorm <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
    rownames(species.CLC.binorm) <- sigma  #ADULTS
    colnames(species.CLC.binorm) <- sigma #JUVENILES
    
    endpoint.CLC.binorm <- c()
    
    for(b in 1:length(sigma)){
      
      for(k in 1:length(sigma)){
        
        
        outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.binorm.clc, resFreq=resFreqMatrix.binorm.clc, iniPA = iniPA, iniPJ = iniPJ, resGen=matrix(c(sigma[b],sigma[k])), 
                                            popSize = popSize, mutProb=mutProb, mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                            threshold = threshold)
        
        
        
        #Filter out similar "species"
        final.data.CLC.binorm <- clc.groups(output = outputCLC)
        
        # Collect Data
        
        species.CLC.binorm[b, k] <- nrow(final.data.CLC.binorm)
        final.data.CLC.binorm$Adult.gen <- c(rep(sigma[b], times = nrow(final.data.CLC.binorm)))
        final.data.CLC.binorm$Juv.gen <- c(rep(sigma[k], times = nrow(final.data.CLC.binorm)))
        final.data.CLC.binorm$ID <- c(rep(id, times = nrow(final.data.CLC.binorm)))
        endpoint.CLC.binorm <- rbind(endpoint.CLC.binorm, final.data.CLC.binorm) 
        id <- id + 1
      }
      
    }
    Total.species.CLC.binorm[[a]] <- species.CLC.binorm
    Total.endpoint.CLC.binorm[[a]] <- endpoint.CLC.binorm
  }
  
  # Calculating mean of 10 runs
  
  # Combine matrices in the list into a 3D array
  array.data.CLC <- array(unlist(Total.species.CLC.binorm), dim = c(dim(Total.species.CLC.binorm[[1]]), length(Total.species.CLC.binorm)))
  
  
  # Calculate mean and standard deviation along the third dimension (across the list)
  Total.mean.CLC.binorm <- apply(array.data.CLC, c(1, 2), mean)
  Total.sd.CLC.binorm <- apply(array.data.CLC, c(1, 2), sd)
  
  
  
  job::export(list(Total.mean.CLC.binorm, Total.sd.CLC.binorm, Total.mean.SLC.binorm, Total.sd.SLC.binorm, Total.endpoint.SLC.binorm, Total.endpoint.CLC.binorm))
}, import = c(resPropMatrix.binorm.clc, resFreqMatrix.binorm.clc, resourceCompetitionCLC, resource.prop.binorm.slc, resource.freq.binorm.slc, resourceCompetitionSLC, clc.groups, slc.groups, sigma,
              popSize, im, fmax, kA, kJ, mutProb, mutVar, time.steps, iniP, iniPA, iniPJ, nmorphs, threshold, maxTr, minTr))



# -------------------------------------------------





# Plotting Mean number of Species ----------------------------------------------

# Name change for plotting reasons ( be sure to save previous run before writing over)

#even <- even.stat
#norm <- norm.stat
#skew <- skew.stat
#binorm <- binorm.stat


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

color_palette <- magma(length(sigma))

even.plot <- ggplot(df.combined, aes(x = Adult.trait, y = Richness, shape = Cycle, color = Juvenile.trait, stroke = 1.7)) +
  geom_point(data = ~filter(.x, Cycle == "Simple"),size = 7, position = position_dodge(0.2), color = "black") +
  geom_point(data = ~filter(.x, Cycle == "Complex"),size = 7, position = position_dodge(0.2)) +
  #geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Number of species") +
  labs(color = "Juvenile \nGeneralism", shape = "Life cycle") +
  ggtitle("Even Resource distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18)) +
  scale_shape_manual(values = c(1,4)) +
  scale_color_manual(values = c(color_palette, "black"))
  
  #+
  #geom_point(df.combined[df.combined$Cycle == "Simple", ], aes(x = Adult.trait, y = Richness, shape = Juvenile.trait, color = Juvenile.trait, stroke = 1.05))


even.plot

# Normal

Total.mean.CLC.norm <- norm$Total.mean.CLC.norm
Total.mean.SLC.norm <- norm$Total.mean.SLC.norm

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


norm.plot <- ggplot(df.combined, aes(x = Adult.trait, y = Richness, shape = Cycle, color = Juvenile.trait, stroke = 1.7)) +
  geom_point(data = ~filter(.x, Cycle == "Simple"),size = 7, position = position_dodge(0.2), color = "black") +
  geom_point(data = ~filter(.x, Cycle == "Complex"),size = 7, position = position_dodge(0.2)) +
  #geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Number of species") +
  labs(color = "Juvenile \nGeneralism", shape = "Life cycle") +
  ggtitle("Normal Resource distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18)) +
  scale_shape_manual(values = c(1,4)) +
  scale_color_manual(values = c(color_palette, "black"))

norm.plot 

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


skew.plot <- ggplot(df.combined, aes(x = Adult.trait, y = Richness, shape = Cycle, color = Juvenile.trait, stroke = 1.7)) +
  geom_point(data = ~filter(.x, Cycle == "Simple"),size = 7, position = position_dodge(0.2), color = "black") +
  geom_point(data = ~filter(.x, Cycle == "Complex"),size = 7, position = position_dodge(0.2)) +
  #geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Number of species") +
  labs(color = "Juvenile \nGeneralism", shape = "Life cycle") +
  ggtitle("Skewed Resource distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18)) +
  scale_shape_manual(values = c(1,4)) +
  scale_color_manual(values = c(color_palette, "black"))

skew.plot

# Binormal

Total.mean.CLC.binorm <- binorm$Total.mean.CLC.binorm
Total.mean.SLC.binorm <- binorm$Total.mean.SLC.binorm

Total.sd.CLC.binorm <-  binorm$Total.sd.CLC.binorm
Total.sd.SLC.binorm <-  binorm$Total.sd.SLC.binorm

x <- as.factor(sigma)



df.CLC <- data.frame(
  Juvenile.trait = rep(x, each = length(x)),
  Adult.trait = rep(x, times = length(x)),
  Richness = as.vector(Total.mean.CLC.binorm),
  sd = as.vector(Total.sd.CLC.binorm),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df.SLC <- data.frame(
  Juvenile.trait = x,
  Adult.trait = x,
  Richness = as.vector(Total.mean.SLC.binorm),
  sd = as.vector(Total.sd.SLC.binorm),
  Cycle = rep("Simple", times = length(x))
)



df.combined <- rbind(df.CLC, df.SLC)


binorm.plot <- ggplot(df.combined, aes(x = Adult.trait, y = Richness, shape = Cycle, color = Juvenile.trait, stroke = 1.7)) +
  geom_point(data = ~filter(.x, Cycle == "Simple"),size = 7, position = position_dodge(0.2), color = "black") +
  geom_point(data = ~filter(.x, Cycle == "Complex"),size = 7, position = position_dodge(0.2)) +
  #geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Number of species") +
  labs(color = "Juvenile \nGeneralism", shape = "Life cycle") +
  ggtitle("Bimodal Normal Resource distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18)) +
  scale_shape_manual(values = c(1,4)) +
  scale_color_manual(values = c(color_palette, "black"))

binorm.plot

# Together:

plot <- list()
plot[[1]] <- even.plot
plot[[2]] <- norm.plot
plot[[3]] <- skew.plot
plot[[4]] <- binorm.plot
wrap_plots(plot)

all.plots <- (even.plot + norm.plot) / (skew.plot + binorm.plot)
all.plots + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A",
  tag_prefix = "(",
  tag_suffix = ")")

#--------------------------------

# Plotting Phenotype Endpoint --------------------

# Even ----------------
# Plotting several runs ---------------------------------------

# Adult = Juvenile sigma

Res <- list()

pdf("plots.even.combined.pdf")

for(s in 1:length(sigma)){
  adu.sigma <- sigma[s]
  
  
  juv.sigma <- adu.sigma
  
  
  last.year.list.even <- data.frame()
  
  for(i in 1:length(even$Total.endpoint.CLC.even)){
    this.run <- even$Total.endpoint.CLC.even[[i]]
    this.run$run <- rep(i, time = nrow(this.run))
    this.run <- this.run[this.run$Adult.gen == adu.sigma, ]
    this.run <- this.run[this.run$Juv.gen == juv.sigma, ]
    last.year.list.even <- rbind(last.year.list.even, this.run)
  }
  
  plot.list.even <- list()
  
  for (i in 1:3){
    
    data <- last.year.list.even[last.year.list.even$run == i, ]
    
    color.palette <- plasma(length(data$Adult_Trait))
    
    plot.list.even[[i]] <- ggplot(data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
      geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
      labs(title = substitute(sigma == value, list(value = adu.sigma)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
      scale_x_continuous(limits = c(-3, 3))+
      scale_y_continuous(limits = c(-3, 3))+
      scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
      theme(aspect.ratio=1) +
      theme_minimal(base_family = "LM Roman 10", base_size = 10)
    
    
  }
  
  
  plots <- wrap_plots(plot.list.even)
  
  
  
  Res[[s]] <- plots + plot_annotation(
    title = 'Even Distribution',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  )+ coord_fixed()
  
  print(Res[[s]])
  

} 

dev.off()

even$Total.endpoint.CLC.even

# Choose Run -------------------------
run <- sample(x = 1:10, size = 1)


last.year.list.even <- even$Total.endpoint.CLC.even[[run]]

#Choose which sigmas

# Same adult sigma 
adu.sigma <- sigma[1]                                    # Choose random: sample(sigma, size = 1)
last.year.list.even.adu <- last.year.list.even[last.year.list.even$Adult.gen == adu.sigma, ]
  
# Same juvenile sigma
juv.sigma <-  sigma[1]                                    #Choose random: sample(sigma, size = 1)
last.year.list.even.juv <- last.year.list.even[last.year.list.even$Juv.gen == juv.sigma, ]

#Plotting different runs

  
# Plotting different sigmas


plot.list.even.adu <- list()
plot.list.even.juv <- list()

A.ids <- unique(last.year.list.even.adu$ID)
J.ids <- unique(last.year.list.even.juv$ID)


for (i in 1:length(sigma)){
  
  adu.data <- last.year.list.even.adu[last.year.list.even.adu$ID == A.ids[i], ]
  juv.data <- last.year.list.even.juv[last.year.list.even.juv$ID == J.ids[i], ]

  
  color.palette <- mako(length(adu.data$Juvenile_Trait))
  
  plot.list.even.adu[[i]] <- ggplot(adu.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  color.palette <- mako(length(juv.data$Adult_Trait))
  
  plot.list.even.juv[[i]] <- ggplot(juv.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
}



combo.plot.list <- list()
midtitle <- textGrob(substitute("Juvenile Generalism" == value, list(value = juv.sigma)), gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.even.adu)*2 + 1)){
  if(i == (length(plot.list.even.adu) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.even.adu) + 1)){
    combo.plot.list[[i]] <- plot.list.even.adu[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.list.even.juv[[i-(length(plot.list.even.juv) + 1)]]
  }
}



layout <- "
ABC
DEF
#G#
HIJ
KLM
"

plots <- wrap_plots(combo.plot.list, design = layout)


plots + plot_annotation(
  title = 'Even Distribution',
  subtitle = substitute("Adult generalism" == value, list(value = adu.sigma)),
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
)+ plot_layout(heights = c(1, 1,  0.4, 1, 1))



# Normal ------------------
# Plotting several runs ---------------------------------------

# Adult = Juvenile sigma

Res <- list()

pdf("plots.norm.combined.pdf")

for(s in 1:length(sigma)){
  adu.sigma <- sigma[s]
  
  juv.sigma <- adu.sigma
  
  
  last.year.list.norm <- data.frame()
  
  for(i in 1:length(norm$Total.endpoint.CLC.norm)){
    this.run <- norm$Total.endpoint.CLC.norm[[i]]
    this.run$run <- rep(i, time = nrow(this.run))
    this.run <- this.run[this.run$Adult.gen == adu.sigma, ]
    this.run <- this.run[this.run$Juv.gen == juv.sigma, ]
    last.year.list.norm <- rbind(last.year.list.norm, this.run)
  }
  
  plot.list.norm <- list()
  
  for (i in 1:3){
    
    data <- last.year.list.norm[last.year.list.norm$run == i, ]
    
    color.palette <- plasma(length(data$Adult_Trait))
    
    plot.list.norm[[i]] <- ggplot(data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
      geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
      labs(title = substitute(sigma == value, list(value = adu.sigma)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
      scale_x_continuous(limits = c(-3, 3))+
      scale_y_continuous(limits = c(-3, 3))+
      scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
      theme_minimal(base_family = "LM Roman 10", base_size = 10)
    
    
  }
  
  
  plots <- wrap_plots(plot.list.norm)
  
  Res[[s]] <- plots + plot_annotation(
    title = 'Normal Distribution',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  )+ coord_fixed()
  print(Res[[s]])
  
  
} 

dev.off()


# Choose Run -------------------------
run <- sample(x = 1:10, size = 1)


last.year.list.norm <- norm$Total.endpoint.CLC.norm[[run]]

#Choose which sigmas

# Same adult sigma 
adu.sigma <- sigma[1]                                    # Choose random: sample(sigma, size = 1)
last.year.list.norm.adu <- last.year.list.norm[last.year.list.norm$Adult.gen == adu.sigma, ]

# Same juvenile sigma
juv.sigma <-  sigma[1]                                    #Choose random: sample(sigma, size = 1)
last.year.list.norm.juv <- last.year.list.norm[last.year.list.norm$Juv.gen == juv.sigma, ]

#Plotting different runs


# Plotting different sigmas


plot.list.norm.adu <- list()
plot.list.norm.juv <- list()

A.ids <- unique(last.year.list.norm.adu$ID)
J.ids <- unique(last.year.list.norm.juv$ID)


for (i in 1:length(sigma)){
  
  adu.data <- last.year.list.norm.adu[last.year.list.norm.adu$ID == A.ids[i], ]
  juv.data <- last.year.list.norm.juv[last.year.list.norm.juv$ID == J.ids[i], ]
  
  
  color.palette <- mako(length(adu.data$Juvenile_Trait))
  
  plot.list.norm.adu[[i]] <- ggplot(adu.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  color.palette <- mako(length(juv.data$Adult_Trait))
  
  plot.list.norm.juv[[i]] <- ggplot(juv.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
}



combo.plot.list <- list()
midtitle <- textGrob(substitute("Juvenile Generalism" == value, list(value = juv.sigma)), gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.norm.adu)*2 + 1)){
  if(i == (length(plot.list.norm.adu) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.norm.adu) + 1)){
    combo.plot.list[[i]] <- plot.list.norm.adu[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.list.norm.juv[[i-(length(plot.list.norm.juv) + 1)]]
  }
}



layout <- "
ABC
DEF
#G#
HIJ
KLM
"


plots <- wrap_plots(combo.plot.list, design = layout)


plots + plot_annotation(
  title = 'Normal Distribution',
  subtitle = substitute("Adult generalism" == value, list(value = adu.sigma)),
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
)+ plot_layout(heights = c(1, 1,  0.4, 1, 1))




# Skewed -------------------
# Plotting several runs ---------------------------------------

# Adult = Juvenile sigma
Res <- list()

pdf("plots.skew.combined.pdf")

for(s in 1:length(sigma)){
  adu.sigma <- sigma[s]
  
  juv.sigma <- adu.sigma
  
  
  last.year.list.skew <- data.frame()
  
  for(i in 1:length(skew$Total.endpoint.CLC.skew)){
    this.run <- skew$Total.endpoint.CLC.skew[[i]]
    this.run$run <- rep(i, time = nrow(this.run))
    this.run <- this.run[this.run$Adult.gen == adu.sigma, ]
    this.run <- this.run[this.run$Juv.gen == juv.sigma, ]
    last.year.list.skew <- rbind(last.year.list.skew, this.run)
  }
  
  plot.list.skew <- list()
  
  for (i in 1:3){
    
    data <- last.year.list.skew[last.year.list.skew$run == i, ]
    
    color.palette <- plasma(length(data$Adult_Trait))
    
    plot.list.skew[[i]] <- ggplot(data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
      geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
      labs(title = substitute(sigma == value, list(value = adu.sigma)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
      scale_x_continuous(limits = c(-3, 3))+
      scale_y_continuous(limits = c(-3, 3))+
      scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
      theme_minimal(base_family = "LM Roman 10", base_size = 10)
    
    
  }
  
  
  plots <- wrap_plots(plot.list.skew)
  
  Res[[s]] <- plots + plot_annotation(
    title = 'Skewed Distribution',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  )+ coord_fixed()
  print(Res[[s]])
  
  
} 

dev.off()



# Choose Run -------------------------
run <- sample(x = 1:10, size = 1)


last.year.list.skew <- skew$Total.endpoint.CLC.skew[[run]]

#Choose which sigmas

# Same adult sigma 
adu.sigma <- sigma[1]                                    # Choose random: sample(sigma, size = 1)
last.year.list.skew.adu <- last.year.list.skew[last.year.list.skew$Adult.gen == adu.sigma, ]

# Same juvenile sigma
juv.sigma <-  sigma[1]                                    #Choose random: sample(sigma, size = 1)
last.year.list.skew.juv <- last.year.list.skew[last.year.list.skew$Juv.gen == juv.sigma, ]

#Plotting different runs


# Plotting different sigmas


plot.list.skew.adu <- list()
plot.list.skew.juv <- list()

A.ids <- unique(last.year.list.skew.adu$ID)
J.ids <- unique(last.year.list.skew.juv$ID)


for (i in 1:length(sigma)){
  
  adu.data <- last.year.list.skew.adu[last.year.list.skew.adu$ID == A.ids[i], ]
  juv.data <- last.year.list.skew.juv[last.year.list.skew.juv$ID == J.ids[i], ]
  
  
  color.palette <- mako(length(adu.data$Juvenile_Trait))
  
  plot.list.skew.adu[[i]] <- ggplot(adu.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  color.palette <- mako(length(juv.data$Adult_Trait))
  
  plot.list.skew.juv[[i]] <- ggplot(juv.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
}



combo.plot.list <- list()
midtitle <- textGrob(substitute("Juvenile Generalism" == value, list(value = juv.sigma)), gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.skew.adu)*2 + 1)){
  if(i == (length(plot.list.skew.adu) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.skew.adu) + 1)){
    combo.plot.list[[i]] <- plot.list.skew.adu[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.list.skew.juv[[i-(length(plot.list.skew.juv) + 1)]]
  }
}



layout <- "
ABC
DEF
#G#
HIJ
KLM
"

plots <- wrap_plots(combo.plot.list, design = layout)


plots + plot_annotation(
  title = 'Skewed Distribution',
  subtitle = substitute("Adult generalism" == value, list(value = adu.sigma)),
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
)+ plot_layout(heights = c(1, 1,  0.4, 1, 1))





# Bi modal Normal -------------------
# Plotting several runs ---------------------------------------

# Adult = Juvenile sigma
Res <- list()

pdf("plots.binorm.combined.pdf")

for(s in 1:length(sigma)){
  adu.sigma <- sigma[s]

  juv.sigma <- adu.sigma
  
  
  last.year.list.binorm <- data.frame()
  
  for(i in 1:length(binorm$Total.endpoint.CLC.binorm)){
    this.run <- binorm$Total.endpoint.CLC.binorm[[i]]
    this.run$run <- rep(i, time = nrow(this.run))
    this.run <- this.run[this.run$Adult.gen == adu.sigma, ]
    this.run <- this.run[this.run$Juv.gen == juv.sigma, ]
    last.year.list.binorm <- rbind(last.year.list.binorm, this.run)
  }
  
  plot.list.binorm <- list()
  
  for (i in 1:3){
    
    data <- last.year.list.binorm[last.year.list.binorm$run == i, ]
    
    color.palette <- plasma(length(data$Adult_Trait))
    
    plot.list.binorm[[i]] <- ggplot(data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
      geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
      labs(title = substitute(sigma == value, list(value = adu.sigma)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
      scale_x_continuous(limits = c(-3, 3))+
      scale_y_continuous(limits = c(-3, 3))+
      scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
      theme_minimal(base_family = "LM Roman 10", base_size = 10)
    
    
  }
  
  
  plots <- wrap_plots(plot.list.binorm)
  
  Res[[s]] <- plots + plot_annotation(
    title = 'Bimodal Normal Distribution',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  )+ coord_fixed()
  print(Res[[s]])
  
  
} 

dev.off()



# Choose Run -------------------------
run <- sample(x = 1:10, size = 1)


last.year.list.binorm <- binorm$Total.endpoint.CLC.binorm[[run]]

#Choose which sigmas

# Same adult sigma 
adu.sigma <- sigma[2]                                    # Choose random: sample(sigma, size = 1)
last.year.list.binorm.adu <- last.year.list.binorm[last.year.list.binorm$Adult.gen == adu.sigma, ]

# Same juvenile sigma
juv.sigma <-  sigma[2]                                    #Choose random: sample(sigma, size = 1)
last.year.list.binorm.juv <- last.year.list.binorm[last.year.list.binorm$Juv.gen == juv.sigma, ]

#Plotting different runs


# Plotting different sigmas


plot.list.binorm.adu <- list()
plot.list.binorm.juv <- list()

A.ids <- unique(last.year.list.binorm.adu$ID)
J.ids <- unique(last.year.list.binorm.juv$ID)


for (i in 1:length(sigma)){
  
  adu.data <- last.year.list.binorm.adu[last.year.list.binorm.adu$ID == A.ids[i], ]
  juv.data <- last.year.list.binorm.juv[last.year.list.binorm.juv$ID == J.ids[i], ]
  
  
  color.palette <- mako(length(adu.data$Juvenile_Trait))
  
  plot.list.binorm.adu[[i]] <- ggplot(adu.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  color.palette <- mako(length(juv.data$Adult_Trait))
  
  plot.list.binorm.juv[[i]] <- ggplot(juv.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
}



combo.plot.list <- list()
midtitle <- textGrob(substitute("Juvenile Generalism" == value, list(value = juv.sigma)), gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.binorm.adu)*2 + 1)){
  if(i == (length(plot.list.binorm.adu) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.binorm.adu) + 1)){
    combo.plot.list[[i]] <- plot.list.binorm.adu[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.list.binorm.juv[[i-(length(plot.list.binorm.juv) + 1)]]
  }
}



layout <- "
ABC
DEF
#G#
HIJ
KLM
"

plots <- wrap_plots(combo.plot.list, design = layout)


plots + plot_annotation(
  title = 'Bimodal Normal Distribution',
  subtitle = substitute("Adult generalism" == value, list(value = adu.sigma)),
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
)+ plot_layout(heights = c(1, 1,  0.4, 1, 1))



#-------------------------------------------------

# Saving data -----------------------------------

save(norm, file = "norm.stat.full.RESULT")
save(even, file = "even.stat.full.RESULT")
save(skew, file = "skew.stat.full.RESULT")
save(binorm, file = "binorm.stat.full.RESULT")





