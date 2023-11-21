
# Script for running simulations



# Initialization ----------------

# Resources

# SLC:
resource.freq <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)           # res. freq. 
resource.prop <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)             # res. property 
abundance <- 20000
resource.abundance <- abundance*resource.freq


# CLC:

resource.frequency <- c(0.1,  0.1,  0.1,  0.1,  0.1, 0.1,  0.1,  0.1,  0.1,  0.1,   # res. freq. of adults (percentage)
                        0.1,  0.1,  0.1,  0.1,  0.1, 0.1,  0.1,  0.1,  0.1,  0.1)   # res. freq. pf juveniles
resource.property  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,                  # res. property of adults
                        1, 2, 3, 4, 5, 6, 7, 8, 9, 10)                   # res. property of juveniles

resource.abundance.adults     <- 20000                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- 20000

resFreqMatrix <- matrix(resource.frequency, nrow=2, ncol=10, byrow = TRUE)
resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

resPropMatrix <- matrix(resource.property, nrow=2, ncol=10, byrow = TRUE) 

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))

rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))

# Normal resources:

m <- 0 
s <- 1
N.resource.frequency <- c()
N.resource.property<- c(seq(from = -2.5, to = 2.5, length.out = 16)) 

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



resource.abundance.adults     <- 20000                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- 20000

# SLC:

resource.prop <- N.resource.property             # res. property 
abundance <- 20000
resource.abundance <- abundance*N.resource.frequency


# CLC:

resFreqMatrix <- matrix(N.resource.frequency, nrow=2, ncol=length(N.resource.frequency), byrow = TRUE)

resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))


resPropMatrix <- matrix(N.resource.property, nrow=2, ncol=length(N.resource.property), byrow = TRUE) 


rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))


# Model runs with varied sigma:

sigma <- c(0.15, 0.3, 0.45, 0.6, 0.75)

Total_species_SLC_single <- c()

Total_species_SLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
rownames(Total_species_SLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
colnames(Total_species_SLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES

Total_species_CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
rownames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
colnames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES


# SLC:

for(i in 1:length(sigma)){
  
  outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
  
  #Filter out similar "species"
  
  final_data_SLC <- slc.groups(output = outputSLC)
  Total_species_SLC_single[i] <- nrow(final_data_SLC)
}


# If adults and juveniles can have different niche width

for(i in 1:length(sigma)){
  print(paste0("loop", i, "started"))
  for(k in 1:length(sigma)){
     
    outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
    
    #Filter out similar "species"
    
    final_data_SLC <- slc.groups(output = outputSLC)
    Total_species_SLC[i, k] <- nrow(final_data_SLC)
  }
 
}




# CLC:


for(i in 1:length(sigma)){
  print(paste0("loop", i, "started"))
  for(k in 1:length(sigma)){
    
    outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
    
    #Filter out similar "species"
    
    final_data_CLC <- clc.groups(output = outputCLC)
    Total_species_CLC[i, k] <- nrow(final_data_CLC)
  }
  
}

# 10 runs: ---------------------------------

# SLC, same sigma for adults and juv

Total_SLC_list <- list()


for(r in 1:10) {
  
  print(paste0("loop ", r, " started"))
  
  Total_species_SLC_single <- c()

  for(i in 1:length(sigma)){
    
    outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
    
    #Filter out similar "species"
    
    final_data_SLC <- slc.groups(output = outputSLC)
    Total_species_SLC_single[i] <- nrow(final_data_SLC)
  }
  
  Total_SLC_list[[r]] <- Total_species_SLC_single
}

# Caluclating mean of 10 runs

Total_mean_SLC <- Reduce(`+`, Total_SLC_list) / length(Total_SLC_list)




# CLC

Total_CLC_list <- list()


for(r in 1:10){
  print(paste0("loop ", r, " started"))
  
  Total_species_CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
  colnames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES
  
  for(i in 1:length(sigma)){
    
    for(k in 1:length(sigma)){
      
      outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
      
      #Filter out similar "species"
      
      final_data_CLC <- clc.groups(output = outputCLC)
      Total_species_CLC[i, k] <- nrow(final_data_CLC)
    }
    
  }
  Total_CLC_list[[r]] <- Total_species_CLC
}

# Caluclating mean of 10 runs

Total_mean_CLC <- Reduce(`+`, Total_CLC_list) / length(Total_CLC_list)

# Saving results

saveRDS(Total_CLC_list, file = "CLC1Alist.RData")

saveRDS(Total_SLC_list, file = "SLC1Alist.RData")
