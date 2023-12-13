# Complex vs Simple number of species sim background


# Results:  Total.mean(number of species) and Total.mean.abund (Total population )


#Functions:

resourceCompetitionCLC <- function(popSize, resProp, resFreq, resGen=matrix(c(0.15,0.15),ncol=1, nrow=2), fmax = 2, 
                                   kA = 0.5, kJ = 0.5,mutProb=0.001, mutVar=0.1, time.steps=200, iniPA=5, iniPJ=5, 
                                   threshold = 0.005, nmorphs = 1, im = 0){
  
  
  pop <- matrix(data = NA, ncol = 4, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniPA
  pop[,3] <- iniPJ
  
  colnames(pop) <- c("Number of indivduals", "Adult trait", "Juvenile trait", "Proxy")
  
  stats         <- matrix(data = c(0, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2]), 
                                   mean(pop[,3]), var(pop[,3])), nrow = 1, ncol = 7)                                                             #Where we will eventually save our stats and phenotypes
  phenotypes <- matrix(data = c(0, popSize, iniPA, iniPJ), nrow = 1, ncol = 4)
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Adult trait", "Juvenile trait")                                        
  
  epsilon <- .Machine$double.eps^10  #Added when some number become zero, very small number
  possAtrait <- seq(from = min(resProp[1,])-1, to = max(resProp[1,])+1, by = mutVar)  # Used when generating immigrants
  possJtrait <- seq(from = min(resProp[2,])-1, to = max(resProp[2,])+1, by = mutVar)
  
  
  for (t in 1:time.steps){
    # Deterministic fecundity proxy alpha ------------------------------------
    
    adults     <- pop                                                            
    alphaA     <- NULL                                                           # Will create a matrix with all alpha values
    alphaSumA  <- NULL
    
    
    resPropAduMatrix <- matrix(data = resProp[1,], ncol = ncol(resProp), nrow = nrow(adults), byrow = T)
    aduTrait <- adults[, 2]
    aduTraitMatrix <- matrix(data = rep(aduTrait, each = ncol(resProp)), ncol = ncol(resProp), nrow = nrow(adults), byrow = T)
    
    alphaA           <- (1/(sqrt(2*pi*resGen[1,1]^2)))*exp(-(((aduTraitMatrix-resPropAduMatrix)^2)/(2*resGen[1,1])^2)) + epsilon                 # Calculation of individual alpha
    adultAbund       <- adults[,1]
    adultAbundMatrix <- matrix(data = rep(adultAbund, each = ncol(resProp)), ncol = ncol(resProp), nrow = nrow(adults), byrow = T)  # Creation of a matrix with population size of each type in the rows
    alphaSumA        <- colSums((alphaA*adultAbundMatrix))                                                                         # Creation of matrix that reflects both the trait but also number of individuals in type
    
    
    RdivAlphaSumA       <- resFreq[1,]/alphaSumA
    RdivAlphaSumATrans  <- matrix(data = RdivAlphaSumA)
    
    Fec <- alphaA%*%RdivAlphaSumATrans
    adults[,4] <- fmax*(Fec/(kA+Fec)) 
    
    # Spawning of offspring -------------------------------------------------
    
    juveniles <- adults                                                         # Create a matrix were we will add juveniles into
    
    juveniles[,1] <- rpois(n = nrow(juveniles), lambda = juveniles[,1]*juveniles[,4]) 
    
    
    # Mutation of offspring -------------------------------------------------
    
    probs <- juveniles[,1]/sum(juveniles[, 1])  # Generates probability of morph being mutated based upon number of individuals.
    N.mut <- as.numeric(rbinom(n = 1, size = sum(juveniles[, 1]), prob = mutProb))
    
    
    if(N.mut > 0){
      random.choice <- c()
      mutation.pos <- c()
      
      for (m in 1:length(N.mut)){
        random.choice <- rmultinom(n = 1 , size = 1, prob = probs)    # Randomly chooses which morphs to mutate based on probs
        mutation.pos[m] <- as.numeric(which(random.choice == 1))
      }
      
      for (i in  1:length(mutation.pos)){
        
        mutChange <- rnorm(n=1, mean=0, sd=mutVar)
        juveniles[mutation.pos[i], 1] <- juveniles[mutation.pos[i], 1] - 1         # Removes the mutated individual from the morph
        
        if(rbinom(n = 1, size = 1, prob = 0.5) == 0){                            # Randomly choose whether adult or juvenile trait gets morphed.
          
          
          new.morph <- matrix(data = c(1, juveniles[mutation.pos[i], 2] + mutChange, #Changes adult trait to a new trait and adds it to the juveniles
                                       juveniles[mutation.pos[i], 3], 0), 
                              ncol = ncol(juveniles), nrow = 1)
          juveniles <- rbind(juveniles, new.morph)
        }
        else {
          
          new.morph <- matrix(data = c(1, juveniles[mutation.pos[i], 2] ,        #Changes juvenile trait to a new trait and adds it to the juveniles
                                       juveniles[mutation.pos[i], 3]+ mutChange, 0), 
                              ncol = ncol(juveniles), nrow = 1)
          juveniles <- rbind(juveniles, new.morph)
        }
        
      }
      
    }
    
    
    
    # Kill off offspring -----------------------------------------------------
    
    alphaSumJ <- NULL
    alphaJ    <- NULL
    # Survival of juveniles also depends on resource availability
    resPropJuvMatrix <- matrix(data = resProp[2,], ncol = ncol(resProp), nrow = nrow(juveniles), byrow = T)
    juvTrait <- juveniles[,3]
    juvTraitMatrix <- matrix(data = rep(juvTrait, each = ncol(resProp)), ncol = ncol(resProp), nrow = nrow(juveniles), byrow = T)
    
    alphaJ            <- (1/(sqrt(2*pi*resGen[2,1]^2)))*exp(-(((juvTraitMatrix-resPropJuvMatrix)^2)/(2*resGen[2,1]^2))) + epsilon
    juvenAbund        <- juveniles[,1]
    juvenAbundMatrix  <- matrix(data = rep(juvenAbund, each = ncol(resProp)), ncol = ncol(resProp), nrow = nrow(juveniles), byrow = T)  # Creation of a matrix with population size of each type in the rows
    alphaSumJ         <- colSums(alphaJ*juvenAbundMatrix)                                                                         # Creation of matrix that reflects both the trait but also number of individuals in type
    
    
    RdivAlphaSumJ     <- resFreq[2,]/alphaSumJ
    RdivAlphaSumJTrans   <- matrix(data = RdivAlphaSumJ)
    
    
    Sur <- alphaJ%*%RdivAlphaSumJTrans
    juveniles[,4] <- (Sur/(kJ+Sur)) 
    
    
    
    
    juveniles[,1] <- rbinom(n = nrow(juveniles) , size = juveniles[,1], prob = juveniles[,4])
    
    pop <- juveniles[juveniles[, 1] != 0, , drop = FALSE]                        # all adults die after reproducing, so the new generation is only juveniles, and all rows with zero individuals are removed.
    
    # Adding immigrants ---------------------------------------------------------------------
    
    if (runif(1) < im){
      
      Atrait <- sample(x = possAtrait, size = 1)
      Jtrait <- sample(x = possJtrait, size = 1)
      
      if(sum(pop[,2] == Atrait & pop[,3] == Jtrait) == 0) {                   # Checks wheter a exact match of immigrant already exists
        rbind(pop, c(1, Atrait, Jtrait, NA))
      } else{
        same <- which(pop[,2] == Atrait & pop[,3] == Jtrait)
        pop[same,1] <- pop[same,1]+1
      }
      
    }
    
    
    
    # extract stats and phenotype ---------------------------------------------
    
    if(nrow(pop) == 0){                                                          # Checks whether population has reached zero, then it breaks the for loop.                                   
      print("Population extinction")
      break
    }
    
    if(sum(which(is.na(pop[,1])) != 0)){                                                          # Checks whether population has reached zero, then it breaks the for loop.                                   
      print("Population extinction")
      break
    }
    
    
    stats <- rbind(stats, c(t, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2]), 
                            mean(pop[,3]), var(pop[,3]))) 
    
    pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2], pop[,3])
    phenotypes <- rbind(phenotypes, pStats)
    
    
  }
  
  # Removing any morphs of very low abundance
  #pop <- pop[pop[, 1] > threshold*stats[nrow(stats), 2], , drop = FALSE] 
  
  #LastStats <- cbind(t, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2]),  mean(pop[,3]), var(pop[,3]))
  #LastPheno <- cbind(rep(time.steps, nrow(pop)), pop[,1], pop[,2], pop[,3])
  
  #colnames(LastStats) <- c("Year", "Population size", "Number of morphs", "mean A trait", "var A", "mean J trait", "var J")
  #colnames(LastPheno) <- c("Year", "Number of indivduals", "Adult Trait", "Juvenile Trait")  
  
  #return output  ------------------------------------------------------------
  colnames(stats) <- c("year", "population size", "Number of morphs", "mean A trait", "var A", "mean J trait", "var J")
  rownames(phenotypes) <- NULL
  
  return(list(stats=stats, phenotypes=phenotypes))  # LastPheno = LastPheno, LastStats = LastStats Add if we want to remove low abundance morphs                               #returns both the stats and the phenotype
  
  
}


resourceCompetitionSLC <- function(popSize, resProp, resFreq, resGen=matrix(c(0.1,0.1),ncol=1, nrow=2), im = 0.001, 
                                   fmax = 2, kA = 0.5, kJ = 0.5, mutProb=0.001, mutVar=0.1, time.steps=200, iniP=5, 
                                   threshold = 0.005, nmorphs = 1){
  
  pop <- matrix(data = NA, ncol = 3, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniP
  
  
  colnames(pop) <- c("Number of indivduals", "Trait", "Proxy")
  
  stats         <- matrix(data = c(0, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2])) 
                          , nrow = 1, ncol = 5)                                                             #Where we will eventually save our stats and phenotypes
  phenotypes <- matrix(data = c(0, popSize, iniP), nrow = 1, ncol = 3)
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Trait")                                        
  
  epsilon <- .Machine$double.eps^10  #Added when some number become zero, very small number
  posstrait <- seq(from = min(resProp)-1, to = max(resProp)+1, by = mutVar)
  
  for (t in 1:time.steps){
    
    # Deterministic fecundity proxy alpha ------------------------------------
    
    adults    <- pop                                                            
    alphaA    <- NULL                                                           # Will create a matrix with all alpha values
    alphaSumA <- NULL
    
    
    resPropAduMatrix <- matrix(data = resProp, ncol = length(resProp), nrow = nrow(adults), byrow = T)
    aduTrait <- adults[, 2]
    aduTraitMatrix <- matrix(data = rep(aduTrait, each = length(resProp)), ncol = length(resProp), nrow = nrow(adults), byrow = T)
    
    alphaA           <- (1/(sqrt(2*pi*resGen[1,1]^2)))*exp(-(((aduTraitMatrix-resPropAduMatrix)^2)/(2*resGen[1,1])^2)) + epsilon                 # Calculation of individual alpha
    adultAbund       <- adults[,1]
    adultAbundMatrix <- matrix(data = rep(adultAbund, each = length(resProp)), ncol = length(resProp), nrow = nrow(adults), byrow = T)  # Creation of a matrix with population size of each type in the rows
    alphaSumA        <- colSums((alphaA*adultAbundMatrix))                                                                         # Creation of matrix that reflects both the trait but also number of individuals in type
    
    
    RdivAlphaSumA       <- resFreq/alphaSumA
    RdivAlphaSumATrans  <- matrix(data = RdivAlphaSumA)
    
    Fec <- alphaA%*%RdivAlphaSumATrans
    adults[,3] <- fmax*(Fec/(kA+Fec)) 
    
    
    
    # Spawning of offspring -------------------------------------------------
    
    juveniles <- adults                                                         # Create a matrix were we will add juveniles into
    
    
    juveniles[,1] <- rpois(n = nrow(juveniles), lambda = juveniles[,1]*juveniles[,3]) 
    
    # Mutation of offspring -------------------------------------------------
    
    probs <- juveniles[,1]/sum(juveniles[, 1])  # Generates probability of morph being mutated based upon number of individuals.
    N.mut <- as.numeric(rbinom(n = 1, size = sum(juveniles[, 1]), prob = mutProb))
    
    
    if(N.mut > 0){
      random.choice <- c()
      mutation.pos <- c()
      
      for (m in 1:length(N.mut)){
        random.choice <- rmultinom(n = 1 , size = 1, prob = probs)              # Randomly chooses which morphs to mutate based on probs
        mutation.pos[m] <- as.numeric(which(random.choice == 1))
      }
      
      for (i in  1:length(mutation.pos)){
        
        mutChange <- rnorm(n=1, mean=0, sd=mutVar)
        juveniles[mutation.pos[i], 1] <- juveniles[mutation.pos[i], 1] - 1       # Removes the mutated individual from the morph
        
        new.morph <- matrix(data = c(1, juveniles[mutation.pos[i], 2] + mutChange, 0),  #Changes trait to a new trait and adds it to the juveniles
                            ncol = ncol(juveniles), nrow = 1)
        juveniles <- rbind(juveniles, new.morph)
        
      }
    }
    
    
    
    # Kill off offspring -----------------------------------------------------
    
    
    alphaSumJ <- NULL
    alphaJ    <- NULL
    
    resPropJuvMatrix <- matrix(data = resProp, ncol = length(resProp), nrow = nrow(juveniles), byrow = T)
    juvTrait <- juveniles[,2]
    juvTraitMatrix <- matrix(data = rep(juvTrait, each = length(resProp)), ncol = length(resProp), nrow = nrow(juveniles), byrow = T)
    
    alphaJ            <- (1/(sqrt(2*pi*resGen[2,1]^2)))*exp(-(((juvTraitMatrix-resPropJuvMatrix)^2)/(2*resGen[2,1]^2))) + epsilon
    juvenAbund        <- juveniles[,1]
    juvenAbundMatrix  <- matrix(data = rep(juvenAbund, each = length(resProp)), ncol = length(resProp), nrow = nrow(juveniles), byrow = T)  # Creation of a matrix with population size of each type in the rows
    alphaSumJ         <- colSums(alphaJ*juvenAbundMatrix)                                                                         # Creation of matrix that reflects both the trait but also number of individuals in type
    
    
    RdivAlphaSumJ     <- resFreq/alphaSumJ
    RdivAlphaSumJTrans   <- matrix(data = RdivAlphaSumJ)
    
    
    Sur <- alphaJ%*%RdivAlphaSumJTrans
    juveniles[,3] <- (Sur/(kJ+Sur)) 
    
    
    
    
    juveniles[,1] <- rbinom(n = nrow(juveniles) , size = juveniles[,1], prob = juveniles[,3])
    
    
    pop <- juveniles[juveniles[, 1] != 0, , drop = FALSE]                        # all adults die after reproducing, so the new generation is only juveniles, and all rows with zero individuals are removed.
    
    # Adding immigrants ---------------------------------------------------------------------
    
    #if (runif(1) < im){
    
    # trait <- sample(x = posstrait, size = 1)
    
    #  if(sum(pop[,2] == trait) == 0) {                   # Checks wheter a exact match of immigrant already exists
    #   rbind(pop, c(1, trait, NA))
    # } else{
    #   same <- which(pop[,2] == trait)
    #   pop[same,1] <- pop[same,1]+1
    # }
    
    # }
    
    
    # extract stats and phenotype ---------------------------------------------
    
    if(nrow(pop) == 0){                                                          # Checks whether population has reached zero, then it breaks the for loop.                                   
      print("Population extinction")
      break
    }
    
    if(sum(which(is.na(pop[,1])) != 0)){                                                          # Checks whether population has reached zero, then it breaks the for loop.                                   
      print("Population extinction")
      break
    }
    
    stats <- rbind(stats, c(t, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2]))) 
    
    pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2])
    phenotypes <- rbind(phenotypes, pStats)
    
    
  }
  
  #return output  ------------------------------------------------------------
  colnames(stats) <- c("year", "population size", "Number of morphs", "mean trait", "var trait")
  rownames(phenotypes) <- NULL
  
  # Removing any morphs of very low abundance
  pop <- pop[pop[, 1] > threshold*stats[nrow(stats), 2], , drop = FALSE] 
  
  LastStats <- cbind(time.steps, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2])) 
  LastPheno <- cbind(rep(time.steps, nrow(pop)), pop[,1], pop[,2])
  
  colnames(LastStats) <- c("year", "population size", "Number of morphs", "mean trait", "var trait")
  colnames(LastPheno) <- c("Year", "Number of indivduals", "Trait")  
  
  
  return(list(stats=stats, phenotypes=phenotypes, LastPheno = LastPheno, LastStats = LastStats))                                 #returns both the stats and the phenotype
  
  
}



slc.groups <- function(output = outputSLC, threshold = 0.2){
  
  phenodataSLC <- data.frame(
    Year = outputSLC$phenotypes[, 1],
    Trait = outputSLC$phenotypes[, 3],
    Num_Individuals = outputSLC$phenotypes[, 2]
  )
  
  last_year_dataSLC <- phenodataSLC[phenodataSLC$Year == max(phenodataSLC$Year), ]
  
  
  last_year_dataS <- subset(last_year_dataSLC, select = -Year)
  last_year_dataS <- subset(last_year_dataS, select = -Num_Individuals)
  rownames(last_year_dataS) <- NULL
  rownames(last_year_dataSLC) <- NULL
  
  
  distance_matrix <- as.matrix(dist(last_year_dataS[, 1, drop = FALSE], method = "euclidean"))
  
  
  distance_matrix[lower.tri(distance_matrix)] <- NA
  
  
  # Find indices of individuals to keep
  
  
  same <- which(distance_matrix < threshold, arr.ind = T)
  same <- same[same[, 1]-same[,2] != 0, , drop = FALSE]
  rownames(same) <- NULL
  
  
  # Initialize an empty list to store groups
  groups <- list()
  
  # Function to find group index for a species
  find_group <- function(species_id) {
    for (g in seq_along(groups)) {
      if (species_id %in% unlist(groups[[g]])) {
        return(g)
      }
    }
    return(0)
  }
  
  # Iterate over rows in the matrix
  for (s in 1:nrow(same)) {
    species1 <- same[s, 1]
    species2 <- same[s, 2]
    
    # Find groups for each species
    group1 <- find_group(species1)
    group2 <- find_group(species2)
    
    if (group1 == 0 & group2 == 0) {
      # Create a new group
      groups <- c(groups, list(c(species1, species2)))
    } else if (group1 == 0) {
      # Add species1 to the group containing species2
      groups[[group2]] <- c(groups[[group2]], species1)
    } else if (group2 == 0) {
      # Add species2 to the group containing species1
      groups[[group1]] <- c(groups[[group1]], species2)
    } else if (group1 != group2) {
      # Merge two groups
      groups[[group1]] <- c(groups[[group1]], groups[[group2]])
      groups <- groups[-group2]
    }
  }
  
  # Filter out duplicate species in each group
  groups <- lapply(groups, function(group) unique(group))
  
  rownames(last_year_dataSLC) <- NULL
  final_data <- last_year_dataSLC         # Place to store filtered data
  total.sub <- c()                     # Place to store subspecies
  
  #Add population count of "subspecies" to main species
  
  for(q in seq_along(groups)){
    combo <- NULL
    combo <- groups[[q]]
    main <- combo[which.max(final_data[combo,3])]
    sub <- combo[-which.max(final_data[combo,3])]
    final_data[main,3] <- final_data[main,3] + sum(final_data[sub,3])
    total.sub <- rbind(c(total.sub, sub))
    
  }
  # Remove subspecies
  final_data <- final_data[-total.sub, ,drop = FALSE]
  return(final_data)
  
  
}



clc.groups <- function(output = outputCLC, threshold = 0.2){
  
  phenodataCLC <- data.frame(
    Year = outputCLC$phenotypes[, 1],
    Adult_Trait = outputCLC$phenotypes[, 3],
    Juvenile_Trait = outputCLC$phenotypes[, 4],
    Num_Individuals = outputCLC$phenotypes[, 2]
  )
  
  last_year_dataCLC <- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  last_year_dataC <- subset(last_year_dataCLC, select = -Year)
  last_year_dataC <- subset(last_year_dataC, select = -Num_Individuals)
  rownames(last_year_dataCLC) <- NULL
  rownames(last_year_dataC) <- NULL
  
  
  distance_matrix_adult <- as.matrix(dist(last_year_dataC[, 1, drop = FALSE], method = "euclidean"))
  distance_matrix_juvenile <- as.matrix(dist(last_year_dataC[, 2, drop = FALSE], method = "euclidean"))
  
  distance_matrix_adult[lower.tri(distance_matrix_adult)] <- NA
  distance_matrix_juvenile[lower.tri(distance_matrix_juvenile)] <- NA
  
  
  # Set a threshold for similarity (adjust as needed)
  threshold <- 0.2
  
  # Find indices of individuals to keep
  
  
  same <- which(distance_matrix_adult < threshold & distance_matrix_juvenile < threshold, arr.ind = T)
  same <- same[same[, 1]-same[,2] != 0, , drop = FALSE]
  rownames(same) <- NULL
  
  
  # Initialize an empty list to store groups
  groups <- list()
  
  # Function to find group index for a species
  find_group <- function(species_id) {
    for (g in seq_along(groups)) {
      if (species_id %in% unlist(groups[[g]])) {
        return(g)
      }
    }
    return(0)
  }
  
  # Iterate over rows in the matrix
  for (s in 1:nrow(same)) {
    species1 <- same[s, 1]
    species2 <- same[s, 2]
    
    # Find groups for each species
    group1 <- find_group(species1)
    group2 <- find_group(species2)
    
    if (group1 == 0 & group2 == 0) {
      # Create a new group
      groups <- c(groups, list(c(species1, species2)))
    } else if (group1 == 0) {
      # Add species1 to the group containing species2
      groups[[group2]] <- c(groups[[group2]], species1)
    } else if (group2 == 0) {
      # Add species2 to the group containing species1
      groups[[group1]] <- c(groups[[group1]], species2)
    } else if (group1 != group2) {
      # Merge two groups
      groups[[group1]] <- c(groups[[group1]], groups[[group2]])
      groups <- groups[-group2]
    }
  }
  
  # Filter out duplicate species in each group
  groups <- lapply(groups, function(group) unique(group))
  
  final_data <- last_year_dataCLC         # Place to store filtered data
  total.sub <- c()                     # Place to store subspecies
  
  #Add population count of "subspecies" to main species
  
  for(q in seq_along(groups)){
    combo <- NULL
    combo <- groups[[q]]
    main <- combo[which.max(final_data[combo,4])]
    sub <- combo[-which.max(final_data[combo,4])]
    final_data[main,4] <- final_data[main,4] + sum(final_data[sub,4])
    total.sub <- rbind(c(total.sub, sub))
    
  }
  # Remove subspecies
  final_data <- final_data[-total.sub, , drop = FALSE]
  
  return(final_data)
}


sigma <- c(0.15, 0.3, 0.45, 0.6, 0.75)


#Even -------------------------------------------------------------------------

# Resources: 

# Evenly distributed Resources

# SLC:
resource.freq <- rep(1/25, times = 25)                                      # res. freq. 
resource.prop <- c(seq(from = -2.5, to = 2.5, length.out = 25))            # res. property 
abundance <- 20000
resource.abundance <- abundance*resource.freq


# CLC:

resource.property<- c(seq(from = -2.5, to = 2.5, length.out = 25)) 

resource.frequency <- rep(1/25, times = 25)

resource.abundance.adults     <- 20000                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- 20000

resFreqMatrix <- matrix(resource.frequency, nrow=2, ncol=length(resource.frequency), byrow = TRUE)
resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

resPropMatrix <- matrix(resource.property, nrow=2, ncol=length(resource.property), byrow = TRUE) 

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))

rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))


#Simming:


Total.species.SLC.single.even <- c()


Total.species.CLC.even <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
rownames(Total.species.CLC.even) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
colnames(Total.species.CLC.even) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES



# SLC

Total.SLC.list.even <- list()
Total.abund.SLC.list.even <- list()

for(r in 1:10) {
  
  print(paste0("loop ", r, " started"))
  
  Total.species.SLC.single.even <- c()
  Abundance.species.SLC.single.even <- c()
  
  for(i in 1:length(sigma)){
    
    outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
    
    #Filter out similar "species"
    
    final.data.SLC.even <- slc.groups(output = outputSLC)
    Total.species.SLC.single.even[i] <- nrow(final.data.SLC.even)
    Abundance.species.SLC.single.even[i] <- sum(final.data.SLC.even[,3])
  }
  
  Total.SLC.list.even[[r]] <- Total.species.SLC.single.even
  Total.abund.SLC.list.even[[r]] <- Abundance.species.SLC.single.even
}

# Caluclating mean of 10 runs

Total.mean.SLC.even <- Reduce(`+`, Total.SLC.list.even) / length(Total.SLC.list.even)
Total.mean.abund.SLC.even <- Reduce(`+`, Total.abund.SLC.list.even) / length(Total.abund.SLC.list.even) 



# CLC

Total.CLC.list.even <- list()
Total.abund.CLC.list.even <- list()

for(r in 1:10){
  print(paste0("loop ", r, " started"))
  
  Total.species.CLC.even <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.even) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
  colnames(Total.species.CLC.even) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES
  Abundance.species.CLC.even <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.even) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
  colnames(Total.species.CLC.even) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES
  
  for(i in 1:length(sigma)){
    
    for(k in 1:length(sigma)){
      
      outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
      
      #Filter out similar "species"
      
      final.data.CLC.even <- clc.groups(output = outputCLC)
      Total.species.CLC.even[i, k] <- nrow(final.data.CLC.even)
      
      Abundance.species.CLC.even[i, k] <- sum(final.data.CLC.even[, 4])
      
    }
    
  }
  Total.CLC.list.even[[r]] <- Total.species.CLC.even
  Total.abund.CLC.list.even[[r]] <- Abundance.species.CLC.even
}

# Calculating mean of 10 runs

Total.mean.CLC.even <- Reduce(`+`, Total.CLC.list.even) / length(Total.CLC.list.even)
Total.mean.abund.CLC.even <- Reduce(`+`, Total.abund.CLC.list.even) / length(Total.abund.CLC.list.even) 



# Normal -----------------------------------------------------------------------



# Normal resources:

m <- 0 
s <- 1
N.resource.frequency <- c()
N.resource.property<- c(seq(from = -2.5, to = 2.5, length.out = 25)) 

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


# Simming

Total.species.SLC.single.normal <- c()


Total.species.CLC.normal <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
rownames(Total.species.CLC.normal) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
colnames(Total.species.CLC.normal) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES



Total.SLC.list.normal <- list()
Total.abund.SLC.list.normal <- list()

for(r in 1:10) {
  
  print(paste0("loop ", r, " started"))
  
  Total.species.SLC.single.normal <- c()
  Abundance.species.SLC.single.normal <- c()
  
  for(i in 1:length(sigma)){
    
    outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
    
    #Filter out similar "species"
    
    final.data.SLC.normal <- slc.groups(output = outputSLC)
    Total.species.SLC.single.normal[i] <- nrow(final.data.SLC.normal)
    Abundance.species.SLC.single.normal[i] <- sum(final.data.SLC.normal[,3])
  }
  
  Total.SLC.list.normal[[r]] <- Total.species.SLC.single.normal
  Total.abund.SLC.list.normal[[r]] <- Abundance.species.SLC.single.normal
}

# Calculating mean of 10 runs

Total.mean.SLC.normal <- Reduce(`+`, Total.SLC.list.normal) / length(Total.SLC.list.normal)
Total.mean.abund.SLC.normal <- Reduce(`+`, Total.abund.SLC.list.normal) / length(Total.abund.SLC.list.normal) 



# CLC

Total.CLC.list.normal <- list()
Total.abund.CLC.list.normal <- list()

for(r in 1:10){
  print(paste0("loop ", r, " started"))
  
  Total.species.CLC.normal <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.normal) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
  colnames(Total.species.CLC.normal) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES
  Abundance.species.CLC.normal <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.normal) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
  colnames(Total.species.CLC.normal) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES
  
  for(i in 1:length(sigma)){
    
    for(k in 1:length(sigma)){
      
      outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
      
      #Filter out similar "species"
      
      final.data.CLC.normal <- clc.groups(output = outputCLC)
      Total.species.CLC.normal[i, k] <- nrow(final.data.CLC.normal)
      
      Abundance.species.CLC.normal[i, k] <- sum(final.data.CLC.normal[, 4])
      
    }
    
  }
  Total.CLC.list.normal[[r]] <- Total.species.CLC.normal
  Total.abund.CLC.list.normal[[r]] <- Abundance.species.CLC.normal
}

# Caluclating mean of 10 runs

Total.mean.CLC.normal <- Reduce(`+`, Total.CLC.list.normal) / length(Total.CLC.list.normal)
Total.mean.abund.CLC.normal <- Reduce(`+`, Total.abund.CLC.list.normal) / length(Total.abund.CLC.list.normal) 


# Skewed ----------------------------------------------------------------------

# Resources

# SLC:

nr.resources <- 25
tot <- (nr.resources*(nr.resources+1))/2

x <- 1/tot
resource.freq <- c()

for (i in 1:nr.resources){ 
  resource.freq[i] <- i*x 
}                                     # res. freq. 

resource.prop <- c(seq(from = -2.5, to = 2.5, length.out = nr.resources))            # res. property 
abundance <- 20000
resource.abundance <- abundance*resource.freq


# CLC:

resource.property<- c(seq(from = -2.5, to = 2.5, length.out = nr.resources)) 


resource.frequency <- c()
for (i in 1:nr.resources){ 
  resource.frequency[i] <- i*x 
}    


resource.abundance.adults     <- 20000                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- 20000

resFreqMatrix <- matrix(resource.frequency, nrow=2, ncol=length(resource.frequency), byrow = TRUE)
resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

resPropMatrix <- matrix(resource.property, nrow=2, ncol=length(resource.property), byrow = TRUE) 

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))

rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))


# Simming


Total.species.SLC.single.skewed <- c()


Total.species.CLC.skewed <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
rownames(Total.species.CLC.skewed) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
colnames(Total.species.CLC.skewed) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES



Total.SLC.list.skewed <- list()
Total.abund.SLC.list.skewed <- list()

for(r in 1:10) {
  
  print(paste0("loop ", r, " started"))
  
  Total.species.SLC.single.skewed <- c()
  Abundance.species.SLC.single.skewed <- c()
  
  for(i in 1:length(sigma)){
    
    outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
    
    #Filter out similar "species"
    
    final.data.SLC.skewed <- slc.groups(output = outputSLC)
    Total.species.SLC.single.skewed[i] <- nrow(final.data.SLC.skewed)
    Abundance.species.SLC.single.skewed[i] <- sum(final.data.SLC.skewed[,3])
  }
  
  Total.SLC.list.skewed[[r]] <- Total.species.SLC.single.skewed
  Total.abund.SLC.list.skewed[[r]] <- Abundance.species.SLC.single.skewed
}

# Caluclating mean of 10 runs

Total.mean.SLC.skewed <- Reduce(`+`, Total.SLC.list.skewed) / length(Total.SLC.list.skewed)
Total.mean.abund.SLC.skewed <- Reduce(`+`, Total.abund.SLC.list.skewed) / length(Total.abund.SLC.list.skewed) 



# CLC

Total.CLC.list.skewed <- list()
Total.abund.CLC.list.skewed <- list()

for(r in 1:10){
  print(paste0("loop ", r, " started"))
  
  Total.species.CLC.skewed <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.skewed) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
  colnames(Total.species.CLC.skewed) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES
  Abundance.species.CLC.skewed <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.skewed) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
  colnames(Total.species.CLC.skewed) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES
  
  for(i in 1:length(sigma)){
    
    for(k in 1:length(sigma)){
      
      outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
      
      #Filter out similar "species"
      
      final.data.CLC.skewed <- clc.groups(output = outputCLC)
      Total.species.CLC.skewed[i, k] <- nrow(final.data.CLC.skewed)
      
      Abundance.species.CLC.skewed[i, k] <- sum(final.data.CLC.skewed[, 4])
      
    }
    
  }
  Total.CLC.list.skewed[[r]] <- Total.species.CLC.skewed
  Total.abund.CLC.list.skewed[[r]] <- Abundance.species.CLC.skewed
}

# Calculating mean of 10 runs

Total.mean.CLC.skewed <- Reduce(`+`, Total.CLC.list.skewed) / length(Total.CLC.list.skewed)
Total.mean.abund.CLC.skewed <- Reduce(`+`, Total.abund.CLC.list.skewed) / length(Total.abund.CLC.list.skewed) 





