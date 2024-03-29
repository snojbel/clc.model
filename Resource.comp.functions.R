
# Resource competition Models

# Simple -------------


resourceCompetitionSLC <- function(popSize, resProp, resFreq, resGen=matrix(c(0.2,0.2),ncol=1, nrow=2), im = 0, 
                                   fmax = 2, kA = 0.5, kJ = 0.5, mutProb=0.0005, mutVar=0.05, time.steps=50000, iniP=0, 
                                   threshold = 0.005, nmorphs = 1, maxTr = 3, minTr = -3){
  
  pop <- matrix(data = NA, ncol = 3, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniP
  
  
  colnames(pop) <- c("Number of indivduals", "Trait", "Proxy")
  
  stats         <- matrix(data = c(0, sum(pop[,1]), 0, nrow(pop), mean(pop[,2]), var(pop[,2])) 
                          , nrow = 1, ncol = 6)                                                             #Where we will eventually save our stats and phenotypes
  phenotypes <- matrix(data = c(0, popSize, iniP), nrow = 1, ncol = 3)
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Trait")                                        
  
  epsilon <- .Machine$double.eps^10  #Added when some number become zero, very small number
  
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
    
    juvenile.pop <- c()
    juvenile.pop <- sum(juveniles[,1])   # To extract number of juveniles
    
    
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
    
    num.of.im <- im*0.05*sum(pop[,1])
    
    for(m in 1:num.of.im){
      trait <- runif(1, min = minTr, max = maxTr)
      
      if(sum(pop[,2] == trait) == 0) {                   # Checks whether a exact match of immigrant already exists
        rbind(pop, c(1, trait, NA))
      } else{
        same <- which(pop[,2] == trait)
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
    
    stats <- rbind(stats, c(t, sum(adults[,1]), juvenile.pop, nrow(pop), mean(pop[,2]), var(pop[,2]))) 
    
    pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2])
    phenotypes <- rbind(phenotypes, pStats)
    
    
  }
  
  #return output  ------------------------------------------------------------
  colnames(stats) <- c("year", "Adult Population size", "Juvenile Population Size", "Number of morphs", "mean trait", "var trait")
  rownames(phenotypes) <- NULL
  
  # Removing any morphs of very low abundance
  #pop <- pop[pop[, 1] > threshold*stats[nrow(stats), 2], , drop = FALSE] 
  
  #LastStats <- cbind(time.steps, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2])) 
  #LastPheno <- cbind(rep(time.steps, nrow(pop)), pop[,1], pop[,2])
  
  #colnames(LastStats) <- c("year", "population size", "Number of morphs", "mean trait", "var trait")
  #colnames(LastPheno) <- c("Year", "Number of indivduals", "Trait")  
  
  
  return(list(stats=stats, phenotypes=phenotypes))                                 #returns both the stats and the phenotype
  
  
}



# Complex -------------------------


resourceCompetitionCLC <- function(popSize, resProp, resFreq, resGen=matrix(c(0.2,0.2),ncol=1, nrow=2), fmax = 2, 
                                   kA = 0.5, kJ = 0.5,mutProb=0.0005, mutVar=0.05, time.steps=50000, iniPA=0, iniPJ=0, 
                                   threshold = 0.005, nmorphs = 1, im = 0, maxTr = 3, minTr = -3){
  
  
  pop <- matrix(data = NA, ncol = 4, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniPA
  pop[,3] <- iniPJ
  
  colnames(pop) <- c("Number of indivduals", "Adult trait", "Juvenile trait", "Proxy")
  
  stats         <- matrix(data = c(0, sum(pop[,1]), 0, nrow(pop), mean(pop[,2]), var(pop[,2]), 
                                   mean(pop[,3]), var(pop[,3])), nrow = 1, ncol = 8)                                                             #Where we will eventually save our stats and phenotypes
  phenotypes <- matrix(data = c(0, popSize, iniPA, iniPJ), nrow = 1, ncol = 4)
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Adult trait", "Juvenile trait")                                        
  
  epsilon <- .Machine$double.eps^10  #Added when some number become zero, very small number

  
  
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
    
    juvenile.pop <- c()
    juvenile.pop <- sum(juveniles[,1])   # To extract number of juveniles
    
    
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
    
    num.of.im <- im*0.05*sum(pop[,1])
    
    for(m in 1:num.of.im) {
    
      Atrait <- trait <- runif(1, min = minTr, max = maxTr)
      Jtrait <- trait <- runif(1, min = minTr, max = maxTr)
      
      if(sum(pop[,2] == Atrait & pop[,3] == Jtrait) == 0) {                   # Checks whether a exact match of immigrant already exists
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
    
    
    stats <- rbind(stats, c(t, sum(adults[,1]), juvenile.pop, nrow(pop), mean(pop[,2]), var(pop[,2]), 
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
  colnames(stats) <- c("year", "Adult population size","Juvenile Population Size", "Number of morphs", "mean A trait", "var A", "mean J trait", "var J")
  rownames(phenotypes) <- NULL
  
  return(list(stats=stats, phenotypes=phenotypes))  # LastPheno = LastPheno, LastStats = LastStats Add if we want to remove low abundance morphs                               #returns both the stats and the phenotype
  
  
}


















