
# Resource competition Models
# This is the the actual simulation algorithm where competition and evolution occurs. 

# Simple life cycle -------------


resourceCompetitionSLC <- function(popSize, resProp, resFreq, resGen=matrix(c(0.2,0.2),ncol=1, nrow=2), im = 0, 
                                   fmax = 2, kA = 0.5, kJ = 0.5, mutProb=0.0005, mutVar=0.05, time.steps=50000, iniP=0, 
                                   threshold = 0.0005, nmorphs = 1, maxTr = 3, minTr = -3){
  
  pop <- matrix(data = NA, ncol = 3, nrow = nmorphs)                            # Each column in this matrix is one species.
  
  pop[,1] <- popSize                                                            # Number of individuals per species is first row.
  pop[,2] <- iniP                                                               # Phenotype for species is second row.
  
  
  colnames(pop) <- c("Number of indivduals", "Trait", "Proxy")                  # Third row is where f and m are stored to calculate fecundity and maturation.
  
  stats         <- matrix(data = c(0, sum(pop[,1]), 0, nrow(pop), mean(pop[,2]), var(pop[,2])) 
                          , nrow = 1, ncol = 6)                                 # Where we will eventually save our data on number of species etc
  phenotypes <- matrix(data = c(0, popSize, iniP), nrow = 1, ncol = 3)          # Where information on each different species is stored
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Trait")                                        
  
  epsilon <- .Machine$double.eps^10                                             # Added when there is risk of r rounding a number down to 0, very small number
  
  for (t in 1:time.steps){
    
    # Deterministic fecundity proxy alpha ------------------------------------
    
    adults    <- pop                                                            # The adults matrix is only created so it easier to mentally seperate adult and juvenile step, also if the adult population wants to be extracted at end of timestep                                                            
    alphaA    <- NULL                                                           # Will create a matrix with all alpha values, is nullified so previous time step is overwritten
    alphaSumA <- NULL
    
    
    resPropAduMatrix <- matrix(data = resProp, ncol = length(resProp),          # Resource property is made into a matrix here so that matrix calculatsion can be used.
                               nrow = nrow(adults), byrow = T)
    aduTrait <- adults[, 2]
    aduTraitMatrix <- matrix(data = rep(aduTrait, each = length(resProp)),      # Made into matrix so matrix calculcations can be used
                             ncol = length(resProp), nrow = nrow(adults), 
                             byrow = T)
    
    alphaA           <- (1/(sqrt(2*pi*resGen[1,1]^2)))*exp(-(((aduTraitMatrix-resPropAduMatrix)^2)/(2*resGen[1,1])^2)) + epsilon                 # Calculation of individual alpha values, equation 2
    adultAbund       <- adults[,1]                                             
    adultAbundMatrix <- matrix(data = rep(adultAbund, each = length(resProp)), ncol = length(resProp), nrow = nrow(adults), byrow = T)  # Creation of a matrix with population size of each type in the rows
    alphaSumA        <- colSums((alphaA*adultAbundMatrix))                      # Creates the denominator for equation 3. Which is each INDIVIDUALS alphaA added together for each resource type, so each column is one resource, each row is one species. 
    
    
    RdivAlphaSumA       <- resFreq/alphaSumA                                    
    RdivAlphaSumATrans  <- matrix(data = RdivAlphaSumA)                         # Is transformed so that matrix caluclations can be done 
    
    Fec <- alphaA%*%RdivAlphaSumATrans                                          # The result is a vector with each species total energy consumption, equation 4. 
    adults[,3] <- fmax*(Fec/(kA+Fec))                                           # Gives f (fecundity) equation 5a
    
    
    
    # Spawning of offspring -------------------------------------------------
    
    juveniles <- adults                                                         # Create a matrix were we will add juveniles into
    
    juveniles[,1] <- rpois(n = nrow(juveniles),                                 # The number of juveniles is drawn from a poisson distribution with mean f x species abundance
                           lambda = juveniles[,1]*juveniles[,3])                # Since each individual of a species has offspring  
    
    juvenile.pop <- c()
    juvenile.pop <- sum(juveniles[,1])                                          # To extract number of juveniles if that wants to be analysed
    
    
    # Mutation of offspring -------------------------------------------------
    
    probs <- juveniles[,1]/sum(juveniles[, 1])    # Generates probability of morph being mutated based upon number of individuals.
    N.mut <- as.numeric(rbinom(n = 1, size = sum(juveniles[, 1]), prob = mutProb)) # Draws the number of mutations this generation
    
    
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
    
    
    
    # Maturation off offspring -----------------------------------------------------
    
    # The calculation of alphaJ is the same as the caluclation of alphaA see above comments for clarifications
    
    alphaSumJ <- NULL
    alphaJ    <- NULL
    
    resPropJuvMatrix <- matrix(data = resProp, ncol = length(resProp), nrow = nrow(juveniles), byrow = T)
    juvTrait <- juveniles[,2]
    juvTraitMatrix <- matrix(data = rep(juvTrait, each = length(resProp)), ncol = length(resProp), nrow = nrow(juveniles), byrow = T)
    
    alphaJ            <- (1/(sqrt(2*pi*resGen[2,1]^2)))*exp(-(((juvTraitMatrix-resPropJuvMatrix)^2)/(2*resGen[2,1]^2))) + epsilon
    juvenAbund        <- juveniles[,1]
    juvenAbundMatrix  <- matrix(data = rep(juvenAbund, each = length(resProp)), ncol = length(resProp), nrow = nrow(juveniles), byrow = T)
    alphaSumJ         <- colSums(alphaJ*juvenAbundMatrix)                                                                         
    
    
    RdivAlphaSumJ     <- resFreq/alphaSumJ
    RdivAlphaSumJTrans   <- matrix(data = RdivAlphaSumJ)
    
    
    Sur <- alphaJ%*%RdivAlphaSumJTrans
    juveniles[,3] <- (Sur/(kJ+Sur))                                             # Gives m, equation 5b
    
    
    
    
    juveniles[,1] <- rbinom(n = nrow(juveniles) , size = juveniles[,1], prob = juveniles[,3])
    
    
    pop <- juveniles[juveniles[, 1] != 0, , drop = FALSE]                        # all adults die after reproducing, so the new generation is only juveniles, and all rows with zero individuals are removed.
    
    # Adding immigrants ---------------------------------------------------------------------
    
    # Only done when im = 1
    if(im == 1) {
      
        trait <- runif(1, min = minTr, max = maxTr)                             # Creates phenotype of immigrant
        
        if(sum(pop[,2] == trait) == 0) {                                        # Checks whether a exact match of immigrant already exists
          pop <- rbind(pop, c(1, trait, NA))                                    # If no duplicate is found, a new species is added to the pop matrix
        } else{
          same <- which(pop[,2] == trait)
          pop[same,1] <- pop[same,1]+1                                          # If so it is just added to that species
        }
      } 

    
    # extract stats and phenotype ---------------------------------------------
    
    if(nrow(pop) == 0){                                                          # Checks whether population has reached zero, then it breaks the for loop.                                   
      print("Population extinction")
      break
    }
    
    if(sum(which(is.na(pop[,1])) != 0)){                                         # Checks whether population has reached zero, then it breaks the for loop.                                   
      print("Population extinction")
      break
    }
    
    # Extracts the results from this generation into the stats and phenotypes matrices. 
    
      stats <- rbind(stats, c(t, sum(adults[,1]), juvenile.pop, nrow(pop), mean(pop[,2]), var(pop[,2]))) 
      
      pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2])
      phenotypes <- rbind(phenotypes, pStats)
 
    
  }
  
  # Removing last time step
  
  stats <- stats[stats[, 1] != time.steps, , drop = FALSE] 
  phenotypes <- phenotypes[phenotypes[,1] != time.steps, , drop = FALSE]
  
  # Removing any morphs of very low abundance for last time step
  pop <- pop[pop[, 1] > threshold*stats[nrow(stats), 2], , drop = FALSE] 
  
  
  # Re-adding modified last time step
  
  stats <- rbind(stats, c(t, sum(adults[,1]), juvenile.pop, nrow(pop), mean(pop[,2]), var(pop[,2]))) 
  
  pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2])
  phenotypes <- rbind(phenotypes, pStats) 
  
  #return output  ------------------------------------------------------------
  colnames(stats) <- c("year", "Adult Population size", "Juvenile Population Size", "Number of morphs", "mean trait", "var trait")
  rownames(phenotypes) <- NULL
  
  
  
  return(list(stats=stats, phenotypes=phenotypes))                                 #returns both the stats and the phenotype
  
  
}



# Complex life cycle -------------------------

# The simple life cycle function is more thoroughly commented, the same methods are used here, please consult above for explanations.
# Only code that differs is explained here. 

resourceCompetitionCLC <- function(popSize, resProp, resFreq, resGen=matrix(c(0.2,0.2),ncol=1, nrow=2), fmax = 2, 
                                   kA = 0.5, kJ = 0.5,mutProb=0.0005, mutVar=0.05, time.steps=50000, iniPA=0, iniPJ=0, 
                                   threshold = 0.0005, nmorphs = 1, im = 0, maxTr = 3, minTr = -3){
  
  
  pop <- matrix(data = NA, ncol = 4, nrow = nmorphs)                            # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniPA                                                              # Adult trait is stored in second row and 
  pop[,3] <- iniPJ                                                              # juvenile trait in third row.
  
  colnames(pop) <- c("Number of indivduals", "Adult trait", "Juvenile trait", "Proxy")
  
  stats         <- matrix(data = c(0, sum(pop[,1]), 0, nrow(pop), mean(pop[,2]), var(pop[,2]), 
                                   mean(pop[,3]), var(pop[,3])), nrow = 1, ncol = 8)                                                             #Where we will eventually save our stats and phenotypes
  phenotypes <- matrix(data = c(0, popSize, iniPA, iniPJ), nrow = 1, ncol = 4)
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Adult trait", "Juvenile trait")                                        
  
  epsilon <- .Machine$double.eps^10                                             # Added when some number become zero, very small number

  
  
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
    
    probs <- juveniles[,1]/sum(juveniles[, 1])    # Generates probability of morph being mutated based upon number of individuals.
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
        juveniles[mutation.pos[i], 1] <- juveniles[mutation.pos[i], 1] - 1      # Removes the mutated individual from the morph
        
        if(rbinom(n = 1, size = 1, prob = 0.5) == 0){                           # Randomly choose whether adult or juvenile trait gets morphed.
          
          
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
    
    
    
    # Maturation off offspring -----------------------------------------------------
    
    alphaSumJ <- NULL
    alphaJ    <- NULL
    # Survival of juveniles also depends on resource availability
    resPropJuvMatrix <- matrix(data = resProp[2,], ncol = ncol(resProp), nrow = nrow(juveniles), byrow = T)
    juvTrait <- juveniles[,3]
    juvTraitMatrix <- matrix(data = rep(juvTrait, each = ncol(resProp)), ncol = ncol(resProp), nrow = nrow(juveniles), byrow = T)
    
    alphaJ            <- (1/(sqrt(2*pi*resGen[2,1]^2)))*exp(-(((juvTraitMatrix-resPropJuvMatrix)^2)/(2*resGen[2,1]^2))) + epsilon
    juvenAbund        <- juveniles[,1]
    juvenAbundMatrix  <- matrix(data = rep(juvenAbund, each = ncol(resProp)), ncol = ncol(resProp), nrow = nrow(juveniles), byrow = T)  
    alphaSumJ         <- colSums(alphaJ*juvenAbundMatrix)                                                                               
    
    
    RdivAlphaSumJ     <- resFreq[2,]/alphaSumJ
    RdivAlphaSumJTrans   <- matrix(data = RdivAlphaSumJ)
    
    
    Sur <- alphaJ%*%RdivAlphaSumJTrans
    juveniles[,4] <- (Sur/(kJ+Sur)) 
    
    
    
    
    juveniles[,1] <- rbinom(n = nrow(juveniles) , size = juveniles[,1], prob = juveniles[,4])
    
    pop <- juveniles[juveniles[, 1] != 0, , drop = FALSE]                        # all adults die after reproducing, so the new generation is only juveniles, and all rows with zero individuals are removed.
    
    # Adding immigrants ---------------------------------------------------------------------
    

    
    if(im == 1) {
      Atrait  <- runif(1, min = minTr, max = maxTr)
      Jtrait  <- runif(1, min = minTr, max = maxTr)
      
      if(sum(pop[,2] == Atrait & pop[,3] == Jtrait) == 0) {                     # Checks whether a exact match of immigrant already exists, both traits are checked in complex
           pop <- rbind(pop, c(1, Atrait, Jtrait, NA))
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
    
    if(sum(which(is.na(pop[,1])) != 0)){                                         # Checks whether population has reached zero, then it breaks the for loop.                                   
      print("Population extinction")
      break
    }
    
    
    stats <- rbind(stats, c(t, sum(adults[,1]), juvenile.pop, nrow(pop), mean(pop[,2]), var(pop[,2]), 
                              mean(pop[,3]), var(pop[,3]))) 
      
    pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2], pop[,3])
    phenotypes <- rbind(phenotypes, pStats)
      
      
    
    
  }
  
  #Removing last time step
  
  stats <- stats[stats[, 1] != time.steps, , drop = FALSE] 
  phenotypes <- phenotypes[phenotypes[, 1] != time.steps, , drop = FALSE]
  
  
  # Removing any morphs of very low abundance for last time step
  pop <- pop[pop[, 1] > threshold*stats[nrow(stats), 2], , drop = FALSE] 
  
  # Readding modified last time step
  stats <- rbind(stats, c(t, sum(adults[,1]), juvenile.pop, nrow(pop), mean(pop[,2]), var(pop[,2]), 
                          mean(pop[,3]), var(pop[,3]))) 
  
  pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2], pop[,3])
  phenotypes <- rbind(phenotypes, pStats)
  
  #return output  ------------------------------------------------------------
  colnames(stats) <- c("year", "Adult population size","Juvenile Population Size", "Number of morphs", "mean A trait", "var A", "mean J trait", "var J")
  rownames(phenotypes) <- NULL
  
  return(list(stats=stats, phenotypes=phenotypes))   #returns both the stats and the phenotype
  
  
}


















