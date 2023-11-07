
# Biggest difference in this version is that there is only one trait and adults and juvenile share a resource



# Full function ----------------------------------------------------------------

resourceCompetitionSLC <- function(popSize, resProp, resFreq, resGen=matrix(c(0.1,0.1),ncol=1, nrow=2), im = 0.001, 
                                fmax = 10, kA = 2, kJ = 0.2, mutProb=0.001, mutVar=0.1, time.steps=200, iniP=4, nmorphs = 1){
  
  pop <- matrix(data = NA, ncol = 3, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniP

  
  colnames(pop) <- c("Number of indivduals", "Trait", "Proxy")
  
  stats         <- matrix(data = c(0, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2])) 
                                  , nrow = 1, ncol = 5)                                                             #Where we will eventually save our stats and phenotypes
  phenotypes <- matrix(data = c(0, popSize, iniP), nrow = 1, ncol = 3)
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Trait")                                        
  
  epsilon <- .Machine$double.eps^10  #Added when some number become zero, very small number
  posstrait <- seq(from = min(resProp), to = max(resProp), by = mutVar)
  
  for (t in 1:time.steps){
    
    # Deterministic fecundity proxy alpha ------------------------------------
    
    adults    <- pop                                                            
    alphaA    <- NULL                                                           # Will create a matrix with all alpha values
    alphaSumA <- NULL
    
    
    resPropAduMatrix <- matrix(data = resProp, ncol = length(resProp), nrow = nrow(adults), byrow = T)
    aduTrait <- adults[, 2]
    aduTraitMatrix <- matrix(data = rep(aduTrait, each = length(resProp)), ncol = length(resProp), nrow = nrow(adults), byrow = T)
    
    alphaA           <- exp(-(((aduTraitMatrix-resPropAduMatrix)^2)/(2*resGen[1,1])^2)) + epsilon                 # Calculation of individual alpha
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
    
    alphaJ            <- exp(-(((juvTraitMatrix-resPropJuvMatrix)^2)/(2*resGen[2,1]^2))) + epsilon
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
    
    stats <- rbind(stats, c(t, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2]))) 
    
    pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2])
    phenotypes <- rbind(phenotypes, pStats)
    
    
  }
  
  #return output  ------------------------------------------------------------
  colnames(stats) <- c("year", "population size", "Number of morphs", "mean trait", "var trait")
  rownames(phenotypes) <- NULL

  
  return(list(stats=stats, phenotypes=phenotypes))                                 #returns both the stats and the phenotype
  
  
}




resource.frequency <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)           # res. freq. 
resource.property <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)             # res. property 
abundance <- 15000
resource.abundance <- abundance*resource.frequency

output <- resourceCompetition(resProp=resource.property, resFreq=resource.abundance, popSize = 10, mutProb=0.001, mutVar=0.05, time.steps = 10000)

stats <- output$stats
phenotypes <- output$phenotypes



# Plotting  -------------------------------------------------------------------

par(mfrow=c(2,1))


plot(x = stats[, 1], y = stats[, 2], xlab = "Year", ylab = "Total population size", type = "l")

plot(x = stats[, 1], y = stats[, 3], xlab = "Year", ylab = "Number of Phenotypes", type = "l")

par(mfrow=c(1,1))

pColors <- rgb(0.7, 0.1, 0.5, alpha = (phenotypes[,2]/(100+phenotypes[,2])))
plot(x = phenotypes[, 1], y = phenotypes[,3], col = pColors, ylab = "Trait", xlab = "Time", pch = 16)


