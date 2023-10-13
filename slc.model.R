
# Biggest difference in this version is that there is only one trait and adults and juvenile share a resource

popSize <- 15                                                                   # As defined below of each morph!
iniP   <- 4
nmorphs <- 1

pop <- matrix(data = NA, ncol = 3, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.

pop[,1] <- popSize
pop[,2] <- iniP

colnames(pop) <- c("Number of indivduals", "Trait", "Proxy")



# Full function ----------------------------------------------------------------

resourceCompetition <- function(popSize, resProp, resFreq, resGen=matrix(c(0.2,0.2),ncol=1, nrow=2), 
                                fmax = 10, kA = 2, kJ = 0.2, mutProb=0.001, mutVar=0.1, time.steps=200, iniP=4, nmorphs = 1){
  
  pop <- matrix(data = NA, ncol = 3, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniP

  
  colnames(pop) <- c("Number of indivduals", "Trait", "Proxy")
  
  stats         <- matrix(data = c(0, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2])) 
                                  , nrow = 1, ncol = 5)                                                             #Where we will eventually save our stats and phenotypes
  phenotypes <- matrix(data = c(0, popSize, iniP), nrow = 1, ncol = 3)
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Trait")                                        
  
  
  for (t in 1:time.steps){
    
    # Deterministic fecundity proxy alpha ------------------------------------
    
    adults    <- pop                                                            
    alphaA    <- NULL                                                           # Will create a matrix with all alpha values
    alphaSumA <- NULL
    
    for(r in 1:length(resProp)){                                                  # Loops through each resource to compute 
      
      rp <- resProp[r]                                                            # Resource property 
      
      alpha     <- exp(-(adults[, 2]-rp)^2/ (2*resGen[1,1])^2 )                 # Compute  alphas for all adults, if different generalities for stages want to be added modify resGen
      alphaA    <- cbind(alphaA, alpha)                                         # Store alphas for all resources for all morphs
      alphaSumA <- c(alphaSumA, sum(alpha*adults[,1]))                          # Store the sum of all alphas for each resource
      
    }
    
    for(i in 1:nrow(adults)){
      
      fec <- 0
      for(r in 1:length(resProp)){                                              # Goes through each resource and calculates source specific  
                                                                                # fecundity and then adds them together to create the total fec..
        R  <- resFreq[r]
        fec <- fec + R * (alphaA[i,r]/alphaSumA[r])                             # Total fecundity is the sum of fecundities gained 
                                                                                # from all the different resources. Which is derived
      }                                                                         # both the ability to consume a resource, alphaA, and the 
      adults[i,3] <- fmax * (fec/(kA+fec))                                      # the total ability of all individuals to consume 
                                                                                # resource, alphaSumA. R = resource available
    }
    
    # Spawning of offspring -------------------------------------------------
    
    juveniles <- adults                                                         # Create a matrix were we will add juveniles into
    
    for(i in 1:nrow(juveniles)) {
      
      juv             <- rpois(n = 1, lambda = adults[i,1] * adults[i,3])       # This calculates the number of juveniles based on a possion distribution. 
      juveniles[i, 1] <- juv                                                    # States number of juveniles in 1 column
      
    }
    
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
    
    for(r in 1:length(resProp)){                                                  # Survival of juveniles also depends on resource availability.
      
      rp <- resProp[r] 
      
      alpha     <- exp(-(juveniles[,2]-rp)^2/ (2*resGen[2,1])^2 )
      alphaJ    <- cbind(alphaJ, alpha)
      alphaSumJ <- c(alphaSumJ, sum(alpha) )
      
    }
    
    for(i in 1:nrow(juveniles)){
      
      sur <- 0
      for(r in 1:length(resProp)){                                                # Goes through each resource and calculates source specific  
                                                                                # survival and then adds them together to create the total sur
        R  <- resFreq[r]
        sur <- sur + R * (alphaJ[i,r]/alphaSumJ[r])                             
        
      }                                                                         
      
      det.sur <- sur/(kJ+sur)                                                   # this is to create a saturating function, survival can never be over 1 higher kJ means flatter curve
      juveniles[i,1] <- rbinom(n = 1 , size = juveniles[i,1], prob = det.sur)   # sees how many survive given their det. survival.
    }
    
    pop <- juveniles[juveniles[, 1] != 0, , drop = FALSE]                        # all adults die after reproducing, so the new generation is only juveniles, and all rows with zero individuals are removed.
    
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




resource.frequency <- c(100, 20, 20, 10, 50)           # res. freq. 
resource.property <- c(1, 3, 4, 5, 10)             # res. property 



output <- resourceCompetition(resProp=resource.property, resFreq=resource.frequency, popSize = 10, mutProb=0.001, mutVar=0.05, time.steps = 500)

stats <- output$stats
phenotypes <- output$phenotypes



# Plotting  -------------------------------------------------------------------

par(mfrow=c(2,1))


plot(x = stats[, 1], y = stats[, 2], xlab = "Year", ylab = "Total population size", type = "l")

plot(x = stats[, 1], y = stats[, 3], xlab = "Year", ylab = "Number of Phenotypes", type = "l")

par(mfrow=c(2,1))

pColors <- rgb(0.7, 0.1, 0.5, alpha = (phenotypes[,2]/(100+phenotypes[,2])))
plot(x = phenotypes[, 1], y = phenotypes[,3], col = pColors, ylab = "Trait", xlab = "Time", pch = 16)


