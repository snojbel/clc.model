



# Second CLC model, optimized, maybe

#Shorthands used: popSize   : Intial population size 
                # resProp   : resource property, a matrix with first row being Adult resources and second row being juveniles resources
                # resFreq   : resource freq, should be relative i.e. add up to one once again first row is adult, second is juvenile. 
                # resGen    : A measure of how generalist the consumers are a matrix with 1 column 2 rows, first row is adults
                # fmax      : A measure of the maximum number of offspring possible
                # kA/kJ     : A measure of how density dependent the fecundity and survival is. Half saturation constant 
                # mutProb   : Mutation probabilty
                # mutVar    : Mutation varitation in the case a mutation occurs
                # time.steps: how long the function runs e.g. years
                # iniP(A/J) : inital phenotype of adult and juvenile(different)
                # nmorphs   : Number of morphs in intial run


resourceCompetition <- function(popSize, resProp, resFreq, resGen=matrix(c(0.1,0.1),ncol=1, nrow=2), fmax = 10, kA = 0.5, kJ = 0.2, mutProb=0.001, mutVar=0.1, time.steps=200, iniPA=4, iniPJ=4, nmorphs = 1){
  
  pop <- matrix(data = NA, ncol = 4, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniPA
  pop[,3] <- iniPJ
  
  colnames(pop) <- c("Number of indivduals", "Adult trait", "Juvenile trait", "Proxy")
  
  stats         <- NULL
  phenotype <- matrix(NA, ncol=4, nrow=0)                                        # Creates a matrix to store all phenotype values of all individuals in.
  colnames(phenotype) <- c("year","Number of indivduals", "Phenotype A", "Phenotype J")
  
  for (t in 1:time.steps){
    
    # Deterministic fecundity proxy alpha ------------------------------------

      adults    <- pop                                                            
      alphaA    <- NULL                                                           # Will create a matrix with all alpha values
      alphaSumA <- NULL
      
      for(r in 1:ncol(resProp)){                                                  # Loops through each resource to compute 
        
        rp <- resProp[1,r]                                                        # Resource property for adult resource [2,]
        
        alpha     <- exp(-(adults[, 2]-rp)^2/ (2*resGen[1,1]) )                   # Compute  alphas for all adults, if different generalities for stages want to be added modify resGen
        alphaA    <- cbind(alphaA, alpha)                                         # Store alphas for all resources for all individuals 
        alphaSumA <- c(alphaSumA, sum(alpha) )                                    # Store the sum of all alphas for each resource
        
      }
      
      for(i in 1:nrow(adults)){
        
        fec <- 0
        for(r in 1:ncol(resProp)){                                                # Goes through each resource and calculates source specific  
                                                                                  # fecundity and then adds them together to create the total fec..
          R  <- resFreq[1,r]
          fec <- fec + R * (alphaA[i,r]/alphaSumA[r])                             # Total fecundity is the sum of fecundities gained 
          # from all the different resources. Which is derived
        }                                                                         # both the ability to consume a resource, alphaA, and the 
        adults[i,4] <- fmax * (fec/(kA+fec))                                       # the total ability of all individuals to consume 
                                                                                  # said resource, alphaSumA. R = resource available
      }
      
      # Spawning of offspring -------------------------------------------------
      
      juveniles <- adults                                                         # Create a matrix were we will add juveniles into
      
      for(i in 1:nrow(juveniles)) {
        
        juv             <- rpois(n = adults[i,1], lambda = adults[i,3])           # This calculates the number of juveniles based on a possion distribution.
        juveniles[i, 1] <- juv                                                    # States number of juveniles in 1 column
                                   
      }
      
      # Mutation of offspring -------------------------------------------------
      
      muts  <- rbinom(n = 1, size = sum(juveniles[, 1]), prob = mutProb)         # Generates number of mutations
      probs <- juveniles[,1]/sum(juveniles[, 1])                                 # Generates probability of morph being mutated based upon number of indivduals. 
      mutation.pos <- rmultinom(n = muts, size = nrow(juveniles), prob = probs)  # Randomly chooses which morphs to mutate based on probs
      
                                      
      for (i in  1:length(mutation.pos)){
        
        mutChange <- rnorm(n=1, mean=0, sd=mutVar)
        
        if(rbinom(n = 1, size = 1, prob = 0.5) == 0){                            # Randomly choose whether adult or juvenile trait gets morphed.
          
          juveniles[mutation.pos[i], 1] <- juveniles[mutation.pos[i], 1] - 1     # Removes the mutated individual from the morph
          new.morph <- matrix(data = c(1, juveniles[mutation.pos[i], 2] + mutChange, #Changes adult trait to a new trait and adds it to the juveniles
                                       1, juveniles[mutation.pos[i], 3], 0), 
                              ncol = ncol(juveniles), nrow = 1)
          juveniles <- rbind(juveniles, new.morph)
        }
        else {
          juveniles[mutation.pos[i], 1] <- juveniles[mutation.pos[i], 1] - 1     # Removes the mutated individual from the morph
          new.morph <- matrix(data = c(1, juveniles[mutation.pos[i], 2] ,        #Changes juvenile trait to a new trait and adds it to the juveniles
                                       1, juveniles[mutation.pos[i], 3]+ mutChange, 0), 
                              ncol = ncol(juveniles), nrow = 1)
          juveniles <- rbind(juveniles, new.morph)
        }
        
      }
      
      # Kill of offspring -----------------------------------------------------
      
      
    
      
  }
  
  
  
}





