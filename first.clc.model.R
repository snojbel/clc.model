
# My first attempt at modifying the Schmid code to my CLC model

#Shorthands used: popSize  : Intial population size 
                # resProp  : resource property, a matrix with first row being juvenile resources and second row being adult resources
                # resFreq  : resource freq, should be relative i.e. add up to one once again first row is juvenile, second is adult. 
                # resGen   : A measure of how generalist the consumers are a matrix with 1 column 2 rows, first row is juvenile
                # fmax     : A measure of the maximum number of offspring possible
                # kA/kJ    : A measure of how density dependent the fecundity and survival is. Half saturation constant 
                # mutProb  : Mutation probabilty
                # mutVar   : Mutation varitation in the case a mutation occurs
                # years    : how long the function runs i.e. time step
                # iniP(A/J): inital phenotype of adult and juvenile(different)
                # iniPvar  : If one wants to add variation in inital phenotype this can be added as an argument
                # dispProb : if one wants to add displacement this can be added.


# Note for myself, there might not be a need for the "life-stage" column since its done in sequence? then again you might want to start 
# with some juveniles as well as adults. 



resourceCompetition <- function(popSize, resProp, resFreq, resGen=matrix(c(0.1,0.1),ncol=1, nrow=2), fmax = 10, kA = 0.5, kJ = 0.2, mutProb=0.001, mutVar=0.1, years=200, iniPA=4, iniPJ=4){
  
  # initialize population ......................................
  pop           <- matrix(NA, ncol=3, nrow=sum(popSize))                        # Creating a matrix where each row is an indivdual. 
  colnames(pop) <- c("phenotypeA", "phenotypeJ", 
                     "fecundity(A)/survival(J)")                                # Their phenotype will be used to determine fecunidty

  pop[,1] <- iniPA                                                              # Sets inital phenotype of all individuals
  pop[,2] <- iniPJ                                                              # rnorm(n=sum(popSize), mean=iniPmean, sd=0.05), inital population can have one phenotype or several
  
  stats         <- NULL
  phenotype <- matrix(NA, ncol=3, nrow=0)                                       # Creates a matrix to store all phenotype values of all individuals in.
  colnames(phenotype) <- c("year", "phenotype A", "Phenotype J")
  
  for(t in 1:years){
    
    # compute fecundity proxy, alpha, for adults ------------------
    adults    <- pop                                                            #This could be used to extracts values in the pop matrix if we distinguish between something.
    alphaSumA <- NULL
    alphaA    <- NULL                                                           # Will create a matrix with all alpha values
    for(r in 1:ncol(resProp)){                                                  # Loops through each resource to compute 
      
      rp <- resProp[2,r]                                                        # Resource property for adult resource [2,]
      
      alpha     <- exp(-(adults[, 1]-rp)^2/ (2*resGen[2,1]) )                   # Compute  alphas for all adults, if different generalities for stages want to be added modify resGen
      alphaA    <- cbind(alphaA, alpha)                                         # Store alphas for all resources for all individuals 
      alphaSumA <- c(alphaSumA, sum(alpha) )                                    # Store the sum of all alphas for each resource
      
    }
    
    # Computing deterministic fecundity (adults) ............................................................
    for(i in 1:nrow(adults)){
      
      fec <- 0
      for(r in 1:ncol(resProp)){                                                # Goes through each resource and calculates source specific  
                                                                                # fecundity and then adds them together to create the total fec..
        R  <- resFreq[2,r]
        fec <- fec + R * (alphaA[i,r]/alphaSumA[r])                             # Total fecundity is the sum of fecundities gained 
                                                                                # from all the different resources. Which is derived
      }                                                                         # both the ability to consume a resource, alphaM, and the 
      adults[i,3] <- fmax* (fec/(kA+fec))                                        # the total ability of all individuals to consume 
                                                                                # said resource, alphaSum. R = resource available
    }                                                                           # the kA business is there to create a saturating function.
    
    
    # Create juveniles based fecundity -------------
    
    juveniles <- matrix(data = NA, nrow = 0, ncol = ncol(pop))                  # Create an empty matrix where I will save the new juveniles
  
    for(i in 1:nrow(adults)) {
      N.juv         <- rpois(n = 1, lambda = adults[i,3])                       # This calculates the number of juveniles based on a possion distribution.
      if(N.juv>0){
      new.juveniles <- matrix(data = c(adults[i,1], adults[i,2], 0), 
                       nrow = N.juv, ncol = ncol(pop), byrow = T)               # Create N.juv rows in an entry with the parents phenotypes
      juveniles     <- rbind(juveniles, new.juveniles)                          # adds together all juveniles created.
      }                                                                          
    }
    
    # Mutate new generation   --------------------
    
    # Mutation of adult trait
    for(i in 1:nrow(juveniles)){
      
      if(runif(n=1, min=0, max=1)<=mutProb){                                    # Currently the chance of mutation is drawn from a random sample
        
        juveniles[i,1] <- juveniles[i,1] + rnorm(n=1, mean=0, sd=mutVar)
        # print("mutation!")
        
      }
    
    
    }
    # Mutation of juvenile trait
    for(i in 1:nrow(juveniles)){
      
      if(runif(n=1, min=0, max=1)<=mutProb){                                    # Currently the chance of mutation is drawn from a random sample
        
        juveniles[i,2] <- juveniles[i,2] + rnorm(n=1, mean=0, sd=mutVar)
        # print("mutation!")
        
      }
      
      
    }
    
    
    # Survival proxy psi?(currently alpha) - juveniles ............................................................
    # juveniles   <- pop
    alphaSumJ <- NULL
    alphaJ    <- NULL
    for(r in 1:ncol(resProp)){                                                  # Survival of juveniles also depends on resource availability.
      
      rp <- resProp[1,r] 
      
      alpha     <- exp(-(juveniles[,2]-rp)^2/ (2*resGen[1,1]) )
      alphaJ    <- cbind(alphaJ, alpha)
      alphaSumJ <- c(alphaSumJ, sum(alpha) )
      
    }
    
    # Survival probability, computed similarly to fecundity .....................................................
    for(i in 1:nrow(juveniles)){
      
      sur <- 0
      for(r in 1:ncol(resProp)){                                                # Goes through each resource and calculates source specific  
                                                                                # survival and then adds them together to create the total sur
        R  <- resFreq[1,r]
        sur <- sur + R * (alphaJ[i,r]/alphaSumJ[r])                             
                                                                                
      }                                                                         
      juveniles[i,3] <- (sur/(kJ+sur))  

    }
    
    # Create next generation of adults based on juvenile survival probability

    for(i in 1:nrow(juveniles)){
      juveniles[i, 3] <- rbinom(n = 1, size = 1, juveniles[i, 3])
    }
    
    pop <- juveniles[juveniles[,3] == 1, ] 
    
    
    # extract stats           ----
    adults.stats       <- adults
    juveniles.stats    <- juveniles
    
    stats <- rbind( stats, c(t,nrow(juveniles.stats), mean(adults.stats[,1]), 
                             var(adults.stats[,1]), mean(juveniles.stats[,2]), 
                             var(juveniles.stats[,2]))) 
    #extract phenotypes of each individual each year
    phenotype.stats <- cbind(rep(t, nrow(juveniles)), juveniles[,1], juveniles[,2])
    phenotype <- rbind(phenotype, phenotype.stats)
    
    if(length(pop) == 0){                                                       # Checks wheter population has reached zero, then it breaks the for loop.                                   
      break
      }
                                                               
  }
  
  # return output stats .............................................
  colnames(stats) <- c("year", "population size", "mean A", "var A", "mean J", "var J")
  colnames(phenotype) <- c("year", "Phenotype A", "Phenotype J")
  return(list(stats=stats, phenotype=phenotype))  #returns both the stats and the phenotype
  
}

# first test run                        ----
resource.frequency <- c(0.2,0.2,0.2,0.2,0.2,  #res. freq. of juveniles
                        0.2,0.2,0.2,0.2,0.2)  #res. freq. pf adults
resource.property<- c(2, 3, 4, 5, 6,  #res. property of juveniles
                      2, 3, 4, 5, 6) #res. property pf adults

resFreqMatrix <- matrix(resource.frequency, nrow=2, ncol=5, byrow = TRUE) 
rownames(resFreqMatrix) <- c("Juvenile", "Adult")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))

resPropMatrix <- matrix(resource.property, nrow=2, ncol=5, byrow = TRUE)          
rownames(resPropMatrix)<-c("Juvenile", "Adult")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))


output <- resourceCompetition(resProp=resPropMatrix, resFreq=resFreqMatrix, popSize = 100, mutProb=0.005, mutVar=0.002, years=100)

stats <- output$stats
phenotypes <- output$phenotypes

