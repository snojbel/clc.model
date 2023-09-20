
# My first attempt at modifying the Schmid code to my "liking"

#Shorthands used: popSize : Intial population size 
                #resProp  : resource property, a matrix with first column being juvenile resources and second column being adult resources
                #resFreq  : resource freq, should be relative i.e. add up to one
                #resGen   : A measure of how generalist the consumers are
                #mutProb  : Mutation probabilty
                #mutVar   : Mutation varitation in the case a mutation occurs
                #years    : how long the function runs i.e. time step
                #iniP     : inital phenotype
                #iniPvar  : If one wants to add variation in inital phenotype this can be added as an argument
                #dispProb : if one wants to add displacement this can be added.


resourceCompetition <- function(popSize,resProp, resFreq, resGen = c(1,1), mutProb=0.001, mutVar=0.1, years=200, iniP=5){
  
  # initialize population ......................................
  pop           <- matrix(NA, ncol=3, nrow=sum(popSize))                        # Creating a matrix where each row is an indivdual. 
  colnames(pop) <- c("Life - stage", "Phenotype", "Fecundity")                  # Their phenotype will be used to determine fecunidty
  
  pop[,1] <- c(rep.int(1,popSize[1]), rep.int(2, popSize[2]))                   # Numbers all the indivduals first column as either juvenile(1) or adult(2)
  pop[,2] <- iniP                                                               # Sets inital phenotype of all individuals
                                                                                # rnorm(n=sum(popSize), mean=iniPmean, sd=0.05), inital population can have one phenotype or several
  
  stats         <- NULL
  phenotype <- matrix(NA, ncol=3, nrow=0)                                       # Creates a matrix to store all phenotype values of all individuals in.
  colnames(phenotype) <- c("year", "Life-stage", "phenotype")
  
  for(t in 1:years){
    
    # compute fecundity proxy, alpha, for adults ------------------
    adults    <- pop[pop[,1]==2,]                                               #This extracts all values in the pop matrix where the life-stage = 2
    alphaSum <- NULL
    alphaM    <- NULL                                                           # Will create a matrix with all alpha values
    for(r in 1:ncol(resProp)){                                                  # Loops through each resource to compute 
      
      rp <- resProp[2,r]                                                        # Resource property for adult resource [2,]
      
      alpha     <- exp(-(adults[,2]-rp)^2/ (2*resGen[2]) )                      # Compute  alphas, if different generalities for stages want to be added modify resGen
      alphaM    <- cbind(alpha, alpha)                                          # Store all alphas
      alphaSum <- c(alphaSum1, sum(alpha) )                                     # Store the sum over all alphas
      
    }
    
    # fec - patch 1 ............................................................
    for(i in 1:nrow(patch1)){
      
      fec <- 0
      for(r in 1:ncol(resProp)){ # r<-1; i<-1
        
        pi  <- resFreq[1,r]
        fec <- fec + pi * (alpha1[i,r]/alphaSum1[r])
        
      }
      patch1[i,3] <- fec
      
    }
    
    # alpha - patch 2 ............................................................
    patch2    <- pop[pop[,1]==2,]
    alphaSum2 <- NULL
    alpha2    <- NULL
    for(r in 1:ncol(resProp)){ # r<-1
      
      alpha     <- exp(-(patch2[,2]-resProp[2,r])^2/ (2*resGen[2,1]) )
      alpha2    <- cbind(alpha2, alpha)
      alphaSum2 <- c(alphaSum2, sum(alpha) )
      
    }
    
    # fec - patch 2 ............................................................
    for(i in 1:nrow(patch2)){
      
      fec <- 0
      for(r in 1:ncol(resProp)){ # r<-1; i<-1
        
        fec <- fec + resFreq[2,r]*alpha2[i,r]/alphaSum2[r]
        
      }
      patch2[i,3] <- fec
      
    }
    
    # create next generation  ----
    patch1new     <- patch1[sample(size=popSize[1], x=1:nrow(patch1),prob=patch1[,3], replace=T),]
    patch1new[,3] <- NA
    
    patch2new     <- patch2[sample(size=popSize[2], x=1:nrow(patch2),prob=patch2[,3], replace=T),]
    patch2new[,3] <- NA
    
    pop <- rbind(patch1new, patch2new)
    
    # mutate new generation   ----
    for(i in 1:nrow(pop)){
      
      if(runif(n=1, min=0, max=1)<=mutProb){
        
        pop[i,2] <- pop[i,2] + rnorm(n=1, mean=0, sd=mutVar)
        # print("mutation!")
        
      }
    }
    
    # dispersal               ----
    dispInds <- runif(n=nrow(pop))<dispProb
    dispInds <- which(dispInds==TRUE)
    if(length(dispInds)>=1){
      for(i in 1:length(dispInds)){
        if(pop[dispInds[i],1]==1){
          pop[dispInds[i],1] <- 2
        } else {
          pop[dispInds[i],1] <- 1
        }
      }
    }
    
    # extract stats           ----
    patch1    <- pop[pop[,1]==1,]
    patch2    <- pop[pop[,1]==2,]
    
    stats <- rbind( stats, c(t, mean(patch1[,2]), var(patch1[,2]), mean(patch2[,2]), var(patch2[,2]))) 
    #extract phenotypes of each individual each year
    phenotypes_patch1 <- cbind(rep(t, nrow(patch1)), rep(1, nrow(patch1)), patch1[,2])
    phenotypes_patch2 <- cbind(rep(t, nrow(patch2)), rep(2, nrow(patch2)), patch2[,2])
    
    phenotype <- rbind(phenotype, phenotypes_patch1, phenotypes_patch2)
  }
  
  # return output stats .............................................
  colnames(stats) <- c("year", "mean1", "var1", "mean2", "var2")
  colnames(phenotype) <- c("year", "patch", "phenotype")
  return(list(stats=stats, phenotype=phenotype))  #returns both the stats and the phenotype
  
}


