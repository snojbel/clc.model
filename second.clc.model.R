



# Second CLC model, optimized, maybe

#Shorthands used: popSize   : Intial population size 
                # resProp   : resource property, a matrix with first row being Adult resources and second row being juveniles resources
                # resFreq   : resource freq, should be relative i.e. add up to one once again first row is adult, second is juvenile. 
                # resGen    : A measure of how generalist the consumers are a matrix with 1 column 2 rows, first row is juvenile
                # fmax      : A measure of the maximum number of offspring possible
                # kA/kJ     : A measure of how density dependent the fecundity and survival is. Half saturation constant 
                # mutProb   : Mutation probabilty
                # mutVar    : Mutation varitation in the case a mutation occurs
                # time.steps: how long the function runs e.g. years
                # iniP(A/J) : inital phenotype of adult and juvenile(different)
                # nmorphs   : Number of morphs in intial run


resourceCompetition <- function(popSize, resProp, resFreq, resGen=matrix(c(0.1,0.1),ncol=1, nrow=2), fmax = 10, kA = 0.5, kJ = 0.2, mutProb=0.001, mutVar=0.1, time.steps=200, iniPA=4, iniPJ=4, nmorphs = 1){
  
  pop <- matrix(data = NA, ncol = 4, nrow = nmorphs)
  
  pop[,1] <- popSize
  pop[,2] <- iniPA
  pop[,3] <- iniPJ
  
  colnames(pop) <- c("Number of indivduals", "Adult trait", "Juvenile trait", "Proxy")
  
  stats         <- NULL
  phenotype <- matrix(NA, ncol=4, nrow=0)                                       # Creates a matrix to store all phenotype values of all individuals in.
  colnames(phenotype) <- c("year","Number of indivduals", "Phenotype A", "Phenotype J")
  
  for (t in 1:time.steps)
  
  
  
}





