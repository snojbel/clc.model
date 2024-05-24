# Parallel code

library(parallel)
library(snow)
library(doParallel)


# Run through of http://www.trutschnig.net/r_parallel.html 





# Exercise 1

adress <- url("http://www.trutschnig.net/RTR2015.RData")

load(adress)
head(RTR2015)

compare2 <- data.frame(funct, runtime)

avg_dl <- mean(RTR2015$rtr_speed_dl)
sd_dl <- sd(RTR2015$rtr_speed_dl)

# lapply
init <- Sys.time()
lapply(RTR2015$rtr_speed_dl, function(x) (x-avg_dl)/sd_dl)
compare2$runtime[1] <- Sys.time()-init


# sapply
init <- Sys.time()
sapply(RTR2015$rtr_speed_dl, function(x) c("init" = x, "std" = (x - avg_dl)/sd_dl))
compare2$runtime[2] <- Sys.time() - init

# Clusters and cores

# Number of workers
no.cores <- detectCores(logical = TRUE) - 4

# Clusters
cl <- makeCluster(no.cores)  # Remeber to release workers bound at a later point with stopCluster(cl)


# parLapply, input can be list of vector, output is list
# parSapply, input can be a list or vector, returns a vector or matrix


# Exercise 2


clusterExport(cl = cl, c("avg_dl", "sd_dl"), envir = .GlobalEnv)

init <- Sys.time()
parLapply(cl, RTR2015$rtr_speed_dl, function(x)  (x - avg_dl)/sd_dl)
compare2$runtime[3] <- Sys.time() - init


init <- Sys.time()
parSapply(cl, RTR2015$rtr_speed_dl, function(x) c("init" = x, "std" = (x - avg_dl)/sd_dl))
compare2$runtime[4] <- Sys.time() - init

stopCluster(cl)

compare2



# Set up for test

mutProbs <- c(0.0005, 0.001, 0.002)

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


# CLC:

resFreqMatrix <- matrix(N.resource.frequency, nrow=2, ncol=length(N.resource.frequency), byrow = TRUE)

resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))


resPropMatrix <- matrix(N.resource.property, nrow=2, ncol=length(N.resource.property), byrow = TRUE) 


rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))

clusterExport(cl = cl, c("resourceCompetitionCLC", "resPropMatrix", "resFreqMatrix", "mutProbs"), envir = .GlobalEnv)


parallel.test <- parSapply(cl, mutProb = mutProbs, FUN = resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, 
                            resGen=matrix(c(0.15,0.15)), popSize = 10, mutVar=0.05, time.steps = 50000))

stopCluster(cl)


