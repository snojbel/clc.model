
# 2-dimensional resource
# The resource has two charactertistics that determine how effective an indivdual is at foraging for it.
# Likewise individuals will have two different traits, but these affect each life stage equally
library(ggplot2) 
library(gridExtra)   
library(extrafont)
library(viridisLite)  # Color things
library(viridis)
library(VGAM)
library(plot3D)

# resGen refers to the two different traits again, but now they affect both lifestages
# cov is the covariance between the two traits feeding efficiency

# Resources ------------------


quantiles <- seq(from = -5, to = 5, length.out = 26)


resource_prop_1 <- c()
for(i in 1:(length(quantiles)-1)){
  resource_prop_1[i] <- (quantiles[i]+quantiles[i+1])/2
}
resource_prop_2 <- resource_prop_1                         # Symmetrical resource properties

resProp1 <- rep(resource_prop_1, each = length(resource_prop_2))  # This is done to make the math in the function nice

resProp2 <- rep(resource_prop_2, times = length(resource_prop_1))

# Remember the number of resources is equal to length resprop1 x length resprop2

# Normal resource:

# Normal resources:

m <- 0 
s <- 3
N.resource.frequency <- c()
N.resource.property<- c(seq(from = -2.5, to = 2.5, length.out = (length(quantiles)-1)^2)) 



N.resource.frequency <- c()
r <- 1

for(i in 1:length(resource_prop_1)){
  for(k in 1:length(resource_prop_2)){
    upper <- pbinorm(quantiles[k+1], quantiles[i+1], mean1 = m, mean2 = m, var1 = s, var2 = s, cov12 = 0)
    rightlower <- pbinorm(quantiles[k+1], quantiles[i], mean1 = m, mean2 = m, var1 = s, var2 = s, cov12 = 0)
    leftlower <- pbinorm(quantiles[k], quantiles[i+1], mean1 = m, mean2 = m, var1 = s, var2 = s, cov12 = 0)
    lowerlower <- pbinorm(quantiles[k], quantiles[i], mean1 = m, mean2 = m, var1 = s, var2 = s, cov12 = 0)
   
    N.resource.frequency[r] <- upper - (rightlower+leftlower) + lowerlower
    r <- r+1
  }
  
}



resource.abundance <- 20000
sum(N.resource.frequency)
resFreq <- N.resource.frequency*resource.abundance 



# Test of resources normality working.

test <- matrix(resFreq, nrow = (length(quantiles)-1), ncol = (length(quantiles)-1), byrow = T)
hist3D(z = test, border = "black")

# Even resources

R <- (length(quantiles)-1)^2
E.resource.frequency <- c(rep(1/R, times = R)) 
resFreqE <- E.resource.frequency*resource.abundance

# Function ------------------

resourceCompetition2dr <- function(popSize, resProp1, resProp2, resFreq, resGen=matrix(c(0.15,0.15),ncol=1, nrow=2), im = 0.01, 
                                   fmax = 2, kA = 0.5, kJ = 0.5, mutProb=0.0005, mutVar=0.05, time.steps=200, iniP1=0, 
                                   iniP2 = 0, threshold = 0.005, nmorphs = 1, cov = 0){
  
  pop <- matrix(data = NA, ncol = 4, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniP1
  pop[,3] <- iniP2
  pop[,4] <- 0
  
  
  colnames(pop) <- c("Number of indivduals", "Trait 1", "Trait 2", "Proxy")
  
  stats         <- matrix(data = c(0, sum(pop[,1]), nrow(pop), mean(pop[,2]), 
                                   var(pop[,2]),mean(pop[,3]), var(pop[,3])) , nrow = 1, ncol = 7)                                                             #Where we will eventually save our stats and phenotypes
  phenotypes <- matrix(data = c(0, popSize, iniP1, iniP2), nrow = 1, ncol = 4)
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Trait 1", "Trait 2")                                        
  
  epsilon <- .Machine$double.eps^10  #Added when some number become zero, very small number
  posstrait1 <- seq(from = min(resProp1)-1, to = max(resProp1)+1, by = mutVar)
  posstrait2 <- seq(from = min(resProp2)-1, to = max(resProp2)+1, by = mutVar)
  
  for (t in 1:time.steps){
    
    # Deterministic fecundity proxy alpha ------------------------------------
    
    adults    <- pop                                                            
    alphaA    <- NULL                                                           # Will create a matrix with all alpha values
    alphaSumA <- NULL
    
    resPropAdu1Matrix <- matrix(data = resProp1, ncol = length(resProp1), nrow = nrow(adults), byrow = T)
    resPropAdu2Matrix <- matrix(data = resProp2, ncol = length(resProp2), nrow = nrow(adults), byrow = T)
    
    
    aduTrait1 <- adults[, 2]
    aduTrait2 <- adults[, 3]
    aduTrait1Matrix <- matrix(data = rep(aduTrait1, each = length(resProp1)), ncol = length(resProp1), nrow = nrow(adults), byrow = T)
    aduTrait2Matrix <- matrix(data = rep(aduTrait2, each = length(resProp2)), ncol = length(resProp2), nrow = nrow(adults), byrow = T)
    
    # The code above makes it so that when I calculate alpha I get a matrix where each row corresponds to a different morph¨
    # and each column a different resource
    
    alphaA           <- (1/(sqrt(((2*pi)^2)*((resGen[1,1]^2)*(resGen[2,1]^2)-cov^2))))*exp(-((((resGen[2,1]^2)*(aduTrait1Matrix-resPropAdu1Matrix)^2)+((resGen[1,1]^2)*(aduTrait2Matrix-resPropAdu2Matrix)^2)-(2*cov*(aduTrait1Matrix-resPropAdu1Matrix)*(aduTrait2Matrix-resPropAdu2Matrix))/(2*(resGen[1,1]^2)*(resGen[2,1]^2)-cov^2)))) + epsilon                 # Calculation of individual alpha
    adultAbund       <- adults[,1]
    adultAbundMatrix <- matrix(data = rep(adultAbund, each = length(resProp1)), ncol = length(resProp1), nrow = nrow(adults), byrow = T)  # Creation of a matrix with population size of each type in the rows
    alphaSumA        <- colSums((alphaA*adultAbundMatrix))                                                                         # Creation of matrix that reflects both the trait but also number of individuals in type
    
    
    RdivAlphaSumA       <- resFreq/alphaSumA
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
        
        if(rbinom(n = 1, size = 1, prob = 0.5) == 0){                            # Randomly choose whether trait 1 or 2 was mutated.
          
          
          new.morph <- matrix(data = c(1, juveniles[mutation.pos[i], 2] + mutChange, #Changes trait 1 to a new trait and adds it to the juveniles
                                       juveniles[mutation.pos[i], 3], 0), 
                              ncol = ncol(juveniles), nrow = 1)
          juveniles <- rbind(juveniles, new.morph)
        }
        else {
          
          new.morph <- matrix(data = c(1, juveniles[mutation.pos[i], 2] ,        #Changes trait 2 to a new trait and adds it to the juveniles
                                       juveniles[mutation.pos[i], 3]+ mutChange, 0), 
                              ncol = ncol(juveniles), nrow = 1)
          juveniles <- rbind(juveniles, new.morph)
        }
        
      }
      
    }
    
    
    
    # Kill off offspring -----------------------------------------------------
    
    
    alphaSumJ <- NULL
    alphaJ    <- NULL
    
    
    resPropJuv1Matrix <- matrix(data = resProp1, ncol = length(resProp1), nrow = nrow(juveniles), byrow = T)
    resPropJuv2Matrix <- matrix(data = resProp2, ncol = length(resProp2), nrow = nrow(juveniles), byrow = T)
    
    
    juvTrait1 <- juveniles[,2]
    juvTrait2 <- juveniles[,3]
    juvTrait1Matrix <- matrix(data = rep(juvTrait1, each = length(resProp1)), ncol = length(resProp1), nrow = nrow(juveniles), byrow = T)
    juvTrait2Matrix <- matrix(data = rep(juvTrait2, each = length(resProp2)), ncol = length(resProp2), nrow = nrow(juveniles), byrow = T)
    
    alphaJ            <- (1/(sqrt(((2*pi)^2)*((resGen[1,1]^2)*(resGen[2,1]^2)-cov^2))))*exp(-((((resGen[2,1]^2)*(juvTrait1Matrix-resPropJuv1Matrix)^2)+((resGen[1,1]^2)*(juvTrait2Matrix-resPropJuv2Matrix)^2)-(2*cov*(juvTrait1Matrix-resPropJuv1Matrix)*(juvTrait2Matrix-resPropJuv2Matrix))/(2*(resGen[1,1]^2)*(resGen[2,1]^2)-cov^2)))) + epsilon
    juvenAbund        <- juveniles[,1]
    juvenAbundMatrix  <- matrix(data = rep(juvenAbund, each = length(resProp1)), ncol = length(resProp1), nrow = nrow(juveniles), byrow = T)  # Creation of a matrix with population size of each type in the rows
    alphaSumJ         <- colSums(alphaJ*juvenAbundMatrix)                                                                         # Creation of matrix that reflects both the trait but also number of individuals in type
    
    
    RdivAlphaSumJ     <- resFreq/alphaSumJ
    RdivAlphaSumJTrans   <- matrix(data = RdivAlphaSumJ)
    
    
    Sur <- alphaJ%*%RdivAlphaSumJTrans
    juveniles[,4] <- (Sur/(kJ+Sur)) 
    
    
    
    
    juveniles[,1] <- rbinom(n = nrow(juveniles) , size = juveniles[,1], prob = juveniles[,4])
    
    
    pop <- juveniles[juveniles[, 1] != 0, , drop = FALSE]                        # all adults die after reproducing, so the new generation is only juveniles, and all rows with zero individuals are removed.
    
    #Adding immigrants ---------------------------------------------------------------------
    
    if (runif(1) < im){
      
      trait1 <- sample(x = posstrait1, size = 1)
      trait2 <- sample(x = posstrait2, size = 1)
      
      if(sum(pop[,2] == trait1 & pop[,3] == trait2) == 0) {                   # Checks whether a exact match of immigrant already exists
        rbind(pop, c(1, trait1, trait2, 0))
      } else{
        same <- which(pop[,2] == trait1 & pop[,3] == trait2)
        pop[same,1] <- pop[same,1] + 1
      }
      
    }
    
    
    # extract stats and phenotype ---------------------------------------------
    
    if(nrow(pop) == 0){                                                          # Checks whether population has reached zero, then it breaks the for loop.                                   
      print("Population extinction")
      break
    }
    
    stats <- rbind(stats, c(t, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2]),mean(pop[,3]), var(pop[,3]))) 
    
    pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2], pop[,3])
    phenotypes <- rbind(phenotypes, pStats)
    
    
  }
  
  #return output  ------------------------------------------------------------
  colnames(stats) <- c("year", "population size", "Number of morphs", "mean trait 1", "var trait 1", "mean trait 2", "var trait 2")
  rownames(phenotypes) <- NULL
  
 
  
  
  return(list(stats=stats, phenotypes=phenotypes))                                 #returns both the stats and the phenotype
  
  
}

# Make it a job to free up R console

job::job(output2dr = {
  output2dr <- resourceCompetition2dr(popSize = 10, iniP1=0, iniP2 = 0, resProp1 = resProp1, resProp2 = resProp2, resFreq = resFreq, 
                                      time.steps = 100000, resGen=matrix(c(0.10,0.10)))

# Control what is returned to the main session
  job::export(output2dr)
}, import = c(resProp1, resProp2, resFreq, resourceCompetition2dr)) #Set was is imported into the job

output2dr <- output2dr$output2dr
stats <- output2dr$stats

# Even job
job::job(output2drEven = {
  output2dr <- resourceCompetition2dr(popSize = 10, iniP1=0, iniP2 = 0, resProp1 = resProp1, resProp2 = resProp2, resFreq = resFreqE, 
                                      time.steps = 100000, resGen=matrix(c(0.10,0.10)))
  
  # Control what is returned to the main session
  job::export(output2drEven)
}, import = c(resProp1, resProp2, resFreqE, resourceCompetition2dr))

output2drEven <- output2drEven$output2dr
stats <- output2dr$stats

# Creating data frame for easy plotting
phenodata2dr <- data.frame(
  Year = output2dr$phenotypes[, 1],
  Trait1 = output2dr$phenotypes[, 3],
  Trait2 = output2dr$phenotypes[, 4],
  Num_Individuals = output2dr$phenotypes[, 2]
)

transparency <- phenodata2dr$Num_Individuals / max(phenodata2dr$Num_Individuals)

# ------------- Trait divergence


evoAdu <- ggplot(phenodata2dr, aes(x=Year, y=Trait1)) + 
  geom_point(size = 2.5, alpha = transparency, color = rgb(0.13, 0.57, 0.55)) +
  xlab("Year") + ylab("Adult Trait") +
  theme_minimal(base_family = "LM Roman 10", base_size = 18)



evoJuv <- ggplot(phenodata2dr, aes(x=Year, y=Trait2)) + 
  geom_point(size = 2.5, alpha = transparency, color = rgb(0.27, 0.001, 0.33)) +
  xlab("Year") + ylab("Juvenile Trait") +
  theme_minimal(base_family = "LM Roman 10", base_size = 18) 


grid.arrange(evoAdu,evoJuv, nrow = 2, widths = c(1))


# -----------------Scatter plot 

last_year_data <- phenodata2dr[phenodata2dr$Year == max(phenodata2dr$Year), ]
color_palette <- mako(length(last_year_data$Trait1))

ggplot(last_year_data, aes(x = Trait2, y = Trait1)) +
  geom_point(aes(size=Num_Individuals), color = color_palette) +                                  # Add points
  labs(x = "Trait 1", y = "Trait 2", size = "Number of individuals") +  
  scale_x_continuous(limits = c(-1.5,1.5)) +
  scale_y_continuous(limits = c(-1.5,1.5)) + 
  theme_minimal(base_family = "LM Roman 10", base_size = 18)










