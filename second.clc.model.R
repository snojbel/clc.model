



# Second CLC model, optimized, maybe

#Shorthands used: popSize   : Initial population size 
                # resProp   : resource property, a matrix with first row being Adult resources and second row being juveniles resources
                # resFreq   : resource freq, should be relative i.e. add up to one once again first row is adult, second is juvenile. 
                # resGen    : A measure of how generalist the consumers are a matrix with 1 column 2 rows, first row is adults
                # fmax      : A measure of the maximum number of offspring possible
                # kA/kJ     : A measure of how density dependent the fecundity and survival is. Half saturation constant 
                # mutProb   : Mutation probability
                # mutVar    : Mutation variation in the case a mutation occurs
                # time.steps: how long the function runs e.g. years
                # iniP(A/J) : initial phenotype of adult and juvenile(different)
                # nmorphs   : Number of morphs in initial run
                # im        : Chance of an immigrant appearing
                # threshold : Sets how many individuals is needed for it to be considered a different morph. Any morph with less than threshold*total pop size at the end of the run will be removed.

# Libraries:
library(extrafont)   #needed to add extra fonts
#font_import()  #Only needed first time in R
#loadfonts()
#fonts() #to check names of fonts
 
library(viridisLite)  # Color things
library(viridis)
library(ggplot2)      # Prettier plots
library(gridExtra)    #For plotting side by side and more in ggplot




# Full function ----------------------------------------------------------------


resourceCompetitionCLC <- function(popSize, resProp, resFreq, resGen=matrix(c(0.15,0.15),ncol=1, nrow=2), fmax = 2, 
                                   kA = 0.5, kJ = 0.5,mutProb=0.001, mutVar=0.1, time.steps=200, iniPA=5, iniPJ=5, 
                                   threshold = 0.005, nmorphs = 1, im = 0){

  
  pop <- matrix(data = NA, ncol = 4, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.
  
  pop[,1] <- popSize
  pop[,2] <- iniPA
  pop[,3] <- iniPJ
  
  colnames(pop) <- c("Number of indivduals", "Adult trait", "Juvenile trait", "Proxy")
 
  stats         <- matrix(data = c(0, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2]), 
                            mean(pop[,3]), var(pop[,3])), nrow = 1, ncol = 7)                                                             #Where we will eventually save our stats and phenotypes
  phenotypes <- matrix(data = c(0, popSize, iniPA, iniPJ), nrow = 1, ncol = 4)
  colnames(phenotypes) <- c("Year", "Number of indivduals", "Adult trait", "Juvenile trait")                                        
  
  epsilon <- .Machine$double.eps^10  #Added when some number become zero, very small number
  possAtrait <- seq(from = min(resProp[1,])-1, to = max(resProp[1,])+1, by = mutVar)  # Used when generating immigrants
  possJtrait <- seq(from = min(resProp[2,])-1, to = max(resProp[2,])+1, by = mutVar)
      
  
  for (t in 1:time.steps){
    # Deterministic fecundity proxy alpha ------------------------------------

      adults     <- pop                                                            
      alphaA     <- NULL                                                           # Will create a matrix with all alpha values
      alphaSumA  <- NULL
      

      resPropAduMatrix <- matrix(data = resProp[1,], ncol = ncol(resProp), nrow = nrow(adults), byrow = T)
      aduTrait <- adults[, 2]
      aduTraitMatrix <- matrix(data = rep(aduTrait, each = ncol(resProp)), ncol = ncol(resProp), nrow = nrow(adults), byrow = T)
      
      alphaA           <- exp(-(((aduTraitMatrix-resPropAduMatrix)^2)/(2*resGen[1,1])^2)) + epsilon                 # Calculation of individual alpha
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
      
      alphaJ            <- exp(-(((juvTraitMatrix-resPropJuvMatrix)^2)/(2*resGen[2,1]^2))) + epsilon
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
      
      #if (runif(1) < im){
        
       # Atrait <- sample(x = possAtrait, size = 1)
       # Jtrait <- sample(x = possJtrait, size = 1)
        
      #  if(sum(pop[,2] == Atrait & pop[,3] == Jtrait) == 0) {                   # Checks wheter a exact match of immigrant already exists
       #   rbind(pop, c(1, Atrait, Jtrait, NA))
       # } else{
       #   same <- which(pop[,2] == Atrait & pop[,3] == Jtrait)
       #   pop[same,1] <- pop[same,1]+1
       # }
        
     # }
      

      
    # extract stats and phenotype ---------------------------------------------
    
   
    
    stats <- rbind(stats, c(t, sum(pop[,1]), nrow(pop), mean(pop[,2]), var(pop[,2]), 
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
  colnames(stats) <- c("year", "population size", "Number of morphs", "mean A trait", "var A", "mean J trait", "var J")
  rownames(phenotypes) <- NULL
  
  return(list(stats=stats, phenotypes=phenotypes))  # LastPheno = LastPheno, LastStats = LastStats Add if we want to remove low abundance morphs                               #returns both the stats and the phenotype
  
  
}

# Appoximation of a normal distribution for resources -------------





# With discrete data
resources <- sort(rpois(10000, 20))
resource.table <- table(resources)
resource.matrix <- matrix(c(as.integer(names(resource.table)), as.numeric(resource.table/sum(resource.table))), nrow = 2, byrow = T)
row.names(resource.matrix) <- c("Resource Characteristic", "Resource Frequency")
hist(resources)

# With cumulative distribution :

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



plot(N.resource.frequency)


resource.abundance.adults     <- 20000                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- 20000

resFreqMatrix <- matrix(N.resource.frequency, nrow=2, ncol=length(N.resource.frequency), byrow = TRUE)

resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))


resPropMatrix <- matrix(N.resource.property, nrow=2, ncol = length(N.resource.frequency), byrow = TRUE) 


rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))

# Initialization ----------------------------------------------------------------
                  
resource.frequency <- c(0.1,  0.1,  0.1,  0.1,  0.1, 0.1,  0.1,  0.1,  0.1,  0.1,   # res. freq. of adults (percentage)
                        0.1,  0.1,  0.1,  0.1,  0.1, 0.1,  0.1,  0.1,  0.1,  0.1)   # res. freq. pf juveniles

resource.property<- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,                  # res. property of adults
                      1, 2, 3, 4, 5, 6, 7, 8, 9, 10)                   # res. property of juveniles

resource.abundance.adults     <- 15000                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- 15000

resFreqMatrix <- matrix(resource.frequency, nrow=2, ncol=10, byrow = TRUE)

resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))


resPropMatrix <- matrix(resource.property, nrow=2, ncol=10, byrow = TRUE) 


rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))



outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix, popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)

statsCLC <- outputCLC$stats
phenotypesCLC <- outputCLC$phenotypes
#LastPhenoCLC <- outputCLC$LastPheno
#LastStatsCLC <- outputCLC$LastStats



# Plotting  -------------------------------------------------------------------


# Population size and number of phenotypes:
par(mfrow=c(2,1))


plot(x = statsCLC[, 1], y = statsCLC[, 2], xlab = "Year", ylab = "Total population size", type = "l")

plot(x = statsCLC[, 1], y = statsCLC[, 3], xlab = "Year", ylab = "Number of Phenotypes", type = "l")

# Trait divergence in adults and juveniles (2 plots)
par(mfrow=c(2,1))

pColors <- rgb(0.5, 0.3, 0.7, alpha = (phenotypesCLC[,2]/(100+phenotypesCLC[,2])))
plot(x = phenotypesCLC[, 1], y = phenotypesCLC[,3], col = pColors, xlab = "Year", ylab = "Adult trait", pch = 16, family = "LM Roman 10", cex.lab = 1.5, cex.axis = 1.3, cex = 1.5)
plot(x = phenotypesCLC[, 1], y = phenotypesCLC[,4], col = pColors, xlab = "Year", ylab = "Juvenile Trait", pch = 16, family = "LM Roman 10", cex.lab = 1.5, cex.axis = 1.3, cex = 1.5)

# On same plot
par(mfrow=c(1,1))
AColors <- rgb(0.5, 0.3, 0.7, alpha = (phenotypesCLC[,2]/(100+phenotypesCLC[,2])))
JColors <- rgb(0.5, 0.7, 0.2, alpha = (phenotypesCLC[,2]/(100+phenotypesCLC[,2])))
plot(x = phenotypesCLC[, 1], y = phenotypesCLC[,3], col = AColors, xlab = "Year", ylab = "Adult trait", pch = 16, cex = 2)
points(x = phenotypesCLC[, 1], y = phenotypesCLC[,4], col = JColors, xlab = "Year", ylab = "Juvenile Trait", pch = 16)

# Prettification of plots -----------------------------------------------------


# Creating data frame for easy plotting
phenodataCLC <- data.frame(
  Year = outputCLC$phenotypes[, 1],
  Adult_Trait = outputCLC$phenotypes[, 3],
  Juvenile_Trait = outputCLC$phenotypes[, 4],
  Num_Individuals = outputCLC$phenotypes[, 2]
)

transparency <- phenodataCLC$Num_Individuals / max(phenodataCLC$Num_Individuals)

# ------------- Trait divergence


evoAdu <- ggplot(phenodataCLC, aes(x=Year, y=Adult_Trait)) + 
                 geom_point(size = 2.5, alpha = transparency, color = rgb(0.13, 0.57, 0.55)) +
                 xlab("Year") + ylab("Adult Trait") +
                 theme_minimal(base_family = "LM Roman 10", base_size = 18)
                 
                 

evoJuv <- ggplot(phenodataCLC, aes(x=Year, y=Juvenile_Trait)) + 
                  geom_point(size = 2.5, alpha = transparency, color = rgb(0.27, 0.001, 0.33)) +
                  xlab("Year") + ylab("Juvenile Trait") +
                  theme_minimal(base_family = "LM Roman 10", base_size = 18) 
                  

grid.arrange(evoAdu,evoJuv, nrow = 2, widths = c(1))


# -----------------Scatter plot 

last_year_data <- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
color_palette <- mako(length(last_year_data$Adult_Trait))

ggplot(last_year_data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color_palette) +                                  # Add points
  labs(x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  theme_minimal(base_family = "LM Roman 10", base_size = 18)


# -----------------------------Making some cuts on the number of species:


last_year_dataCLC <- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
last_year_dataC <- subset(last_year_dataCLC, select = -Year)
last_year_dataC <- subset(last_year_dataC, select = -Num_Individuals)
rownames(last_year_dataC) <- NULL
rownames(last_year_data) <- NULL


distance_matrix_adult <- as.matrix(dist(last_year_dataC[, 1, drop = FALSE], method = "euclidean"))
distance_matrix_juvenile <- as.matrix(dist(last_year_dataC[, 2, drop = FALSE], method = "euclidean"))

distance_matrix_adult[lower.tri(distance_matrix_adult)] <- NA
distance_matrix_juvenile[lower.tri(distance_matrix_juvenile)] <- NA


# Set a threshold for similarity (adjust as needed)
threshold <- 0.2

# Find indices of individuals to keep


same <- which(distance_matrix_adult < threshold & distance_matrix_juvenile < threshold, arr.ind = T)
same <- same[same[, 1]-same[,2] != 0, , drop = FALSE]
rownames(same) <- NULL


# Initialize an empty list to store groups
groups <- list()

# Function to find group index for a species
find_group <- function(species_id) {
  for (i in seq_along(groups)) {
    if (species_id %in% unlist(groups[[i]])) {
      return(i)
    }
  }
  return(0)
}

# Iterate over rows in the matrix
for (i in 1:nrow(same)) {
  species1 <- same[i, 1]
  species2 <- same[i, 2]
  
  # Find groups for each species
  group1 <- find_group(species1)
  group2 <- find_group(species2)
  
  if (group1 == 0 & group2 == 0) {
    # Create a new group
    groups <- c(groups, list(c(species1, species2)))
  } else if (group1 == 0) {
    # Add species1 to the group containing species2
    groups[[group2]] <- c(groups[[group2]], species1)
  } else if (group2 == 0) {
    # Add species2 to the group containing species1
    groups[[group1]] <- c(groups[[group1]], species2)
  } else if (group1 != group2) {
    # Merge two groups
    groups[[group1]] <- c(groups[[group1]], groups[[group2]])
    groups <- groups[-group2]
  }
}

# Filter out duplicate species in each group
groups <- lapply(groups, function(group) unique(group))

rownames(last_year_data) <- NULL
final_data <- last_year_dataCLC         # Place to store filtered data
total.sub <- c()                     # Place to store subspecies

#Add population count of "subspecies" to main species

for(i in seq_along(groups)){
  combo <- NULL
  combo <- groups[[i]]
  main <- combo[which.max(final_data[combo,4])]
  sub <- combo[-which.max(final_data[combo,4])]
  final_data[main,4] <- final_data[main,4] + sum(final_data[sub,4])
  total.sub <- rbind(c(total.sub, sub))
  
}
# Remove subspecies
final_data <- final_data[-total.sub, ]

Total_species <- as.numeric(nrow(final_data))

# Plotting filtered data
color_palette <- mako(length(final_data$Adult_Trait))

ggplot(final_data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color_palette) +                                  # Add points
  labs(x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  theme_minimal(base_family = "LM Roman 10", base_size = 18)




sum(final_data[,4])
sum(last_year_data[,4])

# Little extra check

LastPhenodataCLC <- data.frame(
  Num_Individuals = outputCLC$LastPheno[, 2],
  Adult_Trait = outputCLC$LastPheno[, 3],
  Juvenile_Trait = outputCLC$LastPheno[, 4]
)


color_palette <- mako(length(LastPhenodataCLC$Adult_Trait))

ggplot(LastPhenodataCLC, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color_palette) +                                  # Add points
  labs(x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  theme_minimal(base_family = "LM Roman 10", base_size = 18)






# 3d plot try (not yet successfull) --------------------------------------------



library(plotly)

mycolors <- mako(n = nrow(phenotypes))

plot <- plot_ly( x=~phenotypes[,3], y=~phenotypes[,4], z=~phenotypes[,1], type = 'scatter3d', mode = 'lines',
       
       line = list(width = 4, color = ~1:nrow(phenotypes), colorscale = "Viridis"))

plot <- plot %>% layout(scene = list(zaxis = list(title = "Time"), xaxis = list(title = "Adult Trait"),yaxis = list(title = "Juvenile Trait")))
plot


#Test 2 s text oooh boo

#TESSSSTINGGGGG
#ADDING COMMENTSzzzzzzzzzZ!


#Adding things that might conflict!

#Fixing things

# Add a legend to distinguish between patches
#legend("topright", legend=c("Patch 1", "Patch 2"), col=c(rgb((1-0.9)*0.7,0.59,0.8, alpha =0.6), rgb((2-0.9)*0.7,0.59,0.8, alpha =0.6)), pch=19)

# Add a title to the plot
#title("Individual Phenotypes Over Years for Patch 1 and Patch 2")

# -------------------------------




# Function split into smaller parts(for easier debugging): -----------------------------------------------------------------------------------------

#----------------- Population initialization:
popSize <- 15
iniPA   <- 4
iniPJ   <- 4
nmorphs <- 1

pop <- matrix(data = NA, ncol = 4, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.

pop[,1] <- popSize
pop[,2] <- iniPA
pop[,3] <- iniPJ

colnames(pop) <- c("Number of indivduals", "Adult trait", "Juvenile trait", "Proxy")

# --------------------- Resource initialization (adult)

resFreqA <- c(2,2,5,2,2)              # res. freq. of adults
resPropA <- c(2, 3, 4, 5, 6)          # res. property of adults
resGenA  <- 0.3
kA       <- 0.3
fmax     <- 10 


#----------------- Fecundity  ----------------

birth <- function (pop, resPropA, resFreqA, resGenA, kA, fmax) { 
    
    # Deterministic fecundity proxy alpha ------------------------------------
    
    adults    <- pop                                                            
    alphaA    <- NULL                                                             # Will create a matrix with all alpha values
    alphaSumA <- NULL
    
    for(r in 1:length(resPropA)){                                                 # Loops through each resource to compute 
      
      rp <- resPropA[r]                                                           # Resource property for adult resource [2,]
      
      alpha     <- exp(-(adults[, 2]-rp)^2/ (2*resGenA)^2 )                       # Compute  alphas for all adults, if different generalities for stages want to be added modify resGen
      alphaA    <- cbind(alphaA, alpha)                                           # Store alphas for all resources for all morphs
      alphaSumA <- c(alphaSumA, sum(alpha*adults[,1]))                            # Store the sum of all alphas for each resource
      
    }
    
    for(i in 1:nrow(adults)){
      
      fec <- 0
      for(r in 1:length(resPropA)){                                                  # Goes through each resource and calculates source specific  
        # fecundity and then adds them together to create the total fec..
        R  <- resFreqA[r]
        fec <- fec + R * (alphaA[i,r]/alphaSumA[r])                               # Total fecundity is the sum of fecundities gained 
        # from all the different resources. Which is derived
      }                                                                           # both the ability to consume a resource, alphaA, and the 
      adults[i,4] <- fmax * (fec/(kA+fec))                                        # the total ability of all individuals to consume 
      # said resource, alphaSumA. R = resource available
    }
    
    # Spawning of offspring -------------------------------------------------
    
    juveniles <- adults                                                          # Create a matrix were we will add juveniles into
    
    for(i in 1:nrow(juveniles)) {
      
      juv             <- rpois(n = 1, lambda = adults[i,1] * adults[i,4])        # This calculates the number of juveniles based on a possion distribution. 
      juveniles[i, 1] <- juv                                                     # States number of juveniles in 1 column
      
    }
    
    pop <- juveniles
  
  return(pop)
}

populationSize <- c()
timesteps <- 50

for(t in 1:timesteps){
  
  pop <- birth(pop = pop, resPropA = resPropA , resFreqA = resFreqA, resGenA = resGenA, kA = kA, fmax = 3)
  
  populationSize<- rbind(populationSize, sum(pop[,1])) 
  
}
rownames(populationSize) <- paste0("Year", 1:nrow(populationSize)) 
colnames(populationSize) <- "Number of indivduals"


plot(x = 1:nrow(populationSize), y = populationSize, type = "l", xlab = "Years", ylab = "Abundance")

#----------------- Mutation   -------------------

mutation <- function(pop, mutProb, mutVar){

juveniles <- pop
probs <- juveniles[,1]/sum(juveniles[, 1])  # Generates probability of morph being mutated based upon number of individuals.
N.mut <- as.numeric(rbinom(n = 1, size = sum(juveniles[, 1]), prob = mutProb))


if(N.mut > 0){
  #print("Mutation!")
  random.choice <- c()
  mutation.pos <- c()
  for (m in 1:length(N.mut)){
    random.choice <- rmultinom(n = 1 , size = 1, prob = probs)    # Randomly chooses which morphs to mutate based on probs
    mutation.pos[m] <- as.numeric(which(random.choice == 1))
  }
  
  for (i in  1:length(mutation.pos)){
    
    mutChange <- rnorm(n=1, mean=0, sd=mutVar)
    juveniles[mutation.pos[i], 1] <- juveniles[mutation.pos[i], 1] - 1       # Removes the mutated individual from the morph
    
    if(rbinom(n = 1, size = 1, prob = 0.5) == 0){                            # Randomly choose whether adult or juvenile trait gets morphed.
      
      
      new.morph <- matrix(data = c(1, juveniles[mutation.pos[i], 2] + mutChange, #Changes adult trait to a new trait and adds it to the juveniles
                                   juveniles[mutation.pos[i], 3], 0), 
                          ncol = ncol(juveniles), nrow = 1)
      juveniles <- rbind(juveniles, new.morph)
    }
    else {
      # Removes the mutated individual from the morph
      new.morph <- matrix(data = c(1, juveniles[mutation.pos[i], 2] ,        #Changes juvenile trait to a new trait and adds it to the juveniles
                                   juveniles[mutation.pos[i], 3]+ mutChange, 0), 
                          ncol = ncol(juveniles), nrow = 1)
      juveniles <- rbind(juveniles, new.morph)
    }
    
  }
  
}


pop <- juveniles

return(pop)

}


#----------------- Survival  ------------------------------------

death <- function(pop, resPropJ, resFreqJ, resGenJ, kJ){
  juveniles <- pop
  alphaSumJ <- NULL
  alphaJ    <- NULL
  
  for(r in 1:length(resPropJ)){                                                  # Survival of juveniles also depends on resource availability.
    
    rp <- resPropJ[r] 
    
    alpha     <- exp(-(juveniles[,3]-rp)^2/ (2*resGenJ) )
    alphaJ    <- cbind(alphaJ, alpha)
    alphaSumJ <- c(alphaSumJ, sum(alpha) )
    
  }
  
  for(i in 1:nrow(juveniles)){
    
    sur <- 0
    for(r in 1:length(resPropJ)){                                                # Goes through each resource and calculates source specific  
                                                                                 # survival and then adds them together to create the total sur
      R  <- resFreqJ[r]
      sur <- sur + R * (alphaJ[i,r]/alphaSumJ[r])                             
      
    }                                                                         
    
    det.sur <- sur/(kJ+sur)                                                       # this is to create a saturating function, survival can never be over 1 higher kJ means flatter curve
    juveniles[i,1] <- rbinom(n = 1 , size = juveniles[i,1], prob = det.sur)       # sees how many survive given their det. survival.
  }
  
  pop <- juveniles[juveniles[, 1] != 0, , drop = FALSE]
  
  return(pop)
}

# Test of survival and birth and mutation together ---------------------------------------------------

#----------------- Population initialization
popSize <- 15
iniPA   <- 4
iniPJ   <- 4
nmorphs <- 1

pop <- matrix(data = NA, ncol = 4, nrow = nmorphs)                             # Each column in this matrix is one phenotype combination.

pop[,1] <- popSize
pop[,2] <- iniPA
pop[,3] <- iniPJ

colnames(pop) <- c("Number of indivduals", "Adult trait", "Juvenile trait", "Proxy")

mutProb = 0.005
mutVar = 0.1

populationSize <- c()
numberPhenotypes <- c()
Phenotypes <- matrix(data = c(0, popSize, iniPA, iniPJ), nrow = 1, ncol = 4)
colnames(Phenotypes) <- c("Year", "Number of indivduals", "Adult trait", "Juvenile trait")
timesteps <- 500

# --------------------- Resource initialization (adult)

resFreqA <- c(2,2,5,2,2)              # res. freq. of adults
resPropA <- c(2, 3, 4, 5, 6)          # res. property of adults
resGenA  <- 0.3
kA       <- 0.3
fmax     <- 10 

# --------------------- Resource initialization (juvenile)

resFreqJ <- c(10,3,1,3,10)              # res. freq. of juveniles
resPropJ <- c(1, 3, 4, 5, 7)          # res. property of juveniles
resGenJ  <- 0.2
kJ       <- 0.3


for(t in 1:timesteps){
  
  pop <- birth(pop = pop, resPropA = resPropA , resFreqA = resFreqA, resGenA = resGenA, kA = kA, fmax = 3)

  pop <- mutation(pop=pop, mutProb = mutProb, mutVar = mutVar)
  
  pop <- death(pop = pop, resPropJ = resPropJ, resFreqJ = resFreqJ, resGenJ = resGenJ, kJ = kJ)

  
  
  if (nrow(pop) == 0) {
    print("Extinction!")
    break
  }
  else {
    
   #print(paste0("Loop ", t ," works"))
    
        }
  
  numberPhenotypes <- rbind(numberPhenotypes, nrow(pop))
  pStats <- cbind(rep(t, nrow(pop)), pop[,1], pop[,2], pop[,3])
  Phenotypes <- rbind(Phenotypes, pStats)
  populationSize<- rbind(populationSize, sum(pop[,1])) 
  
}


rownames(populationSize) <- paste0("Year", 1:nrow(populationSize)) 
colnames(populationSize) <- "Number of indivduals"

colnames(numberPhenotypes) <- "# of phenotypes"

rownames(Phenotypes) <- NULL


plot(x = 1:nrow(populationSize), y = populationSize, type = "l", xlab = "Years", ylab = "Abundance")


















