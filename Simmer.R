
# Script for running simulations



# Initialization ----------------

# Resources

# SLC:
resource.freq <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)           # res. freq. 
resource.prop <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)             # res. property 
abundance <- 20000
resource.abundance <- abundance*resource.freq


# CLC:

resource.frequency <- c(0.1,  0.1,  0.1,  0.1,  0.1, 0.1,  0.1,  0.1,  0.1,  0.1,   # res. freq. of adults (percentage)
                        0.1,  0.1,  0.1,  0.1,  0.1, 0.1,  0.1,  0.1,  0.1,  0.1)   # res. freq. pf juveniles
resource.property  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,                  # res. property of adults
                        1, 2, 3, 4, 5, 6, 7, 8, 9, 10)                   # res. property of juveniles

resource.abundance.adults     <- 20000                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- 20000

resFreqMatrix <- matrix(resource.frequency, nrow=2, ncol=10, byrow = TRUE)
resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

resPropMatrix <- matrix(resource.property, nrow=2, ncol=10, byrow = TRUE) 

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))

rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))

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

# SLC:

resource.prop <- N.resource.property             # res. property 
abundance <- 20000
resource.abundance <- abundance*N.resource.frequency


# CLC:

resFreqMatrix <- matrix(N.resource.frequency, nrow=2, ncol=length(N.resource.frequency), byrow = TRUE)

resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))


resPropMatrix <- matrix(N.resource.property, nrow=2, ncol=length(N.resource.property), byrow = TRUE) 


rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))


# Model runs with varied sigma:

sigma <- c(0.15, 0.3, 0.45, 0.6, 0.75)
Total_species_SLC <- c()
Total_species_CLC <- c()

# SLC:


for(i in 1:5){
  print(paste0("loop", i, "started"))
  outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 1000)
  phenotypesSLC <- outputSLC$phenotypes
  
  phenodataSLC <- NULL
  
  phenodataSLC <- data.frame(
    Year = outputSLC$phenotypes[, 1],
    Trait = outputSLC$phenotypes[, 3],
    Num_Individuals = outputSLC$phenotypes[, 2]
  )
  
  last_year_dataSLC <- phenodataSLC[phenodataSLC$Year == max(phenodataSLC$Year), ]
  print(nrow(last_year_dataSLC))
  
  last_year_dataS <- subset(last_year_dataSLC, select = -Year)
  last_year_dataS <- subset(last_year_dataS, select = -Num_Individuals)
  rownames(last_year_dataS) <- NULL
  rownames(last_year_dataSLC) <- NULL
  
  
  distance_matrix <- as.matrix(dist(last_year_dataS[, 1, drop = FALSE], method = "euclidean"))
  
  
  distance_matrix[lower.tri(distance_matrix)] <- NA
  
  
  
  # Set a threshold for similarity (adjust as needed)
  threshold <- 0.2
  
  # Find indices of individuals to keep
  
  
  same <- which(distance_matrix < threshold, arr.ind = T)
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
  
  rownames(last_year_dataSLC) <- NULL
  final_data <- last_year_dataSLC         # Place to store filtered data
  total.sub <- c()                     # Place to store subspecies
  
  #Add population count of "subspecies" to main species
  
  for(i in seq_along(groups)){
    combo <- NULL
    combo <- groups[[i]]
    main <- combo[which.max(final_data[combo,3])]
    sub <- combo[-which.max(final_data[combo,3])]
    final_data[main,3] <- final_data[main,3] + sum(final_data[sub,3])
    total.sub <- rbind(c(total.sub, sub))
    
  }
  # Remove subspecies
  final_data <- final_data[-total.sub, ]
  
  Total_species_SLC[i] <- as.numeric(nrow(final_data))
}





# CLC:

for(i in 1:5){
  
  
}
outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix, popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)

statsCLC <- outputCLC$stats
phenotypesCLC <- outputCLC$phenotypes
LastPhenoCLC <- outputCLC$LastPheno
LastStatsCLC <- outputCLC$LastStats

