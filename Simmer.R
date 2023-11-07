
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

# Model runs:

# SLC:

outputSLC <- resourceCompetitionSLC(resProp=resource.prop, resFreq=resource.abundance, popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)

statsSLC <- outputSLC$stats
phenotypesSLC <- outputSLC$phenotypes

# CLC:

outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)

statsCLC <- outputCLC$stats
phenotypesCLC <- outputCLC$phenotypes

