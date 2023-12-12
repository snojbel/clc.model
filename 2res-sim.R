# Background jobs 2 resources


resource.prop <- c(-1,1)
resource.frequency <- c(0.5, 0.5)
resource.frequency.as <- c(0.2, 0.8)

resource.abundance.adults <- 20000
resource.abundance.juveniles <- 20000

resFreqMatrix <- matrix(resource.frequency, nrow=2, ncol=length(resource.frequency), byrow = TRUE)
resFreqMatrixAs <- matrix(resource.frequency.as, nrow=2, ncol=length(resource.frequency.as), byrow = TRUE)


resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

resFreqMatrixAs[1, ] <- resFreqMatrixAs[1, ]*resource.abundance.adults
resFreqMatrixAs[2, ] <- resFreqMatrixAs[2, ]*resource.abundance.juveniles

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrixAs))

rownames(resFreqMatrixAs) <- c("Adult", "Juvenile")
colnames(resFreqMatrixAs)  <- paste0("Resource ", 1:ncol(resFreqMatrixAs))

resPropMatrix <- matrix(resource.prop, nrow=2, ncol=length(resource.prop), byrow = TRUE) 


rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))


# 2 resource run symmetric -----------------------



sigma <- seq(from = 0.1, to = 0.8, by = 0.05)

last_year_list_2_res <- list()


for(i in 1:length(sigma)){
  
  
  outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
  
  phenodataCLC <- NULL
  
  phenodataCLC <- data.frame(
    Year = outputCLC$phenotypes[, 1],
    Adult_Trait = outputCLC$phenotypes[, 3],
    Juvenile_Trait = outputCLC$phenotypes[, 4],
    Num_Individuals = outputCLC$phenotypes[, 2])
  
  last_year_list_2_res[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  
  
}



# 2 resource run asymmetric 



sigma <- seq(from = 0.1, to = 0.8, by = 0.05)

last_year_list_2_res_as <- list()


for(i in 1:length(sigma)){
  
  
  outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrixAs, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
  
  phenodataCLC <- NULL
  
  phenodataCLC <- data.frame(
    Year = outputCLC$phenotypes[, 1],
    Adult_Trait = outputCLC$phenotypes[, 3],
    Juvenile_Trait = outputCLC$phenotypes[, 4],
    Num_Individuals = outputCLC$phenotypes[, 2])
  
  last_year_list_2_res_as[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  
  
}




