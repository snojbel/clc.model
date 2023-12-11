
# Script for running simulations



# Initialization ----------------

# Evenly distributed Resources

# SLC:
resource.freq <- rep(0.1, times = 16)                                      # res. freq. 
resource.prop <- c(seq(from = -2.5, to = 2.5, length.out = 16))            # res. property 
abundance <- 1000
resource.abundance <- abundance*resource.freq


# CLC:

resource.property<- c(seq(from = -2.5, to = 2.5, length.out = 16)) 

resource.frequency <- rep(0.1, times = 16)

resource.abundance.adults     <- 1000                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- 1000

resFreqMatrix <- matrix(resource.frequency, nrow=2, ncol=length(resource.frequency), byrow = TRUE)
resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

resPropMatrix <- matrix(resource.property, nrow=2, ncol=length(resource.property), byrow = TRUE) 

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


# Skewed resource distribution ---------------------------------------

# SLC:
x <- 1/136
resource.freq <- c()

for (i in 1:16){ 
  resource.freq[i] <- i*x 
  }                                     # res. freq. 

resource.prop <- c(seq(from = -2.5, to = 2.5, length.out = 16))            # res. property 
abundance <- 20000
resource.abundance <- abundance*resource.freq


# CLC:

resource.property<- c(seq(from = -2.5, to = 2.5, length.out = 16)) 


resource.frequency <- c()
for (i in 1:16){ 
  resource.frequency[i] <- i*x 
}    


resource.abundance.adults     <- 20000                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- 20000

resFreqMatrix <- matrix(resource.frequency, nrow=2, ncol=length(resource.frequency), byrow = TRUE)
resFreqMatrix[1, ] <- resFreqMatrix[1, ]*resource.abundance.adults
resFreqMatrix[2, ] <- resFreqMatrix[2, ]*resource.abundance.juveniles

resPropMatrix <- matrix(resource.property, nrow=2, ncol=length(resource.property), byrow = TRUE) 

rownames(resFreqMatrix) <- c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqMatrix))

rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))


# Two resources -------------------------------------------------------

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
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resFreqAsMatrix))

rownames(resFreqMatrixAs) <- c("Adult", "Juvenile")
colnames(resFreqMatrixAs)  <- paste0("Resource ", 1:ncol(resFreqAsMatrix))

resPropMatrix <- matrix(resource.prop, nrow=2, ncol=length(resource.prop), byrow = TRUE) 


rownames(resPropMatrix)<-c("Adult", "Juvenile")
colnames(resFreqMatrix)  <- paste0("Resource ", 1:ncol(resPropMatrix))

# Model runs with varied sigma --------------------------------------------

sigma <- c(0.15, 0.3, 0.45, 0.6, 0.75)

Total_species_SLC_single <- c()

Total_species_SLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
rownames(Total_species_SLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
colnames(Total_species_SLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES

Total_species_CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
rownames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
colnames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES


# SLC:

for(i in 1:length(sigma)){
  
  outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
  
  #Filter out similar "species"
  
  final_data_SLC <- slc.groups(output = outputSLC)
  Total_species_SLC_single[i] <- nrow(final_data_SLC)
}


# If adults and juveniles can have different niche width

for(i in 1:length(sigma)){
  print(paste0("loop", i, "started"))
  for(k in 1:length(sigma)){
     
    outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
    
    #Filter out similar "species"
    
    final_data_SLC <- slc.groups(output = outputSLC)
    Total_species_SLC[i, k] <- nrow(final_data_SLC)
  }
 
}




# CLC:


for(i in 1:length(sigma)){
  print(paste0("loop", i, "started"))
  for(k in 1:length(sigma)){
    
    outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 10000)
    
    #Filter out similar "species"
    
    final_data_CLC <- clc.groups(output = outputCLC)
    Total_species_CLC[i, k] <- nrow(final_data_CLC)
  }
  
}

# 10 runs: ---------------------------------

# SLC, same sigma for adults and juv

Total_SLC_list <- list()
Total_abund_SLC_list <- list()

for(r in 1:10) {
  
  print(paste0("loop ", r, " started"))
  
  Total_species_SLC_single <- c()
  Abundance_species_SLC_single <- c()

  for(i in 1:length(sigma)){
    
    outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
    
    #Filter out similar "species"
    
    final_data_SLC <- slc.groups(output = outputSLC)
    Total_species_SLC_single[i] <- nrow(final_data_SLC)
    Abundance_species_SLC_single[i] <- sum(final_data_SLC[,3])
  }
  
  Total_SLC_list[[r]] <- Total_species_SLC_single
  Total_abund_SLC_list[[r]] <- Abundance_species_SLC_single
}

# Caluclating mean of 10 runs

Total_mean_SLC <- Reduce(`+`, Total_SLC_list) / length(Total_SLC_list)
Total_mean_abund_SLC <- Reduce(`+`, Total_abund_SLC_list) / length(Total_abund_SLC_list) 



# CLC

Total_CLC_list <- list()
Total_abund_CLC_list <- list()

for(r in 1:10){
  print(paste0("loop ", r, " started"))
  
  Total_species_CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
  colnames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES
  Abundance_species_CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #ADULTS
  colnames(Total_species_CLC) <- c(0.15, 0.3, 0.45, 0.6, 0.75) #JUVENILES
  
  for(i in 1:length(sigma)){
    
    for(k in 1:length(sigma)){
      
      outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
      
      #Filter out similar "species"
      
      final_data_CLC <- clc.groups(output = outputCLC)
      Total_species_CLC[i, k] <- nrow(final_data_CLC)
      
      Abundance_species_CLC[i, k] <- sum(final_data_CLC[, 4])
      
    }
    
  }
  Total_CLC_list[[r]] <- Total_species_CLC
  Total_abund_CLC_list[[r]] <- Abundance_species_CLC
}

# Caluclating mean of 10 runs

Total_mean_CLC <- Reduce(`+`, Total_CLC_list) / length(Total_CLC_list)
Total_mean_abund_CLC <- Reduce(`+`, Total_abund_CLC_list) / length(Total_abund_CLC_list) 

# Saving results

saveRDS(Total_CLC_list, file = "CLC1Alist.RData")

saveRDS(Total_SLC_list, file = "SLC1Alist.RData")

# 10 Runs of to see endpoint

last_year_list <- list()
iniPs <- seq(from = -2, to = 2, by = 0.5)

for(i in 1:length(iniPs)){
  print(paste0("loop ", i, " started"))
  outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = iniPs[i], iniPJ = iniPs[i], resGen=matrix(c(0.15, 0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 200000)
  
  
  
  phenodataCLC <- NULL
  
  phenodataCLC <- data.frame(
  Year = outputCLC$phenotypes[, 1],
  Adult_Trait = outputCLC$phenotypes[, 3],
  Juvenile_Trait = outputCLC$phenotypes[, 4],
  Num_Individuals = outputCLC$phenotypes[, 2])
  
  last_year_list[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
}



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



# Different sigma runs endpoint ------------------------

sigmas <- c(0.15, 0.3, 0.45, 0.6, 0.75)

last_year_list_CLC <- list()
last_year_list_SLC <- list()
plot_list_CLC <- list()
plot_list_SLC <- list()

for(i in 1:length(sigmas)){
  print(paste0("loop ", i, " started"))
  outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigmas[i], sigmas[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
  
  
  phenodataCLC <- NULL
  
  phenodataCLC <- data.frame(
    Year = outputCLC$phenotypes[, 1],
    Adult_Trait = outputCLC$phenotypes[, 3],
    Juvenile_Trait = outputCLC$phenotypes[, 4],
    Num_Individuals = outputCLC$phenotypes[, 2])
  
  last_year_list_CLC[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  
  outputSLC <- resourceCompetitionSLC(resProp=resource.prop, iniP = 0, resFreq=resource.abundance, resGen=matrix(c(sigmas[i],sigmas[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
  
  phenodataSLC <- data.frame(
    Year = outputSLC$phenotypes[, 1],
    Trait = outputSLC$phenotypes[, 3],
    Num_Individuals = outputSLC$phenotypes[, 2]
  )
  
  last_year_list_SLC[[i]]<- phenodataSLC[phenodataSLC$Year == max(phenodataSLC$Year), ]
  
  if(i == 1){
    plot_list_CLC[[i]] <- ggplot(last_year_list_CLC[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
        geom_point(colour= "#158FAD", aes(size=Num_Individuals), show.legend = FALSE) +                                  # Add points
        labs(title = substitute(sigma == value, list(value = sigmas[i])), y = "Adult Trait") +                 # Labels for the axes
        scale_x_continuous(limits = c(-3,3)) +
        scale_y_continuous(limits = c(-3,3)) +
        theme_classic(base_family = "LM Roman 10", base_size = 15)+
        theme(axis.text = element_text(family = "LM Roman 10"),
              axis.title = element_text(family = "LM Roman 10", size = 20),
              axis.title.x = element_blank(),
              plot.title = element_text(hjust = 0.5))
    
    plot_list_SLC[[i]] <- ggplot(last_year_list_SLC[[i]], aes(x = Trait, y = 0))+
      geom_point(colour = "#158FAD", aes(size=Num_Individuals), show.legend = FALSE)  +
      geom_segment(data = data.frame(x = c(-2, 0, 2), y = rep(0, 3)),
                   aes(x = x, xend = x, y = -0.1, yend = 0.1), color = "black", size = 0.5) +
      geom_text(data = data.frame(x = c(-2, 0, 2), label = c("-2", "0", "2")),
                aes(x = x, y = -0.2, label = label), vjust = 0.5, hjust = 0.5, size = 7, color = "black", family = "LM Roman 10") +
      annotate("segment",x=-3,xend=3, y=0, yend=0, linewidth=1) +
      annotate("segment",x=-3,xend=-3, y=-0.1,yend=0.1, linewidth=1) +
      annotate("segment",x=3,xend=3, y=-0.1,yend=0.1, linewidth=1) +
      scale_x_continuous(limits = c(-3,3)) +
      scale_y_continuous(limits = c(-1,1)) +
      scale_color_manual(values = unname(colours)) + 
      theme(panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    
  }
  else if(i == 3){
    plot_list_CLC[[i]] <- ggplot(last_year_list_CLC[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
        geom_point(colour= "#158FAD", show.legend = FALSE) +   #aes(size=Num_Individuals),                               # Add points
        labs(title = substitute(sigma == value, list(value = sigmas[i])), x = "Juvenile Trait") +                 # Labels for the axes
        scale_x_continuous(limits = c(-3,3)) +
        scale_y_continuous(limits = c(-3,3))  +
        theme_classic(base_family = "LM Roman 10", base_size = 15)+
        theme(axis.title.y = element_blank(),
              axis.text = element_text(family = "LM Roman 10"),
              axis.title = element_text(family = "LM Roman 10", size = 20),
              plot.title = element_text(hjust = 0.5))
    
    
    
    plot_list_SLC[[i]] <- ggplot(last_year_list_SLC[[i]], aes(x = Trait, y = 0))+
      geom_point(colour = "#158FAD",  show.legend = FALSE)  +        #aes(size=Num_Individuals),
      geom_segment(data = data.frame(x = c(-2, 0, 2), y = rep(0, 3)),
                   aes(x = x, xend = x, y = -0.1, yend = 0.1), color = "black", size = 0.5) +
      geom_text(data = data.frame(x = c(-2, 0, 2), label = c("-2", "0", "2")),
                aes(x = x, y = -0.2, label = label), vjust = 0.5, hjust = 0.5, size = 7, color = "black", family = "LM Roman 10") +
      labs(x = "Species Trait") +
      annotate("segment",x=-3,xend=3, y=0, yend=0, linewidth=1) +
      annotate("segment",x=-3,xend=-3, y=-0.1,yend=0.1, linewidth=1) +
      annotate("segment",x=3,xend=3, y=-0.1,yend=0.1, linewidth=1) +
      scale_x_continuous(limits = c(-3,3)) +
      scale_y_continuous(limits = c(-1,1)) +
      scale_color_manual(values = unname(colours)) +
      theme_tufte(base_family = "LM Roman 10", base_size = 10) +
      theme(axis.title = element_text(family = "LM Roman 10", size = 20),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank())
    
  }
  else{
    plot_list_CLC[[i]] <- ggplot(last_year_list_CLC[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
          geom_point(colour= "#158FAD", show.legend = FALSE) +   #aes(size=Num_Individuals),                               # Add points
          labs(title = substitute(sigma == value, list(value = sigmas[i]))) +                 # Labels for the axes
          scale_x_continuous(limits = c(-3,3)) +
          scale_y_continuous(limits = c(-3,3))  +
          theme_classic(base_family = "LM Roman 10", base_size = 15)+
          theme(axis.text = element_text(family = "LM Roman 10"),
                axis.title = element_text(family = "LM Roman 10", size = 20),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.title = element_text(hjust = 0.5))
    
    plot_list_SLC[[i]] <- ggplot(last_year_list_SLC[[i]], aes(x = Trait, y = 0))+
      geom_point(colour = "#158FAD", show.legend = FALSE)  + #aes(size=Num_Individuals),
      geom_segment(data = data.frame(x = c(-2, 0, 2), y = rep(0, 3)),
                   aes(x = x, xend = x, y = -0.1, yend = 0.1), color = "black", size = 0.5) +
      geom_text(data = data.frame(x = c(-2, 0, 2), label = c("-2", "0", "2")),
                aes(x = x, y = -0.2, label = label), vjust = 0.5, hjust = 0.5, size = 7, color = "black", family = "LM Roman 10") +
      annotate("segment",x=-3,xend=3, y=0, yend=0, linewidth=1) +
      annotate("segment",x=-3,xend=-3, y=-0.1,yend=0.1, linewidth=1) +
      annotate("segment",x=3,xend=3, y=-0.1,yend=0.1, linewidth=1) +
      scale_x_continuous(limits = c(-3,3)) +
      scale_y_continuous(limits = c(-1,1)) +
      scale_color_manual(values = unname(colours)) + 
      theme(panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    
  }
  
  
}


  



grid.arrange(plot_list_CLC[[1]],plot_list_CLC[[2]],plot_list_CLC[[3]],plot_list_CLC[[4]],plot_list_CLC[[5]],
             plot_list_SLC[[1]],plot_list_SLC[[2]],plot_list_SLC[[3]],plot_list_SLC[[4]],plot_list_SLC[[5]],
             nrow = 2, ncol = 5, heights = c(3, 1))




# Different mutationial probabilities (1C) ----------------------------------------

mutations <- c(0.0001, 0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040)

last_year_list_mut <- list()

for (i in 1:length(mutations)) {
  print(paste0("loop ", i, " started"))
  
  outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix, resFreq=resFreqMatrix, iniPA = 0, iniPJ = 0, 
                                    resGen=matrix(c(0.15, 0.15)), popSize = 10, mutProb=mutations[i], mutVar=0.05, time.steps = 100000)
  
  
  phenodataCLC <- NULL
  
  phenodataCLC <- data.frame(
    Year = outputCLC$phenotypes[, 1],
    Adult_Trait = outputCLC$phenotypes[, 3],
    Juvenile_Trait = outputCLC$phenotypes[, 4],
    Num_Individuals = outputCLC$phenotypes[, 2])
  
  last_year_list_mut[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  
  
}


plot_list_mut <- list()

for (i in 1:length(last_year_list_mut)){
  
  color_palette <- mako(length(last_year_list_mut[[i]]$Adult_Trait))
  
  plot_list_mut[[i]] <- ggplot(last_year_list_mut[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color_palette, show.legend = FALSE) +                                  # Add points
    labs(title = substitute(mu == value, list(value = mutations[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3)) +
    scale_y_continuous(limits = c(-3 ,3)) +
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

grid.arrange(grobs = plot_list_2rs_as, ncol = 5, nrow = 3,
             top = text_grob("Effect of mutation chance", size = 10, family = "LM Roman 10"))


