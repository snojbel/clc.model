# Updated simmer with jobs



library(job)
library(gganimate)
library(gifski)
library(png)
library(ggmatplot)
library(tidyverse)
library(viridisLite)  # Color things
library(viridis)
library(ggplot2)      # Prettier plots
library(gridExtra)  
library(ggpubr)
library(cowplot)
library(extrafont) 
library(patchwork)
library(dplyr)
library(grid)

# Resource initializations -------------------


Num.Res <- 16
res.Abund <-  50000

# Evenly distributed Resources

# SLC:
resource.freq.even.slc <- rep(1/Num.Res, times = Num.Res)                                      # res. freq. 
resource.prop.even.slc <- c(seq(from = -2.5, to = 2.5, length.out = Num.Res))            # res. property 
resource.freq.even.slc <- res.Abund*resource.freq.even.slc


# CLC:

resource.property.even.clc <- c(seq(from = -2.5, to = 2.5, length.out = Num.Res)) 

resource.frequency.even.clc <- rep(1/Num.Res, times = Num.Res)

resource.abundance.adults.even.clc     <- res.Abund                              # res. abundance of adults and juveniles
resource.abundance.juveniles.even.clc  <- res.Abund

resFreqMatrix.even.clc <- matrix(resource.frequency.even.clc, nrow=2, ncol=length(resource.frequency.even.clc ), byrow = TRUE)
resFreqMatrix.even.clc [1, ] <- resFreqMatrix.even.clc [1, ]*resource.abundance.adults.even.clc 
resFreqMatrix.even.clc [2, ] <- resFreqMatrix.even.clc [2, ]*resource.abundance.juveniles.even.clc 

resPropMatrix.even.clc <- matrix(resource.property.even.clc, nrow=2, ncol=length(resource.property.even.clc ), byrow = TRUE) 

rownames(resFreqMatrix.even.clc) <- c("Adult", "Juvenile")
colnames(resFreqMatrix.even.clc)  <- paste0("Resource ", 1:ncol(resFreqMatrix.even.clc))

rownames(resPropMatrix.even.clc)<-c("Adult", "Juvenile")
colnames(resPropMatrix.even.clc)  <- paste0("Resource ", 1:ncol(resPropMatrix.even.clc))

# Normal resources:

m <- 0 
s <- 1
N.resource.frequency <- c()
N.resource.property<- c(seq(from = -2.5, to = 2.5, length.out = Num.Res)) 

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



resource.abundance.adults.norm.clc      <- res.Abund                              # res. abundance of adults and juveniles
resource.abundance.juveniles.norm.clc   <- res.Abund

# SLC:

resource.prop.norm.slc  <- N.resource.property             # res. property 
resource.freq.norm.slc  <- res.Abund*N.resource.frequency


# CLC:

resFreqMatrix.norm.clc  <- matrix(N.resource.frequency, nrow=2, ncol=length(N.resource.frequency), byrow = TRUE)

resFreqMatrix.norm.clc [1, ] <- resFreqMatrix.norm.clc [1, ]*resource.abundance.adults.norm.clc 
resFreqMatrix.norm.clc [2, ] <- resFreqMatrix.norm.clc [2, ]*resource.abundance.juveniles.norm.clc 

rownames(resFreqMatrix.norm.clc ) <- c("Adult", "Juvenile")
colnames(resFreqMatrix.norm.clc )  <- paste0("Resource ", 1:ncol(resFreqMatrix.norm.clc ))


resPropMatrix.norm.clc  <- matrix(N.resource.property, nrow=2, ncol=length(N.resource.property), byrow = TRUE) 


rownames(resPropMatrix.norm.clc )<-c("Adult", "Juvenile")
colnames(resPropMatrix.norm.clc )  <- paste0("Resource ", 1:ncol(resPropMatrix.norm.clc))


# Skewed resource distribution

# SLC

tot <- (Num.Res*(Num.Res+1))/2

x <- 1/tot
resource.freq.skew.slc <- c()

for (i in 1:Num.Res){ 
  resource.freq.skew.slc[i] <- i*x 
}                                     # res. freq. 

resource.prop.skew.slc <- c(seq(from = -2.5, to = 2.5, length.out = Num.Res))            # res. property 
resource.freq.skew.slc <- res.Abund*resource.freq.skew.slc


# CLC:

resource.property.skew.clc <- c(seq(from = -2.5, to = 2.5, length.out = Num.Res)) 


resource.frequency.skew.clc <- c()
for (i in 1:Num.Res){ 
  resource.frequency.skew.clc[i] <- i*x 
}    


resource.abundance.adults     <- res.Abund                              # res. abundance of adults and juveniles
resource.abundance.juveniles  <- res.Abund

resFreqMatrix.skew.clc <- matrix(resource.frequency.skew.clc, nrow=2, ncol=length(resource.frequency.skew.clc), byrow = TRUE)
resFreqMatrix.skew.clc[1, ] <- resFreqMatrix.skew.clc[1, ]*resource.abundance.adults
resFreqMatrix.skew.clc[2, ] <- resFreqMatrix.skew.clc[2, ]*resource.abundance.juveniles

resPropMatrix.skew.clc <- matrix(resource.property.skew.clc, nrow=2, ncol=length(resource.property.skew.clc), byrow = TRUE) 

rownames(resFreqMatrix.skew.clc) <- c("Adult", "Juvenile")
colnames(resFreqMatrix.skew.clc)  <- paste0("Resource ", 1:ncol(resFreqMatrix.skew.clc))

rownames(resPropMatrix.skew.clc)<-c("Adult", "Juvenile")
colnames(resPropMatrix.skew.clc)  <- paste0("Resource ", 1:ncol(resPropMatrix.skew.clc))



# Two resources 

resource.prop <- c(-1,1)
resource.frequency <- c(0.5, 0.5)
resource.frequency.as <- c(0.2, 0.8)

resFreqMatrix.2res <- matrix(resource.frequency, nrow=2, ncol=length(resource.frequency), byrow = TRUE)
resFreqMatrixAs.2res <- matrix(resource.frequency.as, nrow=2, ncol=length(resource.frequency.as), byrow = TRUE)


resFreqMatrix.2res[1, ] <- resFreqMatrix.2res[1, ]*res.Abund
resFreqMatrix.2res[2, ] <- resFreqMatrix.2res[2, ]*res.Abund

resFreqMatrixAs.2res[1, ] <- resFreqMatrixAs.2res[1, ]*res.Abund
resFreqMatrixAs.2res[2, ] <- resFreqMatrixAs.2res[2, ]*res.Abund

rownames(resFreqMatrix.2res) <- c("Adult", "Juvenile")
colnames(resFreqMatrix.2res)  <- paste0("Resource ", 1:ncol(resFreqMatrixAs.2res))

rownames(resFreqMatrixAs.2res) <- c("Adult", "Juvenile")
colnames(resFreqMatrixAs.2res)  <- paste0("Resource ", 1:ncol(resFreqMatrixAs.2res))

resPropMatrix.2res <- matrix(resource.prop, nrow=2, ncol=length(resource.prop), byrow = TRUE) 


rownames(resPropMatrix.2res)<-c("Adult", "Juvenile")
colnames(resFreqMatrix.2res)  <- paste0("Resource ", 1:ncol(resPropMatrix.2res))


# ------------------------
# Running simulations -----------------


# Even distribution -----------------------

# SLC

job::job(output.Even.SLC = {
  output.Even.SLC <- resourceCompetitionSLC(resProp=resource.prop.even.slc, iniP = 0, resFreq=resource.freq.even.slc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
  
  # Control what is returned to the main session
  job::export(output.Even.SLC)
}, import = c(resource.prop.even.slc, resource.freq.even.slc, resourceCompetitionSLC))

# CLC

job::job(output.Even.CLC = {
  output.Even.CLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix.even.clc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
  
  # Control what is returned to the main session
  job::export(output.Even.CLC)
}, import = c(resPropMatrix.even.clc, resFreqMatrix.even.clc, resourceCompetitionCLC))



# Normal Distribution --------------------------

# SLC

job::job(output.Norm.SLC = {
  output.Norm.SLC <- resourceCompetitionSLC(resProp=resource.prop.norm.slc, iniP = 0, resFreq=resource.freq.norm.slc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
  
  # Control what is returned to the main session
  job::export(output.Norm.SLC)
}, import = c(resource.prop.norm.slc, resource.freq.norm.slc, resourceCompetitionSLC))

# CLC

job::job(output.Norm.CLC = {
  output.Norm.CLC <- resourceCompetitionCLC(resProp=resPropMatrix.norm.clc, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix.norm.clc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
  
  # Control what is returned to the main session
  job::export(output.Norm.CLC)
}, import = c(resPropMatrix.norm.clc, resFreqMatrix.norm.clc, resourceCompetitionCLC))




# Skewed distribution ------------------------


# SLC

job::job(output.Skew.SLC = {
  output.Skew.SLC <- resourceCompetitionSLC(resProp=resource.prop.skew.slc, iniP = 0, resFreq=resource.freq.skew.slc, 
                                            resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
  
  # Control what is returned to the main session
  job::export(output.Skew.SLC)
}, import = c(resource.prop.skew.slc, resource.freq.skew.slc, resourceCompetitionSLC))

# CLC

job::job(output.Skew.CLC = {
  output.Skew.CLC <- resourceCompetitionCLC(resProp=resPropMatrix.skew.clc, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix.skew.clc, 
                                            resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
  
  # Control what is returned to the main session
  job::export(output.Skew.CLC)
}, import = c(resPropMatrix.skew.clc, resFreqMatrix.skew.clc, resourceCompetitionCLC))


# Results ------------------------------


# Number of species at end of run


species <- matrix(data = NA, nrow = 2, ncol = 3)
rownames(species) <- c("SLC", "CLC")
colnames(species) <- c("Even", "Normal", "Skewed")

unfiltered.species <- matrix(data = NA, nrow = 2, ncol = 3)
rownames(unfiltered.species) <- c("SLC", "CLC")
colnames(unfiltered.species) <- c("Even", "Normal", "Skewed")

# Even

# SLC
output.Even.SLC.result <- output.Even.SLC$output.Even.SLC

stats.even.slc <- output.Even.SLC.result$stats
pheno.even.slc <- output.Even.SLC.result$phenotypes

filt.pheno.even.slc <- slc.groups(output = output.Even.SLC.result)  #Remove too similar "species"
species[1,1] <- as.numeric(nrow(filt.pheno.even.slc))

phenodata.SLC.even <- data.frame(
  Year = output.Even.SLC.result$phenotypes[, 1],
  Trait = output.Even.SLC.result$phenotypes[, 3],
  Num_Individuals = output.Even.SLC.result$phenotypes[, 2]
)
last.year.data.SLC.even <- phenodata.SLC.even[phenodata.SLC.even$Year == max(phenodata.SLC.even$Year), ]
unfiltered.species[1,1] <- as.numeric(nrow(last.year.data.SLC.even))

# CLC

output.Even.CLC.result <- output.Even.CLC$output.Even.CLC

stats.even.clc <- output.Even.CLC.result$stats
pheno.even.clc <- output.Even.CLC.result$phenotypes


filt.pheno.even.clc <- clc.groups(output = output.Even.CLC.result)  #Remove too similar "species"
species[2,1] <- as.numeric(nrow(filt.pheno.even.clc))

phenodata.CLC.even <- data.frame(
  Year = output.Even.CLC.result$phenotypes[, 1],
  Adult_Trait = output.Even.CLC.result$phenotypes[, 3],
  Juvenile_Trait = output.Even.CLC.result$phenotypes[, 4],
  Num_Individuals = output.Even.CLC.result$phenotypes[, 2]
)

last.year.data.CLC.even <- phenodata.CLC.even[phenodata.CLC.even$Year == max(phenodata.CLC.even$Year), ]
unfiltered.species[2,1] <- as.numeric(nrow(last.year.data.CLC.even))

# Normal

# SLC
output.Norm.SLC.result <- output.Norm.SLC$output.Norm.SLC

stats.norm.slc <- output.Norm.SLC.result$stats
pheno.norm.slc <- output.Norm.SLC.result$phenotypes

filt.pheno.norm.slc <- slc.groups(output = output.Norm.SLC.result)  #Remove too similar "species"
species[1,2] <- as.numeric(nrow(filt.pheno.norm.slc))

phenodata.SLC.norm <- data.frame(
  Year = output.Norm.SLC.result$phenotypes[, 1],
  Trait = output.Norm.SLC.result$phenotypes[, 3],
  Num_Individuals = output.Norm.SLC.result$phenotypes[, 2]
)
last.year.data.SLC.norm <- phenodata.SLC.norm[phenodata.SLC.norm$Year == max(phenodata.SLC.norm$Year), ]
unfiltered.species[1,2] <- as.numeric(nrow(last.year.data.SLC.norm))


# CLC

output.Norm.CLC.result <- output.Norm.CLC$output.Norm.CLC

stats.norm.clc <- output.Norm.CLC.result$stats
pheno.norm.clc <- output.Norm.CLC.result$phenotypes

filt.pheno.norm.clc <- clc.groups(output = output.Norm.CLC.result)  #Remove too similar "species"
species[2,2] <- as.numeric(nrow(filt.pheno.norm.clc))

phenodata.CLC.norm <- data.frame(
  Year = output.Norm.CLC.result$phenotypes[, 1],
  Adult_Trait = output.Norm.CLC.result$phenotypes[, 3],
  Juvenile_Trait = output.Norm.CLC.result$phenotypes[, 4],
  Num_Individuals = output.Norm.CLC.result$phenotypes[, 2]
)

last.year.data.CLC.norm <- phenodata.CLC.norm[phenodata.CLC.norm$Year == max(phenodata.CLC.norm$Year), ]
unfiltered.species[2,2] <- as.numeric(nrow(last.year.data.CLC.norm))


# Skewed


# SLC
output.Skew.SLC.result <- output.Skew.SLC$output.Skew.SLC

stats.skew.slc <- output.Skew.SLC.result$stats
pheno.skew.slc <- output.Skew.SLC.result$phenotypes

filt.pheno.skew.slc <- slc.groups(output = output.Skew.SLC.result)  #Remove too similar "species"
species[1,3]<- as.numeric(nrow(filt.pheno.skew.slc))


phenodata.SLC.skew <- data.frame(
  Year = output.Skew.SLC.result$phenotypes[, 1],
  Trait = output.Skew.SLC.result$phenotypes[, 3],
  Num_Individuals = output.Skew.SLC.result$phenotypes[, 2]
)
last.year.data.SLC.skew <- phenodata.SLC.skew[phenodata.SLC.skew$Year == max(phenodata.SLC.skew$Year), ]
unfiltered.species[1,3] <- as.numeric(nrow(last.year.data.SLC.skew))

# CLC

output.Skew.CLC.result <- output.Skew.CLC$output.Skew.CLC

stats.skew.clc <- output.Skew.CLC.result$stats
pheno.skew.clc <- output.Skew.CLC.result$phenotypes

filt.pheno.skew.clc <- clc.groups(output = output.Skew.CLC.result)  #Remove too similar "species"
species[2,3] <- as.numeric(nrow(filt.pheno.skew.clc))

phenodata.CLC.skew <- data.frame(
  Year = output.Skew.CLC.result$phenotypes[, 1],
  Adult_Trait = output.Skew.CLC.result$phenotypes[, 3],
  Juvenile_Trait = output.Skew.CLC.result$phenotypes[, 4],
  Num_Individuals = output.Skew.CLC.result$phenotypes[, 2]
)

last.year.data.CLC.skew <- phenodata.CLC.skew[phenodata.CLC.skew$Year == max(phenodata.CLC.skew$Year), ]
unfiltered.species[2,3] <- as.numeric(nrow(last.year.data.CLC.skew))

# Plotting --------------------------------------


x <- colnames(species)

species.Frame <- data.frame(data.frame(
  SLC = species[1, ],
  CLC = species[2, ],
  dist = c(x)
))


ggplot(data = species.Frame) +
  geom_point(aes(x = dist, y = SLC, shape = "SLC"),color = "lightblue", size = 8) +
  geom_point(aes(x = dist, y = CLC, shape = "CLC"), color = "pink", size = 8) + # Add points
  scale_y_continuous(limits = c(0, 25)) +
  labs(title = " Filtered", x = "Distribution Type", y = "Number of Species", family = "LM Roman 10") +  
  scale_shape_manual(values = c("SLC" = 15, "CLC" = 17), 
                     labels = c("Simple", "Complex")) +  
  theme_classic(base_size = 18)+ 
  theme(axis.text = element_text(family = "LM Roman 10"),
        axis.title = element_text(family = "LM Roman 10", size = 20),
        plot.title = element_text(hjust = 0.5, family = "LM Roman 10"),
        legend.text = element_text(family = "LM Roman 10"),
        legend.title = element_blank()) +
  guides(shape = guide_legend(override.aes = list(color = c("lightblue", "pink"), shape = c(15, 17), size = 7), title = NULL))

# Unfiltered

x <- colnames(species)

unf.species.Frame <- data.frame(data.frame(
  SLC = unfiltered.species[1, ],
  CLC = unfiltered.species[2, ],
  dist = c(x)
))


ggplot(data = unf.species.Frame) +
  geom_point(aes(x = dist, y = SLC, shape = "SLC"),color = "lightblue", size = 8) +
  geom_point(aes(x = dist, y = CLC, shape = "CLC"), color = "pink", size = 8) + # Add points
  scale_y_continuous(limits = c(0, 35)) +
  labs(title = " Unfiltered", x = "Distribution Type", y = "Number of Species", family = "LM Roman 10") +  
  scale_shape_manual(values = c("SLC" = 15, "CLC" = 17), 
                     labels = c("Simple", "Complex")) +  
  theme_classic(base_size = 18)+ 
  theme(axis.text = element_text(family = "LM Roman 10"),
        axis.title = element_text(family = "LM Roman 10", size = 20),
        plot.title = element_text(hjust = 0.5, family = "LM Roman 10"),
        legend.text = element_text(family = "LM Roman 10"),
        legend.title = element_blank()) +
  guides(shape = guide_legend(override.aes = list(color = c("lightblue", "pink"), shape = c(15, 17), size = 7), title = NULL))




# ---------------------
# Running simulations: 9 run to see endpoint -----------------------------

# Normal

job::job(endpoint.normal = {
  
  last.year.list.norm <- list()
  
  for(i in 1:9){
    print(paste0("loop ", i, " started"))
    outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.norm.clc, resFreq=resFreqMatrix.norm.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(0.15, 0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
    
    phenodataCLC <- NULL
    
    phenodataCLC <- data.frame(
      Year = outputCLC$phenotypes[, 1],
      Adult_Trait = outputCLC$phenotypes[, 3],
      Juvenile_Trait = outputCLC$phenotypes[, 4],
      Num_Individuals = outputCLC$phenotypes[, 2])
    
    last.year.list.norm[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  }
  
  
  # Control what is returned to the main session
  job::export(last.year.list.norm)
}, import = c(resPropMatrix.norm.clc, resFreqMatrix.norm.clc, resourceCompetitionCLC))


# Even

job::job(endpoint.even = {
  
  last.year.list.even <- list()
  
  for(i in 1:9){
    print(paste0("loop ", i, " started"))
    outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, resFreq=resFreqMatrix.even.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(0.15, 0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
    
    phenodataCLC <- NULL
    
    phenodataCLC <- data.frame(
      Year = outputCLC$phenotypes[, 1],
      Adult_Trait = outputCLC$phenotypes[, 3],
      Juvenile_Trait = outputCLC$phenotypes[, 4],
      Num_Individuals = outputCLC$phenotypes[, 2])
    
    last.year.list.even[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  }
  
  
  # Control what is returned to the main session
  job::export(last.year.list.even)
}, import = c(resPropMatrix.even.clc, resFreqMatrix.even.clc, resourceCompetitionCLC))

# Skewed

job::job(endpoint.skew = {
  
  last.year.list.skew <- list()
  
  for(i in 1:9){
    print(paste0("loop ", i, " started"))
    outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.skew.clc, resFreq=resFreqMatrix.skew.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(0.15, 0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
    
    phenodataCLC <- NULL
    
    phenodataCLC <- data.frame(
      Year = outputCLC$phenotypes[, 1],
      Adult_Trait = outputCLC$phenotypes[, 3],
      Juvenile_Trait = outputCLC$phenotypes[, 4],
      Num_Individuals = outputCLC$phenotypes[, 2])
    
    last.year.list.skew[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  }
  
  
  # Control what is returned to the main session
  job::export(last.year.list.skew)
}, import = c(resPropMatrix.skew.clc, resFreqMatrix.skew.clc, resourceCompetitionCLC))

# Results


last.year.list.even <- endpoint.even$last.year.list.even
last.year.list.norm <- endpoint.normal$last.year.list.norm
last.year.list.skew <- endpoint.skew$last.year.list.skew




# Plotting 9 runs to see endpoint --------------------------

# Even

plot.list.even <- list()

for (i in 1:9){
  
  color.palette <- mako(length(last.year.list.even[[i]]$Adult_Trait))
  
  plot.list.even[[i]] <- ggplot(last.year.list.even[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) +                                  # Add points
    labs(x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

grid.arrange(grobs = plot.list.even, ncol = 3, nrow = 3, top =text_grob("Even Distribution 50 000 years", size = 10, family = "LM Roman 10"))

# Normal

plot.list.norm <- list()


for (i in 1:3){
  
  color.palette <- mako(length(last.year.list.norm[[i]]$Adult_Trait))
  
  plot.list.norm[[i]] <- ggplot(last.year.list.norm[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) +                                  # Add points
    labs(x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

grid.arrange(grobs = plot.list.norm, ncol = 3, top =text_grob("Normal Distribution 50 000 years", size = 10, family = "LM Roman 10"))

# Skewed

plot.list.skew <- list()

for (i in 1:9){
  
  color.palette <- mako(length(last.year.list.skew[[i]]$Adult_Trait))
  
  plot.list.skew[[i]] <- ggplot(last.year.list.skew[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) +                                  # Add points
    labs(x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

grid.arrange(grobs = plot.list.skew, ncol = 3, nrow = 3, top =text_grob("Skewed Distribution 50 000 years", size = 10, family = "LM Roman 10"))



#-------------------------
# Running simulations to see endpoint and number of species with varied sigma -----------------------------

# Normal

job::job(endpoint.normal.sigma = {
  
  last.year.list.norm <- list()
  filtered.list.norm <- list()
  sigma <- c(0.15, 0.3, 0.45, 0.6, 0.75, 0.90)
  
  for(i in 1:length(sigma)){
    print(paste0("loop ", i, " started"))
    outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.norm.clc, resFreq=resFreqMatrix.norm.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i], sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
    
    phenodataCLC <- NULL
    
    phenodataCLC <- data.frame(
      Year = outputCLC$phenotypes[, 1],
      Adult_Trait = outputCLC$phenotypes[, 3],
      Juvenile_Trait = outputCLC$phenotypes[, 4],
      Num_Individuals = outputCLC$phenotypes[, 2])
    
    last.year.list.norm[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
    filtered.list.norm[[i]] <- clc.groups(output = outputCLC)
  }
  
  
  # Control what is returned to the main session
  job::export(list(last.year.list.norm, filtered.list.norm))
}, import = c(resPropMatrix.norm.clc, resFreqMatrix.norm.clc, resourceCompetitionCLC, clc.groups))


# Even

job::job(endpoint.even.sigma = {
  
  last.year.list.even <- list()
  filtered.list.even <- list()
  sigma <- c(0.15, 0.3, 0.45, 0.6, 0.75, 0.90)
  
  for(i in 1:length(sigma)){
    print(paste0("loop ", i, " started"))
    outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, resFreq=resFreqMatrix.even.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i], sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
    
    phenodataCLC <- NULL
    
    phenodataCLC <- data.frame(
      Year = outputCLC$phenotypes[, 1],
      Adult_Trait = outputCLC$phenotypes[, 3],
      Juvenile_Trait = outputCLC$phenotypes[, 4],
      Num_Individuals = outputCLC$phenotypes[, 2])
    
    last.year.list.even[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
    filtered.list.even[[i]] <- clc.groups(output = outputCLC)
  }
  
  
  # Control what is returned to the main session
  job::export(list(last.year.list.even, filtered.list.even))
}, import = c(resPropMatrix.even.clc, resFreqMatrix.even.clc, resourceCompetitionCLC, clc.groups))

# Skewed

job::job(endpoint.skew.sigma = {
  
  last.year.list.skew <- list()
  filtered.list.skew <- list()
  sigma <- c(0.15, 0.3, 0.45, 0.6, 0.75, 0.90)
  
  for(i in 1:length(sigma)){
    print(paste0("loop ", i, " started"))
    outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.skew.clc, resFreq=resFreqMatrix.skew.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i], sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 50000)
    
    phenodataCLC <- NULL
    
    phenodataCLC <- data.frame(
      Year = outputCLC$phenotypes[, 1],
      Adult_Trait = outputCLC$phenotypes[, 3],
      Juvenile_Trait = outputCLC$phenotypes[, 4],
      Num_Individuals = outputCLC$phenotypes[, 2])
    
    last.year.list.skew[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
    filtered.list.skew[[i]] <- clc.groups(output = outputCLC)
  }
  
  
  # Control what is returned to the main session
  job::export(list(last.year.list.skew, filtered.list.skew))
}, import = c(resPropMatrix.skew.clc, resFreqMatrix.skew.clc, resourceCompetitionCLC, clc.groups))



# Results

last.year.list.even <- endpoint.even.sigma$last.year.list.even
last.year.list.norm <- endpoint.normal.sigma$last.year.list.norm
last.year.list.skew <- endpoint.skew.sigma$last.year.list.skew

filtered.list.even <- endpoint.even.sigma$filtered.list.even
filtered.list.norm <- endpoint.normal.sigma$filtered.list.norm
filtered.list.skew <- endpoint.skew.sigma$filtered.list.skew

sigma <- c(0.15, 0.3, 0.45, 0.6, 0.75, 0.90)

num.of.spe.even <- c()
num.of.spe.norm <- c()
num.of.spe.skew <- c()

for (i in 1:length(sigma)){
  num.of.spe.even[i] <- nrow(filtered.list.even[[i]])
  num.of.spe.norm[i] <- nrow(filtered.list.norm[[i]])
  num.of.spe.skew[i] <- nrow(filtered.list.skew[[i]])
}



# Plotting Abundance vs Trait of Juveniles and Adults ----------------------


combined.even.adult <- list()
combined.norm.adult <- list()
combined.skew.adult <- list()

combined.even.juvenile <- list()
combined.norm.juvenile <- list()
combined.skew.juvenile <- list()

# -------------------------------  Even 


bins <- seq(from = -2.5, to = 2.5, by = 0.5)
mids <- c()
for(i in 1:(length(bins)-1)){
  mids[i] <- (bins[i]+bins[i+1])/2
}



for(i in 1:length(last.year.list.even)){
  last.year.list.even[[i]]$adult.cut <- cut(last.year.list.even[[i]]$Adult_Trait, breaks = bins, labels = F, include.lowest = TRUE)
  last.year.list.even[[i]]$juvenile.cut <- cut(last.year.list.even[[i]]$Juvenile_Trait, breaks = bins, labels = F, include.lowest = TRUE)
  
  # Group by the bin
  grouped.data.adult <- group_by(last.year.list.even[[i]], adult.cut)
  grouped.data.juvenile<- group_by(last.year.list.even[[i]], juvenile.cut)
  
  summarized.data.adult <- c()
  
  # Summarize the abundance
  summarized.data.adult <- summarise(grouped.data.adult, total_abundance = sum(Num_Individuals))
  summarized.data.juvenile <- summarise(grouped.data.juvenile, total_abundance = sum(Num_Individuals))
  
  
  # Convert the summarized data to a dataframe
  
  combined.even.adult[[i]] <- as.data.frame(summarized.data.adult)
  combined.even.juvenile[[i]] <- as.data.frame(summarized.data.juvenile)
}

# Fix first column 

for(i in 1:length(combined.even.adult)){
  for(k in 1:length(combined.even.adult[[i]][,1])){
    combined.even.adult[[i]][k, 1] <- mids[combined.even.adult[[i]][k, 1]]
  }
}

for(i in 1:length(combined.even.juvenile)){
  for(k in 1:length(combined.even.juvenile[[i]][,1])){
    combined.even.juvenile[[i]][k, 1] <- mids[combined.even.juvenile[[i]][k, 1]]
  }
}


plot.list.even.adult <- list()
plot.list.even.juvenile <- list()

for (i in 1:length(sigma)){
  
  
  plot.list.even.adult[[i]] <- ggplot(combined.even.adult[[i]], aes(x = adult.cut, y = total_abundance)) +
    geom_bar(stat = "identity", width = 0.4, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Adult Trait", y = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    #scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  #color.palette <- mako(length(last.year.list.even[[i]]$Juvenile_Trait))
  
  plot.list.even.juvenile[[i]] <- ggplot(combined.even.juvenile[[i]], aes(x = juvenile.cut, y = total_abundance)) +
    geom_bar(stat = "identity", width = 0.4, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    #scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}


layout <- "
ABC
DEF
#G#
HIJ
KLM
"

combo.plot.list <- list()
midtitle <- textGrob("Juvenile Trait", gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.even.adult)*2 + 1)){
  if(i == (length(plot.list.even.adult) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.even.adult) + 1)){
    combo.plot.list[[i]] <- plot.list.even.adult[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.list.even.juvenile[[i-(length(plot.list.even.adult) + 1)]]
  }
}

plots <- wrap_plots(combo.plot.list, design = layout)

plots + plot_annotation(
  title = 'Even Distribution',
  subtitle = 'Adult Trait',
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  #caption = 'Disclaimer: None of these plots are insightful'
) + plot_layout(heights = c(1, 1, 0.4, 1, 1))


# Norm
for(i in 1:length(last.year.list.norm)){
  last.year.list.norm[[i]]$adult.cut <- cut(last.year.list.norm[[i]]$Adult_Trait, breaks = bins, labels = F, include.lowest = TRUE)
  last.year.list.norm[[i]]$juvenile.cut <- cut(last.year.list.norm[[i]]$Juvenile_Trait, breaks = bins, labels = F, include.lowest = TRUE)
  
  # Group by the bin
  grouped.data.adult <- group_by(last.year.list.norm[[i]], adult.cut)
  grouped.data.juvenile<- group_by(last.year.list.norm[[i]], juvenile.cut)
  
  summarized.data.adult <- c()
  
  # Summarize the abundance
  summarized.data.adult <- summarise(grouped.data.adult, total_abundance = sum(Num_Individuals))
  summarized.data.juvenile <- summarise(grouped.data.juvenile, total_abundance = sum(Num_Individuals))
  
  
  # Convert the summarized data to a dataframe
  
  combined.norm.adult[[i]] <- as.data.frame(summarized.data.adult)
  combined.norm.juvenile[[i]] <- as.data.frame(summarized.data.juvenile)
}

# Fix first column 

for(i in 1:length(combined.norm.adult)){
  for(k in 1:length(combined.norm.adult[[i]][,1])){
    combined.norm.adult[[i]][k, 1] <- mids[combined.norm.adult[[i]][k, 1]]
  }
}

for(i in 1:length(combined.norm.juvenile)){
  for(k in 1:length(combined.norm.juvenile[[i]][,1])){
    combined.norm.juvenile[[i]][k, 1] <- mids[combined.norm.juvenile[[i]][k, 1]]
  }
}


plot.list.norm.adult <- list()
plot.list.norm.juvenile <- list()

for (i in 1:length(sigma)){
  
  
  plot.list.norm.adult[[i]] <- ggplot(combined.norm.adult[[i]], aes(x = adult.cut, y = total_abundance)) +
    geom_bar(stat = "identity", width = 0.4, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Adult Trait", y = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    #scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  #color.palette <- mako(length(last.year.list.norm[[i]]$Juvenile_Trait))
  
  plot.list.norm.juvenile[[i]] <- ggplot(combined.norm.juvenile[[i]], aes(x = juvenile.cut, y = total_abundance)) +
    geom_bar(stat = "identity", width = 0.4, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    #scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}


layout <- "
ABC
DEF
#G#
HIJ
KLM
"

combo.plot.list <- list()
midtitle <- textGrob("Juvenile Trait", gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.norm.adult)*2 + 1)){
  if(i == (length(plot.list.norm.adult) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.norm.adult) + 1)){
    combo.plot.list[[i]] <- plot.list.norm.adult[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.list.norm.juvenile[[i-(length(plot.list.norm.adult) + 1)]]
  }
}

plots <- wrap_plots(combo.plot.list, design = layout)

plots + plot_annotation(
  title = 'Normal Distribution',
  subtitle = 'Adult Trait',
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  #caption = 'Disclaimer: None of these plots are insightful'
) + plot_layout(heights = c(1, 1, 0.4, 1, 1))



# Skewed



for(i in 1:length(last.year.list.skew)){
  last.year.list.skew[[i]]$adult.cut <- cut(last.year.list.skew[[i]]$Adult_Trait, breaks = bins, labels = F, include.lowest = TRUE)
  last.year.list.skew[[i]]$juvenile.cut <- cut(last.year.list.skew[[i]]$Juvenile_Trait, breaks = bins, labels = F, include.lowest = TRUE)
  
  # Group by the bin
  grouped.data.adult <- group_by(last.year.list.skew[[i]], adult.cut)
  grouped.data.juvenile<- group_by(last.year.list.skew[[i]], juvenile.cut)
  
  summarized.data.adult <- c()
  
  # Summarize the abundance
  summarized.data.adult <- summarise(grouped.data.adult, total_abundance = sum(Num_Individuals))
  summarized.data.juvenile <- summarise(grouped.data.juvenile, total_abundance = sum(Num_Individuals))
  
  
  # Convert the summarized data to a dataframe
  
  combined.skew.adult[[i]] <- as.data.frame(summarized.data.adult)
  combined.skew.juvenile[[i]] <- as.data.frame(summarized.data.juvenile)
}

# Fix first column 

for(i in 1:length(combined.skew.adult)){
  for(k in 1:length(combined.skew.adult[[i]][,1])){
    combined.skew.adult[[i]][k, 1] <- mids[combined.skew.adult[[i]][k, 1]]
  }
}

for(i in 1:length(combined.skew.juvenile)){
  for(k in 1:length(combined.skew.juvenile[[i]][,1])){
    combined.skew.juvenile[[i]][k, 1] <- mids[combined.skew.juvenile[[i]][k, 1]]
  }
}


plot.list.skew.adult <- list()
plot.list.skew.juvenile <- list()

for (i in 1:length(sigma)){
  
  
  plot.list.skew.adult[[i]] <- ggplot(combined.skew.adult[[i]], aes(x = adult.cut, y = total_abundance)) +
    geom_bar(stat = "identity", width = 0.4, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Adult Trait", y = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    #scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  #color.palette <- mako(length(last.year.list.skew[[i]]$Juvenile_Trait))
  
  plot.list.skew.juvenile[[i]] <- ggplot(combined.skew.juvenile[[i]], aes(x = juvenile.cut, y = total_abundance)) +
    geom_bar(stat = "identity", width = 0.4, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    #scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}


layout <- "
ABC
DEF
#G#
HIJ
KLM
"

combo.plot.list <- list()
midtitle <- textGrob("Juvenile Trait", gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.skew.adult)*2 + 1)){
  if(i == (length(plot.list.skew.adult) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.skew.adult) + 1)){
    combo.plot.list[[i]] <- plot.list.skew.adult[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.list.skew.juvenile[[i-(length(plot.list.skew.adult) + 1)]]
  }
}

plots <- wrap_plots(combo.plot.list, design = layout)

plots + plot_annotation(
  title = 'Skewed Distribution',
  subtitle = 'Adult Trait',
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  #caption = 'Disclaimer: None of these plots are insightful'
) + plot_layout(heights = c(1, 1, 0.4, 1, 1))

# Plotting 9 runs to see endpoint comparison of filtered vs unfiltered --------------------------

# Even

plot.list.even <- list()
plot.filtered.list.even <- list()

for (i in 1:length(sigma)){
  
  color.palette <- mako(length(last.year.list.even[[i]]$Adult_Trait))
  
  plot.list.even[[i]] <- ggplot(last.year.list.even[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-2.5, 2.5))+
    scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  color.palette <- mako(length(filtered.list.even[[i]]$Adult_Trait))
  
  plot.filtered.list.even[[i]] <- ggplot(filtered.list.even[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(color = color.palette, show.legend = FALSE, size = 4) +                                  # Add points
    labs(title = substitute(sigma == value, list(value = sigma[i])), subtitle = substitute("Species" == spec, list(spec = num.of.spe.even[i])),
         x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-2.5, 2.5))+
    scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

 

combo.plot.list <- list()
midtitle <- textGrob("Filtered Endpoint", gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.even)*2 + 1)){
  if(i == (length(plot.list.even) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.even) + 1)){
    combo.plot.list[[i]] <- plot.list.even[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.filtered.list.even[[i-(length(plot.list.even) + 1)]]
  }
}





layout <- "
ABC
DEF
#G#
HIJ
KLM
"

plots <- wrap_plots(combo.plot.list, design = layout)


plots + plot_annotation(
  title = 'Even Distribution',
  subtitle = 'Unfiltered Endpoint',
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  #caption = 'Disclaimer: None of these plots are insightful'
)+ plot_layout(heights = c(1, 1, 0.4, 1, 1))



# Normal

plot.list.norm <- list()
plot.filtered.list.norm <- list()

for (i in 1:length(sigma)){
  
  color.palette <- mako(length(last.year.list.norm[[i]]$Adult_Trait))
  
  plot.list.norm[[i]] <- ggplot(last.year.list.norm[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-2.5, 2.5))+
    scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  color.palette <- mako(length(filtered.list.norm[[i]]$Adult_Trait))
  
  plot.filtered.list.norm[[i]] <- ggplot(filtered.list.norm[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(color = color.palette, show.legend = FALSE, size = 4) +                                  # Add points
    labs(title = substitute(sigma == value, list(value = sigma[i])), subtitle = substitute("Species" == spec, list(spec = num.of.spe.norm[i])),
         x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-2.5, 2.5))+
    scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}



combo.plot.list <- list()
midtitle <- textGrob("Filtered Endpoint", gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.norm)*2 + 1)){
  if(i == (length(plot.list.norm) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.norm) + 1)){
    combo.plot.list[[i]] <- plot.list.norm[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.filtered.list.norm[[i-(length(plot.list.norm) + 1)]]
  }
}





layout <- "
ABC
DEF
#G#
HIJ
KLM
"

plots <- wrap_plots(combo.plot.list, design = layout)


plots + plot_annotation(
  title = 'Normal Distribution',
  subtitle = 'Unfiltered Endpoint',
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  #caption = 'Disclaimer: None of these plots are insightful'
)+ plot_layout(heights = c(1, 1, 0.4, 1, 1))



# Skewed

plot.list.skew <- list()
plot.filtered.list.skew<- list()

for (i in 1:length(sigma)){
  
  color.palette <- mako(length(last.year.list.skew[[i]]$Adult_Trait))
  
  plot.list.skew[[i]] <- ggplot(last.year.list.skew[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = sigma[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-2.5, 2.5))+
    scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
  
  color.palette <- mako(length(filtered.list.skew[[i]]$Adult_Trait))
  
  plot.filtered.list.skew[[i]] <- ggplot(filtered.list.skew[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(color = color.palette, show.legend = FALSE, size = 4) +                                  # Add points
    labs(title = substitute(sigma == value, list(value = sigma[i])), subtitle = substitute("Species" == spec, list(spec = num.of.spe.skew[i])),
         x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-2.5, 2.5))+
    scale_y_continuous(limits = c(-2.5, 2.5))+
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}



combo.plot.list <- list()
midtitle <- textGrob("Filtered Endpoint", gp = gpar(fontsize = 15, fontfamily = "LM Roman 10"),
                     hjust = 0.5)

for(i in 1:(length(plot.list.skew)*2 + 1)){
  if(i == (length(plot.list.skew) + 1)){
    combo.plot.list[[i]] <- midtitle
  }
  else if(i < (length(plot.list.norm) + 1)){
    combo.plot.list[[i]] <- plot.list.skew[[i]]
  }
  else{
    combo.plot.list[[i]] <- plot.filtered.list.skew[[i-(length(plot.list.skew) + 1)]]
  }
}





layout <- "
ABC
DEF
#G#
HIJ
KLM
"

plots <- wrap_plots(combo.plot.list, design = layout)


plots + plot_annotation(
  title = 'Skewed Distribution',
  subtitle = 'Unfiltered Endpoint',
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  #caption = 'Disclaimer: None of these plots are insightful'
)+ plot_layout(heights = c(1, 1, 0.4, 1, 1))

# ------------------------
# Running simulations: compare adult population vs Juvenile population -----------------------------

# Normal

job::job(population.normal = {
  
  #SLC
  sigma <- c(seq(from = 0.05, to = 1.05, by = 0.20))
  
  adult.last.year.norm.SLC <- matrix(data = NA, nrow = 1, ncol = length(sigma))
  colnames(adult.last.year.norm.SLC) <- sigma #JUVENILES
  
  juvenile.last.year.norm.SLC <- matrix(data = NA, nrow = 1, ncol = length(sigma))
  colnames(juvenile.last.year.norm.SLC) <- sigma #JUVENILES
  
  
  for(i in 1:length(sigma)){
    print(paste0("loop ", i, " started"))
    
      outputSLC <- resourceCompetitionSLC(resProp=resource.prop.norm.slc, resFreq=resource.freq.norm.slc, iniP = 0, resGen=matrix(c(sigma[i], sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
      
      
      adult.last.year.norm.SLC[1, i] <- outputSLC$stats[outputSLC$stats[,1] == max(outputSLC$stats[,1]), 2]
      juvenile.last.year.norm.SLC[1, i] <- outputSLC$stats[outputSLC$stats[,1] == max(outputSLC$stats[,1]), 3]
  
    
  }
  # CLC
  
  adult.last.year.norm.CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(adult.last.year.norm.CLC) <- sigma #ADULTS
  colnames(adult.last.year.norm.CLC) <- sigma #JUVENILES
  
  juvenile.last.year.norm.CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(juvenile.last.year.norm.CLC) <- sigma #ADULTS
  colnames(juvenile.last.year.norm.CLC) <- sigma #JUVENILES
  
  
  for(i in 1:length(sigma)){
    print(paste0("loop ", i, " started"))
    for(c in 1:length(sigma)) {
        outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.norm.clc, resFreq=resFreqMatrix.norm.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i], sigma[c])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
        
        adult.last.year.norm.CLC[i, c] <- outputCLC$stats[outputCLC$stats[,1] == max(outputCLC$stats[,1]), 2]
        juvenile.last.year.norm.CLC[i, c] <- outputCLC$stats[outputCLC$stats[,1] == max(outputCLC$stats[,1]), 3]
        
    }
   
  }
  
  
  # Control what is returned to the main session
  job::export(list(adult.last.year.norm.CLC, juvenile.last.year.norm.CLC, adult.last.year.norm.SLC, juvenile.last.year.norm.SLC))
}, import = c(resPropMatrix.norm.clc, resFreqMatrix.norm.clc, resourceCompetitionCLC, resource.prop.norm.slc, resource.freq.norm.slc, resourceCompetitionSLC))


# Even

job::job(population.even = {
  
  #SLC
  sigma <- c(seq(from = 0.05, to = 1.05, by = 0.20))
  
  adult.last.year.even.SLC <- matrix(data = NA, nrow = 1, ncol = length(sigma))
  colnames(adult.last.year.even.SLC) <- sigma #JUVENILES
  
  juvenile.last.year.even.SLC <- matrix(data = NA, nrow = 1, ncol = length(sigma))
  colnames(juvenile.last.year.even.SLC) <- sigma #JUVENILES
  
  
  for(i in 1:length(sigma)){
    print(paste0("loop ", i, " started"))
    
    outputSLC <- resourceCompetitionSLC(resProp=resource.prop.even.slc, resFreq=resource.freq.even.slc, iniP = 0, resGen=matrix(c(sigma[i], sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
    
    
    adult.last.year.even.SLC[1, i] <- outputSLC$stats[outputSLC$stats[,1] == max(outputSLC$stats[,1]), 2]
    juvenile.last.year.even.SLC[1, i] <- outputSLC$stats[outputSLC$stats[,1] == max(outputSLC$stats[,1]), 3]
    
    
  }
  # CLC
  
  adult.last.year.even.CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(adult.last.year.even.CLC) <- sigma #ADULTS
  colnames(adult.last.year.even.CLC) <- sigma #JUVENILES
  
  juvenile.last.year.even.CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(juvenile.last.year.even.CLC) <- sigma #ADULTS
  colnames(juvenile.last.year.even.CLC) <- sigma #JUVENILES
  
  
  for(i in 1:length(sigma)){
    print(paste0("loop ", i, " started"))
    for(c in 1:length(sigma)) {
      outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, resFreq=resFreqMatrix.even.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i], sigma[c])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
      
      adult.last.year.even.CLC[i, c] <- outputCLC$stats[outputCLC$stats[,1] == max(outputCLC$stats[,1]), 2]
      juvenile.last.year.even.CLC[i, c] <- outputCLC$stats[outputCLC$stats[,1] == max(outputCLC$stats[,1]), 3]
      
    
    }
    
  }
  
  
  # Control what is returned to the main session
  job::export(list(adult.last.year.even.CLC, juvenile.last.year.even.CLC, adult.last.year.even.SLC, juvenile.last.year.even.SLC))
}, import = c(resPropMatrix.even.clc, resFreqMatrix.even.clc, resourceCompetitionCLC, resource.prop.even.slc, resource.freq.even.slc, resourceCompetitionSLC))


# Skewed

job::job(population.skew = {
  
  #SLC
  sigma <- c(seq(from = 0.05, to = 1.05, by = 0.20))
  
  adult.last.year.skew.SLC <- matrix(data = NA, nrow = 1, ncol = length(sigma))
  colnames(adult.last.year.skew.SLC) <- sigma #JUVENILES
  
  juvenile.last.year.skew.SLC <- matrix(data = NA, nrow = 1, ncol = length(sigma))
  colnames(juvenile.last.year.skew.SLC) <- sigma #JUVENILES
  
  
  for(i in 1:length(sigma)){
    print(paste0("loop ", i, " started"))
    
    outputSLC <- resourceCompetitionSLC(resProp=resource.prop.skew.slc, resFreq=resource.freq.skew.slc, iniP = 0, resGen=matrix(c(sigma[i], sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
    
    adult.last.year.skew.SLC[1, i] <- outputSLC$stats[outputSLC$stats[,1] == max(outputSLC$stats[,1]), 2]
    juvenile.last.year.skew.SLC[1, i] <- outputSLC$stats[outputSLC$stats[,1] == max(outputSLC$stats[,1]), 3]
    
    
  }
  # CLC
  
  adult.last.year.skew.CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(adult.last.year.skew.CLC) <- sigma #ADULTS
  colnames(adult.last.year.skew.CLC) <- sigma #JUVENILES
  
  juvenile.last.year.skew.CLC <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(juvenile.last.year.skew.CLC) <- sigma #ADULTS
  colnames(juvenile.last.year.skew.CLC) <- sigma #JUVENILES
  
  
  for(i in 1:length(sigma)){
    print(paste0("loop ", i, " started"))
    for(c in 1:length(sigma)) {
      outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.skew.clc, resFreq=resFreqMatrix.skew.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i], sigma[c])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
      
      adult.last.year.skew.CLC[i, c] <- outputCLC$stats[outputCLC$stats[,1] == max(outputCLC$stats[,1]), 2]
      juvenile.last.year.skew.CLC[i, c] <- outputCLC$stats[outputCLC$stats[,1] == max(outputCLC$stats[,1]), 3]
    }
    
  }
  
  
  # Control what is returned to the main session
  job::export(list(adult.last.year.skew.CLC, juvenile.last.year.skew.CLC, adult.last.year.skew.SLC, juvenile.last.year.skew.SLC))
}, import = c(resPropMatrix.skew.clc, resFreqMatrix.skew.clc, resourceCompetitionCLC, resource.prop.skew.slc, resource.freq.skew.slc, resourceCompetitionSLC))



# Results ------------------------


adult.norm.pop.clc <- population.normal$adult.last.year.norm.CLC
juvenile.norm.pop.clc <- population.normal$juvenile.last.year.norm.CLC
adult.norm.pop.slc <- population.normal$adult.last.year.norm.SLC
juvenile.norm.pop.slc <- population.normal$juvenile.last.year.norm.SLC


adult.even.pop.clc <- population.even$adult.last.year.even.CLC
juvenile.even.pop.clc <- population.even$juvenile.last.year.even.CLC
adult.even.pop.slc <- population.even$adult.last.year.even.SLC
juvenile.even.pop.slc <- population.even$juvenile.last.year.even.SLC

adult.skew.pop.clc <- population.skew$adult.last.year.skew.CLC
juvenile.skew.pop.clc <- population.skew$juvenile.last.year.skew.CLC
adult.skew.pop.slc <- population.skew$adult.last.year.skew.SLC
juvenile.skew.pop.slc <- population.skew$juvenile.last.year.skew.SLC

# Plotting abundance of adults and juveniles different sigmas --------------------------

# Normal

x <- rownames(adult.norm.pop.clc)

# Number of species
  

df_adult <- data.frame(
  Juvenile.trait = rep(x, each = 6),
  Abundance = as.vector(adult.norm.pop.clc),
  Adult.trait = rep(x, times = 6),
  Stage = rep("Adult", times = length(x)*length(x)),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df_juvenile <- data.frame(
  Juvenile.trait = rep(x, each = 6),
  Abundance = as.vector(juvenile.norm.pop.clc),
  Adult.trait = rep(x, times = 6),
  Stage = rep("Juvenile", times = length(x)*length(x)),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df_simple_adult <- data.frame(
  Juvenile.trait = x ,
  Abundance = as.vector(adult.norm.pop.slc),
  Adult.trait = x,
  Stage = rep("Adult", times = length(x)),
  Cycle = rep("Simple", times = length(x))
)

df_simple_juvenile <- data.frame(
  Juvenile.trait =x ,
  Abundance = as.vector(juvenile.norm.pop.slc),
  Adult.trait = x,
  Stage = rep("Juvenile", times = length(x)),
  Cycle = rep("Simple", times = length(x))
)



df_combined_complex <- rbind(df_adult, df_juvenile)
df_combined_simple <- rbind(df_simple_adult, df_simple_juvenile)


com <- ggplot(df_combined_complex, aes(x = Adult.trait, y = Abundance, color = Stage, shape = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 7) +
  scale_y_continuous(limits = c(20000, 60000)) +
  xlab("Adult Generalism") +
  ylab("Abundance") +
  ggtitle("Complex Lifecyle") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18),
        legend.position = "none")+
  scale_color_manual(values = c("slateblue", "thistle"))

sim <- ggplot(df_combined_simple, aes(x = Adult.trait, y = Abundance, color = Stage, shape = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 7) +
  scale_y_continuous(limits = c(20000, 60000)) +
  xlab("Adult Generalism") +
  ylab("Abundance") +
  labs(shape = "Juvenile Generalism", color = "Stage") +
  ggtitle("Simple Lifecycle") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))+
  scale_color_manual(values = c("slateblue", "thistle"))+
  guides(shape= guide_legend(override.aes = list(stroke = 1.05)))


layout <- "
AB"

combo.plots <- list(com, sim)

plots <- wrap_plots(combo.plots, design = layout)

plots + plot_annotation(
  title = 'Normal Distribution',
  theme = theme(plot.title = element_text(hjust = 0.35, size = 15, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
)



# Even



df_adult <- data.frame(
  Juvenile.trait = rep(x, each = 6),
  Abundance = as.vector(adult.even.pop.clc),
  Adult.trait = rep(x, times = 6),
  Stage = rep("Adult", times = length(x)*length(x)),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df_juvenile <- data.frame(
  Juvenile.trait = rep(x, each = 6),
  Abundance = as.vector(juvenile.even.pop.clc),
  Adult.trait = rep(x, times = 6),
  Stage = rep("Juvenile", times = length(x)*length(x)),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df_simple_adult <- data.frame(
  Juvenile.trait = x ,
  Abundance = as.vector(adult.even.pop.slc),
  Adult.trait = x,
  Stage = rep("Adult", times = length(x)),
  Cycle = rep("Simple", times = length(x))
)

df_simple_juvenile <- data.frame(
  Juvenile.trait =x ,
  Abundance = as.vector(juvenile.even.pop.slc),
  Adult.trait = x,
  Stage = rep("Juvenile", times = length(x)),
  Cycle = rep("Simple", times = length(x))
)



df_combined_complex <- rbind(df_adult, df_juvenile)
df_combined_simple <- rbind(df_simple_adult, df_simple_juvenile)

com <- ggplot(df_combined_complex, aes(x = Adult.trait, y = Abundance, color = Stage, shape = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 7) +
  scale_y_continuous(limits = c(20000, 60000)) +
  xlab("Adult Generalism") +
  ylab("Abundance") +
  ggtitle("Complex Lifecyle") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18),
        legend.position = "none")+
  scale_color_manual(values = c("slateblue", "thistle"))

sim <- ggplot(df_combined_simple, aes(x = Adult.trait, y = Abundance, color = Stage, shape = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 7) +
  scale_y_continuous(limits = c(20000, 60000)) +
  xlab("Adult Generalism") +
  ylab("Abundance") +
  labs(shape = "Juvenile Generalism", color = "Stage") +
  ggtitle("Simple Lifecycle") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))+
  scale_color_manual(values = c("slateblue", "thistle"))+
  guides(shape= guide_legend(override.aes = list(stroke = 1.05)))


layout <- "
AB"

combo.plots <- list(com, sim)

plots <- wrap_plots(combo.plots, design = layout)

plots + plot_annotation(
  title = 'Even Distribution',
  theme = theme(plot.title = element_text(hjust = 0.35, size = 15, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
)


# Skewed


df_adult <- data.frame(
  Juvenile.trait = rep(x, each = 6),
  Abundance = as.vector(adult.skew.pop.clc),
  Adult.trait = rep(x, times = 6),
  Stage = rep("Adult", times = length(x)*length(x)),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df_juvenile <- data.frame(
  Juvenile.trait = rep(x, each = 6),
  Abundance = as.vector(juvenile.skew.pop.clc),
  Adult.trait = rep(x, times = 6),
  Stage = rep("Juvenile", times = length(x)*length(x)),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df_simple_adult <- data.frame(
  Juvenile.trait = x ,
  Abundance = as.vector(adult.skew.pop.slc),
  Adult.trait = x,
  Stage = rep("Adult", times = length(x)),
  Cycle = rep("Simple", times = length(x))
)

df_simple_juvenile <- data.frame(
  Juvenile.trait =x ,
  Abundance = as.vector(juvenile.skew.pop.slc),
  Adult.trait = x,
  Stage = rep("Juvenile", times = length(x)),
  Cycle = rep("Simple", times = length(x))
)



df_combined_complex <- rbind(df_adult, df_juvenile)
df_combined_simple <- rbind(df_simple_adult, df_simple_juvenile)


com <- ggplot(df_combined_complex, aes(x = Adult.trait, y = Abundance, color = Stage, shape = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 7) +
  scale_y_continuous(limits = c(20000, 60000)) +
  xlab("Adult Generalism") +
  ylab("Abundance") +
  ggtitle("Complex Lifecyle") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18),
        legend.position = "none")+
  scale_color_manual(values = c("slateblue", "thistle"))

sim <- ggplot(df_combined_simple, aes(x = Adult.trait, y = Abundance, color = Stage, shape = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 7) +
  scale_y_continuous(limits = c(20000, 60000)) +
  xlab("Adult Generalism") +
  ylab("Abundance") +
  labs(shape = "Juvenile Generalism", color = "Stage") +
  ggtitle("Simple Lifecycle") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))+
  scale_color_manual(values = c("slateblue", "thistle"))+
  guides(shape= guide_legend(override.aes = list(stroke = 1.05)))


layout <- "
AB"

combo.plots <- list(com, sim)

plots <- wrap_plots(combo.plots, design = layout)

plots + plot_annotation(
  title = 'Skewed Distribution',
  theme = theme(plot.title = element_text(hjust = 0.35, size = 15, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
)




# ----------------------
# Running simulations: 10 runs number of species varied sigmas ----------------

# Even

job::job(ten.run.even = {

sigma <- c(seq(from = 0.05, to = 1.05, by = 0.2))
  
Total.species.SLC.single.even <- c()

Total.species.CLC.even <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
rownames(Total.species.CLC.even) <- sigma  #ADULTS
colnames(Total.species.CLC.even) <- sigma #JUVENILES



# SLC

Total.SLC.list.even <- list()

for(r in 1:10) {
  
  print(paste0("loop ", r, " started"))
  
  Total.species.SLC.single.even <- c()
  
  for(i in 1:length(sigma)){
    

    outputSLC <- resourceCompetitionSLC(resProp=resource.prop.norm.slc, iniP = 0, resFreq=resource.freq.norm.slc, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)

    
    #Filter out similar "species"
    
    final.data.SLC.even <- slc.groups(output = outputSLC)
    Total.species.SLC.single.even[i] <- nrow(final.data.SLC.even)
  }
  
  Total.SLC.list.even[[r]] <- Total.species.SLC.single.even

}


# Caluclating mean and SD of 10 runs



Total.mean.SLC.even <- sapply(1:length(sigma), function(i) mean(sapply(Total.SLC.list.even, function(x) x[i])))

array.data.SLC <- array(unlist(Total.SLC.list.even), dim = c(dim(Total.SLC.list.even[[1]]), length(Total.SLC.list.even)))

Total.sd.SLC.even <- sapply(1:length(sigma), function(i) sd(sapply(Total.SLC.list.even, function(x) x[i])))


# CLC

print("clc start")
Total.CLC.list.even <- list()


for(a in 1:10){
  print(paste0("loop ", a, " started"))
  
  Total.species.CLC.even <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.even) <- sigma #ADULTS
  colnames(Total.species.CLC.even) <- sigma #JUVENILES
  
  for(b in 1:length(sigma)){
    
    for(k in 1:length(sigma)){
      

      outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, resFreq=resFreqMatrix.even.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[b],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)


      
      #Filter out similar "species"
      
      final.data.CLC.even <- clc.groups(output = outputCLC)
      Total.species.CLC.even[b, k] <- nrow(final.data.CLC.even)
      
    }
    
  }
  Total.CLC.list.even[[a]] <- Total.species.CLC.even
}

# Calculating mean of 10 runs

# Combine matrices in the list into a 3D array
array.data.CLC <- array(unlist(Total.CLC.list.even), dim = c(dim(Total.CLC.list.even[[1]]), length(Total.CLC.list.even)))


# Calculate mean and standard deviation along the third dimension (across the list)
Total.mean.CLC.even <- apply(array.data.CLC, c(1, 2), mean)
Total.sd.CLC.even <- apply(array.data.CLC, c(1, 2), sd)



job::export(list(Total.mean.CLC.even, Total.sd.CLC.even, Total.mean.SLC.even, Total.sd.SLC.even))
}, import = c(resPropMatrix.even.clc, resFreqMatrix.even.clc, resourceCompetitionCLC, resource.prop.even.slc, resource.freq.even.slc, resourceCompetitionSLC, clc.groups, slc.groups))


# Normal


job::job(ten.run.norm = {
  
  sigma <- c(seq(from = 0.05, to = 1.05, by = 0.2))
  
  Total.species.SLC.single.norm <- c()
  
  Total.species.CLC.norm <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.norm) <- sigma  #ADULTS
  colnames(Total.species.CLC.norm) <- sigma #JUVENILES
  
  
  
  # SLC
  
  Total.SLC.list.norm <- list()
  
  for(r in 1:10) {
    
    print(paste0("loop ", r, " started"))
    
    Total.species.SLC.single.norm <- c()
    
    for(i in 1:length(sigma)){
      
      outputSLC <- resourceCompetitionSLC(resProp=resource.prop.norm.slc, iniP = 0, resFreq=resource.freq.norm.slc, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
      
      #Filter out similar "species"
      
      final.data.SLC.norm <- slc.groups(output = outputSLC)
      Total.species.SLC.single.norm[i] <- nrow(final.data.SLC.norm)
    }
    
    Total.SLC.list.norm[[r]] <- Total.species.SLC.single.norm
    
  }
  
  
  # Caluclating mean and SD of 10 runs
  
  
  
  Total.mean.SLC.norm <- sapply(1:length(sigma), function(i) mean(sapply(Total.SLC.list.norm, function(x) x[i])))
  
  array.data.SLC <- array(unlist(Total.SLC.list.norm), dim = c(dim(Total.SLC.list.norm[[1]]), length(Total.SLC.list.norm)))
  
  Total.sd.SLC.norm <- sapply(1:length(sigma), function(i) sd(sapply(Total.SLC.list.norm, function(x) x[i])))
  
  
  # CLC
  
  Total.CLC.list.norm <- list()
  
  
  for(r in 1:10){
    print(paste0("loop ", r, " started"))
    
    Total.species.CLC.norm <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
    rownames(Total.species.CLC.norm) <- sigma #ADULTS
    colnames(Total.species.CLC.norm) <- sigma #JUVENILES
    
    for(i in 1:length(sigma)){
      
      for(k in 1:length(sigma)){
        
        outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.norm.clc, resFreq=resFreqMatrix.norm.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
        
        #Filter out similar "species"
        
        final.data.CLC.norm <- clc.groups(output = outputCLC)
        Total.species.CLC.norm[i, k] <- nrow(final.data.CLC.norm)
        
      }
      
    }
    Total.CLC.list.norm[[r]] <- Total.species.CLC.norm
  }
  
  # Calculating mean of 10 runs
  
  # Combine matrices in the list into a 3D array
  array.data.CLC <- array(unlist(Total.CLC.list.norm), dim = c(dim(Total.CLC.list.norm[[1]]), length(Total.CLC.list.norm)))
  
  # Calculate mean and standard deviation along the third dimension (across the list)
  Total.mean.CLC.norm <- apply(array.data.CLC, c(1, 2), mean)
  Total.sd.CLC.norm <- apply(array.data.CLC, c(1, 2), sd)
  
  
  
  job::export(list(Total.mean.CLC.norm, Total.sd.CLC.norm, Total.mean.SLC.norm, Total.sd.SLC.norm))
}, import = c(resPropMatrix.norm.clc, resFreqMatrix.norm.clc, resourceCompetitionCLC, resource.prop.norm.slc, resource.freq.norm.slc, resourceCompetitionSLC, clc.groups, slc.groups))


# Skewed


job::job(ten.run.skew = {
  
  sigma <- c(seq(from = 0.05, to = 1.05, by = 0.2))
  
  Total.species.SLC.single.skew <- c()
  
  Total.species.CLC.skew <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
  rownames(Total.species.CLC.skew) <- sigma  #ADULTS
  colnames(Total.species.CLC.skew) <- sigma #JUVENILES
  
  
  
  # SLC
  
  Total.SLC.list.skew <- list()
  
  for(r in 1:10) {
    
    print(paste0("loop ", r, " started"))
    
    Total.species.SLC.single.skew <- c()
    
    for(i in 1:length(sigma)){
      
      outputSLC <- resourceCompetitionSLC(resProp=resource.prop.skew.slc, iniP = 0, resFreq=resource.freq.skew.slc, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
      
      #Filter out similar "species"
      
      final.data.SLC.skew <- slc.groups(output = outputSLC)
      Total.species.SLC.single.skew[i] <- nrow(final.data.SLC.skew)
    }
    
    Total.SLC.list.skew[[r]] <- Total.species.SLC.single.skew
    
  }
  
  
  # Caluclating mean and SD of 10 runs
  
  
  
  Total.mean.SLC.skew <- sapply(1:length(sigma), function(i) mean(sapply(Total.SLC.list.skew, function(x) x[i])))
  
  array.data.SLC <- array(unlist(Total.SLC.list.skew), dim = c(dim(Total.SLC.list.skew[[1]]), length(Total.SLC.list.skew)))
  
  Total.sd.SLC.skew <- sapply(1:length(sigma), function(i) sd(sapply(Total.SLC.list.skew, function(x) x[i])))
  
  
  # CLC
  
  Total.CLC.list.skew <- list()
  
  
  for(r in 1:10){
    print(paste0("loop ", r, " started"))
    
    Total.species.CLC.skew <- matrix(data = NA, nrow = length(sigma), ncol = length(sigma))
    rownames(Total.species.CLC.skew) <- sigma #ADULTS
    colnames(Total.species.CLC.skew) <- sigma #JUVENILES
    
    for(i in 1:length(sigma)){
      
      for(k in 1:length(sigma)){
        
        outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.skew.clc, resFreq=resFreqMatrix.skew.clc, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[k])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
        
        #Filter out similar "species"
        
        final.data.CLC.skew <- clc.groups(output = outputCLC)
        Total.species.CLC.skew[i, k] <- nrow(final.data.CLC.skew)
        
      }
      
    }
    Total.CLC.list.skew[[r]] <- Total.species.CLC.skew
  }
  
  # Calculating mean of 10 runs
  
  # Combine matrices in the list into a 3D array
  array.data.CLC <- array(unlist(Total.CLC.list.skew), dim = c(dim(Total.CLC.list.skew[[1]]), length(Total.CLC.list.skew)))
  
  # Calculate mean and standard deviation along the third dimension (across the list)
  Total.mean.CLC.skew <- apply(array.data.CLC, c(1, 2), mean)
  Total.sd.CLC.skew <- apply(array.data.CLC, c(1, 2), sd)
  
  
  
  job::export(list(Total.mean.CLC.skew, Total.sd.CLC.skew, Total.mean.SLC.skew, Total.sd.SLC.skew))
}, import = c(resPropMatrix.skew.clc, resFreqMatrix.skew.clc, resourceCompetitionCLC, resource.prop.skew.slc, resource.freq.skew.slc, resourceCompetitionSLC, clc.groups, slc.groups))


# Plotting: 10 runs number of species varied sigmas -----------------------

# Even

Total.mean.CLC.even <- ten.run.even$Total.mean.CLC.even
Total.mean.SLC.even <- ten.run.even$Total.mean.SLC.even

Total.sd.CLC.even <-  ten.run.even$Total.sd.CLC.even
Total.sd.SLC.even <-  ten.run.even$Total.sd.SLC.even

sigma <- c(seq(from = 0.05, to = 1.05, by = 0.2))
x <- as.factor(sigma)



df.CLC <- data.frame(
  Juvenile.trait = rep(x, each = 6),
  Adult.trait = rep(x, times = 6),
  Richness = as.vector(Total.mean.CLC.even),
  sd = as.vector(Total.sd.CLC.even),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df.SLC <- data.frame(
  Juvenile.trait = x,
  Adult.trait = x,
  Richness = as.vector(Total.mean.SLC.even),
  sd = as.vector(Total.sd.SLC.even),
  Cycle = rep("Simple", times = length(x))
)



df.combined <- rbind(df.CLC, df.SLC)


ggplot(df.combined, aes(x = Adult.trait, y = Richness, color = Cycle, shape = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Abundance") +
  labs(shape = "Juvenile Generalism", color = "Life strategy") +
  ggtitle("Even Resource distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))+
  scale_color_manual(values = c("slateblue", "thistle"))


# Normal

Total.mean.CLC.norm <- ten.run.norm$Total.mean.CLC.norm
Total.mean.SLC.norm <- ten.run.norm$Total.mean.SLC.norm

Total.sd.CLC.norm <-  ten.run.norm$Total.sd.CLC.norm
Total.sd.SLC.norm <-  ten.run.norm$Total.sd.SLC.norm

sigma <- c(seq(from = 0.05, to = 1.05, by = 0.2))
x <- as.factor(sigma)




df.CLC <- data.frame(
  Juvenile.trait = rep(x, each = 6),
  Adult.trait = rep(x, times = 6),
  Richness = as.vector(Total.mean.CLC.norm),
  sd = as.vector(Total.sd.CLC.norm),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df.SLC <- data.frame(
  Juvenile.trait = x,
  Adult.trait = x,
  Richness = as.vector(Total.mean.SLC.norm),
  sd = as.vector(Total.sd.SLC.norm),
  Cycle = rep("Simple", times = length(x))
)



df.combined <- rbind(df.CLC, df.SLC)


ggplot(df.combined, aes(x = Adult.trait, y = Richness, color = Cycle, shape = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Abundance") +
  labs(shape = "Juvenile Generalism", color = "Life strategy") +
  ggtitle("Normal Resource Distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))+
  scale_color_manual(values = c("slateblue", "thistle"))


# Skewed





Total.mean.CLC.skew <- ten.run.skew$Total.mean.CLC.skew
Total.mean.SLC.skew <- ten.run.skew$Total.mean.SLC.skew

Total.sd.CLC.skew <-  ten.run.skew$Total.sd.CLC.skew
Total.sd.SLC.skew <-  ten.run.skew$Total.sd.SLC.skew

sigma <- c(seq(from = 0.05, to = 1.05, by = 0.2))
x <- as.factor(sigma)




df.CLC <- data.frame(
  Juvenile.trait = rep(x, each = 6),
  Adult.trait = rep(x, times = 6),
  Richness = as.vector(Total.mean.CLC.skew),
  sd = as.vector(Total.sd.CLC.skew),
  Cycle = rep("Complex", times = length(x)*length(x))
)

df.SLC <- data.frame(
  Juvenile.trait = x,
  Adult.trait = x,
  Richness = as.vector(Total.mean.SLC.skew),
  sd = as.vector(Total.sd.SLC.skew),
  Cycle = rep("Simple", times = length(x))
)



df.combined <- rbind(df.CLC, df.SLC)


ggplot(df.combined, aes(x = Adult.trait, y = Richness, color = Cycle, shape = Juvenile.trait, stroke = 1.05)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin=Richness-sd, ymax=Richness+sd), width=.05) +   #position=position_dodge(.9)
  scale_y_continuous(limits = c(0, 30)) +
  xlab("Adult Generalism") +
  ylab("Abundance") +
  labs(shape = "Juvenile Generalism", color = "Life strategy") +
  ggtitle("Skewed Resource Distribution") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))+
  scale_color_manual(values = c("slateblue", "thistle"))



# -------------------------
# Two resources simulation

job::job(Two.res = {

  last.year.list.2.res <- list()
sigma <- seq(from = 0.25, to = 0.4, by = 0.01)

# Symmetric

for(i in 1:length(sigma)){
  
  
  outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.2res, resFreq=resFreqMatrix.2res, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
  
  phenodataCLC <- NULL
  
  phenodataCLC <- data.frame(
    Year = outputCLC$phenotypes[, 1],
    Adult_Trait = outputCLC$phenotypes[, 3],
    Juvenile_Trait = outputCLC$phenotypes[, 4],
    Num_Individuals = outputCLC$phenotypes[, 2])
  
  last.year.list.2.res[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  
  
}



# Asymmetric 



last.year.list.2.res.as <- list()


for(i in 1:length(sigma)){
  
  
  outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.2res, resFreq=resFreqMatrixAs.2res, iniPA = 0, iniPJ = 0, resGen=matrix(c(sigma[i],sigma[i])), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 20000)
  
  phenodataCLC <- NULL
  
  phenodataCLC <- data.frame(
    Year = outputCLC$phenotypes[, 1],
    Adult_Trait = outputCLC$phenotypes[, 3],
    Juvenile_Trait = outputCLC$phenotypes[, 4],
    Num_Individuals = outputCLC$phenotypes[, 2])
  
  last.year.list.2.res.as[[i]]<- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  
  
}

  job::export(list(last.year.list.2.res, last.year.list.2.res.as))
}, import = c(resPropMatrix.2res, resFreqMatrix.2res, resourceCompetitionCLC, resFreqMatrixAs.2res))



# Results


last.year.list.2.res <- Two.res$last.year.list.2.res
last.year.list.2.res.as <- Two.res$last.year.list.2.res.as


# Two Resources Plotting ----------------------------------------------------

sigmas <- seq(from = 0.25, to = 0.4, by = 0.01)

plot.list.2rs <- list()


for (i in 1:length(last.year.list.2.res)){
  
  color.palette <- mako(length(last.year.list.2.res[[i]]$Adult_Trait))
  
  plot.list.2rs[[i]] <- ggplot(last.year.list.2.res[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) +                                  # Add points
    labs(title = substitute(sigma == value, list(value = sigmas[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-1.5,1.5)) +
    scale_y_continuous(limits = c(-1.5,1.5)) +
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

grid.arrange(grobs = plot.list.2rs, ncol = 4, nrow = 4,
             top = text_grob("Symmetric resources", size = 10, family = "LM Roman 10"))


plot.list.2rs.as <- list()

for (i in 1:length(last.year.list.2.res.as)){
  
  color.palette <- mako(length(last.year.list.2.res.as[[i]]$Adult_Trait))
  
  plot.list.2rs.as[[i]] <- ggplot(last.year.list.2.res.as[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) +                                  # Add points
    labs(title = substitute(sigma == value, list(value = sigmas[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-1.5,1.5)) +
    scale_y_continuous(limits = c(-1.5,1.5)) +
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

grid.arrange(grobs = plot.list.2rs.as, ncol = 4, nrow = 4,
             top = text_grob("Asymmetric resources", size = 10, family = "LM Roman 10"))


