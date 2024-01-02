# Updated simmer with jobs

library(job)
library(ggmatplot)
library(tidyverse)
library(viridisLite)  # Color things
library(viridis)
library(ggplot2)      # Prettier plots
library(gridExtra)  
library(ggpubr)
library(cowplot)
library(extrafont) 

# Resource initializations -------------------

Num.Res <- 25
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


# Running simulations -----------------


# Even distribution -----------------------

# SLC

job::job(output.Even.SLC = {
  output.Even.SLC <- resourceCompetitionSLC(resProp=resource.prop.even.slc, iniP = 0, resFreq=resource.freq.even.slc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Even.SLC)
}, import = c(resource.prop.even.slc, resource.freq.even.slc, resourceCompetitionSLC))

# CLC

job::job(output.Even.CLC = {
  output.Even.CLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix.even.clc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Even.CLC)
}, import = c(resPropMatrix.even.clc, resFreqMatrix.even.clc, resourceCompetitionCLC))



# Normal Distribution --------------------------

# SLC

job::job(output.Norm.SLC = {
  output.Norm.SLC <- resourceCompetitionSLC(resProp=resource.prop.norm.slc, iniP = 0, resFreq=resource.freq.norm.slc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Norm.SLC)
}, import = c(resource.prop.norm.slc, resource.freq.norm.slc, resourceCompetitionSLC))

# CLC

job::job(output.Norm.CLC = {
  output.Norm.CLC <- resourceCompetitionCLC(resProp=resPropMatrix.norm.clc, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix.norm.clc, 
                                           resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Norm.CLC)
}, import = c(resPropMatrix.norm.clc, resFreqMatrix.norm.clc, resourceCompetitionCLC))




# Skewed distribution ------------------------


# SLC

job::job(output.Skew.SLC = {
  output.Skew.SLC <- resourceCompetitionSLC(resProp=resource.prop.skew.slc, iniP = 0, resFreq=resource.freq.skew.slc, 
                                            resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Skew.SLC)
}, import = c(resource.prop.skew.slc, resource.freq.skew.slc, resourceCompetitionSLC))

# CLC

job::job(output.Skew.CLC = {
  output.Skew.CLC <- resourceCompetitionCLC(resProp=resPropMatrix.skew.clc, iniPA = 0, iniPJ = 0, resFreq=resFreqMatrix.skew.clc, 
                                            resGen=matrix(c(0.15,0.15)), popSize = 10, mutProb=0.0005, mutVar=0.05, time.steps = 100000)
  
  # Control what is returned to the main session
  job::export(output.Skew.CLC)
}, import = c(resPropMatrix.skew.clc, resFreqMatrix.skew.clc, resourceCompetitionCLC))


# Results ------------------------------


# Number of species at end of run

species <- matrix(data = NA, nrow = 2, ncol = 3)
rownames(species) <- c("SLC", "CLC")
colnames(species) <- c("Even", "Normal", "Skewed")

# Even

# SLC
output.Even.SLC.result <- output.Even.SLC$output.Even.SLC

stats.even.slc <- output.Even.SLC.result$stats
pheno.even.slc <- output.Even.SLC.result$phenotypes

filt.pheno.even.slc <- slc.groups(output = output.Even.SLC.result)  #Remove too similar "species"
species[1,1] <- as.numeric(nrow(filt.pheno.even.slc))

# CLC

output.Even.CLC.result <- output.Even.CLC$output.Even.CLC

stats.even.clc <- output.Even.CLC.result$stats
pheno.even.clc <- output.Even.CLC.result$phenotypes

filt.pheno.even.clc <- clc.groups(output = output.Even.CLC.result)  #Remove too similar "species"
species[2,1] <- as.numeric(nrow(filt.pheno.even.clc))


# Normal

# SLC
output.Norm.SLC.result <- output.Norm.SLC$output.Norm.SLC

stats.norm.slc <- output.Norm.SLC.result$stats
pheno.norm.slc <- output.Norm.SLC.result$phenotypes

filt.pheno.norm.slc <- slc.groups(output = output.Norm.SLC.result)  #Remove too similar "species"
species[1,2] <- as.numeric(nrow(filt.pheno.norm.slc))

# CLC

output.Norm.CLC.results <- output.Norm.CLC$output.Norm.CLC

stats.norm.clc <- output.Norm.CLC.result$stats
pheno.norm.clc <- output.Norm.CLC.result$phenotypes

filt.pheno.norm.clc <- clc.groups(output = output.Norm.CLC.result)  #Remove too similar "species"
species[2,2] <- as.numeric(nrow(filt.pheno.norm.clc))

# Skewed


# SLC
output.Skew.SLC.result <- output.Skew.SLC$output.Skew.SLC

stats.skew.slc <- output.Skew.SLC.result$stats
pheno.skew.slc <- output.Skew.SLC.result$phenotypes

filt.pheno.skew.slc <- slc.groups(output = output.Skew.SLC.result)  #Remove too similar "species"
species[1,3]<- as.numeric(nrow(filt.pheno.skew.slc))

# CLC

output.Skew.CLC.result <- output.Skew.CLC$output.Skew.CLC

stats.skew.clc <- output.Skew.CLC.result$stats
pheno.skew.clc <- output.Skew.CLC.result$phenotypes

filt.pheno.skew.clc <- clc.groups(output = output.Skew.CLC.result)  #Remove too similar "species"
species[2,3] <- as.numeric(nrow(filt.pheno.skew.clc))

# Plotting


CLC <- species[2,]
SLC <- species[1,]

x <- colnames(species)


SLC <-  data.frame(x = rep(x, length(SLC)), y = SLC)
shapes <- c(rep(x = 8, times = nrow(SLC)))
SLC <- cbind(SLC, shapes)



ggmatplot(x, CLC,
        plot_type = "point",
        xlab = "Distribution Type",
        ylab = "Number of species", size = 8) +
  scale_y_continuous(limits = c(0, 25)) +
  #ggtitle("Number of Species") +
  theme_minimal(base_family = "LM Roman 10", base_size = 15)+
  theme(plot.title = element_text(size = 18)) +                                                  #,panel.grid.major = element_line(colour = "grey", linewidth = 0.3, inherit.blank = FALSE) to add some gridlines
  geom_point(data = SLC, aes(x = x, y = y, group=x), size = 8, shape = shapes)


species.Frame <- data.frame(data.frame(
  SLC = species[1, ],
  CLC = species[2, ],
  dist = c(x)
))


ggplot(data = species.Frame) +
  geom_point(aes(x = dist, y = SLC, shape = "SLC", color = "SLC"), color = "lightblue", fill = "lightblue", size = 8) +
  geom_point(aes(x = dist, y = CLC, shape = "CLC", color = "SLC"), color = "pink", fill = "pink", size = 8) + # Add points
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "Distribution Type", y = "Number of Species", family = "LM Roman 10") +  
  scale_shape_manual(values = c("SLC" = 15, "CLC" = 17), 
                     labels = c("SLC", "CLC")) +  
  theme_classic(base_size = 18)+ 
  theme(axis.text = element_text(family = "LM Roman 10"),
        axis.title = element_text(family = "LM Roman 10", size = 20),
        plot.title = element_text(hjust = 0.5, family = "LM Roman 10"),
        legend.text = element_text(family = "LM Roman 10")) +
  guides(shape = guide_legend(override.aes = list(color = c("lightblue", "pink")), title = NULL))

