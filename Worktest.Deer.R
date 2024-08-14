
# Worktest

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
library(FamilyRank)
hunt_data <- read.csv("Data/hunt_data.csv", header = T, sep = ",", dec = ".")


# Dovhjort ------------------------

# Data prep and processing -------------------------------------------------------------




# Data prep
Year <- unique(hunt_data$Ar)


deer_per_area <- data.frame(    # Data frame for storing processed data
  Year =  unique(hunt_data$Ar)
)

# Area

for(i in 1:length(Year)){
  deer_per_area$Area[i] <- sum(hunt_data$Jaktmarksareal[hunt_data$Ar==Year[i]])
}  # For loop to extract area each year

# Deer shot

for(i in 1:length(Year)){
  deer_per_area$Deer[i] <- sum(hunt_data$Dovhjort[hunt_data$Ar==Year[i]])
}   # Number of deer shot per year



# length(which(hunt_data$Dovhjort[hunt_data$Ar==Year[1]] != 0))  # Checking if numbers where reasonable

# Plotting -------------------

Deer <- ggplot(deer_per_area, aes(x = Year, y = (Deer/Area))) +
  geom_point() +
  geom_smooth(method = lm, color =  "#E16462FF") +
  scale_x_continuous(breaks = c(2011: 2019)) +
  xlab("År") +
  ylab("Avskjutning dovhjort") +
  ggtitle("Omodifierad Data") +
  theme_bw(base_family = "LM Roman 10", base_size = 13)
 

# Optional Modifier of outliers
deer_per_area$ModDeer <- deer_per_area$Deer

deer_per_area$ModDeer[3] <- (deer_per_area$Deer[2] + deer_per_area$Deer[4]) / 2

deer_per_area$ModDeer[6] <- (deer_per_area$Deer[5] + deer_per_area$Deer[7]) / 2


ModDeer <- ggplot(deer_per_area, aes(x = Year, y = ModDeer)) +
  geom_point() +
  geom_smooth(method = lm, color = "#B12A90FF") +
  scale_x_continuous(breaks = c(2011: 2019)) +
  xlab("År") +
  ylab("Avskjutning dovhjort") +
  ggtitle("Bortagning av extremår")+
  theme_bw(base_family = "LM Roman 10", base_size = 13)
  

wrap_plots(Deer, ModDeer) + plot_layout(axis_titles = "collect")



# Vildsvin ------------------------

# Data prep and processing -------------------------------------------------------------




# Data prep
Year <- unique(hunt_data$Ar)


vild_per_area <- data.frame(    # Data frame for storing processed data
  Year =  unique(hunt_data$Ar)
)

# Area

for(i in 1:length(Year)){
  vild_per_area$Area[i] <- sum(hunt_data$Jaktmarksareal[hunt_data$Ar==Year[i]])
}  # For loop to extract area each year

# Animals shot

for(i in 1:length(Year)){
  vild_per_area$vild[i] <- sum(hunt_data$Vildsvin[hunt_data$Ar==Year[i]])
}   # Number of deer shot per year



# Plotting -------------------

Swine <- ggplot(vild_per_area, aes(x = Year, y = (vild/Area))) +
  geom_point() +
  #geom_smooth(method = lm, color =  "#E16462FF") +
  scale_x_continuous(breaks = c(2011: 2019)) +
  xlab("År") +
  ylab("Avskjutning vildsvin") +
  #ggtitle("Omodifierad data") +
  theme_bw(base_family = "LM Roman 10", base_size = 13)



# Stats

# Check for normality
qqnorm(deer_per_area$Deer/deer_per_area$Area)
qqline(deer_per_area$Deer/deer_per_area$Area)

qqnorm(vild_per_area$vild/vild_per_area$Area)
qqline(vild_per_area$vild/vild_per_area$Area)


hjort_modell <- lm((Deer/Area)~Year, data = deer_per_area)

summary(hjort_modell)
par(mfrow = c(2,2))
plot(hjort_modell)

vild_modell <- lm((vild/Area)~Year, data = vild_per_area)

summary(vild_modell)
par(mfrow = c(2,2))
plot(vild_modell)
