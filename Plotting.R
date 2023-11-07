# Plotting 

# Libraries:
library(viridisLite)  # Color things
library(viridis)
library(ggplot2)      # Prettier plots
library(gridExtra)   
library(extrafont)   #needed to add extra fonts
#font_import()  #Only needed first time in R
#loadfonts()
#fonts() #to check names of fonts

# SLC:

phenodataSLC <- data.frame(
  Year = outputSLC$phenotypes[, 1],
  Trait = outputSLC$phenotypes[, 3],
  Num_Individuals = outputSLC$phenotypes[, 2]
)

transparency <- phenodataSLC$Num_Individuals / max(phenodataSLC$Num_Individuals)

# ------------- Trait divergence


evoSLC <- ggplot(phenodataSLC, aes(x=Year, y=Trait)) + 
  geom_point(size = 2.5, alpha = transparency, color = rgb(0.13, 0.57, 0.55)) +
  xlab("Year") + ylab("Adult Trait") +
  theme_minimal(base_family = "LM Roman 10", base_size = 18)


evoSLC



# -----------------Scatter plot 

last_year_data_SLC <- phenodataSLC[phenodataSLC$Year == max(phenodataSLC$Year), ]
color_palette_SLC <- mako(length(last_year_data_SLC$Trait))

ggplot(last_year_data_SLC, aes(x = Trait, y = 1)) +
  geom_point(aes(size=Num_Individuals), color = color_palette_SLC) +                                  # Add points
  labs(x = "Trait", y = " ", size = "Number of individuals") +                 # Labels for the axes
  theme_minimal(base_family = "LM Roman 10", base_size = 18) +
  theme(axis.text.y=element_blank())




# CLC:

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




