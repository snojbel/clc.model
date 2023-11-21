# Plotting 

# Libraries:

library(gganimate)
library(gifski)
library(png)
library(ggmatplot)
library(tidyverse)
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

last_year_dataSLC <- phenodataSLC[phenodataSLC$Year == max(phenodataSLC$Year), ]
color_palette_SLC <- mako(length(last_year_dataSLC$Trait))

ggplot(last_year_dataSLC, aes(x = Trait, y = 1)) +
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


# ------------------ Varied Sigma plot
x <- rownames(Total_mean_CLC)


SLC <-  data.frame(x = rep(x, length(Total_mean_SLC)), y = Total_mean_SLC)
shapes <- c(rep(x = 8, times =nrow(SLC)))
SLC <- cbind(SLC, shapes)


ggmatplot(x, Total_mean_CLC,
          plot_type = "point",
          xlab = "Adult Generalism",
          ylab = "Number of species",
          legend_title = "Juvenile Generalism",
          legend_label = sigma,
          size = 8) +
          theme_minimal(base_family = "LM Roman 10", base_size = 18)+
          geom_point(data = SLC, aes(x = x, y = y, group=x), size = 8, shape = shapes)


# Varied sigma for both CLC and SLC

x <- rownames(Total_species_CLC)


CLCplot <- ggmatplot(x, Total_species_CLC,
          plot_type = "point",
          xlab = "AdultGeneralism",
          ylab = "Number of species",
          legend_title = "Juvenile Generalism",
          legend_label = sigma,
          size = 8,)+
        theme_minimal(base_family = "LM Roman 10", base_size = 18)+
        ggtitle("Complex life cycle")+
        theme(legend.position = "none")
        

SLCplot <- ggmatplot(x, Total_species_SLC,
                     plot_type = "point",
                     xlab = "Adult Generalism",
                     ylab = "Number of species",
                     legend_title = "Juvenile Generalism",
                     legend_label = sigma,
                     size = 8) +
                theme_minimal(base_family = "LM Roman 10", base_size = 18)+
                ggtitle("Simple life cycle")

grid.arrange(CLCplot,SLCplot, ncol = 2, widths = c(1.3,2))

# ----------------- Animated plot

last_year_data <- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
color_palette <- mako(length(last_year_data$Adult_Trait))


plot <- ggplot(phenodataCLC, aes(x = Juvenile_Trait, y = Adult_Trait, size = Num_Individuals)) +
  geom_point(colour = "#158FAD", show.legend = FALSE) +
  scale_x_log10() +
  theme_minimal(base_family = "LM Roman 10", base_size = 18) +
  # gganimate specific bits:
  labs( x = 'Juvenile Trait', y = 'Adult Trait') +
  transition_manual(Year) +
  ease_aes('linear')


anim_save("C:\\Users\\izer4773\\AppData\\Local\\Temp\\RtmpElr2aL\\293c8497101/animated.plot.gif",animation = plot, renderer = gifski_renderer())
