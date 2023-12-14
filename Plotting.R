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
library(ggpubr)
library(cowplot)
library(extrafont) 
#needed to add extra fonts
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
  labs(title = substitute(sigma == value, list(value = sigmas[1])), x = "Juvenile Trait", y = "Adult Trait", family = "LM Roman 10", size = "Number of individuals") +                 # Labels for the axes
  theme_classic(base_size = 18)+ 
   theme(axis.text = element_text(family = "LM Roman 10"),
         axis.title = element_text(family = "LM Roman 10", size = 20),
         axis.title.y = element_blank(),
         plot.title = element_text(hjust = 0.5, family = "LM Roman 10"))



# ------------------ Varied Sigma plot

Total_mean_CLC <- Com.vs.Sim.nr.species.sim_results$Total.mean.CLC.skewed
Total_mean_SLC <- Com.vs.Sim.nr.species.sim_results$Total.mean.SLC.skewed



x <- rownames(Total_mean_CLC)


SLC <-  data.frame(x = rep(x, length(Total_mean_SLC)), y = Total_mean_SLC)
shapes <- c(rep(x = 8, times =nrow(SLC)))
SLC <- cbind(SLC, shapes)

# Number of species

ggmatplot(x, Total_mean_CLC,
          plot_type = "point",
          xlab = "Adult Generalism",
          ylab = "Number of species",
          legend_title = "Juvenile Generalism",
          legend_label = x, size = 8) +
          scale_y_continuous(limits = c(0, 15)) +
          ggtitle("Skewed distribution") +
          theme_minimal(base_family = "LM Roman 10", base_size = 15)+
          theme(plot.title = element_text(size = 18)) +                                                  #,panel.grid.major = element_line(colour = "grey", linewidth = 0.3, inherit.blank = FALSE) to add some gridlines
          geom_point(data = SLC, aes(x = x, y = y, group=x), size = 8, shape = shapes)

# Abundance

SLC.abund <-  data.frame(x = rep(x, length(Total_mean_SLC)), y = Total_mean_abund_SLC)
shapes <- c(rep(x = 8, times =nrow(SLC.abund)))
SLC.abund <- cbind(SLC.abund, shapes)

test <- as.matrix(Total_abund_CLC_list)

ggmatplot(x, test,
          plot_type = "point",
          xlab = "Adult Generalism",
          ylab = "Abundance of total population",
          legend_title = "Juvenile Generalism",
          legend_label = sigma,
          size = 8) +
  theme_minimal(base_family = "LM Roman 10", base_size = 18)+
  geom_point(data = SLC.abund, aes(x = x, y = y, group=x), size = 8, shape = shapes)



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

# Animated plot  -------------------------------

last_year_data <- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
color_palette <- mako(length(last_year_data$Adult_Trait))
#colour = "#158FAD", show.legend = FALSE

phenodataCLC[,1] <- as.factor(phenodataCLC[,1])


plot <- ggplot(phenodataCLC, aes(x = Juvenile_Trait, y = Adult_Trait, size = Num_Individuals)) +
  geom_point(colour = "#158FAD", show.legend = FALSE) +
  theme_minimal(base_family = "LM Roman 10", base_size = 18) +
  # gganimate specific bits:
  labs( x = 'Juvenile Trait', y = 'Adult Trait') +
  transition_time(Year) +
  ease_aes('linear')


anim_save("D:\\Izabel Master thesis\\Code\\clc.model\\plots", animation = plot, renderer = gifski_renderer())


# Many scatter plots -----------------------------------------------

plot_list_10 <- list()

for (i in 1:9){
  
  color_palette <- mako(length(last_year_list[[i]]$Adult_Trait))
  
  plot_list_10[[i]] <- ggplot(last_year_list[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color_palette, show.legend = FALSE) +                                  # Add points
    labs(x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

grid.arrange(grobs = plot_list_10, ncol = 3, nrow = 3, top =textGrob("Different Starting traits", gp = gpar(fontsize = 10, fontfamily = "LM Roman 10")))



# 2 Res Simulation ------------------------------------------------------------

sigmas <- seq(from = 0.1, to = 0.8, by = 0.05)

#1BA
plot_list_2rs <- list()

last_year_list_2_res <- `2res-sim_results`[["last_year_list_2_res"]]

for (i in 1:length(last_year_list_2_res)){
  
  color_palette <- mako(length(last_year_list_2_res[[i]]$Adult_Trait))
  
  plot_list_2rs[[i]] <- ggplot(last_year_list_2_res[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color_palette, show.legend = FALSE) +                                  # Add points
    labs(title = substitute(sigma == value, list(value = sigmas[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-1.5,1.5)) +
    scale_y_continuous(limits = c(-1.5,1.5)) +
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

grid.arrange(grobs = plot_list_2rs, ncol = 5, nrow = 3,
             top = text_grob("Symmetric resources", size = 10, family = "LM Roman 10"))

#1B2

plot_list_2rs_as <- list()
last_year_list_2_res_as <- `2res-sim_results`[["last_year_list_2_res_as"]]

for (i in 1:length(last_year_list_2_res_as)){
  
  color_palette <- mako(length(last_year_list_2_res_as[[i]]$Adult_Trait))
  
  plot_list_2rs_as[[i]] <- ggplot(last_year_list_2_res_as[[i]], aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color_palette, show.legend = FALSE) +                                  # Add points
    labs(title = substitute(sigma == value, list(value = sigmas[i])), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-1.5,1.5)) +
    scale_y_continuous(limits = c(-1.5,1.5)) +
    theme_minimal(base_family = "LM Roman 10", base_size = 10)
}

grid.arrange(grobs = plot_list_2rs_as, ncol = 5, nrow = 3,
             top = text_grob("Asymmetric resources", size = 10, family = "LM Roman 10"))
