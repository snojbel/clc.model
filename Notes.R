


outputCLC <- output.Even.CLC.0$output.Even.CLC


phenodataCLC <- data.frame(
  Year = outputCLC$phenotypes[, 1],
  Adult_Trait = outputCLC$phenotypes[, 3],
  Juvenile_Trait = outputCLC$phenotypes[, 4],
  Num_Individuals = outputCLC$phenotypes[, 2]
)

transparency <- phenodataCLC$Num_Individuals / max(phenodataCLC$Num_Individuals)

# ------------- Trait divergence


evoAdu <- ggplot(phenodataCLC, aes(x=Year, y=Adult_Trait)) + 
  geom_point(size = 2.5, alpha = transparency, color = "#3B0F70FF")  +
  xlab("Year") + ylab("Adult Trait") +
  theme_minimal(base_family = "LM Roman 10", base_size = 18)



evoJuv <- ggplot(phenodataCLC, aes(x=Year, y=Juvenile_Trait)) + 
  geom_point(size = 2.5, alpha = transparency, color = "#FE9F6DFF") +
  xlab("Year") + ylab("Juvenile Trait") +
  theme_minimal(base_family = "LM Roman 10", base_size = 18) 


grid.arrange(evoAdu,evoJuv, nrow = 2, widths = c(1))




last_year_data <- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
color_palette <- mako(length(last_year_data$Adult_Trait))
#colour = "#158FAD", show.legend = FALSE

phenodataCLC[,1] <- as.factor(phenodataCLC[,1])


plot <- ggplot(phenodataCLC, aes(x = Juvenile_Trait, y = Adult_Trait, size = Num_Individuals)) +
  geom_point(colour = "#FE9F6DFF", show.legend = FALSE) +
  theme_minimal(base_family = "LM Roman 10", base_size = 18) +
  # gganimate specific bits:
  labs( x = 'Juvenile Trait', y = 'Adult Trait') +
  transition_time(Year) +
  ease_aes('linear')


anim_save("D:\\Izabel Master thesis\\Code\\clc.model\\plots", animation = plot, renderer = gifski_renderer())





# Abundance adult vs juvenile vs number of traits


statdataCLC <- data.frame(
  Year = outputCLC$stats[, 1],
  Adult = outputCLC$stats[, 2],
  Juvenile = outputCLC$stats[, 3],
  Morphs = outputCLC$stats[, 4]
)



Abund <- ggplot(statdataCLC, aes(x=Year)) + 
  geom_line(aes(y = Adult), color = "#721F81FF") +
  xlab("Time") + ylab("Abundance") +
  theme_minimal(base_family = "LM Roman 10", base_size = 18) 

Species <- ggplot(statdataCLC, aes(x=Year)) + 
  geom_line(aes(y = Morphs), color = "#FE9F6DFF") +
  xlab("Time") + ylab("Number of Species") +
  theme_minimal(base_family = "LM Roman 10", base_size = 18)

all.plots <- (Abund / Species)

#Checking number of species SLC

outputSLC <- output.Even.SLC.IMI$output.Even.SLC

phenodataSLC <- data.frame(
  Year = outputSLC$phenotypes[, 1],
  Trait = outputSLC$phenotypes[, 3],
  Num_Individuals = outputSLC$phenotypes[, 2]
)

last.year.data.SLC <- phenodataSLC[phenodataSLC$Year == max(phenodataSLC$Year), ]

nrow(last.year.data.SLC)


# Regular plot

#Set sigma 
sigma = 0.05

#plasma(5)
#"#0D0887FF" "#7E03A8FF" "#CC4678FF" "#F89441FF" "#F0F921FF"

outputCLC <- output.Even.CLC.mutprob0.005.200$output.Even.CLC



phenodataCLC <- data.frame(
  Year = outputCLC$phenotypes[, 1],
  Adult_Trait = outputCLC$phenotypes[, 3],
  Juvenile_Trait = outputCLC$phenotypes[, 4],
  Num_Individuals = outputCLC$phenotypes[, 2]
)

last.year.data.CLC <- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
nrow(last.year.data.CLC)

ggplot(last.year.data.CLC, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = "#CC4678FF", show.legend = FALSE) + 
  labs(title = substitute(sigma == value, list(value = 0.2)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
  coord_fixed() +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)


# Over time plot


outputCLC <- output.Even.CLC.IMI$output.Even.CLC

phenodataCLC <- data.frame(
  Year = outputCLC$phenotypes[, 1],
  Adult_Trait = outputCLC$phenotypes[, 3],
  Juvenile_Trait = outputCLC$phenotypes[, 4],
  Num_Individuals = outputCLC$phenotypes[, 2]
)



ggplot(phenodataCLC, aes(x = Adult_Trait, y = Juvenile_Trait)) +
  geom_point(aes(size=Num_Individuals, color = Year)) +                                  # Add points
  geom_point(data = ~filter(.x, Year == max(phenodataCLC$Year)), color = "#7700b3", shape = 4, size = 3) + 
  labs(x = "Adult Trait", y = " Juvenile Trait ", size = "Number of individuals", color = " Time") +# Labels for the axes
  scale_color_viridis(option = "A") +
  guides(size = "none") +                  
  theme_minimal(base_family = "LM Roman 10", base_size = 18)+
  scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
  coord_fixed()
  


