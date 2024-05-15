
# Resource distributions
  
  
# Resource initializations -------------------


Num.Res <- 16
res.Abund <-  50000

# Evenly distributed Resources


resource.freq.even.slc <- rep(1/Num.Res, times = Num.Res)                                      # res. freq. 
resource.prop.even.slc <- c(seq(from = -2.5, to = 2.5, length.out = Num.Res))            # res. property 
resource.freq.even.slc <- res.Abund*resource.freq.even.slc

even.res <- data.frame(
  Frequency = resource.freq.even.slc,
  Property = resource.prop.even.slc 
)


even.res.plot <- ggplot(data = even.res, aes(x = Property, y = Frequency)) +
  geom_bar(color = "#C33D80FF" , fill = "#C33D80FF", stat = "identity", width=0.25) + 
  labs(title = "Even Distribution", x = "Resource characteristic", y = "Resource abundance") +                 # Labels for the axes
  scale_x_continuous(limits = c(-2.7, 2.7))+
  scale_y_continuous(limits = c(0, 7000))+
  theme_minimal(base_family = "LM Roman 10", base_size = 10)+
  theme(plot.title = element_text(size = 18)) 





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


# SLC:

resource.prop.norm.slc  <- N.resource.property             # res. property 
resource.freq.norm.slc  <- res.Abund*N.resource.frequency

norm.res <- data.frame(
  Frequency = resource.freq.norm.slc,
  Property = resource.prop.norm.slc 
)

norm.res.plot <- ggplot(data = norm.res, aes(x = Property, y = Frequency)) +
  geom_bar(color = "#D35171FF" , fill = "#D35171FF", stat = "identity", width=0.25) + 
  labs(title = "Normal Distribution", x = "Resource characteristic", y = "Resource abundance") +                 # Labels for the axes
  scale_x_continuous(limits = c(-2.7, 2.7)) +
  scale_y_continuous(limits = c(0, 7000)) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10) +
  theme(plot.title = element_text(size = 18)) 


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

skew.res <- data.frame(
  Frequency = resource.freq.skew.slc,
  Property = resource.prop.skew.slc 
)

skew.res.plot <- ggplot(data = skew.res, aes(x = Property, y = Frequency)) +
  geom_bar(color = "#E16462FF" , fill = "#E16462FF", stat = "identity", width=0.25) + 
  labs(title = "Skewed  Distribution", x = "Resource characteristic", y = "Resource abundance") +                 # Labels for the axes
  scale_x_continuous(limits = c(-2.7, 2.7))+
  scale_y_continuous(limits = c(0, 7000))+
  theme_minimal(base_family = "LM Roman 10", base_size = 10) +
  theme(plot.title = element_text(size = 18))





# Bimodal resources


m1 <- -1.25
m2 <-  1.25
s <- 0.5
Bi.resource.frequency <- c()
Bi.resource.property<- c(seq(from = -2.5, to = 2.5, length.out = Num.Res)) 



mid.add <- c()
midpoint <- c()

for(i in 1:(length(Bi.resource.property))){
  mid.add <- (Bi.resource.property[i+1]-Bi.resource.property[i])/2
  high.midpoint <- Bi.resource.property[(i)]+mid.add
  low.midpoint <- Bi.resource.property[(i)]-mid.add
  if(i == 1){
    Bi.resource.frequency[i] <- pnorm(high.midpoint, mean = m1, sd = s)/2 
  }else if(i == length(Bi.resource.property)){
    low.midpoint <- Bi.resource.property[(i-1)] + (Bi.resource.property[i]-Bi.resource.property[i-1])/2
    Bi.resource.frequency[i] <- pnorm(low.midpoint, mean = m2, sd = s, lower.tail = FALSE)
  }else if (Bi.resource.property[i]<0) {
    Bi.resource.frequency[i] <- (pnorm(high.midpoint, mean = m1, sd = s) - pnorm(low.midpoint, mean = m1, sd = s))/2
  }else{
    Bi.resource.frequency[i] <- (pnorm(high.midpoint, mean = m2, sd = s) - pnorm(low.midpoint, mean = m2, sd = s))/2
  }
}


# SLC:

resource.prop.binorm.slc  <- Bi.resource.property             # res. property 
resource.freq.binorm.slc  <- res.Abund*Bi.resource.frequency



binorm.res <- data.frame(
  Frequency = resource.freq.binorm.slc,
  Property = resource.prop.binorm.slc 
)

binorm.res.plot <- ggplot(data = binorm.res, aes(x = Property, y = Frequency)) +
  geom_bar(color = "#ED7953FF" , fill = "#ED7953FF", stat = "identity", width=0.25) + 
  labs(title = "Bimodal Normal Distribution", x = "Resource characteristic", y = "Resource abundance") +                 # Labels for the axes
  scale_x_continuous(limits = c(-2.7, 2.7))+
  scale_y_continuous(limits = c(0, 7000))+
  theme_minimal(base_family = "LM Roman 10", base_size = 10) +
  theme(plot.title = element_text( size = 18))




res.plots <- (even.res.plot + norm.res.plot + skew.res.plot + binorm.res.plot) 

res.plots + plot_layout(ncol = 4, nrow = 1) 


# Community shape plots

#Three communities: Single Axes, Double axes, and Mixed

mutational.probability <- round(seq(0.000001, 0.00001 , length.out = 10), digits = 8)

comm.data <- matrix(data = NA, nrow = 3, ncol = 10)

colnames(comm.data) <- mutational.probability
rownames(comm.data) <- c("Single Axis", "Mixed", "Double axes")

comm.data[,1] <- c(10, 0, 0)
comm.data[,2] <- c(9, 2, 0)
comm.data[,3] <- c(4, 5, 1)
comm.data[,4] <- c(1, 8, 1)
comm.data[,5] <- c(2, 4, 4)
comm.data[,6] <- c(2, 7, 2)
comm.data[,7] <- c(0, 7, 3)
comm.data[,8] <- c(0, 5, 5)
comm.data[,9] <- c(1, 2, 7)
comm.data[,10] <- c(0, 4, 6)



Community <- data.frame(
  Mutational.Probability = rep(mutational.probability, times = 3),
  Community = rep(c("Single Axis", "Mixed", "Double Axes"), each = 10),
  Times =c(comm.data[1,],comm.data[2,], comm.data[3,] )
)

color_palette <- magma(6)
"#000004FF" "#3B0F70FF" "#8C2981FF" "#DE4968FF" "#FE9F6DFF" "#FCFDBFFF"


qqnorm(Community$Times)
qqline(Community$Times)

summary(lm(Community$Times~Community$Mutational.Probability*Community$Community))



ggplot(data = Community, aes(x = Mutational.Probability, y = Times, color = Community)) +
  geom_point()+
  geom_smooth(method = "lm") +
  xlab("Mutational Probability") +
  ylab("Community Occurance") +
  labs(color = "Community") +
  scale_color_manual(values = c("#FE9F6DFF", "#DE4968FF", "#8C2981FF"))+
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))




# Choose plots for data lets goo

single <- even$Total.endpoint.CLC.even[[5]][[1]]
mix <- even$Total.endpoint.CLC.even[[10]][[3]]
double <- even$Total.endpoint.CLC.even[[7]][[7]]
    
    
color.palette <- magma(length(single$Adult_Trait))
    
single.plot <- ggplot(data = single, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(mu == value, list(value = 0.000001)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
    theme_minimal(base_family = "LM Roman 10", base_size = 10)+
    coord_fixed()

color.palette <- magma(length(mix$Adult_Trait))
   
mix.plot <- ggplot(data = mix, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
  labs(title = substitute(mu == value, list(value = 0.000004)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)+
  coord_fixed()    

color.palette <- magma(length(double$Adult_Trait))

double.plot <- ggplot(data = double, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
  labs(title = substitute(mu == value, list(value = 0.000009)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)+
  coord_fixed()  

  
plots <- wrap_plots(list(single.plot, mix.plot, double.plot))
  
plots + plot_annotation(
  tag_levels = "A",
  tag_prefix = "(",
  tag_suffix = ")")


# Several sigmas and both CLC plus SLC endpoint same plot

# SLC

sigma <- seq(from = 0.05, to = 0.8, length.out = 6)
plot.list.even.SLC <- list()

for(s in 1:length(sigma)){
  adu.sigma <- sigma[s]
  
  this.run <- even$Total.endpoint.SLC.even[[1]][[s]]
  
  data <- this.run
  data$Y <- 0
  
  color.palette <- magma(length(data$Trait))
  
  plot.list.even.SLC[[s]] <- ggplot(data, aes(x =Trait, y = Y)) +
    geom_hline(yintercept = 0, color = "grey")+ 
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    scale_x_continuous(limits = c(-3, 3))+
    ylim(0,0)+
    scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
    theme_minimal(base_family = "LM Roman 10", base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 14),
          panel.grid = element_blank())
  
  
}
plot.list.even.SLC[[1]]

# CLC

sigma <- seq(from = 0.05, to = 0.8, length.out = 6)
plot.list.even <- list()

for(s in 1:length(sigma)){
  adu.sigma <- sigma[s]
  
  juv.sigma <- adu.sigma
  
  this.run <- even$Total.endpoint.CLC.even[[1]]
  this.run <- this.run[this.run$Adult.gen == adu.sigma, ]
  this.run <- this.run[this.run$Juv.gen == juv.sigma, ]

  data <- this.run
    
  color.palette <- magma(length(data$Adult_Trait))
    
  plot.list.even[[s]] <- ggplot(data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
    geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
    labs(title = substitute(sigma == value, list(value = adu.sigma)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
    theme_minimal(base_family = "LM Roman 10", base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_blank(),
          axis.text = element_text(size = 14)) +
    coord_fixed()
    
    
  }

plot.list.even[[1]] + plot.list.even[[2]] + plot.list.even[[3]] + plot.list.even[[4]] + plot.list.even[[5]] + plot.list.even[[6]] +
  plot.list.even.SLC[[1]] + plot.list.even.SLC [[2]] + plot.list.even.SLC [[3]] + plot.list.even.SLC [[4]] + plot.list.even.SLC [[5]] + plot.list.even.SLC [[6]] + 
    plot_layout(ncol = 6, heights = c(3, 1) ) %>%
  annotate_figure(top=text_grob("h"))


# Four distribution phase plane plot
  
even.data <- even$Total.endpoint.CLC.even[[1]]
even.data <- even.data[even.data$Adult.gen == 0.2, ]
even.data <- even.data[even.data$Juv.gen == 0.2, ]

norm.data <- norm$Total.endpoint.CLC.norm[[1]]
norm.data <- norm.data[norm.data$Adult.gen == 0.2, ]
norm.data <- norm.data[norm.data$Juv.gen == 0.2, ]

skew.data <- skew$Total.endpoint.CLC.skew[[1]]
skew.data <- skew.data[skew.data$Adult.gen == 0.2, ]
skew.data <- skew.data[skew.data$Juv.gen == 0.2, ]

binorm.data <- binorm$Total.endpoint.CLC.binorm[[1]]
binorm.data <- binorm.data[binorm.data$Adult.gen == 0.2, ]
binorm.data <- binorm.data[binorm.data$Juv.gen == 0.2, ]

color.palette <- magma(length(even.data$Adult_Trait))

even.plot <- ggplot(data = even.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
  labs(title = substitute(sigma == value, list(value = 0.2)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)+
  coord_fixed()

color.palette <- magma(length(norm.data$Adult_Trait))

norm.plot <- ggplot(data = norm.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
  labs(title = substitute(sigma == value, list(value = 0.2)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)+
  coord_fixed()   

color.palette <- magma(length(skew.data$Adult_Trait))

skew.plot <- ggplot(data = skew.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
  labs(title = substitute(sigma == value, list(value = 0.2)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)+
  coord_fixed() 

color.palette <- magma(length(binorm.data$Adult_Trait))

binorm.plot <- ggplot(data = binorm.data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
  labs(title = substitute(sigma == value, list(value = 0.2)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)+
  coord_fixed()


even.plot + norm.plot + skew.plot + binorm.plot + plot_layout()

plots + plot_annotation(
  tag_levels = "A",
  tag_prefix = "(",
  tag_suffix = ")")

 