
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

mutational.probability <- round(seq(0.000001, 0.00001 , length.out = 8), digits = 8)

comm.data <- matrix(data = NA, nrow = 3, ncol = 8)

colnames(comm.data) <- mutational.probability
rownames(comm.data) <- c("Single Axis", "Mixed", "Double axes")

comm.data[,1] <- c(5, 1, 0)
comm.data[,2] <- c(0, 3, 3)
comm.data[,3] <- c(0, 1, 5)
comm.data[,4] <- c(0, 1, 5)
comm.data[,5] <- c(0, 1, 5)
comm.data[,6] <- c(0, 2, 4)
comm.data[,7] <- c(0, 1, 5)
comm.data[,8] <- c(0, 0, 6)
comm.data[,9] <- c()
comm.data[,10] <- c()



Community <- data.frame(
  Mutational.Probability = rep(mutational.probability, times = 3),
  Community = rep(c("Single Axis", "Mixed", "Double Axes"), each = 8),
  Times =c(comm.data[1,],comm.data[2,], comm.data[3,] )
)

color_palette <- magma(6)
"#000004FF" "#3B0F70FF" "#8C2981FF" "#DE4968FF" "#FE9F6DFF" "#FCFDBFFF"


ggplot(data = Community, aes(x = Mutational.Probability, y = Times, color = Community)) +
  geom_line(size = 1,) +
  xlab("Mutational Probability") +
  ylab("Community Occurance") +
  labs(color = "Community") +
  scale_color_manual(values = c("#FE9F6DFF", "#DE4968FF", "#8C2981FF"))+
  theme_minimal(base_family = "LM Roman 10", base_size = 15) +
  theme(plot.title = element_text(size = 18))


 