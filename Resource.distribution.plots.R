
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
  labs(title = "Even Resource Distribution", x = "Resource characteristic", y = "Resource abundance") +                 # Labels for the axes
  scale_x_continuous(limits = c(-2.7, 2.7))+
  scale_y_continuous(limits = c(0, 7000))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)





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
  labs(title = "Normal Resource Distribution", x = "Resource characteristic", y = "Resource abundance") +                 # Labels for the axes
  scale_x_continuous(limits = c(-2.7, 2.7))+
  scale_y_continuous(limits = c(0, 7000))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)


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
  labs(title = "Skewed Resource Distribution", x = "Resource characteristic", y = "Resource abundance") +                 # Labels for the axes
  scale_x_continuous(limits = c(-2.7, 2.7))+
  scale_y_continuous(limits = c(0, 7000))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)





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
  labs(title = "Bimodal Normal Resource Distribution", x = "Resource characteristic", y = "Resource abundance") +                 # Labels for the axes
  scale_x_continuous(limits = c(-2.7, 2.7))+
  scale_y_continuous(limits = c(0, 7000))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal(base_family = "LM Roman 10", base_size = 10)




res.plots <- (even.res.plot + norm.res.plot + skew.res.plot + binorm.res.plot) #& xlab(NULL) & ylab(NULL) &theme(plot.margin = margin(5.5, 5.5, 0, 0))

res.plots + plot_layout(ncol = 4, nrow = 1) #+
  #labs(tag = c("Resource characteristic")) +
  #theme(
  #  plot.tag = element_text(size = rel(1)),
   # plot.tag.position = c("bottom")
  #)

