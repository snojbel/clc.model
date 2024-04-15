



job::job(even = {
  
  
  rep <- 3
  
  Total.species.CLC.single.even <- c()

  
  
  
  # SLC
  
  Total.species.CLC.even <- list()
  Total.endpoint.CLC.even <- list()
  
  for(r in 1:rep) {
    
    
    id <- 1
    
    print(paste0("loop ", r, " started"))
    
    Number.species.CLC.even <- c()
    endpoint.CLC.even <- list()
    
    for(i in 1:length(mutProb)){
      
      
      outputCLC <- resourceCompetitionCLC(resProp=resPropMatrix.even.clc, resFreq=resFreqMatrix.even.clc, iniPA = iniPA, iniPJ = iniPJ, resGen=matrix(c(sigma,sigma)), 
                                          popSize = popSize, mutProb=mutProb[i], mutVar=mutVar, time.steps = time.steps, im = im, fmax = fmax, kA = kA, nmorphs = nmorphs,
                                          threshold = threshold)
      
      
      #Filter out similar "species" and collect number of species data
      
      final.data.CLC.even <- clc.groups(output = outputCLC)
      Number.species.CLC.even[i] <- nrow(final.data.CLC.even)
      
      #Collect endpoint data
      final.data.CLC.even$ID <- c(rep(id, times = nrow(final.data.CLC.even)))
      
      
      endpoint.CLC.even[[i]] <- final.data.CLC.even 
      id <- id + 1 
    }
    
    Total.species.CLC.even[[r]] <- Number.species.CLC.even
    Total.endpoint.CLC.even[[r]] <- endpoint.CLC.even
  }
  
  
  
  
  
  job::export(list(Total.endpoint.CLC.even))
}, import = c(resPropMatrix.even.clc, resFreqMatrix.even.clc, resourceCompetitionCLC, resource.prop.even.slc, resource.freq.even.slc, resourceCompetitionSLC, 
              clc.groups, slc.groups, sigma, popSize, im, fmax, kA, kJ, mutProb, mutVar, time.steps, iniP, iniPA, iniPJ, nmorphs, threshold, maxTr, minTr))

even$Total.endpoint.CLC.even[[1]]



Res <- list()

pdf("plots.even.combined.pdf")

for(s in 1:length(mutProb)){
  adu.mutProb <- mutProb[s]
  
  last.year.list.even <- data.frame()
  
  for(i in 1:length(even$Total.endpoint.CLC.even)){
    this.run <- even$Total.endpoint.CLC.even[[i]][[s]]
    this.run$run <- i
    last.year.list.even <- rbind(last.year.list.even, this.run)
  }
  
  plot.list.even <- list()
  
  for (i in 1:3){
    
    data <- last.year.list.even[last.year.list.even$run == i, ]
    
    color.palette <- plasma(length(data$Adult_Trait))
    
    plot.list.even[[i]] <- ggplot(data, aes(x = Juvenile_Trait, y = Adult_Trait)) +
      geom_point(aes(size=Num_Individuals), color = color.palette, show.legend = FALSE) + 
      labs(title = substitute(sigma == value, list(value = adu.mutProb)), x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
      scale_x_continuous(limits = c(-3, 3))+
      scale_y_continuous(limits = c(-3, 3))+
      scale_size_continuous(limits=c(1,40000),breaks=c(seq(from = 0, to = 40000, by = 5000))) +
      theme(aspect.ratio=1) +
      theme_minimal(base_family = "LM Roman 10", base_size = 10)
    
    
  }
  
  
  plots <- wrap_plots(plot.list.even)
  
  
  
  Res[[s]] <- plots + plot_annotation(
    title = 'Even Distribution',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  )+ coord_fixed()
  
  print(Res[[s]])
  
  
} 

dev.off()




Res <- list()

pdf("plots.even.combined.mutprob.pdf")

for(i in 1:length(even$Total.endpoint.CLC.even)){
  
  
  placeholder <- even$Total.endpoint.CLC.even[[i]][[1]]
  placeholder$run <- i
  
  
  
  
  
  Res[[s]] <- plots + plot_annotation(
    title = 'Even Distribution',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"), plot.subtitle = element_text(hjust = 0.5, size = 15, family = "LM Roman 10"))
  )+ coord_fixed()
  
  print(Res[[s]])
  
  
  
}



dev.off()
