# Statistics test

species <- rbind(two.res.sym$Total.species.CLC.even, two.res.assym$Total.species.CLC.even)
rownames <- c("Symmetric", "Asymmteric")

chisq.test(species)



