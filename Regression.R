
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

single <- filter(Community, Community == "Single Axis")
mixed <- filter(Community, Community == "Mixed")
double <- filter(Community, Community == "Double Axes")

single.test <- summary(lm(single$Times~single$Mutational.Probability))
mixed.test <- summary(lm(mixed$Times~mixed$Mutational.Probability))
double.test <- summary(lm(double$Times~double$Mutational.Probability))


