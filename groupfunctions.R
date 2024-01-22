# Grouping functions


# SLC -------------------------

slc.groups <- function(output = outputSLC, threshold = 0.2){
  outputSLC <- output
  
  phenodataSLC <- data.frame(
  Year = outputSLC$phenotypes[, 1],
  Trait = outputSLC$phenotypes[, 3],
  Num_Individuals = outputSLC$phenotypes[, 2]
      )

  last_year_dataSLC <- phenodataSLC[phenodataSLC$Year == max(phenodataSLC$Year), ]

  
  last_year_dataS <- subset(last_year_dataSLC, select = -Year)
  last_year_dataS <- subset(last_year_dataS, select = -Num_Individuals)
  rownames(last_year_dataS) <- NULL
  rownames(last_year_dataSLC) <- NULL
  
  
  distance_matrix <- as.matrix(dist(last_year_dataS[, 1, drop = FALSE], method = "euclidean"))
  
  
  distance_matrix[lower.tri(distance_matrix)] <- NA
  
  
  # Find indices of individuals to keep
  
<<<<<<< HEAD
  if(sum(which(distance_matrix  < threshold, arr.ind = T)) == 0){
=======
  if(sum(which(distance_matrix < threshold, arr.ind = T)) == 0){
>>>>>>> 6df208ae489e2818008299a76854cffae3373ce0
    return(last_year_dataSLC)
  } # Checks if there are zero individuals who are alike.
  
  same <- which(distance_matrix < threshold, arr.ind = T)
  same <- same[same[, 1]-same[,2] != 0, , drop = FALSE]
  rownames(same) <- NULL
  
  
  # Initialize an empty list to store groups
  groups <- list()
  
  # Function to find group index for a species
  find_group <- function(species_id) {
    for (g in seq_along(groups)) {
      if (species_id %in% unlist(groups[[g]])) {
        return(g)
      }
    }
    return(0)
  }
  
  # Iterate over rows in the matrix
  for (s in 1:nrow(same)) {
    species1 <- same[s, 1]
    species2 <- same[s, 2]
    
    # Find groups for each species
    group1 <- find_group(species1)
    group2 <- find_group(species2)
    
    if (group1 == 0 & group2 == 0) {
      # Create a new group
      groups <- c(groups, list(c(species1, species2)))
    } else if (group1 == 0) {
      # Add species1 to the group containing species2
      groups[[group2]] <- c(groups[[group2]], species1)
    } else if (group2 == 0) {
      # Add species2 to the group containing species1
      groups[[group1]] <- c(groups[[group1]], species2)
    } else if (group1 != group2) {
      # Merge two groups
      groups[[group1]] <- c(groups[[group1]], groups[[group2]])
      groups <- groups[-group2]
    }
  }
  
  # Filter out duplicate species in each group
  groups <- lapply(groups, function(group) unique(group))
  
  rownames(last_year_dataSLC) <- NULL
  final_data <- last_year_dataSLC         # Place to store filtered data
  total.sub <- c()                     # Place to store subspecies
  
  #Add population count of "subspecies" to main species
  
  for(q in seq_along(groups)){
    combo <- NULL
    combo <- groups[[q]]
    main <- combo[which.max(final_data[combo,3])]
    sub <- combo[-which.max(final_data[combo,3])]
    final_data[main,3] <- final_data[main,3] + sum(final_data[sub,3])
    total.sub <- rbind(c(total.sub, sub))
    
  }
  # Remove subspecies
  final_data <- final_data[-total.sub, ,drop = FALSE]
  return(final_data)

  
}

# CLC --------------------------------------------------------


clc.groups <- function(output = outputCLC, threshold = 0.2){
  
  outputCLC <- output
  phenodataCLC <- data.frame(
    Year = outputCLC$phenotypes[, 1],
    Adult_Trait = outputCLC$phenotypes[, 3],
    Juvenile_Trait = outputCLC$phenotypes[, 4],
    Num_Individuals = outputCLC$phenotypes[, 2]
  )
  
  last_year_dataCLC <- phenodataCLC[phenodataCLC$Year == max(phenodataCLC$Year), ]
  last_year_dataC <- subset(last_year_dataCLC, select = -Year)
  last_year_dataC <- subset(last_year_dataC, select = -Num_Individuals)
  rownames(last_year_dataCLC) <- NULL
  rownames(last_year_dataC) <- NULL
  
  
  distance_matrix_adult <- as.matrix(dist(last_year_dataC[, 1, drop = FALSE], method = "euclidean"))
  distance_matrix_juvenile <- as.matrix(dist(last_year_dataC[, 2, drop = FALSE], method = "euclidean"))
  
  distance_matrix_adult[lower.tri(distance_matrix_adult)] <- NA
  distance_matrix_juvenile[lower.tri(distance_matrix_juvenile)] <- NA
  
  
  # Set a threshold for similarity (adjust as needed)
  threshold <- 0.2
  
  # Find indices of individuals to keep
  
  if(sum(which(distance_matrix_adult < threshold & distance_matrix_juvenile < threshold, arr.ind = T)) == 0){
    return(last_year_dataCLC)
  }   # Checks if there are zero individuals who are alike.
  
  same <- which(distance_matrix_adult < threshold & distance_matrix_juvenile < threshold, arr.ind = T)
  same <- same[same[, 1]-same[,2] != 0, , drop = FALSE]
  rownames(same) <- NULL
  
  
  # Initialize an empty list to store groups
  groups <- list()
  
  # Function to find group index for a species
  find_group <- function(species_id) {
    for (g in seq_along(groups)) {
      if (species_id %in% unlist(groups[[g]])) {
        return(g)
      }
    }
    return(0)
  }
  
  # Iterate over rows in the matrix
  for (s in 1:nrow(same)) {
    species1 <- same[s, 1]
    species2 <- same[s, 2]
    
    # Find groups for each species
    group1 <- find_group(species1)
    group2 <- find_group(species2)
    
    if (group1 == 0 & group2 == 0) {
      # Create a new group
      groups <- c(groups, list(c(species1, species2)))
    } else if (group1 == 0) {
      # Add species1 to the group containing species2
      groups[[group2]] <- c(groups[[group2]], species1)
    } else if (group2 == 0) {
      # Add species2 to the group containing species1
      groups[[group1]] <- c(groups[[group1]], species2)
    } else if (group1 != group2) {
      # Merge two groups
      groups[[group1]] <- c(groups[[group1]], groups[[group2]])
      groups <- groups[-group2]
    }
  }
  
  # Filter out duplicate species in each group
  groups <- lapply(groups, function(group) unique(group))
  
  final_data <- last_year_dataCLC         # Place to store filtered data
  total.sub <- c()                     # Place to store subspecies
  
  #Add population count of "subspecies" to main species
  
  for(q in seq_along(groups)){
    combo <- NULL
    combo <- groups[[q]]
    main <- combo[which.max(final_data[combo,4])]
    sub <- combo[-which.max(final_data[combo,4])]
    final_data[main,4] <- final_data[main,4] + sum(final_data[sub,4])
    total.sub <- rbind(c(total.sub, sub))
    
  }
  # Remove subspecies
  final_data <- final_data[-total.sub, , drop = FALSE]
  
  return(final_data)
}


final_data_CLC <- clc.groups()

color_palette <- mako(length(final_data_CLC$Adult_Trait))

ggplot(final_data_CLC, aes(x = Juvenile_Trait, y = Adult_Trait)) +
  geom_point(aes(size=Num_Individuals), color = color_palette) +                                  # Add points
  labs(x = "Juvenile Trait", y = "Adult Trait", size = "Number of individuals") +                 # Labels for the axes
  theme_minimal(base_family = "LM Roman 10", base_size = 18)

nrow(final_data_CLC)
