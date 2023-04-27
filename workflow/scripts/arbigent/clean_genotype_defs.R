add_gts_revisited <- function(data, bias_factor, bias_add_factor, cutoff) {
  # Initialize genotype column with "UNK"
  data$genotype <- "UNK"
  
  # Assign genotype based on confidence levels
  high_confidence <- data$confidence_nobias_over_hard >= bias_factor
  low_confidence <- data$confidence_nobias_over_hard < bias_factor | data$confidence_nobias_over_hard < bias_add_factor
  second_confidence <- data$confidence_hard_over_second < cutoff
  
  data$genotype[high_confidence] <- data$pred_nobias[high_confidence]
  data$genotype[high_confidence & second_confidence] <- paste0(data$pred_nobias[high_confidence & second_confidence])
  data$genotype[low_confidence] <- data$pred_hard[low_confidence]
  data$genotype[low_confidence & second_confidence] <- data$pred_hard[low_confidence & second_confidence]
  
  # Assign "noreads" genotype to rows with zero confidence
  no_reads <- data$confidence_hard_over_second == 0 & data$confidence_nobias_over_hard == 0
  data$genotype[no_reads] <- "noreads"
  
  # Return modified data
  return(data)
}

# Function to add genotypes to a table based on certain confidence thresholds
# Arguments:
# - tab: a data frame containing prediction data
# - bias_factor: a threshold for bias confidence
# - bias_add_factor: a threshold for adding bias confidence
# - cutoff: a threshold for second best confidence
# Returns: a modified data frame with genotypes added based on confidence thresholds
add_gts_revisited_lowconf <- function(tab, bias_factor, bias_add_factor, cutoff) {
  
  # Initialize genotype column as "UNK" (unknown)
  tab$GT <- "UNK"
  
  # Print status update
  cat("Adding genotypes based on confidence thresholds...\n")
  
  # Add genotypes for high-confidence predictions
  high_confidence <- tab$confidence_nobias_over_hard >= bias_factor
  tab$GT[high_confidence] <- tab$pred_nobias[high_confidence]
  
  # Add genotypes for low-confidence predictions
  low_confidence <- high_confidence & tab$confidence_hard_over_second < cutoff
  tab$GT[low_confidence] <- paste0(tab$pred_nobias[low_confidence])
  
  # Add genotypes for bias and bias-add confidence
  bias_confidence <- tab$confidence_nobias_over_hard < bias_factor | tab$confidence_nobias_over_hard < bias_add_factor
  tab$GT[bias_confidence] <- tab$pred_hard[bias_confidence]
  
  # Add genotypes for low-confidence bias and bias-add predictions
  low_bias_confidence <- bias_confidence & tab$confidence_hard_over_second < cutoff
  tab$GT[low_bias_confidence] <- paste0(tab$pred_hard[low_bias_confidence], "_lowconf")
  
  # Add "noreads" genotype for cases where no reads were detected
  tab$GT[tab$confidence_hard_over_second == 0 & tab$confidence_nobias_over_hard == 0] <- "noreads"
  
  # Print status update
  cat("Genotypes added successfully.\n")
  
  # Return modified data frame
  return(tab)
}

# This function takes in a data frame `tab`, a bias factor `bf`, a bias add factor `baf`, and a cutoff value `cutoff`.
# It creates a new column called "GTL" and populates it based on the values in other columns of the data frame.
# Rows are classified into one of four categories based on the values in the "confidence_hard_over_second" and "confidence_nobias_over_hard" columns.
# The "GTL" column is then populated with a string that is dependent on which category the row belongs to.
add_long_gts_revisited <- function(tab, bf, baf, cutoff) {
  
  # Set the "GTL" column of `tab` to "UNK" for all rows.
  tab$GTL <- "UNK"
  
  # Round the "confidence_hard_over_second" and "confidence_nobias_over_hard" columns to one decimal place
  # and store the result in new columns called "cohos" and "conhoh", respectively.
  tab$cohos <- round(tab$confidence_hard_over_second, 1)
  tab$conhoh <- round(tab$confidence_nobias_over_hard, 1)
  
  # Create logical vectors `gt1` and `gt2` to identify rows that meet certain conditions.
  gt1 <- (tab$conhoh >= bf)
  gt2 <- (gt1 & (tab$cohos < cutoff))
  
  # Create logical vectors `lt1` and `lt2` to identify rows that meet certain conditions.
  lt1 <- ((tab$conhoh < bf) | (tab$conhoh < baf))
  lt2 <- (lt1 & (tab$cohos < cutoff))
  
  # Set the "GTL" column of rows that meet the condition `gt1` to a string made up of
  # the "pred_nobias", "pred_hard", "cohos", and "conhoh" columns, separated by colons.
  tab[gt1, "GTL"] <- paste(tab[gt1, "pred_nobias"], "T", tab[gt1, "pred_hard"], tab[gt1, "cohos"] + tab[gt1, "conhoh"], tab[gt1, "cohos"], sep = ":")
  
  # Set the "GTL" column of rows that meet the condition `gt2` to the same string as `gt1`.
  tab[gt2, "GTL"] <- paste(tab[gt2, "pred_nobias"], "T", tab[gt2, "pred_hard"], tab[gt2, "cohos"] + tab[gt2, "conhoh"], tab[gt2, "cohos"], sep = ":")
  
  # Set the "GTL" column of rows that meet the condition `lt1` to a string made up of
  # the "pred_hard", "pred_nobias", "cohos", and "conhoh" columns, separated by colons.
  tab[lt1, "GTL"] <- paste(tab[lt1, "pred_hard"], "T", tab[lt1, "pred_nobias"], tab[lt1, "cohos"], tab[lt1, "cohos"] + tab[lt1, "conhoh"], sep = ":")
  
  # Set the "GTL" column of rows that meet the condition `lt2` to a string made up of
  # the "pred_hard", "F", "pred_nobias", "cohos", and "conhoh" columns, separated by colons.
  tab[lt2, "GTL"] <- paste(tab[lt2, "pred_hard"], "F", tab[lt2, "pred_nobias"], tab[lt2, "cohos"], tab[lt2, "cohos"] + tab[lt2, "conhoh"], sep = ":")

  tab[(tab$cohos == 0) & (tab$conhoh == 0), "GTL"] <- "noreads"
  return(tab)
}

add_simple_gts <- function(data, bias_factor, bias_add_factor, cutoff) {
  # Initialize genotype column with "UNK"
  data$GTs <- "UNK"
  
  # Assign genotypes based on confidence levels
  high_confidence <- data$confidence_nobias_over_hard >= bias_factor
  low_confidence <- data$confidence_nobias_over_hard < bias_factor | data$confidence_nobias_over_hard < bias_add_factor
  second_confidence <- data$confidence_hard_over_second < cutoff
  
  data$GTs[high_confidence & second_confidence] <- "complex"
  data$GTs[low_confidence & second_confidence] <- "simple_lowconf"
  data$GTs[high_confidence | low_confidence] <- "complex"
  data$GTs[!(high_confidence | low_confidence | second_confidence)] <- "simple"
  
  # Assign "noreads" genotype to rows with zero confidence
  no_reads <- data$confidence_hard_over_second == 0 & data$confidence_nobias_over_hard == 0
  data$GTs[no_reads] <- "noreads"
  
  # Return modified data
  return(data)
}