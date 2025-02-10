### Supplementary Table RNA-- add distance to closest UGVs

sigs_05_f<-read.csv("temp/comparisons/shared_genes_FIRSTann.csv")
rna_list<-read.csv("data/rnaseq_results_batch_sigs0.1_edited.csv") 
## will eventually need to re-generate this file, but for efficiency I just took it from the previous run

# Function to search for the closest significant gene variants

find_it <- function(data, chrom_value, lowerPOS, upperPOS) {
  # Filter data for the specified CHROM value
  filtered_data <- data[data$CHROM == chrom_value, ]
  
  # Filter for POS values within the specified range
  range <- filtered_data[filtered_data$POS >= lowerPOS & filtered_data$POS <= upperPOS, ]
  
  # If there are any rows within the range, return those
  if (nrow(range) > 0) {
    return(paste("Within Gene:", paste(range$POS, collapse = ", ")))
  } else {
    # If no values fall within the range, find closest smaller and larger POS values
    smaller_pos <- filtered_data[filtered_data$POS < lowerPOS, ]
    larger_pos <- filtered_data[filtered_data$POS > upperPOS, ]
    
    # Find the closest POS value smaller than lowerPOS
    closest_smaller <- if (nrow(smaller_pos) > 0) {
      max(smaller_pos$POS)
    } else {
      NA  # No smaller value found
    }
    
    # Find the closest POS value larger than upperPOS
    closest_larger <- if (nrow(larger_pos) > 0) {
      min(larger_pos$POS)
    } else {
      NA  # No larger value found
    }
    
    # Combine results into a string
    return(paste(
      if (!is.na(closest_smaller)) paste(closest_smaller) else "",
      if (!is.na(closest_larger)) paste(closest_larger) else "",
      sep = " , "
    ))
  }
}

find_sigs <- function(data, range_table) {
  # Add a new column to the range_table with results from find_it
  range_table$result <- apply(range_table, 1, function(row) {
    chrom_value <- row["CHROM"]
    lowerPOS <- as.numeric(row["lowerPOS"])
    upperPOS <- as.numeric(row["upperPOS"])
    
    # Call the find_it function
    find_it(data, chrom_value, lowerPOS, upperPOS)
  })
  
  # Return the updated range_table
  return(range_table)
}


# Run the function
result <- find_sigs(sigs_05_f, rna_list)

# View the result
print(result)
View(result)

write.csv(as.data.frame(result), file = "temp_tables/closest_variants_to_DGEs.csv", row.names = FALSE, quote=FALSE) 


### slightly different version of table that separates it out

find_it <- function(data, chrom_value, lowerPOS, upperPOS) {
  # Filter data for the specified CHROM value
  filtered_data <- data[data$CHROM == chrom_value, ]
  
  # Filter for POS values within the specified range
  range <- filtered_data[filtered_data$POS >= lowerPOS & filtered_data$POS <= upperPOS, ]
  
  # Initialize results
  in_range <- NA
  closest_smaller <- NA
  closest_larger <- NA
  
  # If there are any rows within the range, set `in_range` column
  if (nrow(range) > 0) {
    in_range <- paste(range$POS, collapse = "|")
  }
  
  # If no values fall within the range, find closest smaller and larger POS values
  if (nrow(range) == 0) {
    smaller_pos <- filtered_data[filtered_data$POS < lowerPOS, ]
    larger_pos <- filtered_data[filtered_data$POS > upperPOS, ]
    
    # Find the closest POS value smaller than lowerPOS
    if (nrow(smaller_pos) > 0) {
      closest_smaller <- max(smaller_pos$POS)
    }
    
    # Find the closest POS value larger than upperPOS
    if (nrow(larger_pos) > 0) {
      closest_larger <- min(larger_pos$POS)
    }
  }
  
  # Return results as a list
  return(list(in_range = in_range, closest_smaller = closest_smaller, closest_larger = closest_larger))
}

find_sigs <- function(data, range_table) {
  # Apply the `find_it` function to each row of `range_table`
  results <- apply(range_table, 1, function(row) {
    chrom_value <- row["CHROM"]
    lowerPOS <- as.numeric(row["lowerPOS"])
    upperPOS <- as.numeric(row["upperPOS"])
    
    # Call the `find_it` function
    find_it(data, chrom_value, lowerPOS, upperPOS)
  })
  
  # Convert the results (list of lists) into a data frame
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  
  # Combine the original range_table with the results data frame
  final_table <- cbind(range_table, results_df)
  
  # Return the updated table
  return(final_table)
}

# Run the function
resultb <- find_sigs(sigs_05_f, rna_list)

# View the result

View(resultb)

write.csv(as.data.frame(resultb), file = "temp_tables/closest_variants_to_DGEs_v2.csv", row.names = FALSE, quote=FALSE) 

View(sigs_05_f)

### Now collect the p-values associated with each of the positions identified in the previous function

find_it <- function(data, chrom_value, lowerPOS, upperPOS) {
  # Filter data for the specified CHROM value
  filtered_data <- data[data$CHROM == chrom_value, ]
  
  # Filter for POS values within the specified range
  range <- filtered_data[filtered_data$POS >= lowerPOS & filtered_data$POS <= upperPOS, ]
  
  # Initialize results
  in_range <- NA
  in_range_pval <- NA
  closest_smaller <- NA
  closest_smaller_pval <- NA
  closest_larger <- NA
  closest_larger_pval <- NA
  
  # If there are any rows within the range, set `in_range` and its p-values
  if (nrow(range) > 0) {
    in_range <- paste(range$POS, collapse = "|")
    in_range_pval <- paste(range$pval, collapse = "|")
  }
  
  # If no values fall within the range, find closest smaller and larger POS values
  if (nrow(range) == 0) {
    smaller_pos <- filtered_data[filtered_data$POS < lowerPOS, ]
    larger_pos <- filtered_data[filtered_data$POS > upperPOS, ]
    
    # Find the closest POS value smaller than lowerPOS and its p-value
    if (nrow(smaller_pos) > 0) {
      closest_smaller <- max(smaller_pos$POS)
      closest_smaller_pval <- smaller_pos$pval[which.max(smaller_pos$POS)]
    }
    
    # Find the closest POS value larger than upperPOS and its p-value
    if (nrow(larger_pos) > 0) {
      closest_larger <- min(larger_pos$POS)
      closest_larger_pval <- larger_pos$pval[which.min(larger_pos$POS)]
    }
  }
  
  # Return results as a list
  return(list(
    in_range = in_range,
    in_range_pval = in_range_pval,
    closest_smaller = closest_smaller,
    closest_smaller_pval = closest_smaller_pval,
    closest_larger = closest_larger,
    closest_larger_pval = closest_larger_pval
  ))
}

find_sigs <- function(data, range_table) {
  # Apply the `find_it` function to each row of `range_table`
  results <- apply(range_table, 1, function(row) {
    chrom_value <- row["CHROM"]
    lowerPOS <- as.numeric(row["lowerPOS"])
    upperPOS <- as.numeric(row["upperPOS"])
    
    # Call the `find_it` function
    find_it(data, chrom_value, lowerPOS, upperPOS)
  })
  
  # Convert the results (list of lists) into a data frame
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  
  # Combine the original range_table with the results data frame
  final_table <- cbind(range_table, results_df)
  
  # Return the updated table
  return(final_table)
}

# Run the function
resultc <- find_sigs(sigs_05_f, rna_list)

# View the result
write.csv(as.data.frame(resultc), file = "temp_tables/closest_variants_to_DGEs_v3.csv", row.names = FALSE, quote=FALSE) 
