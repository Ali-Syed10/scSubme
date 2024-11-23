## run_find_markers_subsample function is where seurat object, target group which is ident.1 normally in seurat find markers and comparison_group which is ident.2 and n_interations is 
#how many samples you would like to take. 


run_find_markers_subsample <- function(seurat_obj, target_group, comparison_group, n_iterations) {
  results_list <- list() # Create List
  
  # Get the cells in the target group (mEndoA2)
  target_cells <- WhichCells(seurat_obj, idents = target_group) # Take cells from your group of cells under investigation. 
  n_target <- length(target_cells)  # Number of cells in the small group
  
  # Loop for the number of iterations
  for (i in seq_len(n_iterations)) {
    # Get all cells in the comparison group (all other groups)
    comparison_cells <- WhichCells(seurat_obj, idents = comparison_group)
    
    # Randomly sample from the comparison group
    sampled_comparison <- sample(comparison_cells, n_target, replace = FALSE)
    
    # Create a new Seurat object with the sampled cells
    combined_cells <- c(target_cells, sampled_comparison)
    combined_data <- subset(seurat_obj, cells = combined_cells)
    
    # Perform differential expression analysis
    markers <- FindMarkers(combined_data, ident.1 = target_group, ident.2 = comparison_group)
    
    # Store results with iteration number
    markers$iteration <- i
    results_list[[i]] <- markers
  }
  
  return(results_list)
}


### Combin e the results as a list is returned. 

combine_results <- function(results_list) {
  # Combine all results into a single data frame
  combined_df <- do.call(rbind, lapply(results_list, function(x) {
    x$gene <- rownames(x)  # Add gene names as a column
    return(x)
  }))

## Average the results out. 
  
  # Group by gene name and summarize p-values and other metrics
  aggregated_results <- combined_df %>%
    group_by(gene) %>%
    summarize(mean_pval = mean(p_val, na.rm = TRUE), 
              mean_log2FC = mean(avg_log2FC, na.rm = TRUE),
              p_val_adj = mean(p_val_adj, na.rm = TRUE),
              mean_pct_1 = mean(pct.1, na.rm = TRUE),
              mean_pct_2 = mean(pct.2, na.rm = TRUE),
              .groups = 'drop') %>%
    arrange(mean_pval)
  
  return(aggregated_results)
}
