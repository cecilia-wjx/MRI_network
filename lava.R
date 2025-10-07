library(LAVA)
library(readxl)
library(tidyverse)
library(writexl)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)

# Read in locus info file
loci = read.loci("/path")

process_locus_complete <- function(i, loci, input) {
  tryCatch({
    locus <- LAVA::process.locus(loci[i,], input)
    
    if(!is.null(locus)) {
      result <- LAVA::run.univ.bivar(locus)
      
      if(!is.null(result)) {
        if(!is.null(result$univ)) {
          result$univ$locus_id <- i
          result$univ$chr <- loci$CHR[i]
          result$univ$start <- loci$START[i]
          result$univ$stop <- loci$STOP[i]
        }
        
        if(!is.null(result$bivar)) {
          result$bivar$locus_id <- i
          result$bivar$chr <- loci$CHR[i]
          result$bivar$start <- loci$START[i]
          result$bivar$stop <- loci$STOP[i]
        }
        
        return(result)
      }
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

# Batch specific parameters
START_IDX <- BATCH_START
END_IDX <- BATCH_END
info_files <- paste0("/path/clique_", START_IDX:END_IDX, ".txt")
output_biv <- "/path/LAVA/result_biv"
output_univ <- "/path/LAVA/result_univ"

dir.create(output_biv, showWarnings = FALSE, recursive = TRUE)
dir.create(output_univ, showWarnings = FALSE, recursive = TRUE)

setwd("/path/LAVA/")

for(file_idx_local in 1:length(info_files)) {
  file_idx_global <- START_IDX + file_idx_local - 1
  
  # Check if file exists
  if(!file.exists(info_files[file_idx_local])) {
    cat("File not found, skipping:", info_files[file_idx_local], "\\n")
    next
  }
  
  tryCatch({
    current_info <- read.table(info_files[file_idx_local], header = TRUE, stringsAsFactors = FALSE)
    if(nrow(current_info) == 0) {
      cat("Empty file, skipping:", info_files[file_idx_local], "\\n")
      next
    }
    
    cat("Processing file:", info_files[file_idx_local], "\\n")
    
    input <- process.input(input.info.file = info_files[file_idx_local],
                           sample.overlap.file = NULL,
                           ref.prefix = "lava-ukb-v1.1",
                           phenos = current_info$phenotype)
    
    n_cores <- 10
    cat("Using", n_cores, "cores for parallel computation\\n")
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    all_results <- foreach(i = 1:nrow(loci), .packages = c("LAVA")) %dopar% {
      process_locus_complete(i, loci, input)
    }
    
    stopCluster(cl)
    
    # Filter out NULL results
    all_results <- all_results[!sapply(all_results, is.null)]
    
    if(length(all_results) > 0) {
      # Extract univariate results
      univ_results_list <- lapply(all_results, function(x) x$univ)
      biv_results_list <- lapply(all_results, function(x) x$bivar)
      
      # Filter out NULL results
      univ_results_list <- univ_results_list[!sapply(univ_results_list, is.null)]
      biv_results_list <- biv_results_list[!sapply(biv_results_list, is.null)]
      
      # Save univariate results
      if(length(univ_results_list) > 0) {
        univ_df <- do.call(rbind, univ_results_list)
        univ_df$file_id <- file_idx_global
        univ_df$file_name <- basename(info_files[file_idx_local])
        
        output_file_univ <- file.path(output_univ, paste0("univ_clique_", file_idx_global, ".csv"))
        fwrite(univ_df, output_file_univ)
        cat("Saved univariate results to:", output_file_univ, "\\n")
      }
      
      # Save bivariate results
      if(length(biv_results_list) > 0) {
        biv_df <- do.call(rbind, biv_results_list)
        biv_df$file_id <- file_idx_global
        biv_df$file_name <- basename(info_files[file_idx_local])
        
        output_file_biv <- file.path(output_biv, paste0("biv_clique_", file_idx_global, ".csv"))
        fwrite(biv_df, output_file_biv)
        cat("Saved bivariate results to:", output_file_biv, "\\n")
      }
    } else {
      cat("No valid results for:", info_files[file_idx_local], "\\n")
    }
    
  }, error = function(e) {
    cat("Error processing file:", info_files[file_idx_local], "Error message:", e$message, "\\n")
    # Continue to next iteration without using next
  })
}

cat("Batch", BATCH_NUM, "completed: files", START_IDX, "to", END_IDX, "\\n")