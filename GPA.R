library(tidyverse)
library("GPA")
library(readxl)
library(data.table)
library(parallel)
library(foreach)
library(doParallel)
library(writexl)

# Set up parallel processing
num_cores <- 10
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Read data
cmd_pair <- read_excel("/path/cmd_heart_brain.xlsx")

# Initialize result matrix
result <- foreach(i = 1:nrow(cmd_pair), .combine = rbind) %dopar% {
  library(tidyverse)
  library(GPA)
  library(data.table)
  
  result_row <- rep(NA, 12)
  result_row[1:2] <- unlist(cmd_pair[i, 1:2])
  
  # Read and process data
  brain <- fread(paste0("/path/extra_brain/", cmd_pair[i,4]))
  heart <- fread(paste0("/path/extra_heart/", cmd_pair[i,3]))
  
  brain <- brain[,c(1,2,4,5,9)]
  heart <- heart[,c(1,2,4,5,10)]
  colnames(heart) <- c("chr","rsid","A1","A2","p")
  
  # Merge and filter data
  data <- heart %>%
    inner_join(brain, by = c("chr", "rsid")) %>%
    filter(
      (A1.x == A1.y & A2.x == A2.y) |
        (A1.x == A2.y & A2.x == A1.y)
    )
  
  # GPA analysis
  fit.GPA.noAnn <- GPA(data[,c(5,8)], NULL)
  fit.GPA.pleiotropy.H0 <- GPA(data[,c(5,8)], NULL, pleiotropyH0=TRUE)
  test.GPA.pleiotropy <- pTest(fit.GPA.noAnn, fit.GPA.pleiotropy.H0)
  
  # Store results
  result_row[3:12] <- c(
    unname(test.GPA.pleiotropy$pi),
    unname(test.GPA.pleiotropy$piSE),
    unname(test.GPA.pleiotropy$statistics),
    unname(test.GPA.pleiotropy$pvalue)
  )
  
  result_row
}

# Format results
result <- as.data.frame(result)
colnames(result) <- c("heart","brain","p00","p10","p01","p11",
                      "se00","se10","se01","se11","stat","p")

# Save results
setwd("/path/result")
write_xlsx(result, "res_heart_brain.xlsx")

# Close cluster
stopCluster(cl)
