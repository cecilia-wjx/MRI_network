library(tidyverse)
library(data.table)
library(dplyr)
library(mice)
library(survival)
library(bigreadr)

#read into files
setwd("/path")

# Logistic Regression Functions
# Function for brain imaging traits (includes scanner position covariates)
log_reg_brain <- function(phenotype, inputdata){
  logres <- as.data.frame(matrix(NA, ncol=6, nrow=length(phenotype)))
  
  for(var in phenotype){
    # Model 1: Basic adjustment 
    formula1 <- formula(paste("status ~", var, "+age + sex + ethnic + BMI + site + scan_x_de + scan_y_de + scan_z_de + scan_table_de"))
    log_model1 <- glm(formula1, family = binomial(link = "logit"), data = inputdata)
    
    # Model 2: + Socioeconomic factors
    formula2 <- formula(paste("status ~", var, "+  age + sex + ethnic + BMI + TDI + edu + site + scan_x_de + scan_y_de + scan_z_de + scan_table_de"))
    log_model2 <- glm(formula2, family = binomial(link = "logit"), data = inputdata)
    
    # Model 3: + Lifestyle factors
    formula3 <- formula(paste("status ~", var, "+  age + sex + ethnic + BMI + TDI + edu + smoke + alcohol + METmin + site + scan_x_de + scan_y_de + scan_z_de + scan_table_de + group_time"))
    log_model3 <- glm(formula3, family = binomial(link = "logit"), data = inputdata)
    
    # Calculate OR and 95% CI for each model
    OR_results1 <- as.data.frame(exp(cbind("OR"=summary(log_model1)$coefficients[2,1],
                                           "LL"=summary(log_model1)$coefficients[2,1]-1.96*summary(log_model1)$coefficients[2,2],
                                           "UL"=summary(log_model1)$coefficients[2,1]+1.96*summary(log_model1)$coefficients[2,2])))
    
    OR_results2 <- as.data.frame(exp(cbind("OR"=summary(log_model2)$coefficients[2,1],
                                           "LL"=summary(log_model2)$coefficients[2,1]-1.96*summary(log_model2)$coefficients[2,2],
                                           "UL"=summary(log_model2)$coefficients[2,1]+1.96*summary(log_model2)$coefficients[2,2])))
    
    OR_results3 <- as.data.frame(exp(cbind("OR"=summary(log_model3)$coefficients[2,1],
                                           "LL"=summary(log_model3)$coefficients[2,1]-1.96*summary(log_model3)$coefficients[2,2],
                                           "UL"=summary(log_model3)$coefficients[2,1]+1.96*summary(log_model3)$coefficients[2,2])))
    
    # Format results
    result = data.table(
      or1=paste0(round(OR_results1$OR,2)," (",
                 round(OR_results1$LL,2),"-",
                 round(OR_results1$UL,2),")"),
      p1=summary(log_model1)$coefficients[2,"Pr(>|z|)"],
      or2=paste0(round(OR_results2$OR,2)," (",
                 round(OR_results2$LL,2),"-",
                 round(OR_results2$UL,2),")"),
      p2=summary(log_model2)$coefficients[2,"Pr(>|z|)"],
      or3=paste0(round(OR_results3$OR,2)," (",
                 round(OR_results3$LL,2),"-",
                 round(OR_results3$UL,2),")"),
      p3=summary(log_model3)$coefficients[2,"Pr(>|z|)"]
    )
    
    logres[which(phenotype == var),] <- result
  }
  
  colnames(logres) <- c("OR1","p1","OR2","p2","OR3","p3")
  logres <- cbind(phenotype, logres)
  return(logres)
}

# Function for abdominal imaging traits
log_reg <- function(phenotype, inputdata){
  logres <- as.data.frame(matrix(NA, ncol=6, nrow=length(phenotype)))
  
  for(var in phenotype){
    # Model 1: Basic adjustment
    formula1 <- formula(paste("status ~", var, "+age + sex + ethnic + BMI + site"))
    log_model1 <- glm(formula1, family = binomial(link = "logit"), data = inputdata)
    
    # Model 2: + Socioeconomic factors
    formula2 <- formula(paste("status ~", var, "+  age + sex + ethnic + BMI + TDI + edu + site"))
    log_model2 <- glm(formula2, family = binomial(link = "logit"), data = inputdata)
    
    # Model 3: + Lifestyle factors
    formula3 <- formula(paste("status ~", var, "+  age + sex + ethnic + BMI + TDI + edu + smoke + alcohol + METmin + site + group_time"))
    log_model3 <- glm(formula3, family = binomial(link = "logit"), data = inputdata)
    
    # Calculate OR and 95% CI
    OR_results1 <- as.data.frame(exp(cbind("OR"=summary(log_model1)$coefficients[2,1],
                                           "LL"=summary(log_model1)$coefficients[2,1]-1.96*summary(log_model1)$coefficients[2,2],
                                           "UL"=summary(log_model1)$coefficients[2,1]+1.96*summary(log_model1)$coefficients[2,2])))
    
    OR_results2 <- as.data.frame(exp(cbind("OR"=summary(log_model2)$coefficients[2,1],
                                           "LL"=summary(log_model2)$coefficients[2,1]-1.96*summary(log_model2)$coefficients[2,2],
                                           "UL"=summary(log_model2)$coefficients[2,1]+1.96*summary(log_model2)$coefficients[2,2])))
    
    OR_results3 <- as.data.frame(exp(cbind("OR"=summary(log_model3)$coefficients[2,1],
                                           "LL"=summary(log_model3)$coefficients[2,1]-1.96*summary(log_model3)$coefficients[2,2],
                                           "UL"=summary(log_model3)$coefficients[2,1]+1.96*summary(log_model3)$coefficients[2,2])))
    
    # Format results
    result = data.table(
      or1=paste0(round(OR_results1$OR,2)," (",
                 round(OR_results1$LL,2),"-",
                 round(OR_results1$UL,2),")"),
      p1=summary(log_model1)$coefficients[2,"Pr(>|z|)"],
      or2=paste0(round(OR_results2$OR,2)," (",
                 round(OR_results2$LL,2),"-",
                 round(OR_results2$UL,2),")"),
      p2=summary(log_model2)$coefficients[2,"Pr(>|z|)"],
      or3=paste0(round(OR_results3$OR,2)," (",
                 round(OR_results3$LL,2),"-",
                 round(OR_results3$UL,2),")"),
      p3=summary(log_model3)$coefficients[2,"Pr(>|z|)"]
    )
    
    logres[which(phenotype == var),] <- result
  }
  
  colnames(logres) <- c("OR1","p1","OR2","p2","OR3","p3")
  logres <- cbind(phenotype, logres)
  return(logres)
}

# Function for cardiac imaging traits 
log_reg_heart <- function(phenotype, inputdata){
  logres <- as.data.frame(matrix(NA, ncol=6, nrow=length(phenotype)))
  
  for(var in phenotype){
    # Model 1: Basic adjustment
    formula1 <- formula(paste("status ~", var, "+age + sex + ethnic + BMI + site"))
    log_model1 <- glm(formula1, family = binomial(link = "logit"), data = inputdata)
    
    # Model 2: + Socioeconomic factors
    formula2 <- formula(paste("status ~", var, "+  age + sex + ethnic + BMI + TDI + edu + site"))
    log_model2 <- glm(formula2, family = binomial(link = "logit"), data = inputdata)
    
    # Model 3: + Lifestyle factors + body surface area
    formula3 <- formula(paste("status ~", var, "+  age + sex + ethnic + BMI + TDI + edu + smoke + alcohol + METmin + site + bodysurface + group_time"))
    log_model3 <- glm(formula3, family = binomial(link = "logit"), data = inputdata)
    
    # Calculate OR and 95% CI
    OR_results1 <- as.data.frame(exp(cbind("OR"=summary(log_model1)$coefficients[2,1],
                                           "LL"=summary(log_model1)$coefficients[2,1]-1.96*summary(log_model1)$coefficients[2,2],
                                           "UL"=summary(log_model1)$coefficients[2,1]+1.96*summary(log_model1)$coefficients[2,2])))
    
    OR_results2 <- as.data.frame(exp(cbind("OR"=summary(log_model2)$coefficients[2,1],
                                           "LL"=summary(log_model2)$coefficients[2,1]-1.96*summary(log_model2)$coefficients[2,2],
                                           "UL"=summary(log_model2)$coefficients[2,1]+1.96*summary(log_model2)$coefficients[2,2])))
    
    OR_results3 <- as.data.frame(exp(cbind("OR"=summary(log_model3)$coefficients[2,1],
                                           "LL"=summary(log_model3)$coefficients[2,1]-1.96*summary(log_model3)$coefficients[2,2],
                                           "UL"=summary(log_model3)$coefficients[2,1]+1.96*summary(log_model3)$coefficients[2,2])))
    
    # Format results
    result = data.table(
      or1=paste0(round(OR_results1$OR,2)," (",
                 round(OR_results1$LL,2),"-",
                 round(OR_results1$UL,2),")"),
      p1=summary(log_model1)$coefficients[2,"Pr(>|z|)"],
      or2=paste0(round(OR_results2$OR,2)," (",
                 round(OR_results2$LL,2),"-",
                 round(OR_results2$UL,2),")"),
      p2=summary(log_model2)$coefficients[2,"Pr(>|z|)"],
      or3=paste0(round(OR_results3$OR,2)," (",
                 round(OR_results3$LL,2),"-",
                 round(OR_results3$UL,2),")"),
      p3=summary(log_model3)$coefficients[2,"Pr(>|z|)"]
    )
    
    logres[which(phenotype == var),] <- result
  }
  
  colnames(logres) <- c("OR1","p1","OR2","p2","OR3","p3")
  logres <- cbind(phenotype, logres)
  return(logres)
}

# ==============================================================================
# Run Association Analysis for Each Organ
# ==============================================================================

# Brain imaging traits
logres_cmd_dep_brain <- log_reg_brain(pheno_brain, final_cmd_dep_brain)

# Cardiac imaging traits
logres_cmd_dep_heart <- log_reg_heart(pheno_heart, final_cmd_dep_heart)

# Abdominal imaging traits
logres_cmd_dep_ab <- log_reg(pheno_ab, final_cmd_dep_ab)

# ==============================================================================
# FDR Correction and Selection of Significant Traits
# ==============================================================================

# Apply FDR correction to p-values from Model 3
logres_cmd_dep_brain$fdr_p <- p.adjust(logres_cmd_dep_brain$p3, method="fdr")
sig_fdr_cmd_dep_brain <- logres_cmd_dep_brain[logres_cmd_dep_brain$fdr_p<0.05,]

logres_cmd_dep_heart$fdr_p <- p.adjust(logres_cmd_dep_heart$p3, method="fdr")
sig_fdr_cmd_dep_heart <- logres_cmd_dep_heart[logres_cmd_dep_heart$fdr_p<0.05,]

logres_cmd_dep_ab$fdr_p <- p.adjust(logres_cmd_dep_ab$p3, method="fdr")
sig_fdr_cmd_dep_ab <- logres_cmd_dep_ab[logres_cmd_dep_ab$fdr_p<0.05,]

# Combine significant traits from all organs
sig_fdr_cmd_dep <- rbind(sig_fdr_cmd_dep_ab, sig_fdr_cmd_dep_heart, sig_fdr_cmd_dep_brain)
