library(dplyr)
library(readxl)
library(dplyr)
library(survival)
library(glmnet)
library(parallel)
library(doParallel)
library(caret)
library(openxlsx)
library(CsChange)
library(data.table)
library(mice)
library(bigreadr)
library(Hmisc)
library(survival)       
library(prodlim)
library(pec)
library(tidyverse)
library(survIDINRI)

#---------------------read files-----------------------------------------------#
setwd("/path")

disease_code<-fread("DiseaseCMDandDep.csv") 
colnames(disease_code)[1]<-"eid"

#abdominal 
#liver PDFF/liver iron/liver volume/pancreas PDFF/pancreas volume/Subcutaneous fat volume/Visceral fat volume
comp_abmri_out<-fread("comp_abmri_out.csv") %>% 
  dplyr::select(.,-1) %>% 
  dplyr::select(.,c("eid","liver_iron","liver_pdff","n_21080_2_0","n_21090_2_0","n_21087_2_0",
                    "n_21086_2_0","n_21085_2_0"))
#heart
comp_heartmri_out<-fread("comp_heartmri_out.csv") %>% dplyr::select(.,-1)

#brain
#DTI-FA skeleton
comp_DTI_out <- fread("comp_DTI_out.csv") %>% dplyr::select(.,-1)
#subcortical volumetric segmentation
comp_volume_out <- fread("comp_volume_out.csv") %>% dplyr::select(.,-1)

brainlist <- list(comp_DTI_out,comp_volume_out)
brain <- brainlist %>% reduce(inner_join, by = "eid") 
#other covariates
covariate<-fread("COV2.csv")
colnames(covariate)[1]<-"eid"
#biochemical model
bioche<-read.xlsx("Bioche.xlsx")
colnames(bioche)[1]<-"eid"
dput(names(bioche))
#---------------------------select traits--------------------------------------#
data <- read.xlsx("surface_traits.xlsx")
dput(data$Code)
liver <- comp_abmri_out %>% dplyr::select(.,"eid","n_21080_2_0")
heart <- comp_heartmri_out %>% dplyr::select(.,"eid","n_24181_2_0","n_24134_2_0","n_24128_2_0","n_24129_2_0")
brain <- brain %>% dplyr::select(.,"eid","n_25059_2_0", "n_25082_2_0", "n_25085_2_0", "n_25367_2_0", 
                                 "n_25347_2_0", "n_25346_2_0", "n_25371_2_0", "n_25370_2_0", "n_25373_2_0", 
                                 "n_25372_2_0", "n_25364_2_0", "n_25375_2_0", "n_25369_2_0", "n_25368_2_0", 
                                 "n_25126_2_0", "n_25005_2_0")

trait_cmd_list <- list(liver,heart,brain,covariate,bioche,disease_code)
trait_cmd <- trait_cmd_list %>% reduce(inner_join, by = "eid") 

vars <- c("Total_chol", "CRP", "HbA1c", "HDL", "LDL", "Triglycerides", "SBP")
data.frame(
  variable = vars,
  missing = sapply(vars, function(x) sum(is.na(trait_cmd[[x]]))),
  total = nrow(trait_cmd),
  missing_rate = sapply(vars, function(x) round(sum(is.na(trait_cmd[[x]]))/nrow(trait_cmd)*100, 2))
)

data_to_impute <- trait_cmd[, ..vars]
imputed_data <- mice(data_to_impute, method = 'pmm', m = 10, seed = 123)
completed_data <- complete(imputed_data)
trait_cmd[, vars] <- completed_data

#---------------------------imputation functions--------------------------------------#

imputation_cov<-function(merge_cov_out){
  cov_set <- c("age","sex","BMI","TDI","edu","smoke","alcohol","METmin","Ethnicity")
  
  merge_cov_out$sex<-merge_cov_out$n_31_0_0
  merge_cov_out$age<-merge_cov_out$n_21003_2_0
  merge_cov_out$BMI<-merge_cov_out$n_21001_2_0
  merge_cov_out$TDI<-merge_cov_out$n_22189_0_0
  merge_cov_out$edu<-case_when(merge_cov_out$n_6138_0_0 == 1 | merge_cov_out$n_6138_0_0 == 6 ~ 1,
                               merge_cov_out$n_6138_0_0 == 2 | merge_cov_out$n_6138_0_0 == 3 | merge_cov_out$n_6138_0_0 == 4 ~ 2,
                               merge_cov_out$n_6138_0_0 == 5 ~ 3,
                               merge_cov_out$n_6138_0_0 == -7 ~ 4,
                               merge_cov_out$n_6138_0_0 == -3 ~ NA)
  merge_cov_out <- merge_cov_out %>%
    mutate(Ethnicity = case_when(
      n_21000_0_0 >= 1001 & n_21000_0_0 <= 1003 ~ "1",
      n_21000_0_0 %in% c(-1, -3) ~ NA_character_,
      !is.na(n_21000_0_0) ~ "2"
    ))
  merge_cov_out$smoke<-ifelse(merge_cov_out$n_20116_2_0==-3,NA,merge_cov_out$n_20116_2_0)
  merge_cov_out$alcohol<-ifelse(merge_cov_out$n_20117_2_0==-3,NA,merge_cov_out$n_20117_2_0)
  
  merge_cov_out$MET1<-merge_cov_out$n_864_2_0
  merge_cov_out$MET2<-merge_cov_out$n_874_2_0
  merge_cov_out$MET3<-merge_cov_out$n_884_2_0
  merge_cov_out$MET4<-merge_cov_out$n_894_2_0
  merge_cov_out$MET5<-merge_cov_out$n_904_2_0
  merge_cov_out$MET6<-merge_cov_out$n_914_2_0
  merge_cov_out$MET1[merge_cov_out$MET1==-1|merge_cov_out$MET1==-2|merge_cov_out$MET1==-3]=0
  merge_cov_out$MET2[merge_cov_out$MET2==-1|merge_cov_out$MET2==-3|is.na(merge_cov_out$MET2)]=0
  merge_cov_out$MET3[merge_cov_out$MET3==-1|merge_cov_out$MET3==-2|merge_cov_out$MET3==-3]=0
  merge_cov_out$MET4[merge_cov_out$MET4==-1|merge_cov_out$MET4==-3|is.na(merge_cov_out$MET4)]=0
  merge_cov_out$MET5[merge_cov_out$MET5==-1|merge_cov_out$MET5==-2|merge_cov_out$MET5==-3]=0
  merge_cov_out$MET6[merge_cov_out$MET6==-1|merge_cov_out$MET6==-3|is.na(merge_cov_out$MET6)]=0
  merge_cov_out$METmin=merge_cov_out$MET1*merge_cov_out$MET2+merge_cov_out$MET3*merge_cov_out$MET4*3+merge_cov_out$MET5*merge_cov_out$MET6*6
  merge_cov_out$METmin<-ifelse(merge_cov_out$METmin %in% c(0), NA, merge_cov_out$METmin) 
  
  merge_cov_out$sex <- factor(merge_cov_out$sex)
  merge_cov_out$edu <- factor(merge_cov_out$edu)
  merge_cov_out$smoke <- factor(merge_cov_out$smoke)
  merge_cov_out$alcohol <- factor(merge_cov_out$alcohol)
  merge_cov_out$Ethnicity <- factor(merge_cov_out$Ethnicity)
  
  data_to_impute <- merge_cov_out[, ..cov_set]
  imputed_data <- mice(data_to_impute, method = 'pmm', m = 10, seed = 123)
  completed_data <- complete(imputed_data)
  merge_cov_out[, cov_set] <- completed_data
  
  cov_set <- c("age","sex","BMI","TDI","edu","smoke","alcohol","METmin","Ethnicity",
               "Total_chol", "CRP", "HbA1c", "HDL", "LDL", "Triglycerides", "SBP")
  Mdata<-merge_cov_out %>% dplyr::select(.,c("eid",cov_set))
  return(Mdata)
}

Mdata_cmd<-imputation_cov(trait_cmd)

#--------------------------------Multimorbidity progression-------------------------#

#Depression--MRI--CMDs

cox_cmd_dep_1_1 <- disease_code %>% 
  filter(.,time_CMD1 > 0 & time_DEP <= 0) %>% 
  mutate(.,status=ifelse((event_CMD1==1 & event_DEP==1),1,0)) %>%
  as.data.frame(.) 

cox_cmd_dep_1_1$time <- apply(cox_cmd_dep_1_1,1,function(x){
  max(x[3],x[21])
})

final_cmd_dep_1_1list <- list(Mdata_cmd,liver,heart,brain,cox_cmd_dep_1_1)
final_cmd_dep_1_1 <- final_cmd_dep_1_1list %>% reduce(inner_join, by = "eid")

#CMDs--MRI--Depression
cox_cmd_dep_1_2 <- disease_code %>% 
  filter(.,time_CMD1 <=0  & time_DEP >0) %>% 
  mutate(.,status=ifelse((event_CMD1==1 & event_DEP==1),1,0)) %>%
  as.data.frame(.) 

cox_cmd_dep_1_2$time <- apply(cox_cmd_dep_1_2,1,function(x){
  max(x[3],x[21])
})

final_cmd_dep_1_2list <- list(Mdata_cmd,liver,heart,brain,cox_cmd_dep_1_2)
final_cmd_dep_1_2 <- final_cmd_dep_1_2list %>% reduce(inner_join, by = "eid")

#-----------------------------Multimorbidity development---------------------------#

# free of both conditions--MRI--Multimorbidity
cox_cmd_dep_2 <- disease_code %>% 
  filter(.,time_CMD1 > 0 & time_DEP>0) %>% 
  mutate(.,status=ifelse((event_CMD1==1 & event_DEP==1),1,0)) %>%
  as.data.frame(.) 

cox_cmd_dep_2$time <- apply(cox_cmd_dep_2,1,function(x){
  max(x[3],x[21])
})

final_cmd_dep_2list <- list(Mdata_cmd,liver,heart,brain,cox_cmd_dep_2)
final_cmd_dep_2 <- final_cmd_dep_2list %>% reduce(inner_join, by = "eid")


#------------------------------Model evaluation-----------------------------------#

library(glmnet)
library(survival)
library(caret)
library(doParallel)
library(foreach)
library(openxlsx)
library(rms)
library(pec)

# construct lifestyle score, biochemical score, and MRI score
build_scores_cv <- function(data) {
  data[,c(11:38)] <-  as.data.frame(scale(data[, 11:38]))
  set.seed(123)
  folds <- createFolds(1:nrow(data), k = 10, list = TRUE, returnTrain = FALSE)
  score_vars <- list(
    Lifestyle_score = c("BMI", "smoke", "alcohol", "METmin"),
    Biochemical_score = c("Total_chol", "CRP", "HbA1c", "HDL", "LDL", "Triglycerides", "SBP"),
    MRI_score = paste0("n_", c(21080, 24181, 24134, 24128, 24129, 25059, 25082, 25085, 
                               25367, 25347, 25346, 25371, 25370, 25373, 25372, 25364, 
                               25375, 25369, 25368, 25126, 25005), "_2_0")
  )
  
  scores_df <- data.frame(eid = data$eid)
  coef_results <- list()
  
  for (score_name in names(score_vars)) {
    cat("construct", score_name, "...\n")
    
    current_vars <- score_vars[[score_name]]
    available_vars <- current_vars[current_vars %in% colnames(data)]
    
    if (length(available_vars) == 0) {
      warning("scores ", score_name, " NA")
      scores_df[[score_name]] <- 0
      next
    }
    
    cv_coefs <- list()
    cv_scores <- matrix(NA, nrow = nrow(data), ncol = 10)
    
    for (fold in 1:10) {
      test_idx <- folds[[fold]]
      train_idx <- setdiff(1:nrow(data), test_idx)
      
      train_data <- as.data.frame(data[train_idx, ])
      test_data <- as.data.frame(data[test_idx, ])
      
      x_train_raw <- train_data[, available_vars, drop = FALSE]
      x_test_raw <- test_data[, available_vars, drop = FALSE]
      x_train_scaled <- x_train_raw
      x_test_scaled <- x_test_raw
      
      x_train <- model.matrix(~ . - 1, data = x_train_scaled)
      x_test <- model.matrix(~ . - 1, data = x_test_scaled)
      
      y_train <- Surv(train_data$time, train_data$status)
      
      cv_fit <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0.5, nfolds = 10)
      model <- glmnet(x_train, y_train, family = "cox", alpha = 0.5, lambda = cv_fit$lambda.min)
      coefs <- as.vector(coef(model))
      
      scores <- as.vector(x_test %*% coefs)
      
      cv_coefs[[fold]] <- coefs
      cv_scores[test_idx, fold] <- scores
    }
    
    if (length(cv_coefs) > 0) {
      max_length <- max(sapply(cv_coefs, length))
      coef_matrix <- matrix(NA, nrow = max_length, ncol = 10)
      
      for (fold in 1:10) {
        if (length(cv_coefs[[fold]]) > 0) {
          coef_matrix[1:length(cv_coefs[[fold]]), fold] <- cv_coefs[[fold]]
        }
      }
      
      avg_coefs <- rowMeans(coef_matrix, na.rm = TRUE)
      
      if (exists("x_train") && length(avg_coefs) == ncol(x_train)) {
        names(avg_coefs) <- colnames(x_train)
      } else {
        names(avg_coefs) <- paste0("Var_", 1:length(avg_coefs))
      }
    } else {
      avg_coefs <- numeric(0)
    }
    coef_results[[score_name]] <- avg_coefs
    
    final_scores <- rowMeans(cv_scores, na.rm = TRUE)
    scores_df[[score_name]] <- final_scores
  }
  
  return(list(scores = scores_df, coefficients = coef_results))
}

FUNChange <- function(model0, model1, data) {
  modelcompare <- CsChange(model0, model1, data = data, nb = 200)
  change <- round(modelcompare[[1]]$change, 3)
  change.lwr <- round(modelcompare[[1]]$low, 3)
  change.upr <- round(modelcompare[[1]]$up, 3)
  paste0(change, " (", change.lwr, ", ", change.upr, ")")
}

get_royston_d <- function(model) {
  d = royston(model)[1]
  se_d = royston(model)[2]
  d.lwr = d - 1.96 * se_d
  d.upr = d + 1.96 * se_d
  paste0(round(d, 3), " (", round(d.lwr, 3), ", ", round(d.upr, 3), ")")
}
calculate_nri_idi <- function(data, base_model, new_model, time_point = 8, 
                              risk_thresholds = c(0, 0.1, 0.2, 1),
                              n_boot = 50, continuous_nri = TRUE) {
  
  # ---- 辅助函数：计算个体风险 ----
  get_risk <- function(model, data, time_point) {
    tryCatch({
      lp <- predict(model, type = "lp", newdata = data)
      base_haz <- basehaz(model, centered = FALSE)
      
      if (time_point <= max(base_haz$time)) {
        closest_idx <- which.min(abs(base_haz$time - time_point))
        base_cumhaz <- base_haz$hazard[closest_idx]
      } else {
        base_cumhaz <- base_haz$hazard[nrow(base_haz)]
      }
      
      risk <- 1 - exp(-base_cumhaz * exp(lp))
      
      risk <- pmin(pmax(risk, 0), 1)
      
      return(risk)
      
    }, error = function(e) {
      warning("Using simplified risk calculation due to error: ", e$message)
      lp <- predict(model, type = "lp", newdata = data)
      risk <- 1 - exp(-exp(lp) * time_point)
      return(pmin(pmax(risk, 0), 1))
    })
  }
  
  # ---- sampling ----#
  core_calc <- function(index = 1:nrow(data)) {
    d <- data[index, ]
    
    base_risk <- get_risk(base_model, d, time_point)
    new_risk  <- get_risk(new_model, d, time_point)
    
    events <- ifelse(d$time <= time_point & d$status == 1, 1, 0)
    event_indices <- which(events == 1)
    nonevent_indices <- which(events == 0)
    
    if (length(event_indices) == 0 || length(nonevent_indices) == 0) {
      return(c(nri = NA, idi = NA))
    }
    
    # -------- NRI calculation -------- #
    if (continuous_nri) {
      event_improve <- mean(new_risk[event_indices] > base_risk[event_indices], na.rm = TRUE)
      event_worsen <- mean(new_risk[event_indices] < base_risk[event_indices], na.rm = TRUE)
      
      nonevent_improve <- mean(new_risk[nonevent_indices] < base_risk[nonevent_indices], na.rm = TRUE)
      nonevent_worsen <- mean(new_risk[nonevent_indices] > base_risk[nonevent_indices], na.rm = TRUE)
      
      event_diff <- event_improve - event_worsen
      nonevent_diff <- nonevent_improve - nonevent_worsen
      nri_total <- event_diff + nonevent_diff
      
    } else {
      # NRI classification
      base_cat <- cut(base_risk, breaks = risk_thresholds, include.lowest = TRUE, labels = FALSE)
      new_cat  <- cut(new_risk, breaks = risk_thresholds, include.lowest = TRUE, labels = FALSE)
      
      # Check if there is a valid classification
      if (all(is.na(base_cat)) || all(is.na(new_cat))) {
        return(c(nri = NA, idi = NA))
        
      }
      
      event_up    <- sum(events == 1 & new_cat > base_cat, na.rm = TRUE)
      event_down  <- sum(events == 1 & new_cat < base_cat, na.rm = TRUE)
      nonevent_up   <- sum(events == 0 & new_cat > base_cat, na.rm = TRUE)
      nonevent_down <- sum(events == 0 & new_cat < base_cat, na.rm = TRUE)
      
      # components of NRI
      if (length(event_indices) > 0) {
        nri_event <- (event_up - event_down) / length(event_indices)
      } else {
        nri_event <- 0
      }
      
      if (length(nonevent_indices) > 0) {
        nri_nonevent <- (nonevent_down - nonevent_up) / length(nonevent_indices)
      } else {
        nri_nonevent <- 0
      }
      
      nri_total <- nri_event + nri_nonevent
    }
    
    # -------- IDI calculation -------- #
    # IS: Mean prediction probability difference in the event group
    is_diff <- mean(new_risk[event_indices], na.rm = TRUE) - mean(base_risk[event_indices], na.rm = TRUE)
    
    # IP: Difference in mean model-predicted probabilities for the non-event group
    ip_diff <- mean(new_risk[nonevent_indices], na.rm = TRUE) - mean(base_risk[nonevent_indices], na.rm = TRUE)
    
    # IDI = IS - IP
    idi <- is_diff - ip_diff
    
    return(c(nri = nri_total, idi = idi))
  }
  
  # ---- Original estimation ---- #
  tryCatch({
    est <- core_calc()
    
    # In case of estimation failure, the output defaults to NA
    if (any(is.na(est))) {
      return(list(
        nri = c(NA, NA, NA),
        idi = c(NA, NA, NA)
      ))
    }
    
    # ---- Bootstrap confidence intervals ----
    boot_results <- replicate(n_boot, {
      idx <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
      core_calc(idx)
    }, simplify = FALSE)
    
    # extract bootstrap results
    boot_matrix <- do.call(rbind, boot_results)
    valid_boots <- complete.cases(boot_matrix)
    
    if (sum(valid_boots) < 5) {
      warning("Insufficient valid bootstrap samples for confidence intervals")
      return(list(
        nri = c(est["nri"], NA, NA),
        idi = c(est["idi"], NA, NA)
      ))
    }
    # Confidence intervals
    valid_boot_matrix <- boot_matrix[valid_boots, , drop = FALSE]
    nri_ci <- quantile(valid_boot_matrix[, "nri"], c(0.025, 0.975), na.rm = TRUE)
    idi_ci <- quantile(valid_boot_matrix[, "idi"], c(0.025, 0.975), na.rm = TRUE)
    
    return(list(
      nri = c(est["nri"], nri_ci),
      idi = c(est["idi"], idi_ci)
    ))
    
  }, error = function(e) {
    warning("NRI/IDI calculation failed: ", e$message)
    return(list(
      nri = c(NA, NA, NA),
      idi = c(NA, NA, NA)
    ))
  })
}

# Functions for model evaluation

evaluate_prediction_models <- function(data_with_scores) {
  model_formulas <- list(
    model1 = "Surv(time, status) ~ age + sex + Ethnicity",
    model2 = "Surv(time, status) ~ age + sex + Ethnicity + TDI + edu",
    model3 = "Surv(time, status) ~ age + sex + Ethnicity + TDI + edu + Lifestyle_score",
    model4 = "Surv(time, status) ~ age + sex + Ethnicity + TDI + edu + Biochemical_score",
    model5 = "Surv(time, status) ~ age + sex + Ethnicity + TDI + edu + MRI_score",
    model6 = "Surv(time, status) ~ age + sex + Ethnicity + TDI + edu + Lifestyle_score + MRI_score",
    model7 = "Surv(time, status) ~ age + sex + Ethnicity + TDI + edu + Biochemical_score + MRI_score"
  )
  
  models <- list()
  for (i in 1:7) {
    models[[paste0("model", i)]] <- coxph(as.formula(model_formulas[[i]]), data = data_with_scores, x=TRUE, y=TRUE)
  }
  
  results <- list()
  
  # ---- 1. C-index ----
  c_indices <- sapply(models, function(m) {
    c_obj <- concordance(m)
    c_val <- c_obj$concordance
    c_se <- sqrt(c_obj$var)
    c_lwr <- c_val - 1.96 * c_se
    c_upr <- c_val + 1.96 * c_se
    paste0(round(c_val, 3), " (", round(c_lwr, 3), ", ", round(c_upr, 3), ")")
  })

  # ---- 2. Royston’s D ----
  royston_d_results <- sapply(models, function(m) {
    tryCatch({
      get_royston_d(m)
    }, error = function(e) "Failed")
  })
  
  # ---- 3. Brier score ----
  tryCatch({
    pred_obj <- pec(object = models,
                    formula = Surv(time, status) ~ .,
                    data = data_with_scores,
                    times = 8)
    
    brier_scores <- tryCatch({
      AppErr <- pred_obj$AppErr
      if ("8" %in% rownames(AppErr)) {
        scores <- AppErr["8", ]
      } else {
        scores <- apply(AppErr, 2, function(x) mean(x, na.rm = TRUE))  # 平均Brier score
      }
      round(as.numeric(scores), 3)
    }, error = function(e) {
      warning("Brier score calculation failed: ", e$message)
      rep(NA, length(models))
    })
    names(brier_scores) <- names(models)
  })
  
  # ---- 4. Model comparison ----
  model_comparisons <- list()
  for (i in 1:6) {
    for (j in (i+1):7) {
      model_i_name <- paste0("model", i)
      model_j_name <- paste0("model", j)
      tryCatch({
        comparison_name <- paste0(model_i_name, "_vs_", model_j_name)
        model_comparisons[[comparison_name]] <- FUNChange(models[[model_i_name]],
                                                          models[[model_j_name]],
                                                          data_with_scores)
      }, error = function(e) {
        comparison_name <- paste0(model_i_name, "_vs_", model_j_name)
        model_comparisons[[comparison_name]] <- "Failed"
      })
    }
  }

  ##---- 5. NRI and IDI ----
  nri_idi_results <- list()
  base_model <- models$model1
  for (i in 2:7) {
    model_name <- paste0("model", i)
    nri_idi_result <- calculate_nri_idi(data_with_scores,
                                        base_model = base_model,
                                        new_model = models[[model_name]],
                                        time_point = 8,
                                        continuous_nri = TRUE)
    if (!is.na(nri_idi_result$nri[1])) {
      nri_formatted <- paste0(round(nri_idi_result$nri[1], 3), " (",
                              round(nri_idi_result$nri[2], 3), ", ",
                              round(nri_idi_result$nri[3], 3), ")")
      idi_formatted <- paste0(round(nri_idi_result$idi[1], 3), " (",
                              round(nri_idi_result$idi[2], 3), ", ",
                              round(nri_idi_result$idi[3], 3), ")")
    } else {
      nri_formatted <- "Failed"
      idi_formatted <- "Failed"
    }
    nri_idi_results[[paste0(model_name, "_NRI")]] <- nri_formatted
    nri_idi_results[[paste0(model_name, "_IDI")]] <- idi_formatted
  }

  results$c_indices <- c_indices
  results$royston_d <- royston_d_results
  results$brier_scores <- brier_scores
  results$model_comparisons <- model_comparisons
  results$nri_idi <- nri_idi_results
  results$models <- models
  
  return(results)
}

run_complete_analysis <- function(datasets) {
  all_results <- list()
  
  for (dataset_name in names(datasets)) {
    cat("Dataset:", dataset_name, "\n")
    # 1. construct scores
    score_results <- build_scores_cv(datasets[[dataset_name]])
    # 2. combine dataset
    data_with_scores <- merge(as.data.frame(datasets[[dataset_name]]), 
                              score_results$scores, by = "eid", all.x = TRUE)
    # 3. model evaluation
    model_results <- evaluate_prediction_models(data_with_scores)
    # ---- 4. Output calibration curves as PDF  ----
    pdf(paste0("Calibration_", dataset_name, ".pdf"), width = 12, height = 9)
    par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
    for (i in 1:7) {
      m_name <- paste0("model", i)
      tryCatch({
        # Use existing model for direct calibration calculation
        model <- model_results$models[[m_name]]
        # Calculate predicted risk
        pred_lp <- predict(model, type = "lp")
        pred_risk <- 1 - exp(-exp(pred_lp) * 8)  # 8-year risk
        # Observed outcomes
        observed <- ifelse(data_with_scores$time <= 8 & data_with_scores$status == 1, 1, 0)
        # Group by risk for calibration
        risk_groups <- cut(pred_risk, breaks = quantile(pred_risk, probs = seq(0,1,0.1)),
                           include.lowest = TRUE)
        expected_risk <- tapply(pred_risk, risk_groups, mean, na.rm = TRUE)
        observed_risk <- tapply(observed, risk_groups, mean, na.rm = TRUE)
        valid_points <- !is.na(expected_risk) & !is.na(observed_risk)
        # Plot calibration curve
        plot(expected_risk[valid_points], observed_risk[valid_points],
             xlim = c(0,1), ylim = c(0,1),
             xlab = "Predicted 8-year risk", ylab = "Observed 8-year risk",
             main = paste0(dataset_name, " - ", m_name),
             pch = 16, col = "blue")
        lines(expected_risk[valid_points], observed_risk[valid_points], col="blue", lwd=2)
        # Ideal calibration line
        abline(0, 1, col = "red", lty = 2, lwd = 2)
        # Add fitted line
        if (sum(valid_points) > 2) {
          fit <- lm(observed_risk[valid_points] ~ expected_risk[valid_points])
          abline(fit, col = "darkgreen", lwd = 2)
        }
      }, error = function(e) {
        plot(0, 0, type = "n", axes = FALSE, ann = FALSE,
             main = paste0(dataset_name, " - ", m_name))
        text(0, 0, "Calibration failed", cex = 1.2)
      })
    }
    dev.off()
    all_results[[dataset_name]] <- list(
      score_coefficients = score_results$coefficients,
      individual_scores = score_results$scores,
      model_performance = model_results,
      data_with_scores = data_with_scores
    )
  }
  return(all_results)
}
# Function to save complete results
save_complete_results <- function(all_results, filename) {
  wb <- createWorkbook()
  # 1. Score coefficients table
  addWorksheet(wb, "Score_Coefficients")
  coef_summary <- data.frame()
  for (dataset_name in names(all_results)) {
    for (score_name in names(all_results[[dataset_name]]$score_coefficients)) {
      coefs <- all_results[[dataset_name]]$score_coefficients[[score_name]]
      temp_df <- data.frame(
        Dataset = dataset_name,
        Score = score_name,
        Variable = names(coefs),
        Coefficient = as.numeric(coefs),
        stringsAsFactors = FALSE
      )
      coef_summary <- rbind(coef_summary, temp_df)
    }
  }
  writeData(wb, "Score_Coefficients", coef_summary)
  # 2. Model C-index table
  addWorksheet(wb, "Model_C_Index")
  c_index_df <- data.frame(
    Dataset = names(all_results),
    stringsAsFactors = FALSE
  )
  for (i in 1:7) {
    model_name <- paste0("Model", i)
    c_index_df[[model_name]] <- sapply(all_results, function(x) x$model_performance$c_indices[[paste0("model", i)]])
  }
  writeData(wb, "Model_C_Index", c_index_df)
  # 3. Royston's D statistics table
  addWorksheet(wb, "Royston_D")
  royston_d_df <- data.frame(
    Dataset = names(all_results),
    stringsAsFactors = FALSE
  )
  for (i in 1:7) {
    model_name <- paste0("Model", i)
    royston_d_df[[model_name]] <- sapply(all_results, function(x) x$model_performance$royston_d[[paste0("model", i)]])
  }
  writeData(wb, "Royston_D", royston_d_df)
  # 4. Brier Scores table
  addWorksheet(wb, "Brier_Scores")
  brier_df <- data.frame(
    Dataset = names(all_results),
    stringsAsFactors = FALSE
  )
  for (i in 1:7) {
    model_name <- paste0("Model", i)
    brier_df[[model_name]] <- sapply(all_results, function(x) x$model_performance$brier_scores[[paste0("model", i)]])
  }
  writeData(wb, "Brier_Scores", brier_df)
  # 5. Model comparisons table
  addWorksheet(wb, "Model_Comparisons")
  comparison_names <- names(all_results[[1]]$model_performance$model_comparisons)
  comparison_df <- data.frame(
    Dataset = names(all_results),
    stringsAsFactors = FALSE
  )
  for (comp_name in comparison_names) {
    comparison_df[[comp_name]] <- sapply(all_results, function(x) x$model_performance$model_comparisons[[comp_name]])
  }
  writeData(wb, "Model_Comparisons", comparison_df)
  # 6. NRI and IDI results table
  addWorksheet(wb, "NRI_IDI_Results")
  nri_idi_names <- names(all_results[[1]]$model_performance$nri_idi)
  nri_idi_df <- data.frame(
    Dataset = names(all_results),
    stringsAsFactors = FALSE
  )
  for (nri_idi_name in nri_idi_names) {
    nri_idi_df[[nri_idi_name]] <- sapply(all_results, function(x) x$model_performance$nri_idi[[nri_idi_name]])
  }
  writeData(wb, "NRI_IDI_Results", nri_idi_df)
  # 7. Individual score data
  for (dataset_name in names(all_results)) {
    sheet_name <- paste0(dataset_name, "_Scores")
    if (nchar(sheet_name) > 31) {
      sheet_name <- substr(sheet_name, 1, 31)
    }
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, all_results[[dataset_name]]$individual_scores)
  }
  saveWorkbook(wb, filename, overwrite = TRUE)
}
# Execute complete analysis
datasets <- list(
  "CMD_DEP_1_1" = final_cmd_dep_1_1,
  "CMD_DEP_1_2" = final_cmd_dep_1_2,  
  "CMD_DEP_2" = final_cmd_dep_2
)
# Run analysis
complete_results <- run_complete_analysis(datasets)
# Save results
save_complete_results(complete_results, "Complete_CMD_Analysis_Results_4.xlsx")
data_score <- read.xlsx("Complete_CMD_Analysis_Results_4.xlsx",sheet = 9)
final_cmd_dep_analysis1 <- merge(final_cmd_dep_2,data_score,by="eid")
model_formulas <- list(
  model1 = "Surv(time, status) ~ age + sex",
  model2 = "Surv(time, status) ~ age + sex + TDI + edu",
  model3 = "Surv(time, status) ~ age + sex + TDI + edu + Lifestyle_score",
  model4 = "Surv(time, status) ~ age + sex + TDI + edu + Biochemical_score",
  model5 = "Surv(time, status) ~ age + sex + TDI + edu + MRI_score",
  model6 = "Surv(time, status) ~ age + sex + TDI + edu + Lifestyle_score + MRI_score",
  model7 = "Surv(time, status) ~ age + sex + TDI + edu + Biochemical_score + MRI_score"
)
# Fit all models
models <- list()
for (i in 1:7) {
  models[[paste0("model", i)]] <- coxph(as.formula(model_formulas[[i]]), data = final_cmd_dep_analysis1, x=TRUE, y=TRUE)
}
pred_obj <- pec(
  object = models,
  formula = Surv(time, status) ~ 1,   # Only specify survival outcome, no additional variables
  data = final_cmd_dep_analysis1,
  cens.model = "cox",
  splitMethod = "boot",
  B = 200,  # bootstrap iterations, recommended 100-200
  times = seq(4, max(final_cmd_dep_analysis1$time), by = 1) # multiple time points
)
boot_results <- pred_obj$AppErr
brier_with_ci <- lapply(names(models), function(model_name) {
  boot_mat <- boot_results[[model_name]]
  if (is.null(boot_mat)) {
    return(NULL)  # Return NULL if no results
  }
  # If vector, convert to 1-row matrix
  if (is.vector(boot_mat)) {
    boot_mat <- matrix(boot_mat, nrow = 1)
    rownames(boot_mat) <- names(boot_results[[model_name]])
  }
  # Ensure matrix before apply
  res <- apply(boot_mat, 1, function(x) {
    m <- mean(x, na.rm = TRUE)
    l <- quantile(x, 0.025, na.rm = TRUE)
    u <- quantile(x, 0.975, na.rm = TRUE)
    sprintf("%.3e (%.3e, %.3e)", m, l, u)
  })
  return(res)
})
brier_ci_tables <- lapply(seq_along(brier_with_ci), function(i) {
  df <- as.data.frame(t(brier_with_ci[[i]]))
  df$time <- as.numeric(rownames(df))
  df$model <- names(models)[i]
  colnames(df)[1] <- "Brier_score"
  df
})
brier_ci_df <- do.call(rbind, brier_ci_tables)
brier_ci_df <- brier_ci_df[, c("model", "time", "Brier_score")]
head(brier_ci_df)
write.xlsx(brier_ci_df, "Brier_scores_with_CI_2.xlsx", overwrite = TRUE)
