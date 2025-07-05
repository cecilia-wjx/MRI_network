library(tidyverse)
library(bigreadr)
library(writexl)
library(readxl)
library(stringr)

source("/FunctionSet.R")

for(i in 1:nrow(all)){
  name_heart<-all$heart_file[i]
  name_ab<-all$ab_file[i]
  name_brain<-all$brain_file[i]
  heart_set<-heart_list[[name_heart]]
  brain_set<-brain_list[[name_brain]]
  ab_set<-ab_list[[name_ab]]
  
  pairs3_list<-list(heart_set,ab_set,brain_set)
  pairs3<-pairs3_list %>% reduce(inner_join,by=c("rsid"))
  
  pairs3$beta_ab <- ifelse((pairs3$A1.x == pairs3$A1.y & pairs3$A2.x == pairs3$A2.y), 
                           pairs3$beta_ab, pairs3$beta_ab* (-1))
  pairs3$beta_brain <- ifelse((pairs3$A1.x == pairs3$A1 & pairs3$A2.x == pairs3$A2), 
                              pairs3$beta_brain, pairs3$beta_brain * (-1))
  
  
  # Remove SNPs with Z-score larger than 1.96 or smaller than 1.96
  pairs3$z1<-pairs3$beta_ab/pairs3$se_ab
  pairs3$z2<-pairs3$beta_heart/pairs3$se_heart
  pairs3$z3<-pairs3$beta_brain/pairs3$se_brain
  
  Index<-which(abs(pairs3$z1)>1.96 | abs(pairs3$z2)>1.96 
               |abs(pairs3$z3)>1.96)
  
  pairs3_forcorr<-pairs3[-Index,26:28]
  
  # Correlation matrix
  corMatrix<-cor(pairs3_forcorr)
  
  #read into org file
  
  ab_org<-fread2("/path")
  heart_org<-fread2("/path")
  brain_org<-fread2("/path")
  
  # CPASSOC analysis
  ab_CPASSOC<-ab_org %>% 
    rename(
      rsid = variant_id,
      chr_ab = chromosome,
      beta_ab=beta,
      se_ab=standard_error,
      pval_ab=p_value,
      A1=effect_allele,
      A2=other_allele
    ) %>%
    mutate(id_ab=name_ab,N_ab=N_ab[N_ab$trait==name_ab,2]) %>% 
    dplyr::select(.,rsid,chr_ab,A1,A2,beta_ab,se_ab,pval_ab,id_ab,N_ab)
  
  heart_CPASSOC<-heart_org %>% 
    rename(
      rsid = SNP,
      chr_heart = CHR,
      beta_heart=BETA,
      se_heart=SE,
      N_heart=N,
      pval_heart=P
    ) %>%
    mutate(id_heart=name_heart) %>% 
    dplyr::select(.,rsid,chr_heart,A1,A2,beta_heart,se_heart,pval_heart,id_heart,N_heart)
  
  brain_CPASSOC<-brain_org %>% 
    rename(
      rsid = rsid,
      chr_brain = chr,
      beta_brain=beta,
      se_brain=se,
      pval_brain=p,
    ) %>%
    mutate(id_brain=name_brain,N_brain=N_brain[N_brain$traits==name_brain,1]) %>% 
    dplyr::select(.,rsid,chr_brain,A1,A2,beta_brain,se_brain,pval_brain,id_brain,N_brain)
  
  # merge
  merge_3 <- list(heart_CPASSOC,ab_CPASSOC,brain_CPASSOC)
  merge_3_match <- merge_3 %>% reduce(inner_join, by = c("rsid"))
  
  # rm palindromic variants
  var_logic<-c()
  var_logic<-apply(merge_3_match[,c(3:4,11:12,19:20)],1,
                   function(row){any(nchar(row)>1)})
  
  merge_3_match2<-merge_3_match[!var_logic,]
  
  # allel--A1/A2
  merge_3_match2$beta_ab <-ifelse((merge_3_match2$A1.x == merge_3_match2$A1.y & merge_3_match2$A2.x == merge_3_match2$A2.y), 
                                  merge_3_match2$beta_ab, merge_3_match2$beta_ab* (-1))
  merge_3_match2$beta_brain <- ifelse((merge_3_match2$A1.x == merge_3_match2$A1 & merge_3_match2$A2.x == merge_3_match2$A2), 
                                      merge_3_match2$beta_brain, merge_3_match2$beta_brain * (-1))
  
  # calculate Z1/Z2
  merge_3_match2<-merge_3_match2 %>% 
    mutate(.,z1=beta_ab/se_ab,
           z2=beta_heart/se_heart,
           z3=beta_brain/se_brain)
  
  # extract X matrix
  X<-merge_3_match2 %>% dplyr::select(.,z1,z2,z3)
  
  size_ab<-unique(pairs3$N_ab)
  size_heart<-median(pairs3$N_heart)
  size_brain<-unique(pairs3$N_brain)
  
  S<-c(size_ab,size_heart,size_brain)
  
  ### Shom calculation
  S_hom <- SHom(X = X, SampleSize = S, CorrMatrix = corMatrix)
  p_Shom <- pchisq(S_hom, df = 1, ncp = 0, lower.tail = F)
  
  sig_hom<-merge_3_match2[which(p_Shom<5*10^-8),]
  
  #find p
  p_index_hom <- which(p_Shom<5*10^-8)
  sig_hom2 <- cbind(sig_hom,p_index_hom)
  
  sig_final_hom_3<- sig_hom2[(sig_hom2$pval_heart<10^-3 & sig_hom2$pval_ab<10^-3 & sig_hom2$pval_brain<10^-3),]
  
  p_hom_3 <- as.data.frame(cbind(sig_final_hom_3$rsid,p_Shom[sig_final_hom_3$p_index_hom]))
  
  colnames(p_hom_3)<-c("rsid","p")
  
  ### Shet calculation
  ### the main important (longest) calculation:
  
  para <- EstimateGamma(N = 1E4, SampleSize = S, CorrMatrix = corMatrix,
                        correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = T)
  j <- SHet(X = X, SampleSize = S, CorrMatrix = corMatrix)
  p_Shet <- pgamma(q = j-para[3], shape = para[1], scale = para[2], lower.tail = F)
  
  sig_het<-merge_3_match2[which(p_Shet<5*10^-8),]
  
  #find p
  p_index <- which(p_Shet<5*10^-8)
  sig_het2 <- cbind(sig_het,p_index)
  
  sig_final_het_3 <- sig_het2[(sig_het2$pval_heart<10^-3 & sig_het2$pval_ab<10^-3 & sig_het2$pval_brain<10^-3),]
  
  p_het_3 <- as.data.frame(cbind(sig_final_het_3$rsid,p_Shet[sig_final_het_3$p_index]))
  
  colnames(p_het_3)<-c("rsid","p")
  
  # write files
  setwd("/path")
  
  write_xlsx(list(res_shom_3 = sig_final_hom_3,
                  res_shet_3 = sig_final_het_3,
                  p_hom_3=p_hom_3,
                  p_het_3=p_het_3),
             path = paste0(name_ab,"_",name_heart,"_",name_brain,".xlsx"))
}














