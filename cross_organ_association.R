library("dplyr")
library("purrr")
library("tidyverse")
library("bigreadr")
library("writexl")
library("stringr")
library("readxl")


assoc_3<-function(pheno_ab,pheno_heart,pheno_brain,all_data,write_file){
  
  res_heart_liver1<-as.data.frame(matrix(NA,ncol=6,nrow=length(pheno_ab)*length(pheno_heart)))
  res_heart_brain1<-as.data.frame(matrix(NA,ncol=6,nrow=length(pheno_heart)*length(pheno_brain)))
  res_liver_brain1<-as.data.frame(matrix(NA,ncol=6,nrow=length(pheno_ab)*length(pheno_brain)))
  
  res_heart_liver2<-as.data.frame(matrix(NA,ncol=6,nrow=length(pheno_ab)*length(pheno_heart)))
  res_heart_brain2<-as.data.frame(matrix(NA,ncol=6,nrow=length(pheno_heart)*length(pheno_brain)))
  res_liver_brain2<-as.data.frame(matrix(NA,ncol=6,nrow=length(pheno_ab)*length(pheno_brain)))
  
  res_heart_liver3<-as.data.frame(matrix(NA,ncol=6,nrow=length(pheno_ab)*length(pheno_heart)))
  res_heart_brain3<-as.data.frame(matrix(NA,ncol=6,nrow=length(pheno_heart)*length(pheno_brain)))
  res_liver_brain3<-as.data.frame(matrix(NA,ncol=6,nrow=length(pheno_ab)*length(pheno_brain)))
  
  
  row1<-1
  for(i in pheno_heart){
    for(j in pheno_ab){
      #model1
      formula1 <- formula(paste(i,"~",j, "+ age + sex + ethnic + TDI + edu + site"))
      model1 <- lm(formula1, data = all_data)
      res_heart_liver1[row1,1:2]<-c(i,j)
      res_heart_liver1[row1,3:6]<-summary(model1)$coefficients[2,]
      
      #model2
      formula2 <- formula(paste(i,"~",j, "+ age + sex + ethnic + TDI + edu + smoke + alcohol + METmin + BMI + site"))
      model2 <- lm(formula2, data = all_data)
      res_heart_liver2[row1,1:2]<-c(i,j)
      res_heart_liver2[row1,3:6]<-summary(model2)$coefficients[2,]
      
      #model3
      formula3 <- formula(paste(i,"~",j, "+ age + sex + ethnic + TDI + edu + smoke + alcohol + METmin + BMI + site + his_hyplip + his_CVD + his_t2d + his_hypten + his_dep + bodysurface"))
      model3 <- lm(formula3, data = all_data)
      res_heart_liver3[row1,1:2]<-c(i,j)
      res_heart_liver3[row1,3:6]<-summary(model3)$coefficients[2,]
      
      row1<-row1+1
    }
  }
  
  
  row2<-1
  for(i in pheno_brain){
    for(j in pheno_heart){
      #model1
      formula1 <- formula(paste(i,"~",j, "+ age + sex + ethnic + TDI + edu + site + scan_x_de + scan_y_de + scan_z_de + scan_table_de"))
      model1 <- lm(formula1, data = all_data)
      res_heart_brain1[row2,1:2]<-c(i,j)
      res_heart_brain1[row2,3:6]<-summary(model1)$coefficients[2,]
      
      #model2
      formula2 <- formula(paste(i,"~",j, "+ age + sex + ethnic + TDI + edu + smoke + alcohol + METmin + BMI + site + scan_x_de + scan_y_de + scan_z_de + scan_table_de"))
      model2 <- lm(formula2, data = all_data)
      res_heart_brain2[row2,1:2]<-c(i,j)
      res_heart_brain2[row2,3:6]<-summary(model2)$coefficients[2,]
      
      #model3
      formula3 <- formula(paste(i,"~",j, "+ age + sex + ethnic + TDI + edu + smoke + alcohol + METmin + BMI + site + his_hyplip + his_CVD + his_t2d + his_hypten + his_dep + scan_x_de + scan_y_de + scan_z_de + scan_table_de + bodysurface"))
      model3 <- lm(formula3, data = all_data)
      res_heart_brain3[row2,1:2]<-c(i,j)
      res_heart_brain3[row2,3:6]<-summary(model3)$coefficients[2,]
      
      row2<-row2+1
    }
  }
  
  row3<-1
  for(i in pheno_brain){
    for(j in pheno_ab){
      #model1
      formula1 <- formula(paste(i,"~",j, "+ age + sex + ethnic + TDI + edu + site + scan_x + scan_y + scan_z + scan_table"))
      model1 <- lm(formula1, data = all_data)
      res_liver_brain1[row3,1:2]<-c(i,j)
      res_liver_brain1[row3,3:6]<-summary(model1)$coefficients[2,]
      
      #model2
      formula2 <- formula(paste(i,"~",j, "+ age + sex + ethnic + TDI + edu + smoke + alcohol + METmin + BMI + site + scan_x_de + scan_y_de + scan_z_de + scan_table_de"))
      model2 <- lm(formula2, data = all_data)
      res_liver_brain2[row3,1:2]<-c(i,j)
      res_liver_brain2[row3,3:6]<-summary(model2)$coefficients[2,]
      
      #model3
      formula3 <- formula(paste(i,"~",j, "+ age + sex + ethnic + TDI + edu + smoke + alcohol + METmin + BMI + site + his_hyplip + his_CVD + his_t2d + his_hypten + his_dep + scan_x_de + scan_y_de + scan_z_de + scan_table_de"))
      model3 <- lm(formula3, data = all_data)
      res_liver_brain3[row3,1:2]<-c(i,j)
      res_liver_brain3[row3,3:6]<-summary(model3)$coefficients[2,]
      
      
      row3<-row3+1
    }
  }
  
  colnames(res_heart_liver1)<-c("phenotype1","phenotype2","Estimate","Std. Error","t value","p value")
  colnames(res_heart_brain1)<-c("phenotype1","phenotype2","Estimate","Std. Error","t value","p value")
  colnames(res_liver_brain1)<-c("phenotype1","phenotype2","Estimate","Std. Error","t value","p value")
  
  colnames(res_heart_liver2)<-c("phenotype1","phenotype2","Estimate","Std. Error","t value","p value")
  colnames(res_heart_brain2)<-c("phenotype1","phenotype2","Estimate","Std. Error","t value","p value")
  colnames(res_liver_brain2)<-c("phenotype1","phenotype2","Estimate","Std. Error","t value","p value")
  
  
  colnames(res_heart_liver3)<-c("phenotype1","phenotype2","Estimate","Std. Error","t value","p value")
  colnames(res_heart_brain3)<-c("phenotype1","phenotype2","Estimate","Std. Error","t value","p value")
  colnames(res_liver_brain3)<-c("phenotype1","phenotype2","Estimate","Std. Error","t value","p value")
  
  setwd("/path")
  write_xlsx(list(res_heart_liver3 = res_heart_liver3,
                  res_heart_brain3 = res_heart_brain3,
                  res_liver_brain3 = res_liver_brain3
  ), path = write_file)
  
}

assoc_3(sel_ab_cmd,sel_heart_cmd,sel_brain_cmd,all_data,"res.xlsx")