install.packages("haven")
library(haven)
library(survival)

aarp_data <- read_sas("../AARP-data/Chuck Matthews 26Jul2023/chuck_matthews_24jul2023.sas7bdat")
head(aarp_data)
dim(aarp_data)[1] # 163781
# P4: excluded 1660 individuals (1.0%) whose risk factor questionnaire was completed by a proxy 
# and 2157 individuals (1.3%) with incomplete physical activity data, 
# resulting in 159 937 former smokers (97.7%)

sum(is.na(aarp_data$FL_PROXY))
163781-sum(aarp_data$FL_FUQ)
table(aarp_data$FL_FUQ)
table(aarp_data$RF_SUBCOHORT)
sum(is.na(aarp_data$FL_PROXY))
table(aarp_data$RF_PHYS_LIGHT_CURR)

# loop for all the columns
col_names <- colnames(aarp_data)
freq_col <- list()
for (i in 1:length(col_names)){
  if(length(table(aarp_data[col_names[i]]))<25)
    print(table(aarp_data[col_names[i]]))
}
table(aarp_data[col_names[10]])
length(table(aarp_data[col_names[i]]))

mode(aarp_data[col_names[i]])

summary(aarp_data)
class(aarp_data[col_names[i]])
class(aarp_data$CANCER)

(aarp_data[col_names[10]])
freq_col[10]

# Table 4:
# age, sex, race, and ethnicity, educational level, perceived general health, 
# time since quitting and smoking intensity
# and individual lifestyle recommendation??
aarp_data$SEX
aarp_data$RACEI
aarp_data$EDUCM
aarp_data$HEALTH
aarp_data$ENTRY_AGE
aarp_data$EXIT_AGE
aarp_data$RF_ENTRY_AGE
aarp_data$FUQ_SMOKE_STOP
aarp_data$SMOKE_DOSE

# RF_PHYS_MODVIG_CURR 
# + ENTRY_AGE + SEX + RACEI + EDUCM + HEALTH
# + BMI_CUR + HEI2015_TOTAL_SCORE + MPED_A_BEV_NOFOOD
# + FUQ_SMOKE_STOP + SMOKE_DOSE
# 1. check NA
summary(aarp_data$BMI_CUR)
sum(table(aarp_data$BMI_CUR1))
trans_bmi <- function(x){
  if (is.na(x)) { 
    return(9)
  } else if (x>=30 |x<18.5) {
    return(0)
  } else if  (x>=25) {
    return(1)
  } else if  (x>=18.5) {
    return(2)
  } else {
    return(9)
  }
}
aarp_data$BMI_CUR1 <- unlist(lapply(aarp_data$BMI_CUR, trans_bmi))
q25 <- quantile(aarp_data$HEI2015_TOTAL_SCORE)[2]
q50 <- quantile(aarp_data$HEI2015_TOTAL_SCORE)[3]
q75 <- quantile(aarp_data$HEI2015_TOTAL_SCORE)[4]
trans_hei <- function(x){
  if (x<q25) { 
    return(0)
  } else if (x<q50) {
    return(1)
  } else if (x<q75) {
    return(2)
  } else {
    return(3)
  }
}
aarp_data$HEI2015_TOTAL_SCORE1 <- unlist(lapply(aarp_data$HEI2015_TOTAL_SCORE, trans_hei))
sum(table(aarp_data$HEI2015_TOTAL_SCORE1))
summary(aarp_data$FUQ_SMOKE_STOP)
trans_age <- function(x){
  if (x<55) { 
    return(0)
  } else if (x<=59) {
    return(1)
  } else if  (x<=64) {
    return(2)
  } else if  (x<=69) {
    return(3)
  } else {
    return(4)
  }
}
aarp_data$ENTRY_AGE1 <- unlist(lapply(aarp_data$ENTRY_AGE, trans_age))
summary(aarp_data$MPED_A_BEV_NOFOOD)
sum(table(aarp_data$BMI_CUR1))

summary(aarp_data$RF_PHYS_MODVIG_CURR)

trans_alcohol <- function(x,y){
  if (x==1 & y<=1) { 
    return(1)
  } else if (x==1 & y>1) {
    return(0)
  } else if (x==0 & y<=2) {
    return(1)
  } else {
    return(0)
  }
}
aarp_data$MPED_A_BEV_NOFOOD1 <- mapply(trans_alcohol, aarp_data$SEX, aarp_data$MPED_A_BEV_NOFOOD)
sum(table(aarp_data$MPED_A_BEV_NOFOOD1))

sum(table(aarp_data$SMOKE_DOSE))
summary(aarp_data$FUQ_SMOKE_STOP)
table(aarp_data$FUQ_SMOKE_STOP)



# body weight
aarp_data$HT_CUR
aarp_data$WT_CUR
aarp_data$BMI_CUR
# diet
# alcohol intake
aarp_data$HEI2015_TOTAL_SCORE #diet
aarp_data$MPED_A_BEV_NOFOOD

# treatment
aarp_data$RF_PHYS_LIGHT_CURR
aarp_data$RF_PHYS_MODVIG_CURR

sum(table(aarp_data$RF_PHYS_LIGHT_CURR))
sum(table(aarp_data$RF_PHYS_MODVIG_CURR))

# response
table(aarp_data$LUNG_MORT)
table(aarp_data$CANCER_MORT)
table(aarp_data$CANCER)
hist(aarp_data$PERSONYRS)
aarp_data$EXIT_AGE
aarp_data$PERSONYRS # Person time from baseline to exit for mortality


table(aarp_data$NDI_ICD10_RECODE_113)
table(aarp_data$NDI_ICD_VERSION)

sum(ifelse(aarp_data$NDI_ICD9_RECODE_72%in%c("510","520","530","540","550",
                                         "560","570","580"),1,0))

sum(ifelse(aarp_data$NDI_ICD10_RECODE_113%in%c("076","077","078","082","083",
                                             "084","085","086"),1,0))
trans_respiratory <- function(x,y){
  if (x%in%c("510","520","530","540","550",
             "560","570","580") | y%in%c("076","077","078","082","083",
                                         "084","085","086")) {
    return(1)
  } else {
    return(0)
  }
}
m1 <- ifelse(aarp_data$NDI_ICD9_RECODE_72%in%c("510","520","530","540","550",
                                         "560","570","580"),1,0)
m2 <- ifelse(aarp_data$NDI_ICD10_RECODE_113%in%c("076","077","078","082","083",
                                                 "084","085","086"),1,0)
aarp_data$respiratory_mort <- ifelse(m1+m2>=1,1,0)
 mapply(trans_respiratory, 
                                    aarp_data$NDI_ICD9_RECODE_72, 
                                     aarp_data$NDI_ICD10_RECODE_113)

sum(aarp_data$respiratory_mort)

# cox PH model
# try to repeat the all causes deaths first
# no cause death?
aarp_data$LUNG_MORT
aarp_data$LUNG_MORT1 <- ifelse(aarp_data$LUNG_MORT<1,0,1)
aarp_data$RF_PHYS_MODVIG_CURR1 <- ifelse(aarp_data$RF_PHYS_MODVIG_CURR<2,0,
                                        ifelse(aarp_data$RF_PHYS_MODVIG_CURR<=3,1,
                                               ifelse(aarp_data$RF_PHYS_MODVIG_CURR<=5,2)))
trans_pa <- function(x){
  if (x<=2) { 
    return(0)
  } else if (x<=3) {
    return(1)
  } else if  (x<=5) {
    return(2)
  } else {
    return(x)
  }
}

trans_pa(3)
aarp_data$RF_PHYS_MODVIG_CURR1 <- unlist(lapply(aarp_data$RF_PHYS_MODVIG_CURR, trans_pa))

sum(is.na(aarp_data$RF_PHYS_MODVIG_CURR))
# include all possible treatment values
fit.naive <- coxph(Surv(PERSONYRS, LUNG_MORT) ~ RF_PHYS_MODVIG_CURR 
                   + ENTRY_AGE + SEX + RACEI + EDUCM + HEALTH
                   + BMI_CUR + HEI2015_TOTAL_SCORE + MPED_A_BEV_NOFOOD
                   + FUQ_SMOKE_STOP + SMOKE_DOSE, data=aarp_data[aarp_data$A!=2,])
fit.naive <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ factor(RF_PHYS_MODVIG_CURR1) 
                   + factor(ENTRY_AGE1) + SEX + RACEI + EDUCM + HEALTH
                   + BMI_CUR + HEI2015_TOTAL_SCORE + MPED_A_BEV_NOFOOD
                   + FUQ_SMOKE_STOP + SMOKE_DOSE, data=aarp_data[aarp_data$RF_PHYS_MODVIG_CURR1!=9,])
fit.naive <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ factor(RF_PHYS_MODVIG_CURR1) 
                   + factor(ENTRY_AGE1) + SEX + RACEI + EDUCM + HEALTH
                   + BMI_CUR + HEI2015_TOTAL_SCORE + MPED_A_BEV_NOFOOD
                   + factor(SMOKE_DOSE), data=aarp_data[aarp_data$RF_PHYS_MODVIG_CURR1!=9,])

fit.naive <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ factor(RF_PHYS_MODVIG_CURR1) 
                   + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                   + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                   + factor(SMOKE_DOSE), data=aarp_data[aarp_data$RF_PHYS_MODVIG_CURR1!=9,])

fit.naive1 <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ factor(RF_PHYS_MODVIG_CURR1) 
                   + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                   + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                   + factor(FUQ_SMOKE_STOP) + factor(SMOKE_DOSE), data=data2)

fit.naive2 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ factor(RF_PHYS_MODVIG_CURR1) 
                    + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                    + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                    + factor(FUQ_SMOKE_STOP) + factor(SMOKE_DOSE), data=data2)
fit.naive2 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ factor(RF_PHYS_MODVIG_CURR1) 
                    + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                    + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                    + factor(SMOKE_DOSE), data=data1)

data1 <- aarp_data[aarp_data$RF_PHYS_MODVIG_CURR1!=9,]
data2 <- data1[!is.na(data1$FUQ_SMOKE_STOP),]

# naive causal hazard ratio
# include all possible treatment values
fit.naive <- coxph(Surv(Y, D) ~ A + sex + age + health
                   + years_quit + cpd, data=data)
exp(fit.naive$coefficients[2])
# 1 vs 0 - 0.88
# 2 vs 0 - 0.83
fit.naive.1 <- coxph(Surv(Y, D) ~ A + sex + age + health
                     + years_quit + cpd, data = data[data$A!=2,])
exp(fit.naive.1$coefficients[1:2])
fit.naive.2 <- coxph(Surv(Y, D) ~ A + sex + age + health
                     + years_quit + cpd, data = data[data$A!=1,])
exp(fit.naive.2$coefficients[1:2])

  