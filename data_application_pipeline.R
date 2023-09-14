rm(list=ls())

# install packages
library(survival)
library(survSens)
library(simex)

######################
# data pre processing
#####################
# ENTRY_AGE
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
# BMI_CUR
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
# HEI2015_TOTAL_SCORE
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
# MPED_A_BEV_NOFOOD
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
# RF_PHYS_MODVIG_CURR
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
aarp_data$RF_PHYS_MODVIG_CURR1 <- unlist(lapply(aarp_data$RF_PHYS_MODVIG_CURR, trans_pa))
# LUNG_MORT
aarp_data$LUNG_MORT1 <- ifelse(aarp_data$LUNG_MORT<1,0,1)
# all_mort
code10_all <- paste("0", as.character(seq(19,43,1)), sep="")
m1 <- ifelse(aarp_data$NDI_ICD9_RECODE_72%in%c("150","160","170","180","190",
                                               "200","210","220","230",'240'),1,0)
m2 <- ifelse(aarp_data$NDI_ICD10_RECODE_113%in%code10_all,1,0)
aarp_data$all_mort <- ifelse(m1+m2>=1,1,0)

# respiratory_mort
m1 <- ifelse(aarp_data$NDI_ICD9_RECODE_72%in%c("510","520","530","540","550",
                                               "560","570","580"),1,0)
m2 <- ifelse(aarp_data$NDI_ICD10_RECODE_113%in%c("076","077","078","082","083",
                                                 "084","085","086"),1,0)
aarp_data$respiratory_mort <- ifelse(m1+m2>=1,1,0)

######################################################
## STEP 0: Obtain the naive estimate in the Cox model
#####################################################
table(aarp_data$RF_PHYS_MODVIG_CURR1)
# 0     1     2     9 
# 39524 40481 81522  2254 

# naive model - lung_mort
data1 <- aarp_data[aarp_data$RF_PHYS_MODVIG_CURR1!=9,]
data2 <- data1[!is.na(data1$FUQ_SMOKE_STOP),]

fit.naive.lungcancer.1 <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ factor(RF_PHYS_MODVIG_CURR1) 
                   + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                   + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                   + factor(SMOKE_DOSE), data=data1)

fit.naive.lungcancer.2 <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ factor(RF_PHYS_MODVIG_CURR1) 
                    + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                    + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                    + factor(FUQ_SMOKE_STOP) + factor(SMOKE_DOSE), data=data2)

# naive model - all_mort
3657 + 8807 + 12403
sum(aarp_data$all_mort)
fit.naive.allcancer.1 <- coxph(Surv(PERSONYRS, all_mort) ~ factor(RF_PHYS_MODVIG_CURR1) 
                                + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                                + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                                + factor(SMOKE_DOSE), data=data1)

fit.naive.allcancer.2 <- coxph(Surv(PERSONYRS, all_mort) ~ factor(RF_PHYS_MODVIG_CURR1) 
                                + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                                + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                                + factor(FUQ_SMOKE_STOP) + factor(SMOKE_DOSE), data=data2)

# naive model - respiratory_mort
1465 + 2634 + 2972
sum(aarp_data$respiratory_mort)
fit.naive.RD.1 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ factor(RF_PHYS_MODVIG_CURR1) 
                        + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                        + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                        + factor(SMOKE_DOSE), data=data1)

fit.naive.RD.2 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ factor(RF_PHYS_MODVIG_CURR1) 
                    + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                    + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                    + factor(FUQ_SMOKE_STOP) + factor(SMOKE_DOSE), data=data2)


######################################################
## STEP 1: Suppose that there is an unobserved confounder
## keep in mind that we want to know at what level of alpha and beta 
#####################################################
# use fit.naive.RD.2 as a start point
# define the range of alpha and beta (the key point is how to ??)
# stochastic EM with regression
data3 <- data2[data2$RF_PHYS_MODVIG_CURR1%in%c(0,1),]
X <- cbind(factor(data3$ENTRY_AGE1), factor(data3$SEX), factor(data3$RACEI), factor(data3$EDUCM), factor(data3$HEALTH),
           factor(data3$BMI_CUR1), factor(data3$HEI2015_TOTAL_SCORE1), factor(data3$MPED_A_BEV_NOFOOD1),
           factor(data3$FUQ_SMOKE_STOP), factor(data3$SMOKE_DOSE))
beta <- 0
alpha <- 0
fit.naive.RD.3 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ factor(RF_PHYS_MODVIG_CURR1) 
                        + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI) + factor(EDUCM) + factor(HEALTH)
                        + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                        + factor(FUQ_SMOKE_STOP) + factor(SMOKE_DOSE), data=data3)
fit.stoEM_reg <- survSensitivity(data3$PERSONYRS, data3$respiratory_mort, 
                                 data3$RF_PHYS_MODVIG_CURR1, X, "stoEM_reg",
                                 zetaT = beta, zetaZ = alpha,
                                 B = 20)
theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se

######################################################
## STEP 2: Suppose that there is misclassification error in SMOKE_DOSE
#####################################################

# the key is what should be a proper misclassification matrix?
# mcsimex
data3$SMOKE_DOSE_f <- factor(data3$SMOKE_DOSE)
cols <- c("RF_PHYS_MODVIG_CURR1", "ENTRY_AGE1", "SEX", "RACEI", "EDUCM", "HEALTH",
          "BMI_CUR1", "HEI2015_TOTAL_SCORE1", "MPED_A_BEV_NOFOOD1",
          "FUQ_SMOKE_STOP", "SMOKE_DOSE")
data3[cols] <- lapply(data3[cols], factor)
p_ij <- matrix(c(0.8, 0.1, 0.01, 0.01, 0.01, 0.01,
                 0.16, 0.8, 0.1, 0.01, 0.01, 0.01,
                 0.01, 0.07, 0.8, 0.1, 0.01, 0.01,
                 0.01, 0.01, 0.07, 0.8, 0.07, 0.01,
                 0.01, 0.01, 0.01, 0.07, 0.8, 0.16,
                 0.01, 0.01, 0.01, 0.01, 0.1, 0.8), nrow = 6, byrow = TRUE)
dimnames(p_ij) <- list(levels(data3$SMOKE_DOSE), levels(data3$SMOKE_DOSE))
fit.naive.RD.3.t <- coxph(Surv(PERSONYRS, respiratory_mort) ~ RF_PHYS_MODVIG_CURR1
                        + ENTRY_AGE1 + SEX +RACEI + EDUCM + HEALTH
                        + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                        + FUQ_SMOKE_STOP + SMOKE_DOSE, data=data3, model=TRUE)
#####################################
# this one took really a long time to run
####################################
fit.mcsimex.1 <- mcsimex(fit.naive.RD.3.t, mc.matrix = list(SMOKE_DOSE=p_ij), 
                         SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE)
theta.hat.mcsimex.1 <- fit.mcsimex.1$coefficients[1]
theta.hat.mcsimex.se.1 <- sqrt(fit.mcsimex.1$variance.jackknife[1,1])

###########################################
# the above two steps can be done separated
###########################################

######################################################
## STEP 3: Suppose that there is an unobserved confounder as well as a misclassification error in SMOKE_DOSE
#####################################################
# adjusted both
data3$Usim <- SimulateU_surv(data3$PERSONYRS, data3$respiratory_mort, 
                                as.numeric(data3$RF_PHYS_MODVIG_CURR1)-1, X, 
                                zetat = beta, zetaz = alpha, theta = 0.5, offset = TRUE)$U 
# SimulateU_surv include iter=20 but not B=5
naive.model.U <- coxph(Surv(PERSONYRS, respiratory_mort) ~ RF_PHYS_MODVIG_CURR1
                       + ENTRY_AGE1 + SEX +RACEI + EDUCM + HEALTH
                       + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                       + FUQ_SMOKE_STOP + SMOKE_DOSE + Usim, data=data3, model=TRUE) # simex cannot take offset

fit.mcsimex.1.U <- mcsimex(naive.model.U, mc.matrix = list(SMOKE_DOSE=p_ij), 
                         SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE)
theta.hat.simex.U.3 <- fit.mcsimex.1.U$coefficients[1]
theta.hat.simex.se.U.3 <- sqrt(fit.mcsimex.1.U$variance.jackknife[1,1])
