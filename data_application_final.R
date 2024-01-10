# rm(list=ls())

# install packages
library(survival)
library(survSens)
library(simex)
library(ggplot2)
library(haven)
library(gridExtra)
library(tidyverse)

aarp_data1 <- read_sas("../AARP-data/Chuck Matthews 11Oct2023/chuck_matthews_09oct2023.sas7bdat")
head(aarp_data1)
dim(aarp_data1)[1]

##################
# data exclusion
##################
table(aarp_data1$FL_RF_PROXY) # 1 - 1660
aarp_data2 <- aarp_data1[aarp_data1$FL_RF_PROXY !=1,]
table(aarp_data2$RF_PHYS_MODVIG_CURR) # 9 - 2157
aarp_data <- aarp_data2[aarp_data2$RF_PHYS_MODVIG_CURR !=9,]
dim(aarp_data)[1] 
# TODO: didn't match with 159964

# 159937/163781
# 163781-1660-2157
# 163781-1660-159937
# 159964/163781

#####################
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

# RACEI
trans_race <- function(x){
  if (x%in%c(4,5,6)) { 
    return(4)
  } else {
    return(x)
  }
}
aarp_data$RACEI1 <- unlist(lapply(aarp_data$RACEI, trans_race))

# EDUCM
trans_educ <- function(x){
  if (x<=2) { 
    return(1)
  } else if (x==3) {
    return(2)
  } else if  (x==4) {
    return(3)
  } else if  (x==5) {
    return(4)
  } else {
    return(x)
  }
}
aarp_data$EDUCM1 <- unlist(lapply(aarp_data$EDUCM, trans_educ))

# HEALTH
trans_health <- function(x){
  if (x<=2) { 
    return(1)
  } else if (x<=3) {
    return(2)
  } else if  (x<=5) {
    return(3)
  } else {
    return(4)
  }
}
aarp_data$HEALTH1 <- unlist(lapply(aarp_data$HEALTH, trans_health))

# smoking behavior
# SMOKE_DOSE
summary(aarp_data$SMOKE_DOSE) # no NA
table(aarp_data$SMOKE_DOSE)
# 1     2     3     4     5     6 
# 40439 49256 31351 21025 13907  3986
sum(table(aarp_data$SMOKE_DOSE))

# SMOKE_QUIT
# transfer smoke_quit_detailed to smoke_quit
summary(aarp_data$SMOKE_QUIT_DETAILED) # no NA
table(aarp_data$SMOKE_QUIT_DETAILED)
trans_smokequit <- function(x){
  if (x == 0) { 
    return(0)
  } else if (x>=1 & x<=6) {
    return(1)
  } else if (x>=7 & x<=12) {
    return(2)
  } else if (x>=13 & x<=18) {
    return(3)
  } else if (x>=19 & x<=24) {
    return(4)
  } else if (x>=25 & x<=30) {
    return(5)
  } else {
    return(9)
  }
}
aarp_data$SMOKE_QUIT1 <- unlist(lapply(aarp_data$SMOKE_QUIT_DETAILED, trans_smokequit))
sum(table(aarp_data$SMOKE_QUIT1))

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
# q25 <- 62.5
# q50 <- 69.3
# q75 <- 75.2


trans_hei <- function(x){
  if (x<=q25) { 
    return(0)
  } else if (x<=q50) {
    return(1)
  } else if (x<=q75) {
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
  if (x<2) { 
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
# binned 0 and 1
aarp_data$RF_PHYS_MODVIG_CURR2 <- ifelse(aarp_data$RF_PHYS_MODVIG_CURR1<=1,0,1)
########################################
# binning 2
# 0-2 <1hr; 3-5 >1hr
########################################
# RF_PHYS_MODVIG_CURR
# > table(aarp_data$RF_PHYS_MODVIG_CURR)
# 
# 0     1     2     3     4     5 
# 4875 17651 16495 40118 41965 38860 
trans_pa_b2 <- function(x){
  if (x<3) { 
    return(0)
  } else {
    return(1)
  }
}
aarp_data$RF_PHYS_MODVIG_CURR3 <- unlist(lapply(aarp_data$RF_PHYS_MODVIG_CURR, trans_pa_b2))
table(aarp_data$RF_PHYS_MODVIG_CURR3)

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

# use functions in data_application_pipeline.R to pre-process the variables
# table 1
table(aarp_data$BMI_CUR1)
summary(aarp_data$BMI_CUR)
# TODO: 2728 NA
table(aarp_data$RF_PHYS_MODVIG_CURR1)
# 0     1     2 
# 22526 56613 80825
table(aarp_data$HEI2015_TOTAL_SCORE1)
# 0     1     2     3 
# 39991 39991 39991 39991
table(aarp_data$MPED_A_BEV_NOFOOD1)
# 0      1 
# 27897 132067

# generate BMI
summary(aarp_data$HT_CUR)
summary(aarp_data$WT_CUR)
gen_bmi <- function(x, y){
  if (is.na(x) | (is.na(y))) {
    return(NA)
  } else {
    return(y/(x^2))
  }
}
aarp_data$bmi <- mapply(gen_bmi, aarp_data$HT_CUR, aarp_data$WT_CUR)
summary(aarp_data$bmi)

# table 2
# The distribution of all the other baseline confounders looks reasonable.
# Reasonable means the percentages falls in the range shown in Table 2.
table(aarp_data$SEX)/159964
table(aarp_data$ENTRY_AGE1)/159964
table(aarp_data$RACEI1)/159964
table(aarp_data$EDUCM1)/159964
table(aarp_data$HEALTH1)/159964
table(aarp_data$SMOKE_DOSE)/159964
table(aarp_data$SMOKE_QUIT1)/159964

# with a mean (SD) follow-up of 18.9 (6.3) years
summary(aarp_data$PERSONYRS)
mean(aarp_data$PERSONYRS) # 18.97
sd(aarp_data$PERSONYRS) # 6.24

aarp_data$NDI_ICD9_RECODE_72
##########################################

#########################################
# Descriptive statistics for result part
########################################
mean(aarp_data$ENTRY_AGE)
sd(aarp_data$ENTRY_AGE)
#####################
# Generate Table 1
####################
sum(table(aarp_data$RF_PHYS_MODVIG_CURR2))
PA2_0_cnt <- table(aarp_data$RF_PHYS_MODVIG_CURR2)[1]
PA2_1_cnt <- table(aarp_data$RF_PHYS_MODVIG_CURR2)[2]
total_cnt <- sum(table(aarp_data$RF_PHYS_MODVIG_CURR2))
results <- aarp_data.1 %>%
  group_by(RF_PHYS_MODVIG_CURR2) %>%
  summarise_at(c("ENTRY_AGE1", "SEX"), ~n())

aarp_data_long <- aarp_data %>%
  select(cols) %>%
  pivot_longer(cols = -"RF_PHYS_MODVIG_CURR2", names_to = "variable", values_to = "value")

counts <- aarp_data_long %>%
  group_by(RF_PHYS_MODVIG_CURR2, variable, value) %>%
  summarise(count = n(), .groups = 'drop') 

table1_pa0 <- counts %>%
  filter(RF_PHYS_MODVIG_CURR2==0) %>%
  mutate(prop=round(count/PA2_0_cnt*100, 1))

table1_pa1 <- counts %>%
  filter(RF_PHYS_MODVIG_CURR2==1) %>%
  mutate(prop=round(count/PA2_1_cnt*100, 1))

table1 <- left_join(table1_pa0, table1_pa1, by = c("variable", "value"))

table1 <- table1 %>% 
  mutate(total = count.x + count.y) %>%
  mutate(prop.total = round(total/total_cnt*100,1))

# paste point estimate and confidence interval
paste.output <- function(x, y){
  return(paste(as.character(x), " (", as.character(y), ")", sep=""))
}

table1$output0 <- mapply(paste.output, table1$count.x, table1$prop.x)
table1$output1 <- mapply(paste.output, table1$count.y, table1$prop.y)
table1$output.total <- mapply(paste.output, table1$total, table1$prop.total)

write.csv(table1, file="table1.csv")
##########################
# naive estimator
##########################
fit.naive.RD.1 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ factor(RF_PHYS_MODVIG_CURR2) 
                        + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
                        + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                        + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data)

fit.naive.all.1 <- coxph(Surv(PERSONYRS, all_mort) ~ factor(RF_PHYS_MODVIG_CURR2) 
                         + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
                         + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                         + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data)

fit.naive.lung.1 <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ factor(RF_PHYS_MODVIG_CURR2) 
                          + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
                          + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                          + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data)

est <- exp(fit.naive.RD.1$coefficients[1])
se <- sqrt(fit.naive.RD.1$var[1,1])
exp.se <- est*se
c(est-qnorm(1-0.05/2,0,1)*exp.se, est+qnorm(1-0.05/2,0,1)*exp.se)

est <- exp(fit.naive.all.1$coefficients[1])
se <- sqrt(fit.naive.all.1$var[1,1])
exp.se <- est*se
c(est-qnorm(1-0.05/2,0,1)*exp.se, est+qnorm(1-0.05/2,0,1)*exp.se)

est <- exp(fit.naive.lung.1$coefficients[1])
se <- sqrt(fit.naive.lung.1$var[1,1])
exp.se <- est*se
c(est-qnorm(1-0.05/2,0,1)*exp.se, est+qnorm(1-0.05/2,0,1)*exp.se)

##########################
# unobserved confounding
#########################
X <- cbind(factor(aarp_data$ENTRY_AGE1), factor(aarp_data$SEX), factor(aarp_data$RACEI1), factor(aarp_data$EDUCM1), factor(aarp_data$HEALTH1),
           factor(aarp_data$BMI_CUR1), factor(aarp_data$HEI2015_TOTAL_SCORE1), factor(aarp_data$MPED_A_BEV_NOFOOD1),
           factor(aarp_data$SMOKE_QUIT1), factor(aarp_data$SMOKE_DOSE))

# RD + binning 1
# 24 in total around 2hrs
beta <- c(0, 0.5, 1, 1.5, 2, 2.5)
alpha <- c(-1.5, -1, -0.5, 0)
system.time(fit.stoEM_reg.b1.RD <- survSensitivity(aarp_data$PERSONYRS, aarp_data$respiratory_mort, 
                                                   aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                   zetaT = beta, zetaZ = alpha,
                                                   B = 5))
save(fit.stoEM_reg.b1.RD, file="fit.stoEM_reg.b1.RD.RData")



# contour plot
rst4unobs.RD <- fit.stoEM_reg.b1.RD
rst4unobs.RD$tau1 <- exp(rst4unobs.RD$tau1)
# p <- plotsens(rst4unobs.RD, coeff0 = round(exp(fit.naive.RD.1$coefficients[1]),4))
# p1 <- p + ggtitle("Outcome: respiratory disease") + theme(legend.position="bottom")
# save(rst4unobs.RD, file="rst4unobs.RD.RData")

# all + binning 1
system.time(fit.stoEM_reg.b1.all_mort <- survSensitivity(aarp_data$PERSONYRS, 
                                                           aarp_data$all_mort, 
                                                           aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                           zetaT = beta, zetaZ = alpha,
                                                           B = 5))
save(fit.stoEM_reg.b1.all_mort, file="fit.stoEM_reg.b1.all_mort.RData")

# contour plot
rst4unobs.all <- fit.stoEM_reg.b1.all_mort
rst4unobs.all$tau1 <- exp(rst4unobs.all$tau1)
# p <- plotsens(rst4unobs.all, coeff0 = round(exp(fit.naive.all.1$coefficients[1]),4))
# p2 <- p + ggtitle("Outcome: all cancer mortality") + theme(legend.position="bottom")
# save(rst4unobs.all, file="rst4unobs.all.RData")

# lung + binning 1
system.time(fit.stoEM_reg.b1.lung_mort <- survSensitivity(aarp_data$PERSONYRS, 
                                                            aarp_data$LUNG_MORT1, 
                                                            aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                            zetaT = beta, zetaZ = alpha,
                                                            B = 5))
save(fit.stoEM_reg.b1.lung_mort, file="fit.stoEM_reg.b1.lung_mort.RData")

# contour plot
rst4unobs.lung <- fit.stoEM_reg.b1.lung_mort
rst4unobs.lung$tau1 <- exp(rst4unobs.lung$tau1)
# p <- plotsens(rst4unobs.lung, coeff0 = round(exp(fit.naive.lung.1$coefficients[1]),4))
# p3 <- p + ggtitle("Outcome: lung cancer mortality") + theme(legend.position="bottom")
# save(rst4unobs.lung, file="rst4unobs.lung.RData")
# 
# 
# grid.arrange(p1, p2, p3, ncol = 3)

head(rst4unobs.RD)
# 1. exp(tau1)
# 2. transform se
# 3. ll and ul
# 4. confidence interval and output

# add outcome category
# put contour plot figures together
rst4unobs.RD$outcome <- "Respiratory Disease"
rst4unobs.lung$outcome <- "Lung Cancer"
rst4unobs.all$outcome <- "Cancer"

rst4unobs.plot <- bind_rows(rst4unobs.RD, 
                            rst4unobs.lung,
                            rst4unobs.all)

rst4unobs.plot$outcome <- factor(rst4unobs.plot$outcome,
                                 levels = c("Respiratory Disease", 
                                            "Lung Cancer",
                                            "Cancer"))

rst4unobs.plot %>%
  ggplot(aes(x = zetaz, y = zetat1)) +
  stat_contour(aes(z = tau1)) + 
  stat_contour(aes(z = t), colour = "red", 
               linetype="dotdash", breaks = c(-1.96,1.96)) +
  annotate("text", x = 0, y = 0, label = 0.67) +
  metR::geom_text_contour(aes(z = tau1), skip = 0, min.size = 0) +
  facet_grid(~ outcome, scales = "free") +
  labs(x = "association between U and A", y = "association between U and T") +
  theme_bw()

# ggsave(file="Figures/Fig20_coverage_Sim1.eps", width = 290,
#        height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig2a.plot) #saves g


# transform estimate and standard error
rst4unobs.RD$exp.se <- rst4unobs.RD$tau1*rst4unobs.RD$tau1.se
rst4unobs.lung$exp.se <- rst4unobs.lung$tau1*rst4unobs.lung$tau1.se
rst4unobs.all$exp.se <- rst4unobs.all$tau1*rst4unobs.all$tau1.se

rst4unobs.lung$ll <- rst4unobs.lung$tau1-qnorm(1-0.05/2,0,1)*rst4unobs.lung$exp.se
rst4unobs.lung$ul <- rst4unobs.lung$tau1+qnorm(1-0.05/2,0,1)*rst4unobs.lung$exp.se
rst4unobs.all$ll <- rst4unobs.all$tau1-qnorm(1-0.05/2,0,1)*rst4unobs.all$exp.se
rst4unobs.all$ul <- rst4unobs.all$tau1+qnorm(1-0.05/2,0,1)*rst4unobs.all$exp.se
rst4unobs.RD$ll <- rst4unobs.RD$tau1-qnorm(1-0.05/2,0,1)*rst4unobs.RD$exp.se
rst4unobs.RD$ul <- rst4unobs.RD$tau1+qnorm(1-0.05/2,0,1)*rst4unobs.RD$exp.se

# paste lower bound and upper bound
paste.ci <- function(x, y){
  ci <- paste(as.character(round(x,2)),
              as.character(round(y,2)),
              sep="-")
  return(paste("(", ci, ")", sep=""))
}

rst4unobs.RD$CI <- mapply(paste.ci, rst4unobs.RD$ll, rst4unobs.RD$ul)
rst4unobs.lung$CI <- mapply(paste.ci, rst4unobs.lung$ll, rst4unobs.lung$ul)
rst4unobs.all$CI <- mapply(paste.ci, rst4unobs.all$ll, rst4unobs.all$ul)

# paste point estimate and confidence interval
paste.rst <- function(x, y){
  return(paste(as.character(round(x,2)), y, sep=" "))
}

rst4unobs.RD$output <- mapply(paste.rst, rst4unobs.RD$tau1, rst4unobs.RD$CI)
rst4unobs.lung$output <- mapply(paste.rst, rst4unobs.lung$tau1, rst4unobs.lung$CI)
rst4unobs.all$output <- mapply(paste.rst, rst4unobs.all$tau1, rst4unobs.all$CI)

# organize tables for Table 3
rst4unobs.RD.24 <- rst4unobs.RD %>% 
    select(c("zetat1", "zetaz", "output")) %>% 
    spread(zetaz, value = output) %>%
    select(c(1,5,4,3,2))
write.csv(rst4unobs.RD.24, file="rst4unobs.RD.24.csv")

rst4unobs.lung.24 <- rst4unobs.lung %>% 
    select(c("zetat1", "zetaz", "output")) %>% 
    spread(zetaz, value = output) %>%
    select(c(1,5,4,3,2))
write.csv(rst4unobs.lung.24, file="rst4unobs.lung.24.csv")

rst4unobs.all.24 <- rst4unobs.all %>% 
  select(c("zetat1", "zetaz", "output")) %>% 
  spread(zetaz, value = output) %>%
  select(c(1,5,4,3,2))
write.csv(rst4unobs.all.24, file="rst4unobs.all.24.csv")

############################
# check the range of the coefficients of the observed covariates
###########################
coef_YX.RD <- as.data.frame(summary(fit.naive.RD.1)$coefficients[,1:3])
summary(fit.naive.RD.1)$coefficients
coef_YX.RD %>% 
  arrange(desc(abs(coef)))
# beta can be as large as 2 (but the first couples are age and health status)
# then follow by smoking quit and intensity - which are around 1.27
# it depends on whether in literature we can find other risk factors leading to RD that are stronger than 
# smoking

coef_YX.all <- as.data.frame(summary(fit.naive.all.1)$coefficients[,1:3])
coef_YX.all %>% 
  arrange(desc(abs(coef)))
# beta can be as large as 1.35 (but the first couples are age and health status)
# then follow by smoking quit and intensity - which are around 0.55
# but 0.1 can be a very reasonable value for coefficient b/c observed covariates such as race, sex and diet even educ have coef around 0.1.

coef_YX.lung <- as.data.frame(summary(fit.naive.lung.1)$coefficients[,1:3])
coef_YX.lung %>% 
  arrange(desc(abs(coef)))
# beta can be as large as 1.66 (but the first is age)
# then follow by smoking quit and intensity - which are around 1.4
# but 0.1 can be a very reasonable value for coefficient b/c observed covariates such as race, sex and diet even educ have coef around 0.1.

fit.PA.1 <- glm(factor(RF_PHYS_MODVIG_CURR2) ~ factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
    + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
    + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data, family = "binomial")
summary.fit.PA.1 <- summary(fit.PA.1)$coefficients
coef_AX <- as.data.frame(summary.fit.PA.1[, c("Estimate", "Std. Error")])
coef_AX %>% 
  arrange(desc(abs(Estimate)))

coef_AX %>%
  rename("Est" = "Estimate",
         "Se" = "Std. Error") %>%
  mutate(ll = Est-1.96*Se,
         ul = Est+1.96*Se) 
round(coef_AX, 2) %>%
# alpha can be as large as 0.9

# check simulated U for standard error

0.5*(1-0.5)
0.1*0.9

##########################
# residual confounding
#########################
# RD + binning 1
#################
cols <- c("RF_PHYS_MODVIG_CURR2", "RF_PHYS_MODVIG_CURR3", "ENTRY_AGE1", "SEX", "RACEI1", "EDUCM1", "HEALTH1",
          "BMI_CUR1", "HEI2015_TOTAL_SCORE1", "MPED_A_BEV_NOFOOD1",
          "SMOKE_QUIT1", "SMOKE_DOSE")
aarp_data.1 <- aarp_data
aarp_data.1[cols] <- lapply(aarp_data.1[cols], factor)

fit.naive.RD.t <- coxph(Surv(PERSONYRS, respiratory_mort) ~ RF_PHYS_MODVIG_CURR2
                        + ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
                        + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                        + SMOKE_QUIT1 + SMOKE_DOSE, data=aarp_data.1, model = TRUE)

# sym 1 mild
system.time(fit.mcsimex.1.1 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.1), 
                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.1.1, file="fit.mcsimex.1.1.RData")

# sym 1 moderate
system.time(fit.mcsimex.1.2 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.2), 
                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.1.2, file="fit.mcsimex.1.2.RData")

# sym 1 extreme
system.time(fit.mcsimex.1.3 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.3), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.1.3, file="fit.mcsimex.1.3.RData")

# asym 1 mild
system.time(fit.mcsimex.3.1 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.1), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.3.1, file="fit.mcsimex.3.1.RData")

# asym 1 moderate
system.time(fit.mcsimex.3.2 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.2), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.3.2, file="fit.mcsimex.3.2.RData")

# asym 1 extreme
system.time(fit.mcsimex.3.3 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.3), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.3.3, file="fit.mcsimex.3.3.RData")

# asym 2 mild
system.time(fit.mcsimex.4.1 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.1), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.4.1, file="fit.mcsimex.4.1.RData")

# asym 2 moderate
system.time(fit.mcsimex.4.2 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.2), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.4.2, file="fit.mcsimex.4.2.RData")

# asym 2 extreme
system.time(fit.mcsimex.4.3 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.3), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.4.3, file="fit.mcsimex.4.3.RData")
######################

# all + binning 1
#################
fit.naive.all.t <- coxph(Surv(PERSONYRS, all_mort) ~ RF_PHYS_MODVIG_CURR2
                         + ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
                         + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                         + SMOKE_QUIT1 + SMOKE_DOSE, data=aarp_data.1, model = TRUE)

# sym 1 mild
system.time(fit.mcsimex.all.1.1 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.1), 
                               SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.all.1.1, file="fit.mcsimex.all.1.1.RData")

# sym 1 moderate
system.time(fit.mcsimex.all.1.2 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.2), 
                               SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.all.1.2, file="fit.mcsimex.all.1.2.RData")

# sym 1 extreme
system.time(fit.mcsimex.all.1.3 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.3), 
                                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.all.1.3, file="fit.mcsimex.all.1.3.RData")

# asym 1 mild
system.time(fit.mcsimex.all.3.1 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.1), 
                                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.all.3.1, file="fit.mcsimex.all.3.1.RData")

# asym 1 moderate
system.time(fit.mcsimex.all.3.2 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.2), 
                                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.all.3.2, file="fit.mcsimex.all.3.2.RData")

# asym 1 extreme
system.time(fit.mcsimex.all.3.3 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.3), 
                                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.all.3.3, file="fit.mcsimex.all.3.3.RData")

# asym 2 mild
system.time(fit.mcsimex.all.4.1 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.1), 
                                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.all.4.1, file="fit.mcsimex.all.4.1.RData")

# asym 2 moderate
system.time(fit.mcsimex.all.4.2 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.2), 
                                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.all.4.2, file="fit.mcsimex.all.4.2.RData")

# asym 2 extreme
system.time(fit.mcsimex.all.4.3 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.3), 
                                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.all.4.3, file="fit.mcsimex.all.4.3.RData")
#####################

# lung + binning 1
###################
fit.naive.lung.t <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ RF_PHYS_MODVIG_CURR2
                          + ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
                          + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                          + SMOKE_QUIT1 + SMOKE_DOSE, data=aarp_data.1, model = TRUE)
# sym 1 mild
system.time(fit.mcsimex.lung.1.1 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.1), 
                                SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.lung.1.1, file="fit.mcsimex.lung.1.1.RData")

# sym 1 moderate
system.time(fit.mcsimex.lung.1.2 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.2), 
                                SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.lung.1.2, file="fit.mcsimex.lung.1.2.RData")

# sym 1 extreme
system.time(fit.mcsimex.lung.1.3 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.3), 
                                            SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.lung.1.3, file="fit.mcsimex.lung.1.3.RData")

# asym 1 mild
system.time(fit.mcsimex.lung.3.1 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.1), 
                                            SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.lung.3.1, file="fit.mcsimex.lung.3.1.RData")

# asym 1 moderate
system.time(fit.mcsimex.lung.3.2 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.2), 
                                            SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.lung.3.2, file="fit.mcsimex.lung.3.2.RData")

# asym 1 extreme
system.time(fit.mcsimex.lung.3.3 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.3), 
                                            SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.lung.3.3, file="fit.mcsimex.lung.3.3.RData")

# asym 2 mild
system.time(fit.mcsimex.lung.4.1 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.1), 
                                            SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.lung.4.1, file="fit.mcsimex.lung.4.1.RData")

# asym 2 moderate
system.time(fit.mcsimex.lung.4.2 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.2), 
                                            SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.lung.4.2, file="fit.mcsimex.lung.4.2.RData")

# asym 2 extreme
system.time(fit.mcsimex.lung.4.3 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.3), 
                                            SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.lung.4.3, file="fit.mcsimex.lung.4.3.RData")

################
# combine 9*3=27 results
outcomes <- c("RD", "lung", "all")
id1 <- c(1,3,4)
id2 <- c(1,2,3)
id.comb <- tidyr::crossing(id1, id2)
id.list <- apply(id.comb, 1, paste , collapse = "." )
rst4resid <- data.frame(matrix(nrow = 0, ncol = 4)) 
colnames(rst4resid) <- c("outcome", "type", "est", "se")
for (oc in outcomes){
  for (id in id.list){
    data_name = paste("fit.mcsimex.", oc, ".", id, sep="")
    rst <- get(data_name)
    est <- exp(rst$coefficients[1])
    se <- sqrt(rst$variance.jackknife[1,1])
    rst <- data.frame(outcome=oc, type=id, est=est, se=se)
    rst4resid <- rbind(rst4resid, rst)
  }
}

# output to csv

split.1 <- function(x){
  id.1 <- strsplit(x, "\\.")[[1]][1]
  if(id.1=="1"){
    return("symmetric")
  } else if (id.1=="3"){
    return("asymmetric")
  } else {
    return("asymmetric_slow")
  }
}
split.2 <- function(x){
  id.2 <- strsplit(x, "\\.")[[1]][2]
  if(id.2=="1"){
    return("mild")
  } else if (id.2=="2"){
    return("moderate")
  } else {
    return("extreme")
  }
}

rst4resid$matrix.type <- unlist(lapply(rst4resid$type, split.1))
rst4resid$mis.level <- unlist(lapply(rst4resid$type, split.2))

rst4resid$exp.se <- rst4resid$est*rst4resid$se
rst4resid$ll <- rst4resid$est-qnorm(1-0.05/2,0,1)*rst4resid$exp.se
rst4resid$ul <- rst4resid$est+qnorm(1-0.05/2,0,1)*rst4resid$exp.se
row.names(rst4resid) <- NULL

rst4resid$CI <- mapply(paste.ci, rst4resid$ll, rst4resid$ul)
rst4resid$output <- mapply(paste.rst, rst4resid$est, rst4resid$CI)

rst4resid <- rst4resid[with(rst4resid, order(outcome, mis.level)), ]

write.csv(rst4resid, file="rst4resid.csv")
########
# TODO
#######
# TODO: confirm the meaning of B=5
# TODO: binning 2
# TODO: change mcsimex to mcsimex.RD

###########
# Appendix
###########
# generate 9 misclassification matrix
########################
# sym 1 mild
true_prob <- 0.9
false_prob <- 1 - true_prob
p_ij.1.1 <- matrix(c(true_prob, false_prob/2, 0, 0, 0, 0,
                     false_prob, true_prob, false_prob/2, 0, 0, 0,
                     0, false_prob/2, true_prob, false_prob/2, 0, 0,
                     0, 0, false_prob/2, true_prob, false_prob/2, 0,
                     0, 0, 0, false_prob/2, true_prob, false_prob,
                     0, 0, 0, 0, false_prob/2, true_prob), nrow = 6, byrow = TRUE)
dimnames(p_ij.1.1) <- list(levels(aarp_data$SMOKE_DOSE), levels(aarp_data$SMOKE_DOSE))
p_ij.1.1 <- build.mc.matrix(p_ij.1.1, method = "jlt") # do it twice
check.mc.matrix(list(p_ij.1.1))

# moderate
true_prob <- 0.65
false_prob <- 1 - true_prob
p_ij.1.2 <- matrix(c(true_prob, false_prob/2, 0, 0, 0, 0,
                     false_prob, true_prob, false_prob/2, 0, 0, 0,
                     0, false_prob/2, true_prob, false_prob/2, 0, 0,
                     0, 0, false_prob/2, true_prob, false_prob/2, 0,
                     0, 0, 0, false_prob/2, true_prob, false_prob,
                     0, 0, 0, 0, false_prob/2, true_prob), nrow = 6, byrow = TRUE)
dimnames(p_ij.1.2) <- list(levels(aarp_data$SMOKE_DOSE), levels(aarp_data$SMOKE_DOSE))
p_ij.1.2 <- build.mc.matrix(p_ij.1.2, method = "jlt") # not random did twice
check.mc.matrix(list(p_ij.1.2))

# extreme
true_prob <- 0.35
false_prob <- 1 - true_prob
p_ij.1.3 <- matrix(c(true_prob, false_prob/2, 0, 0, 0, 0,
                     false_prob, true_prob, false_prob/2, 0, 0, 0,
                     0, false_prob/2, true_prob, false_prob/2, 0, 0,
                     0, 0, false_prob/2, true_prob, false_prob/2, 0,
                     0, 0, 0, false_prob/2, true_prob, false_prob,
                     0, 0, 0, 0, false_prob/2, true_prob), nrow = 6, byrow = TRUE)
dimnames(p_ij.1.3) <- list(levels(aarp_data$SMOKE_DOSE), levels(aarp_data$SMOKE_DOSE))
p_ij.1.3 <- build.mc.matrix(p_ij.1.3, method = "jlt") # not random did twice
check.mc.matrix(list(p_ij.1.3))

# asym 1 mild
true_prob <- 0.9
false_prob <- 1 - true_prob
prop <- 1 + 1/2 + 1/4 + 1/8 + 1/16
row <- c(true_prob, false_prob*(1/prop), false_prob*(1/(2*prop)),
         false_prob*(1/(4*prop)), false_prob*(1/(8*prop)), false_prob*(1/(16*prop)))

p_ij.3.1 <- matrix(c(row[1],row[2:6],
                     row[2:1],row[3:6],
                     row[3:1],row[4:6],
                     row[4:1],row[5:6],
                     row[5:1],row[6],
                     row[6:1]), nrow = 6, byrow = FALSE)
check.mc.matrix(list(p_ij.3.1))
dimnames(p_ij.3.1) <- list(levels(aarp_data$SMOKE_DOSE), levels(aarp_data$SMOKE_DOSE))

# asym 1 moderate
true_prob <- 0.7
false_prob <- 1 - true_prob
row <- c(true_prob, false_prob*(1/prop), false_prob*(1/(2*prop)),
         false_prob*(1/(4*prop)), false_prob*(1/(8*prop)), false_prob*(1/(16*prop)))
p_ij.3.2 <- matrix(c(row[1],row[2:6],
                     row[2:1],row[3:6],
                     row[3:1],row[4:6],
                     row[4:1],row[5:6],
                     row[5:1],row[6],
                     row[6:1]), nrow = 6, byrow = FALSE)
check.mc.matrix(list(p_ij.3.2))
dimnames(p_ij.3.2) <- list(levels(aarp_data$SMOKE_DOSE), levels(aarp_data$SMOKE_DOSE))

# asym 1 extreme
true_prob <- 0.5
false_prob <- 1 - true_prob
row <- c(true_prob, false_prob*(1/prop), false_prob*(1/(2*prop)),
         false_prob*(1/(4*prop)), false_prob*(1/(8*prop)), false_prob*(1/(16*prop)))
p_ij.3.3 <- matrix(c(row[1],row[2:6],
                     row[2:1],row[3:6],
                     row[3:1],row[4:6],
                     row[4:1],row[5:6],
                     row[5:1],row[6],
                     row[6:1]), nrow = 6, byrow = FALSE)
check.mc.matrix(list(p_ij.3.3))
dimnames(p_ij.3.3) <- list(levels(aarp_data$SMOKE_DOSE), levels(aarp_data$SMOKE_DOSE))

# asym 2 mild
true_prob <- 0.9
false_prob <- 1 - true_prob
prop <- 1 + 0.75 + (0.75)^2 + (0.75)^3 + (0.75)^4
row <- c(true_prob, false_prob*(1/prop), false_prob*(0.75/prop),
         false_prob*(0.75^2/prop), false_prob*(0.75^3/prop), false_prob*(0.75^4/prop))
p_ij.4.1 <- matrix(c(row[1],row[2:6],
                     row[2:1],row[3:6],
                     row[3:1],row[4:6],
                     row[4:1],row[5:6],
                     row[5:1],row[6],
                     row[6:1]), nrow = 6, byrow = FALSE)
check.mc.matrix(list(p_ij.4.1))
dimnames(p_ij.4.1) <- list(levels(aarp_data$SMOKE_DOSE), levels(aarp_data$SMOKE_DOSE))

# asym 2 moderate
true_prob <- 0.7
false_prob <- 1 - true_prob
row <- c(true_prob, false_prob*(1/prop), false_prob*(0.75/prop),
         false_prob*(0.75^2/prop), false_prob*(0.75^3/prop), false_prob*(0.75^4/prop))
p_ij.4.2 <- matrix(c(row[1],row[2:6],
                     row[2:1],row[3:6],
                     row[3:1],row[4:6],
                     row[4:1],row[5:6],
                     row[5:1],row[6],
                     row[6:1]), nrow = 6, byrow = FALSE)
check.mc.matrix(list(p_ij.4.2))
dimnames(p_ij.4.2) <- list(levels(aarp_data$SMOKE_DOSE), levels(aarp_data$SMOKE_DOSE))

# asym 2 extreme
true_prob <- 0.5
false_prob <- 1 - true_prob
row <- c(true_prob, false_prob*(1/prop), false_prob*(0.75/prop),
         false_prob*(0.75^2/prop), false_prob*(0.75^3/prop), false_prob*(0.75^4/prop))
p_ij.4.3 <- matrix(c(row[1],row[2:6],
                     row[2:1],row[3:6],
                     row[3:1],row[4:6],
                     row[4:1],row[5:6],
                     row[5:1],row[6],
                     row[6:1]), nrow = 6, byrow = FALSE)
check.mc.matrix(list(p_ij.4.3))
dimnames(p_ij.4.3) <- list(levels(aarp_data$SMOKE_DOSE), levels(aarp_data$SMOKE_DOSE))
