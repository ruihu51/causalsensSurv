rm(list=ls())

library(haven)

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

######################
# data pre processing
#####################

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

# 
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

############
#
############

fit.naive.RD.1 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ factor(RF_PHYS_MODVIG_CURR2) 
                        + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
                        + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                        + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data)
X <- cbind(factor(aarp_data$ENTRY_AGE1), factor(aarp_data$SEX), factor(aarp_data$RACEI1), factor(aarp_data$EDUCM1), factor(aarp_data$HEALTH1),
           factor(aarp_data$BMI_CUR1), factor(aarp_data$HEI2015_TOTAL_SCORE1), factor(aarp_data$MPED_A_BEV_NOFOOD1),
           factor(aarp_data$SMOKE_QUIT1), factor(aarp_data$SMOKE_DOSE))
beta <- c(0.5, 1, 1.5, 2, 2.5)
alpha <- c(-1.5, -1, -0.5)
beta <- c(2, 2.5)
alpha <- c(-0.5)
beta <- c(0.5, 1, 1.5)
alpha <- c(-1.5, -1, -0.5)
system.time(fit.stoEM_reg.b1.RD <- survSensitivity(aarp_data$PERSONYRS, aarp_data$respiratory_mort, 
                                 aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                 zetaT = beta, zetaZ = alpha,
                                 B = 5))
theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se
save(fit.stoEM_reg.b1.RD, file="fit.stoEM_reg.b1.RD.RData")

fit.stoEM_reg.b1.RD.1 <- fit.stoEM_reg

exp(theta.hat.stoEM_reg)


exp(theta.hat.stoEM_reg.se)

fit.stoEM_reg.b1.RD.1
fit.stoEM_reg.b1.RD.2
colnames(fit.stoEM_reg.b1.RD.3)[1:4]

fit.stoEM_reg.b1.RD
rst4unobs <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(rst4unobs) <- colnames(fit.stoEM_reg.b1.RD.3)[1:4]
rst4unobs[1,] <- c(-0.5, 0.5, 0.8556484, 0.0248)
rst4unobs[2,] <- c(-0.5, 1, 0.94, 0.025)

rst4unobs[3,] <- c(-0.5, 1.5, 1.014, 0.027)

rst4unobs[4,] <- c(-1, 0.5, 0.9386, 0.0248)
rst4unobs[5,] <- c(-1, 1, 1.12, 0.025)
rst4unobs[6,] <- c(-1, 1.5, 1.31, 0.025)

rst4unobs[7,] <- c(-1.5, 0.5, 1.01605, 0.0248)
rst4unobs[8,] <- c(-1.5, 1, 1.31, 0.025)
rst4unobs[9,] <- c(-1.5, 1.5, 1.65, 0.027)

rst4unobs <- rbind(rst4unobs, fit.stoEM_reg.b1.RD.1)

rst4unobs <- rbind(rst4unobs, fit.stoEM_reg.b1.RD.2)
rst4unobs <- rbind(rst4unobs, fit.stoEM_reg.b1.RD.3)
fit.stoEM_reg.b1.RD.2$tau1 <- exp(fit.stoEM_reg.b1.RD.2$tau1)
fit.stoEM_reg.b1.RD.3$tau1 <- exp(fit.stoEM_reg.b1.RD.3$tau1)
fit.stoEM_reg.b1.RD.1[,1:4]
fit.stoEM_reg.b1.RD.1$tau1 <- exp(fit.stoEM_reg.b1.RD.1$tau1)
fit.stoEM_reg.b1.RD$tau1 <- exp(fit.stoEM_reg.b1.RD$tau1)
rst4unobs <- fit.stoEM_reg.b1.RD

rst4unobs.RD <- rst4unobs

rst4unobs
plotsens(rst4unobs, coeff0 = 0.79)
##########
aarp_data$all_mort
fit.naive.all.1 <- coxph(Surv(PERSONYRS, all_mort) ~ factor(RF_PHYS_MODVIG_CURR2) 
                        + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
                        + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                        + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data)
fit.naive.all.1$coefficients
beta <- c(0.5, 1, 1.5, 2, 2.5)
alpha <- c(-1.5, -1, -0.5)
system.time(fit.stoEM_reg.b1.all_mort <- survSensitivity(aarp_data$PERSONYRS, 
                                             aarp_data$all_mort, 
                                             aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                             zetaT = beta, zetaZ = alpha,
                                             B = 5))
plotsens(fit.stoEM_reg.b1.all_mort, coeff0 = -0.05)
beta <- c(0)
alpha <- c(-1.5, -1, -0.5, 0)
system.time(fit.stoEM_reg.b1.all_mort.1 <- survSensitivity(aarp_data$PERSONYRS, 
                                                         aarp_data$all_mort, 
                                                         aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                         zetaT = beta, zetaZ = alpha,
                                                         B = 5))
plotsens(fit.stoEM_reg.b1.all_mort, coeff0 = -0.05)
save(fit.stoEM_reg.b1.all_mort.1, file="fit.stoEM_reg.b1.all_mort.1.RData")
beta <- c(0.5, 1, 1.5, 2, 2.5)
alpha <- c(0)
system.time(fit.stoEM_reg.b1.all_mort.2 <- survSensitivity(aarp_data$PERSONYRS, 
                                                           aarp_data$all_mort, 
                                                           aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                           zetaT = beta, zetaZ = alpha,
                                                           B = 5))
plotsens(fit.stoEM_reg.b1.all_mort, coeff0 = -0.05)
save(fit.stoEM_reg.b1.all_mort.2, file="fit.stoEM_reg.b1.all_mort.2.RData")

# TODO: add  0 in alpha and beta


fit.stoEM_reg.b1.all_mort
fit.stoEM_reg.b1.all_mort.1
fit.stoEM_reg.b1.all_mort.2
rst4unobs.all <- rbind(fit.stoEM_reg.b1.all_mort, fit.stoEM_reg.b1.all_mort.1)
rst4unobs.all <- rbind(rst4unobs.all, fit.stoEM_reg.b1.all_mort.2)
rst4unobs.all$tau1 <- exp(rst4unobs.all$tau1)
plotsens(rst4unobs.all, coeff0 = round(exp(fit.naive.all.1$coefficients[1]),4))
save(rst4unobs.all, file="rst4unobs.all.RData")

#######
aarp_data$LUNG_MORT1
fit.naive.lung.1 <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ factor(RF_PHYS_MODVIG_CURR2) 
                         + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
                         + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                         + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data)
fit.naive.lung.1$coefficients
beta <- c(0, 0.5, 1, 1.5)
alpha <- c(-1.5, -1, -0.5, 0)
system.time(fit.stoEM_reg.b1.lung_mort <- survSensitivity(aarp_data$PERSONYRS, 
                                                         aarp_data$LUNG_MORT1, 
                                                         aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                         zetaT = beta, zetaZ = alpha,
                                                         B = 5))
save(fit.stoEM_reg.b1.lung_mort, file="fit.stoEM_reg.b1.lung_mort.RData")
beta <- c(2, 2.5)
alpha <- c(-1.5, -1, -0.5, 0)
system.time(fit.stoEM_reg.b1.lung_mort.1 <- survSensitivity(aarp_data$PERSONYRS, 
                                                          aarp_data$LUNG_MORT1, 
                                                          aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                          zetaT = beta, zetaZ = alpha,
                                                          B = 5))
save(fit.stoEM_reg.b1.lung_mort.1, file="fit.stoEM_reg.b1.lung_mort.1.RData")

rst4unobs.lung <- rbind(fit.stoEM_reg.b1.lung_mort, fit.stoEM_reg.b1.lung_mort.1)
rst4unobs.lung$tau1 <- exp(rst4unobs.lung$tau1)
plotsens(rst4unobs.lung, coeff0 = round(exp(fit.naive.lung.1$coefficients[1]),4))
save(rst4unobs.lung, file="rst4unobs.lung.RData")

contour.data <- rbind(fit.stoEM_reg.b1.lung_mort, fit.stoEM_reg.b1.lung_mort.1)
contour.data$tau1 <- exp(contour.data$tau1)
plotsens(contour.data, coeff0 = 0.9355)

fit.stoEM_reg

ggplot(data = rst4unobs, aes(x = zetaz, y = zetat1)) +
  stat_contour(aes(z = tau1, colour = after_stat(NULL))) +
  # stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
  xlab("Coef. on U in model for treatment") +
  ylab("Coef. on U in model for response") +
  theme_bw() +
  annotate("text", x = 0, y = 0, label = 0.79) +
  metR::geom_text_contour(aes(z = tau1), skip = 0, min.size = 0)

install.packages("metR")
library(metR)

plotsens(contour.data, coeff0 = 0.79)

contour.data <- fit.stoEM_reg
contour.data$tau1 <- exp(fit.stoEM_reg$tau1)


#######################
#
#######################
cols <- c("RF_PHYS_MODVIG_CURR2", "RF_PHYS_MODVIG_CURR3", "ENTRY_AGE1", "SEX", "RACEI1", "EDUCM1", "HEALTH1",
          "BMI_CUR1", "HEI2015_TOTAL_SCORE1", "MPED_A_BEV_NOFOOD1",
          "SMOKE_QUIT1", "SMOKE_DOSE")
aarp_data.1 <- aarp_data
aarp_data.1[cols] <- lapply(aarp_data.1[cols], factor)
fit.naive.RD.t <- coxph(Surv(PERSONYRS, respiratory_mort) ~ RF_PHYS_MODVIG_CURR2
                            + ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
                            + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                            + SMOKE_QUIT1 + SMOKE_DOSE, data=aarp_data.1, model = TRUE)
dimnames(p_ij.1.1) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
fit.mcsimex.1.1 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.1), 
                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE)
save(fit.mcsimex.1.1, file="fit.mcsimex.1.1.RData")
save(p_ij.1.1, file="p_ij.1.1.RData")
names(p_ij.1.1)
exp(fit.mcsimex.1.1$coefficients[1])
sqrt(fit.mcsimex.1.1$variance.jackknife[1,1])
dimnames(p_ij.1.2) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
save(p_ij.1.2, file="p_ij.1.2.RData")
fit.mcsimex.1.2 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.2), 
                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE)
save(fit.mcsimex.1.2, file="fit.mcsimex.1.2.RData")

dimnames(p_ij.1.3) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
save(p_ij.1.3, file="p_ij.1.3.RData")
system.time(fit.mcsimex.1.3 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.3), 
                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.1.3, file="fit.mcsimex.1.3.RData")

# asym 1 mild
dimnames(p_ij.3.1) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
save(p_ij.3.1, file="p_ij.3.1.RData")
system.time(fit.mcsimex.3.1 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.1), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.3.1, file="fit.mcsimex.3.1.RData")

# asym 1 moderate
dimnames(p_ij.3.2) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
save(p_ij.3.2, file="p_ij.3.2.RData")
system.time(fit.mcsimex.3.2 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.2), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.3.2, file="fit.mcsimex.3.2.RData")

# asym 1 extreme
dimnames(p_ij.3.3) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
save(p_ij.3.3, file="p_ij.3.3.RData")
system.time(fit.mcsimex.3.3 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.3.3), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.3.3, file="fit.mcsimex.3.3.RData")

# asym 2 mild
dimnames(p_ij.4.1) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
save(p_ij.4.1, file="p_ij.4.1.RData")
system.time(fit.mcsimex.4.1 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.1), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.4.1, file="fit.mcsimex.4.1.RData")

# asym 2 moderate
dimnames(p_ij.4.2) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
save(p_ij.4.2, file="p_ij.4.2.RData")
system.time(fit.mcsimex.4.2 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.2), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.4.2, file="fit.mcsimex.4.2.RData")

# asym 2 extreme
dimnames(p_ij.4.3) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
save(p_ij.4.3, file="p_ij.4.3.RData")
system.time(fit.mcsimex.4.3 <- mcsimex(fit.naive.RD.t, mc.matrix = list(SMOKE_DOSE=p_ij.4.3), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.4.3, file="fit.mcsimex.4.3.RData")

######
aarp_data.1$LUNG_MORT1
fit.naive.lung.t <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ RF_PHYS_MODVIG_CURR2
                        + ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
                        + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                        + SMOKE_QUIT1 + SMOKE_DOSE, data=aarp_data.1, model = TRUE)
fit.mcsimex.lung.1.1 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.1), 
                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE)
save(fit.mcsimex.lung.1.1, file="fit.mcsimex.lung.1.1.RData")
fit.mcsimex.lung.1.2 <- mcsimex(fit.naive.lung.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.2), 
                           SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE)
save(fit.mcsimex.lung.1.2, file="fit.mcsimex.lung.1.2.RData")
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

######################
# all cancer
aarp_data.1$all_mort
fit.naive.all.t <- coxph(Surv(PERSONYRS, all_mort) ~ RF_PHYS_MODVIG_CURR2
                          + ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
                          + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                          + SMOKE_QUIT1 + SMOKE_DOSE, data=aarp_data.1, model = TRUE)

# sym 1 mild
fit.mcsimex.all.1.1 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.1), 
                                SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE)
save(fit.mcsimex.all.1.1, file="fit.mcsimex.all.1.1.RData")

# sym 1 moderate
fit.mcsimex.all.1.2 <- mcsimex(fit.naive.all.t, mc.matrix = list(SMOKE_DOSE=p_ij.1.2), 
                                SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE)
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
##############


exp(fit.mcsimex.all.4.2$coefficients[1])
load(file = file_name)
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
paste("fit.mcsimex.", "RD", ".", "1.1", sep="")
save(rst4resid, file="rst4resid.RData")
fit.mcsimex.RD.4.3 <- fit.mcsimex.4.3
rst4resid$est <- round(rst4resid$est, 4)

######################
# output check
######################

# 
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

# sqrt(0.00105)*1.14

# check I: randomly select
rst1 <- survSensitivity(aarp_data$PERSONYRS, aarp_data$respiratory_mort, 
                aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                zetaT = 1, zetaZ = -1.5,
                B = 5)
exp(rst1$tau1)

rst2 <- survSensitivity(aarp_data$PERSONYRS, aarp_data$all_mort, 
                        aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                        zetaT = 0.5, zetaZ = -1,
                        B = 5)
exp(rst2$tau1)

rst3 <- survSensitivity(aarp_data$PERSONYRS, aarp_data$LUNG_MORT1, 
                        aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                        zetaT = 2, zetaZ = -1,
                        B = 5)
exp(rst3$tau1)

# check II: why ests are the same when zetaT=0
rst4unobs.RD[rst4unobs.RD$zetat1==0,]
rst4unobs.RD[rst4unobs.RD$zetaz==0,]

# zero situation won't violate the unobserved confounding assumption no need to report only needed for plot

# naive
naive <- rst4unobs.RD[(rst4unobs.RD$zetat1==0) & (rst4unobs.RD$zetaz==0),]
est <- naive$tau1
se <- naive$tau1.se
exp.se <- est*se
c(est-qnorm(1-0.05/2,0,1)*exp.se, est+qnorm(1-0.05/2,0,1)*exp.se)

naive <- rst4unobs.all[(rst4unobs.all$zetat1==0) & (rst4unobs.all$zetaz==0),]
est <- naive$tau1
se <- naive$tau1.se
exp.se <- est*se
c(est-qnorm(1-0.05/2,0,1)*exp.se, est+qnorm(1-0.05/2,0,1)*exp.se)

naive <- rst4unobs.lung[(rst4unobs.lung$zetat1==0) & (rst4unobs.lung$zetaz==0),]
est <- naive$tau1
se <- naive$tau1.se
exp.se <- est*se
c(est-qnorm(1-0.05/2,0,1)*exp.se, est+qnorm(1-0.05/2,0,1)*exp.se)

# add data
beta <- c(0, 0.5, 1, 1.5, 2, 2.5)
alpha <- c(0.5, 1, 1.5)
system.time(fit.stoEM_reg.b1.RD.4 <- survSensitivity(aarp_data$PERSONYRS, aarp_data$respiratory_mort, 
                                                     aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                     zetaT = beta, zetaZ = alpha,
                                                     B = 5))
save(fit.stoEM_reg.b1.RD.4, file="fit.stoEM_reg.b1.RD.4.RData")

load("rst4unobs.RD.RData")
rst4unobs.RD


d1 <- fit.stoEM_reg.b1.RD.4
d1$tau1 <- exp(d1$tau1)
d1$exp.se <- d1$tau1*d1$tau1.se
d1$ll <- d1$tau1-qnorm(1-0.05/2,0,1)*d1$exp.se
d1$ul <- d1$tau1+qnorm(1-0.05/2,0,1)*d1$exp.se
d1$CI <- mapply(paste.ci, d1$ll, d1$ul)
d1$output <- mapply(paste.rst, d1$tau1, d1$CI)

rst4unobs.RD <- rbind(rst4unobs.RD, d1)

