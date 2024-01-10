beta <- c(0, 0.5, 1, 1.5, 2, 2.5)
alpha <- c(0.5, 1, 1.5)
system.time(fit.stoEM_reg.b1.RD.4 <- survSensitivity(aarp_data$PERSONYRS, aarp_data$respiratory_mort, 
                                                   aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                   zetaT = beta, zetaZ = alpha,
                                                   B = 5))
save(fit.stoEM_reg.b1.RD.4, file="fit.stoEM_reg.b1.RD.4.RData")

beta <- c(0, 0.1, 0.5)
alpha <- c(-0.1, 0.1, 0.5)
system.time(fit.stoEM_reg.b1.lung_mort.5 <- survSensitivity(aarp_data$PERSONYRS, 
                                                            aarp_data$LUNG_MORT1, 
                                                            aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                            zetaT = beta, zetaZ = alpha,
                                                            B = 5))
save(fit.stoEM_reg.b1.lung_mort.5, file="fit.stoEM_reg.b1.lung_mort.5.RData")

beta <- 0.1
alpha <- -0.5

beta <- c(0, 0.1, 0.5)
alpha <- c(-0.1, 0.1, 0.5)
system.time(fit.stoEM_reg.b1.all_mort.5 <- survSensitivity(aarp_data$PERSONYRS, 
                                                            aarp_data$all_mort, 
                                                            aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                            zetaT = beta, zetaZ = alpha,
                                                            B = 5))
save(fit.stoEM_reg.b1.all_mort.5, file="fit.stoEM_reg.b1.all_mort.5.RData")

# add data for appendix
beta <- c(-1, -0.5)
alpha <- c(-1, -0.5, 0, 0.5, 1)
system.time(fit.stoEM_reg.b1.RD.6 <- survSensitivity(aarp_data$PERSONYRS, aarp_data$respiratory_mort, 
                                                     aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                     zetaT = beta, zetaZ = alpha,
                                                     B = 5))
save(fit.stoEM_reg.b1.RD.6, file="fit.stoEM_reg.b1.RD.6.RData")

beta <- c(-0.5, -0.1)
alpha <- c(-0.5, -0.1, 0, 0.1, 0.5)
system.time(fit.stoEM_reg.b1.lung_mort.6 <- survSensitivity(aarp_data$PERSONYRS, 
                                                            aarp_data$LUNG_MORT1, 
                                                            aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                            zetaT = beta, zetaZ = alpha,
                                                            B = 5))
save(fit.stoEM_reg.b1.lung_mort.6, file="fit.stoEM_reg.b1.lung_mort.6.RData")

beta <- c(-0.5, -0.1)
alpha <- c(-0.5, -0.1, 0, 0.1, 0.5)
system.time(fit.stoEM_reg.b1.all_mort.6 <- survSensitivity(aarp_data$PERSONYRS, 
                                                           aarp_data$all_mort, 
                                                           aarp_data$RF_PHYS_MODVIG_CURR2, X, "stoEM_reg",
                                                           zetaT = beta, zetaZ = alpha,
                                                           B = 5))
save(fit.stoEM_reg.b1.all_mort.6, file="fit.stoEM_reg.b1.all_mort.6.RData")

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

table(aarp_data$RF_PHYS_MODVIG_CURR3)
table(aarp_data$RF_PHYS_MODVIG_CURR3, aarp_data$respiratory_mort)
table(aarp_data$RF_PHYS_MODVIG_CURR3, aarp_data$LUNG_MORT1)
table(aarp_data$RF_PHYS_MODVIG_CURR3, aarp_data$all_mort)

fit.naive.RD.2 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ factor(RF_PHYS_MODVIG_CURR3) 
                        + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
                        + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                        + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data)

fit.naive.all.2 <- coxph(Surv(PERSONYRS, all_mort) ~ factor(RF_PHYS_MODVIG_CURR3) 
                         + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
                         + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                         + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data)

fit.naive.lung.2 <- coxph(Surv(PERSONYRS, LUNG_MORT1) ~ factor(RF_PHYS_MODVIG_CURR3) 
                          + factor(ENTRY_AGE1) + factor(SEX) + factor(RACEI1) + factor(EDUCM1) + factor(HEALTH1)
                          + factor(BMI_CUR1) + factor(HEI2015_TOTAL_SCORE1) + factor(MPED_A_BEV_NOFOOD1)
                          + factor(SMOKE_QUIT1) + factor(SMOKE_DOSE), data=aarp_data)

est <- exp(fit.naive.RD.2$coefficients[1])
se <- sqrt(fit.naive.RD.2$var[1,1])
exp.se <- est*se
c(est-qnorm(1-0.05/2,0,1)*exp.se, est+qnorm(1-0.05/2,0,1)*exp.se)

est <- exp(fit.naive.all.2$coefficients[1])
se <- sqrt(fit.naive.all.2$var[1,1])
exp.se <- est*se
c(est-qnorm(1-0.05/2,0,1)*exp.se, est+qnorm(1-0.05/2,0,1)*exp.se)

est <- exp(fit.naive.lung.2$coefficients[1])
se <- sqrt(fit.naive.lung.2$var[1,1])
exp.se <- est*se
c(est-qnorm(1-0.05/2,0,1)*exp.se, est+qnorm(1-0.05/2,0,1)*exp.se)

X <- cbind(factor(aarp_data$ENTRY_AGE1), factor(aarp_data$SEX), factor(aarp_data$RACEI1), factor(aarp_data$EDUCM1), factor(aarp_data$HEALTH1),
           factor(aarp_data$BMI_CUR1), factor(aarp_data$HEI2015_TOTAL_SCORE1), factor(aarp_data$MPED_A_BEV_NOFOOD1),
           factor(aarp_data$SMOKE_QUIT1), factor(aarp_data$SMOKE_DOSE))

beta <- c(0, 0.5, 1)
alpha <- c(-1, -0.5, 0, 0.5, 1)
system.time(fit.stoEM_reg.b2.RD <- survSensitivity(aarp_data$PERSONYRS, aarp_data$respiratory_mort, 
                                                   aarp_data$RF_PHYS_MODVIG_CURR3, X, "stoEM_reg",
                                                   zetaT = beta, zetaZ = alpha,
                                                   B = 5))
save(fit.stoEM_reg.b2.RD, file="fit.stoEM_reg.b2.RD.RData")

# all + binning 2
beta <- c(0, 0.1, 0.5)
alpha <- c(-0.5, -0.1, 0, 0.1, 0.5)
system.time(fit.stoEM_reg.b2.all_mort <- survSensitivity(aarp_data$PERSONYRS, 
                                                         aarp_data$all_mort, 
                                                         aarp_data$RF_PHYS_MODVIG_CURR3, X, "stoEM_reg",
                                                         zetaT = beta, zetaZ = alpha,
                                                         B = 5))
save(fit.stoEM_reg.b2.all_mort, file="fit.stoEM_reg.b2.all_mort.RData")

# lung + binning 2
system.time(fit.stoEM_reg.b2.lung_mort <- survSensitivity(aarp_data$PERSONYRS, 
                                                          aarp_data$LUNG_MORT1, 
                                                          aarp_data$RF_PHYS_MODVIG_CURR3, X, "stoEM_reg",
                                                          zetaT = beta, zetaZ = alpha,
                                                          B = 5))
save(fit.stoEM_reg.b2.lung_mort, file="fit.stoEM_reg.b2.lung_mort.RData")

#####################################
# organize
#####################################
d1 <- fit.stoEM_reg.b2.lung_mort
d1$tau1 <- exp(d1$tau1)
d1$exp.se <- d1$tau1*d1$tau1.se
d1$ll <- d1$tau1-qnorm(1-0.05/2,0,1)*d1$exp.se
d1$ul <- d1$tau1+qnorm(1-0.05/2,0,1)*d1$exp.se
d1$CI <- mapply(paste.ci, d1$ll, d1$ul)
d1$output <- mapply(paste.rst, d1$tau1, d1$CI)

rst4unobs.b2.lung_mort <- d1 %>% 
  select(c("zetat1", "zetaz", "output")) %>% 
  spread(zetaz, value = output)

write.csv(rst4unobs.b2.lung_mort, file="../output.data/rst4unobs.b2.lung_mort.csv")

rst4unobs.b2.all_mort
################################3

load("rst4unobs.RD.RData")
rst4unobs.RD



load("fit.stoEM_reg.b1.lung_mort.6.RData")

d1 <- fit.stoEM_reg.b1.RD.6
d1$tau1 <- exp(d1$tau1)
d1$exp.se <- d1$tau1*d1$tau1.se
d1$ll <- d1$tau1-qnorm(1-0.05/2,0,1)*d1$exp.se
d1$ul <- d1$tau1+qnorm(1-0.05/2,0,1)*d1$exp.se
d1$CI <- mapply(paste.ci, d1$ll, d1$ul)
d1$output <- mapply(paste.rst, d1$tau1, d1$CI)

d1
rst4unobs.RD <- rbind(rst4unobs.RD, d1)

rst4unobs.RD %>% 
  select(c("zetat1", "zetaz", "output")) %>% 
  spread(zetaz, value = output)

d1 <- fit.stoEM_reg.b1.lung_mort.5
d1$tau1 <- exp(d1$tau1)
d1$exp.se <- d1$tau1*d1$tau1.se
d1$ll <- d1$tau1-qnorm(1-0.05/2,0,1)*d1$exp.se
d1$ul <- d1$tau1+qnorm(1-0.05/2,0,1)*d1$exp.se
d1$CI <- mapply(paste.ci, d1$ll, d1$ul)
d1$output <- mapply(paste.rst, d1$tau1, d1$CI)

d1
load("rst4unobs.lung.RData")
rst4unobs.RD <- rbind(rst4unobs.RD, d1)
rst4unobs.lung

rst4unobs.lung <- rbind(rst4unobs.lung, d1)

save(rst4unobs.lung, file="rst4unobs.lung.RData")

d1 <- fit.stoEM_reg.b1.all_mort.4
d1$tau1 <- exp(d1$tau1)
d1$exp.se <- d1$tau1*d1$tau1.se
d1$ll <- d1$tau1-qnorm(1-0.05/2,0,1)*d1$exp.se
d1$ul <- d1$tau1+qnorm(1-0.05/2,0,1)*d1$exp.se
d1$CI <- mapply(paste.ci, d1$ll, d1$ul)
d1$output <- mapply(paste.rst, d1$tau1, d1$CI)

load("rst4unobs.all.RData")
d1
rst4unobs.all[1:24,]
rst4unobs.all <- rbind(rst4unobs.all[1:24,], d1)

save(rst4unobs.all, file="rst4unobs.all.RData")



rst4unobs.RD.42 <- rst4unobs.RD %>% 
  select(c("zetat1", "zetaz", "output")) %>% 
  spread(zetaz, value = output)

write.csv(rst4unobs.RD.42, file="rst4unobs.RD.42.csv")



rst4unobs.lung.33 <- rst4unobs.lung %>% 
  select(c("zetat1", "zetaz", "output")) %>% 
  spread(zetaz, value = output)


write.csv(rst4unobs.lung.33, file="rst4unobs.lung.33.csv")

rst4unobs.all.33 <- rst4unobs.all %>% 
  select(c("zetat1", "zetaz", "output")) %>% 
  spread(zetaz, value = output)


write.csv(rst4unobs.all.33, file="rst4unobs.all.33.csv")

rst4unobs.RD$t

plotsens(rst4unobs.lung, coeff0 = round(exp(fit.naive.lung.1$coefficients[1]),4))
ggplot(data = tau.res, aes(x = zetaz, y = zetat1)) +
  stat_contour(aes(z = tau1, colour = after_stat(level))) +
  stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
  xlab("Coef. on U in model for treatment") +
  ylab("Coef. on U in model for response") +
  theme_bw() +
  annotate("text", x = 0, y = 0, label = coeff0) +
  geom_text_contour(aes(z = tau1), skip = 0, min.size = 0)

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

fit.naive.RD.t <- coxph(Surv(PERSONYRS, respiratory_mort) ~ RF_PHYS_MODVIG_CURR2
                        + ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
                        + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                        + SMOKE_QUIT1 + SMOKE_DOSE, data=aarp_data.1, model = TRUE)

fit.naive.RD.t.1 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ RF_PHYS_MODVIG_CURR2
                        + ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
                        + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                        + SMOKE_QUIT1, data=aarp_data.1, model = TRUE)

fit.naive.RD.t.2 <- coxph(Surv(PERSONYRS, respiratory_mort) ~ RF_PHYS_MODVIG_CURR2
                          + ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
                          + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
                          + SMOKE_DOSE, data=aarp_data.1, model = TRUE)

table(aarp_data.1$SMOKE_DOSE)

load("../output.data/p_ij.4.3.RData")

p_ij.4.3

summary(fit.mcsimex.4.3)


# A ~ X
glm_AX <- glm(RF_PHYS_MODVIG_CURR2 ~ ENTRY_AGE1 + SEX + RACEI1 + EDUCM1 + HEALTH1
    + BMI_CUR1 + HEI2015_TOTAL_SCORE1 + MPED_A_BEV_NOFOOD1
    + SMOKE_QUIT1 + SMOKE_DOSE, data = aarp_data.1, family = "binomial")
summary(glm_AX)

prop.table(table(aarp_data$RF_PHYS_MODVIG_CURR2, aarp_data$SMOKE_DOSE), margin=1)
prop.table(table(aarp_data$respiratory_mort, aarp_data$SMOKE_DOSE), margin=1)
?prop.table

p <- c()
n <- 100
for (i in 1:6){
  p <- append(p, sample(1:n, 1))
  n <- n - p[length(p)]
}
p2 <- p
p <- append(p, sample(1:100, 1))
p
p[length(p)]
100 - p[-1]
n
p_ij.random <- matrix(c(p1,p2,p3,p4,p5,p6)/100, nrow = 6, byrow = FALSE)
check.mc.matrix(list(p_ij.random))
p_ij.random <- build.mc.matrix(p_ij.random, method = "jlt") # not random did twice
dimnames(p_ij.random) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
system.time(fit.mcsimex.random <- mcsimex(fit.naive.RD.t, 
                                          mc.matrix = list(SMOKE_DOSE=p_ij.random), 
                                       SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.random, file="fit.mcsimex.random.RData")
summary(fit.mcsimex.random)
exp(-0.21388)

save(aarp_data.1, file="aarp_data.1.RData")
save(aarp_data, file="aarp_data.RData")

p_ij.random <- matrix(c(p1[c(6,1,3,2,5,4)],
                        p2[c(5,3,2,1,4,6)],
                        p3,
                        p4,
                        p5[c(3,1,2,5,4,6)],
                        p6[c(6,3,4,2,5,1)])/100, nrow = 6, byrow = FALSE)

check.mc.matrix(list(p_ij.random))
p_ij.random <- build.mc.matrix(p_ij.random, method = "jlt") # not random did twice
p_ij.random.1 <- p_ij.random

dimnames(p_ij.random.1) <- list(levels(aarp_data.1$SMOKE_DOSE), levels(aarp_data.1$SMOKE_DOSE))
system.time(fit.mcsimex.random.1 <- mcsimex(fit.naive.RD.t, 
                                          mc.matrix = list(SMOKE_DOSE=p_ij.random.1), 
                                          SIMEXvariable = c("SMOKE_DOSE"), asymptotic = FALSE))
save(fit.mcsimex.random.1, file="fit.mcsimex.random.1.RData")

summary(fit.mcsimex.random.1)
exp(-0.21489)

######################
#
generate.data <- function(n, alpha=0, beta=0, theta = c(0,1), 
                          eta_a = list(eta_a1=c(0,-0.1),
                                       eta_a2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_t3=c(0,-0.2,-0.4,-0.6,-0.8,-1))){
  
  # covariates
  X1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
  X2 <- factor(sample(0:4, size = n, replace = TRUE, 
                      prob = c(0.12, 0.2, 0.29, 0.35, 0.04)))
  X3 <- factor(sample(0:5, size = n, replace = TRUE, 
                      prob = c(0.21, 0.30, 0.21, 0.15, 0.10, 0.03)))
  # U: unmeasured confounder
  U <- rbinom(n, 1, p=0.5) # dist-1
  # U <- rnorm(n, 0, 1) # dist-2
  # U <- runif(n, 0, 1) # dist-3
  sim.data <- data.frame(X1=X1, X2=X2, X3=X3, U=U)
  
  # treatment 
  glm_a <- with(sim.data, eta_a$eta_a1[X1] + eta_a$eta_a2[X2] 
                + eta_a$eta_a3[X3] + alpha*U)
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  rate <- exp(with(sim.data, theta[A] + eta_t$eta_t1[X1] + eta_t$eta_t2[X2] 
                   + eta_t$eta_t3[X3] + beta*U))
  # outcome
  T <- rexp(n, rate)
  C <- runif(n, 0.5, 1)
  
  # observed outcome
  Y <- pmin(T, C)
  D <- (T<=C)
  sim.data$Y <- Y
  sim.data$D <- D
  
  return(sim.data)
}


sim.data <- generate.data(n=160000, alpha=alpha, beta=beta, theta = c(0,0))
fit.true <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3 + U, data=sim.data)
fit.true
sqrt(0.25/1000)
sqrt(0.25/160000)

aarp_data %>%
  select(NDI_ICD9_RECODE_72, NDI_ICD10_RECODE_113, NDI_ICD_VERSION) %>%
  filter((is.na(NDI_ICD_VERSION)) & (NDI_ICD9_RECODE_72!=""))

summary(aarp_data$NDI_ICD_VERSION)
sum(table(aarp_data$NDI_ICD_VERSION))

aarp_data %>%
  select(NDI_ICD9_RECODE_72, NDI_ICD10_RECODE_113, NDI_ICD_VERSION) %>%
  filter(!is.na(NDI_ICD_VERSION)) %>%
  group_by(NDI_ICD10_RECODE_113) %>%
  summarise(n=n()) %>% 
  View

aarp_data %>%
  select(NDI_ICD9_RECODE_72, NDI_ICD10_RECODE_113, NDI_ICD_VERSION) %>%
  filter((NDI_ICD10_RECODE_113==".N"))

sum(table(aarp_data$EDUCM1))

coef_YX.RD <- as.data.frame(summary(fit.naive.RD.1)$coefficients[,1:3])

coef_YX.lung <- as.data.frame(summary(fit.naive.lung.1)$coefficients[,1:3])
coef_YX.lung.output <- round(coef_YX.lung %>%
  rename("exp.est"="exp(coef)",
         "se"="se(coef)") %>%
  mutate(exp.se = exp.est*se) %>%
  mutate(ll = exp.est-1.96*exp.se,
         ul = exp.est+1.96*exp.se) %>%
  select(c("exp.est", "ll", "ul")), 2)

coef_YX.all <- as.data.frame(summary(fit.naive.all.1)$coefficients[,1:3])
coef_YX.all.output <- round(coef_YX.all %>%
  rename("exp.est"="exp(coef)",
        "se"="se(coef)") %>%
  mutate(exp.se = exp.est*se) %>%
  mutate(ll = exp.est-1.96*exp.se,
        ul = exp.est+1.96*exp.se) %>%
  select(c("exp.est", "ll", "ul")), 2)
paste.ci <- function(x, y){
  ci <- paste(as.character(round(x,2)),
              as.character(round(y,2)),
              sep="-")
  return(paste("(", ci, ")", sep=""))
}

coef_YX.lung.output$CI <- mapply(paste.ci, 
                                 coef_YX.lung.output$ll, 
                                 coef_YX.lung.output$ul)


# paste point estimate and confidence interval
paste.rst <- function(x, y){
  return(paste(as.character(round(x,2)), y, sep=" "))
}

coef_YX.lung.output$output <- mapply(paste.rst, 
                                     coef_YX.lung.output$exp.est, 
                                     coef_YX.lung.output$CI)

coef_YX.lung.output
write.csv(coef_YX.lung.output, file="coef_YX.lung.output.csv")

coef_YX.all.output$CI <- mapply(paste.ci, 
                                 coef_YX.all.output$ll, 
                                 coef_YX.all.output$ul)
coef_YX.all.output$output <- mapply(paste.rst, 
                                     coef_YX.all.output$exp.est, 
                                     coef_YX.all.output$CI)
write.csv(coef_YX.all.output, file="coef_YX.all.output.csv")


load("rst4unobs.RD.RData")
load("rst4unobs.lung.RData")


rst4unobs.RD$outcome <- "Respiratory Disease"
rst4unobs.lung$outcome <- "Lung Cancer"

rst4unobs.lung.tmp <- rst4unobs.lung %>%
  filter((zetaz<0.5) & (zetaz>=-0.5) & (abs(zetat1)<1))

rst4unobs.lung.tmp[12,1:2] <- c(0.0, 0.1)
rst4unobs.lung.tmp[12,-c(1:2)] <- rst4unobs.lung.tmp[8,-c(1:2)]

rst4unobs.plot <- bind_rows(rst4unobs.RD %>%
                              filter((abs(zetaz)<=1) & (abs(zetat1)<=1)), 
                            rst4unobs.lung.tmp)

rst4unobs.plot$outcome <- factor(rst4unobs.plot$outcome,
                                 levels = c("Respiratory Disease", 
                                            "Lung Cancer"))



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
rst4unobs.RD %>%
  filter((abs(zetaz)<=1) & (abs(zetat1)<=1)) %>%
  ggplot(aes(x = zetaz, y = zetat1)) +
  stat_contour(aes(z = tau1)) + 
  stat_contour(aes(z = t), colour = "red", 
               linetype="dotdash", breaks = c(-1.96,1.96)) +
  metR::geom_text_contour(aes(z = tau1), skip = 0, min.size = 0) + 
  annotate("text", x = 0, y = 0, label = 0.67)

rst4unobs.lung.tmp %>%
  ggplot(aes(x = zetaz, y = zetat1)) +
    stat_contour(aes(z = tau1)) + 
    stat_contour(aes(z = t), colour = "red", 
                 linetype="dotdash", breaks = c(-1.96,1.96)) +
    annotate("text", x = 0, y = 0, label = 0.67) +
    metR::geom_text_contour(aes(z = tau1), skip = 0, min.size = 0)

rst4unobs.lung.tmp %>%
  ggplot(aes(x = zetaz, y = zetat1)) +
  stat_contour(aes(z = tau1)) +
  metR::geom_text_contour(aes(z = tau1), skip = 0, min.size = 0)

p <- rst4unobs.lung.tmp %>%
  ggplot(aes(x = zetaz, y = zetat1, z = tau1)) +
  geom_contour() +
  geom_dl(aes(label = ..level..), method = "top.bumptwice")

ggplot(rst4unobs.RD %>%
         filter((abs(zetaz)<=1) & (abs(zetat1)<=1)), aes(x = zetaz, y = zetat1, z = tau1)) +
  geom_contour() +
  geom_dl(aes(label = format(tau1, digits = 2)), method = "smart.grid")

ggplot(rst4unobs.lung.tmp, aes(x = zetaz, y = zetat1, z = tau1)) +
  geom_contour() +
  geom_dl(aes(label = format(tau1, digits = 2)), method = "smart.grid")

direct.label(p, aes(label = tau1))
library(lattice)
library(directlabels)
install.packages("directlabels")
dat <- melt(volcano)
brks <- c(100, 120, 140, 160)
g <- ggplot(dat, aes(x = Var1, y = Var2, z = value)) +
  geom_contour(colour='black', breaks = brks)+
  geom_dl(aes(label=..level..), method="bottom.pieces", 
          stat="contour",breaks = brks)
g
