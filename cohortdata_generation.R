rm(list=ls())

# install packages
library(survival)

# generate data following the data structure of the cohort data
# sample sizes
n <- 160000

# treatment - physical activity
A <- sample(0:2, size = n, replace = TRUE, prob = c(0.15, 0.35, 0.5))
cat("distribution of the treatment: ", prop.table(table(A)))
data <- data.frame(sort(A))
colnames(data) <- c("A")

# covariates
# binary - eg. gender
# 0-male # 1-female
sex  <- c()
prob_list <- matrix(c(0.72, 0.28,
                      0.71, 0.29,
                      0.59, 0.41), nrow = 3, byrow = TRUE)
size_list <- table(A)
for (i in 1:length(size_list)){
  sex  <- append(sex, 
                   sample(0:1, size = size_list[i], replace = TRUE, 
                          prob = prob_list[i,]))
}
data$sex <- factor(sex)

# ordinal - eg. age
# sample age based on A (but here assume age is independent of A)
age <- c()
for (size in table(A)){
  age <- append(age, 
               sample(0:4, size = size, replace = TRUE, 
                      prob = c(0.12, 0.2, 0.29, 0.35, 0.04)))
}
data$age <- factor(age)

# categorical - eg. health status
# sample health status based on A
# better health status indicates higher physical activity score
health <- c()
prob_list <- matrix(c(0.24, 0.42, 0.33, 0.01,
                      0.15, 0.38, 0.46, 0.01,
                      0.06, 0.25, 0.68, 0.01), nrow = 3, byrow = TRUE)
for (i in 1:length(size_list)){
  health <- append(health, 
                sample(0:3, size = size_list[i], replace = TRUE, 
                       prob = prob_list[i,])) # - as.factor
}
data$health <- factor(health)

# smoking behavior
# years since cessation
years_quit <- c()
prob_list <- matrix(c(0.66, 0.18, 0.11, 0.05,
                      0.71, 0.15, 0.10, 0.4,
                      0.82, 0.11, 0.05, 0.02), nrow = 3, byrow = TRUE)
for (i in 1:length(size_list)){
  years_quit <- append(years_quit, 
                   sample(0:3, size = size_list[i], replace = TRUE, 
                          prob = prob_list[i,])) # - as.factor
}
data$years_quit <- factor(years_quit)

# CPDs
cpd<- c()
prob_list <- matrix(c(0.18, 0.26, 0.21, 0.16, 0.14, 0.05,
                      0.21, 0.30, 0.21, 0.15, 0.10, 0.03,
                      0.33, 0.33, 0.17, 0.10, 0.05, 0.02), nrow = 3, byrow = TRUE)
for (i in 1:length(size_list)){
  cpd <- append(cpd, 
             sample(0:5, size = size_list[i], replace = TRUE, 
                    prob = prob_list[i,])) # - as.factor
}
data$cpd <- factor(cpd)

data$A <- factor(data$A)

# try to mimic Table 2
prop.table(table(data$A, data$sex), margin=1)
prop.table(table(data$A, data$age), margin=1)
prop.table(table(data$A, data$health), margin=1)
prop.table(table(data$A, data$years_quit), margin=1)
prop.table(table(data$A, data$cpd), margin=1)

# sample T conditional on A and X
# small trick: multiply factor variables
rate <- (1/30)*exp(with(data, c(0, -0.12, -0.19)[A] + c(0, -0.05)[sex] 
                        + c(0, 0.01, 0.01, 0.01, 0.01)[age] + c(0, -0.03, -0.01, 0)[health] 
                        + c(0, 0.03, 0.05, 0.07)[years_quit] + c(0, 0.02, 0.03, 0.05, 0.07, 0.09)[cpd]))
# outcome
# survival time T
# censoring time C
T <- rexp(n, rate)
C <- runif(n, 24, 26)

# observed outcome
# observed event time
# censoring indicator
Y <- pmin(T, C)
cat("follow-up years: ", mean(Y), sd(Y)) # mean - 18.9, sd - 6.3
D <- (T<=C)
cat("percentage of deaths: ", mean(D)) # around 0.5
data$Y <- Y
data$D <- D

# naive causal hazard ratio
# include all possible treatment values
fit.naive <- coxph(Surv(Y, D) ~ A + sex + age + health
                   + years_quit + cpd, data=data)
exp(fit.naive$coefficients[1:2])
# 1 vs 0 - 0.88
# 2 vs 0 - 0.83
fit.naive.1 <- coxph(Surv(Y, D) ~ A + sex + age + health
                   + years_quit + cpd, data = data[data$A!=2,])
exp(fit.naive.1$coefficients[1:2])
fit.naive.2 <- coxph(Surv(Y, D) ~ A + sex + age + health
                     + years_quit + cpd, data = data[data$A!=1,])
exp(fit.naive.2$coefficients[1:2])
# coefficient slightly different
# TODO: can we do separate models?


# possible measurement error (misclassification in this case because the variables are categorical)
# TODO: measurement error model for ordinal variables

# possible unmeasurement variable
U <- rbinom(n, 1, 0.5) # the simplest setup
# generate continuous variable and transform it into categorical variables
U <- runif(n, 0, 1)
U <- rnorm(n, 0, 1)

# add real data interpretation
# data (Y,D,A,W,X.star) A-0,1,2
# naive model
# 1 vs 0 - 0.88
# 2 vs 0 - 0.83
# take 1 vs 0 as the example
theta.hat.naive <- -0.12
# confidence intervals - (ll,ul), 0 - not include zero

# unmeasured confounding
# (alpha, beta) = (1,1), (2,2)
fit.stoEM_reg <- survSensitivity(sim.data$Y, sim.data$D, A, X, "stoEM_reg",
                                 zetaT = beta, zetaZ = alpha,
                                 B = 5)
theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se

# result table & interpretation
# alpha | beta | theta.hat.stoEM | CI
# -------------------------------------
# 0 | 0 | -0.12 | (ll,ul), 0
# 1 | 1 | -0.07 | (ll,ul), 0
# 2 | 2 | -0.01 | (ll,0,ul)
# problematic when the true effect is zero and the naive estimate suggests significantly different from zero
# add citation: a strong unmeasured confounder would be usually known

# okay if the previous result is conservative, i.e., theta.hat.naive underestimate
# say theta.hat.stoEM_reg=-0.25
# naive estimator tend to overestimate

# TODO: may be add simulation at the null

# residual confounding
sim.data$X3.star.cont <- as.numeric(sim.data$X3.star)-1
naive.model.continuous <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont, data=sim.data, model = TRUE)
var.true <- 0.5 # 0.8, 1, 1.5 or related to var(Y) 
# TODO: add simulation
fit.simex.1 <- simex(naive.model.continuous, measurement.error = var.true, 
                     SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
theta.hat.simex.1 <- fit.simex.1$coefficients[1]
theta.hat.simex.se.1 <- sqrt(fit.simex.1$variance.jackknife[1,1])

# need to find theoretical evidence saying that simex is not sensitive to misspecification?
# result table & interpretation
# sigma^2 | theta.hat.simex | CI
# -------------------------------------
# 0.4 | 0.0119 |
# 0.5 | 0.0087 | 
# 0.6 | 0.00731 |
# 0.8 | -0.0104 |
# when the true effect is zero, theta.hat.naive = 0.0186, theta.hat.simex = 0.0087


# TODO: check sec 6


#######################################################################################
# TODO: check the simulation result - unmeasured confounding
load(file="../sim.data/sim.1.1.1.1.2500.pa.RData")
sim.1.1 <- sim.1.1.1.1.2500.pa %>%
  gather(key='type', value='est', c("naive", "stoEM_reg", "naive.se", "stoEM_reg.se")) %>%
  filter(type %in% c("naive","stoEM_reg")) %>%
  mutate(bias=est-(1),
         setting="Moderate unmeasured confounder")

sim.1.1 %>%
  ggplot(aes(x=est, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 1, colour="blue", linetype = "longdash")
############

# TODO: check simulation result when true effect is zero what theta.hat.naive looks like
load(file = "../sim.data/sim.2.surv.4.4.1.2500.pa.RData")
colnames <- colnames(sim.2.surv.4.4.1.2500.pa)
ret.1.est <- sim.2.surv.4.4.1.2500.pa %>%
  gather(key='type', value='est', colnames[1:22]) %>%
  filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4",
                     "mcsimex.1", "mcsimex.2", "mcsimex.3", "mcsimex.4")) %>%
  mutate(bias=est-(0))

ret.1.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))
ret.1.est %>%
  filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4")) %>%
  ggplot(aes(x=est, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash")

# when the true effect is not zero
# naive estimate tend to overestimate
load(file = "../sim.data/sim.2.surv.4.4.3.2500.pa.RData")
colnames <- colnames(sim.2.surv.4.4.3.2500.pa)
ret.3.est <- sim.2.surv.4.4.3.2500.pa %>%
  gather(key='type', value='est', colnames[1:22]) %>%
  filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4",
                     "mcsimex.1", "mcsimex.2", "mcsimex.3", "mcsimex.4")) %>%
  mutate(bias=est-(0.5))

ret.3.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))
ret.3.est %>%
  filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4")) %>%
  ggplot(aes(x=est, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash")
















