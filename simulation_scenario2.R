rm(list=ls())

# install packages
library(survival)
library(simex)
library(mecor)
library(tidyr)
library(dplyr)
library(ggplot2)

# Scenario II
# Aim:
# 1) survival data - SIMEX estimator property
# 2) gaussian outcome - comparison between SIMEX and regression calibration


# 1)
generate.data <- function(n, theta = c(0,1), p_ij=p_ij,
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
  X3.star <- factor(unlist(lapply(X3, gen_chd_given_parents, category=levels(X3), prob_matrix=p_ij)))
  sim.data <- data.frame(X1=X1, X2=X2, X3=X3, X3.star=X3.star)
  
  # treatment 
  glm_a <- with(sim.data, eta_a$eta_a1[X1] + eta_a$eta_a2[X2] 
                + eta_a$eta_a3[X3])
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  rate <- exp(with(sim.data, theta[A] + eta_t$eta_t1[X1] + eta_t$eta_t2[X2] 
                   + eta_t$eta_t3[X3]))
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
# true misclassification matrix
p_ij <- matrix(c(0.8, 0.1, 0.01, 0.01, 0.01, 0.01,
                  0.16, 0.8, 0.1, 0.01, 0.01, 0.01,
                  0.01, 0.07, 0.8, 0.1, 0.01, 0.01,
                  0.01, 0.01, 0.07, 0.8, 0.07, 0.01,
                  0.01, 0.01, 0.01, 0.07, 0.8, 0.16,
                  0.01, 0.01, 0.01, 0.01, 0.1, 0.8), nrow = 6, byrow = TRUE)
p_ij.2 <- matrix(c(0.8, 0.07, 0.01, 0.01, 0.01, 0.01,
                   0.1, 0.8, 0.07, 0.01, 0.01, 0.01,
                   0.07, 0.1, 0.8, 0.07, 0.01, 0.01,
                   0.01, 0.01, 0.1, 0.8, 0.1, 0.07,
                   0.01, 0.01, 0.01, 0.1, 0.8, 0.1,
                   0.01, 0.01, 0.01, 0.01, 0.07, 0.8), nrow = 6, byrow = TRUE)
p_ij.3 <- matrix(c(0.75, 0.04, 0.02, 0.02, 0.02, 0.02,
                   0.15, 0.75, 0.04, 0.02, 0.02, 0.02,
                   0.04, 0.15, 0.75, 0.04, 0.02, 0.02,
                   0.02, 0.02, 0.15, 0.75, 0.15, 0.04,
                   0.02, 0.02, 0.02, 0.15, 0.75, 0.15,
                   0.02, 0.02, 0.02, 0.02, 0.04, 0.75), nrow = 6, byrow = TRUE) # why 0.7,0.2 doesn't work
p_ij.4 <- matrix(c(0.85, 0.02, 0.01, 0.01, 0.01, 0.01,
                   0.1, 0.85, 0.02, 0.01, 0.01, 0.01,
                   0.02, 0.1, 0.85, 0.02, 0.01, 0.01,
                   0.01, 0.01, 0.1, 0.85, 0.1, 0.02,
                   0.01, 0.01, 0.01, 0.1, 0.85, 0.1,
                   0.01, 0.01, 0.01, 0.01, 0.02, 0.85), nrow = 6, byrow = TRUE)
n <- 2500
sim.data <- generate.data(n=n, theta = c(0,0.5), p_ij = p_ij,
                          eta_a = list(eta_a1=c(0,-0.1),
                                       eta_a2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_t3=c(0,-0.2,-0.4,-0.6,-0.8,-1)))

dimnames(p_ij) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))
dimnames(p_ij.2) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))
dimnames(p_ij.3) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))
dimnames(p_ij.4) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))

seed <- ret$seed[183]
for (i in 183){
  seed <- ret$seed[i]
  set.seed(seed)
  cat(i, seed, "\n")
  sim.data <- generate.data(n=n, theta = c(0,0.5), p_ij = p_ij,
                            eta_a = list(eta_a1=c(0,-0.1),
                                         eta_a2=c(0,-0.01,-0.02,-0.03,-0.04),
                                         eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                            eta_t = list(eta_t1=c(0,-0.1),
                                         eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                         eta_t3=c(0,-0.2,-0.4,-0.6,-0.8,-1)))
  fit.true <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3, data=sim.data)
  fit.obs <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data)
  fit.crude <- coxph(Surv(Y, D) ~ A + X1 + X2, data=sim.data)
  
  naive.model <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data, model = TRUE)
  sim.data$X3.star.cont <- as.numeric(sim.data$X3.star)
  naive.model.continuous <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont, data=sim.data, model = TRUE)
}
# true
fit.true <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3, data=sim.data)
theta.hat.true <- fit.true$coefficients[1]
theta.hat.true.se <- sqrt(fit.true$var[1,1])

# observed adjusted causal hazard ratio
fit.obs <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data)
theta.hat.obs <- fit.obs$coefficients[1]
theta.hat.obs.se <- sqrt(fit.obs$var[1,1])

# crude effect
fit.crude <- coxph(Surv(Y, D) ~ A + X1 + X2, data=sim.data)
theta.hat.crude <- fit.crude$coefficients[1]
theta.hat.crude.se <- sqrt(fit.crude$var[1,1])

# mcsimex
naive.model <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data, model = TRUE)
fit.mcsimex.1 <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij), 
                         SIMEXvariable = c("X3.star"), asymptotic = FALSE)
theta.hat.mcsimex.1 <- fit.mcsimex.1$coefficients[1]
theta.hat.mcsimex.se.1 <- sqrt(fit.mcsimex.1$variance.jackknife[1,1])

fit.mcsimex.2 <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.2), 
                         SIMEXvariable = c("X3.star"), asymptotic = FALSE)
theta.hat.mcsimex.2 <- fit.mcsimex.2$coefficients[1]
theta.hat.mcsimex.se.2 <- sqrt(fit.mcsimex.2$variance.jackknife[1,1])

fit.mcsimex.3 <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.3), 
                         SIMEXvariable = c("X3.star"), asymptotic = FALSE)
theta.hat.mcsimex.3 <- fit.mcsimex.3$coefficients[1]
theta.hat.mcsimex.se.3 <- sqrt(fit.mcsimex.3$variance.jackknife[1,1])

fit.mcsimex.4 <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.4), 
                         SIMEXvariable = c("X3.star"), asymptotic = FALSE)
theta.hat.mcsimex.4 <- fit.mcsimex.4$coefficients[1]
theta.hat.mcsimex.se.4 <- sqrt(fit.mcsimex.4$variance.jackknife[1,1])

# simex
sim.data$X3.star.cont <- as.numeric(sim.data$X3.star)-1
naive.model.continuous <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont, data=sim.data, model = TRUE)
var.true <- var(as.numeric(sim.data$X3)-as.numeric(sim.data$X3.star))
fit.simex.1 <- simex(naive.model.continuous, measurement.error = var.true, 
                     SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
theta.hat.simex.1 <- fit.simex.1$coefficients[1]
theta.hat.simex.se.1 <- sqrt(fit.simex.1$variance.jackknife[1,1])

fit.simex.2 <- simex(naive.model.continuous, measurement.error = 0.5, 
                     SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
theta.hat.simex.2 <- fit.simex.2$coefficients[1]
theta.hat.simex.se.2 <- sqrt(fit.simex.2$variance.jackknife[1,1])

fit.simex.3 <- simex(naive.model.continuous, measurement.error = 0.8, 
                     SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
theta.hat.simex.3 <- fit.simex.3$coefficients[1]
theta.hat.simex.se.3 <- sqrt(fit.simex.3$variance.jackknife[1,1])

fit.simex.4 <- simex(naive.model.continuous, measurement.error = 0.4, 
                     SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
theta.hat.simex.4 <- fit.simex.4$coefficients[1]
theta.hat.simex.se.4 <- sqrt(fit.simex.4$variance.jackknife[1,1])


  
ret <- data.frame()
for (n in 2500) {
  for (j in 110:500) {
    
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)
    cat(n, j, seed, '\n')
    
    # generate data
    sim.data <- generate.data(n=n, theta = c(0,0.1), p_ij = p_ij,
                              eta_a = list(eta_a1=c(0,-0.1),
                                           eta_a2=c(0,-0.01,-0.02,-0.03,-0.04),
                                           eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                              eta_t = list(eta_t1=c(0,-0.1),
                                           eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                           eta_t3=c(0,-0.2,-0.4,-0.6,-0.8,-1)))
    
    try({
      # true
      fit.true <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3, data=sim.data)
      theta.hat.true <- fit.true$coefficients[1]
      theta.hat.true.se <- sqrt(fit.true$var[1,1])
      
      # observed adjusted causal hazard ratio
      fit.obs <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data)
      theta.hat.obs <- fit.obs$coefficients[1]
      theta.hat.obs.se <- sqrt(fit.obs$var[1,1])
      
      # crude effect
      fit.crude <- coxph(Surv(Y, D) ~ A + X1 + X2, data=sim.data)
      theta.hat.crude <- fit.crude$coefficients[1]
      theta.hat.crude.se <- sqrt(fit.crude$var[1,1])
      
      # mcsimex
      naive.model <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data, model = TRUE)
      fit.mcsimex.1 <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij), 
                             SIMEXvariable = c("X3.star"), asymptotic = FALSE)
      theta.hat.mcsimex.1 <- fit.mcsimex.1$coefficients[1]
      theta.hat.mcsimex.se.1 <- sqrt(fit.mcsimex.1$variance.jackknife[1,1])
      
      fit.mcsimex.2 <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.2), 
                               SIMEXvariable = c("X3.star"), asymptotic = FALSE)
      theta.hat.mcsimex.2 <- fit.mcsimex.2$coefficients[1]
      theta.hat.mcsimex.se.2 <- sqrt(fit.mcsimex.2$variance.jackknife[1,1])
      
      fit.mcsimex.3 <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.3), 
                               SIMEXvariable = c("X3.star"), asymptotic = FALSE)
      theta.hat.mcsimex.3 <- fit.mcsimex.3$coefficients[1]
      theta.hat.mcsimex.se.3 <- sqrt(fit.mcsimex.3$variance.jackknife[1,1])
      
      fit.mcsimex.4 <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.4), 
                               SIMEXvariable = c("X3.star"), asymptotic = FALSE)
      theta.hat.mcsimex.4 <- fit.mcsimex.4$coefficients[1]
      theta.hat.mcsimex.se.4 <- sqrt(fit.mcsimex.4$variance.jackknife[1,1])
      
      # simex
      sim.data$X3.star.cont <- as.numeric(sim.data$X3.star)-1
      naive.model.continuous <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont, data=sim.data, model = TRUE)
      var.true <- var(as.numeric(sim.data$X3)-as.numeric(sim.data$X3.star))
      fit.simex.1 <- simex(naive.model.continuous, measurement.error = var.true, 
                             SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.1 <- fit.simex.1$coefficients[1]
      theta.hat.simex.se.1 <- sqrt(fit.simex.1$variance.jackknife[1,1])
      
      fit.simex.2 <- simex(naive.model.continuous, measurement.error = 0.5, 
                           SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.2 <- fit.simex.2$coefficients[1]
      theta.hat.simex.se.2 <- sqrt(fit.simex.2$variance.jackknife[1,1])
      
      fit.simex.3 <- simex(naive.model.continuous, measurement.error = 0.8, 
                           SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.3 <- fit.simex.3$coefficients[1]
      theta.hat.simex.se.3 <- sqrt(fit.simex.3$variance.jackknife[1,1])
      
      fit.simex.4 <- simex(naive.model.continuous, measurement.error = 0.4, 
                           SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.4 <- fit.simex.4$coefficients[1]
      theta.hat.simex.se.4 <- sqrt(fit.simex.4$variance.jackknife[1,1])
      
      if (length(ret)>0){
        ret.tmp <- data.frame(true = theta.hat.true, true.se = theta.hat.true.se,
                              observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                              crude = theta.hat.crude, crude.se = theta.hat.crude.se, 
                              mcsimex.1 = theta.hat.mcsimex.1, mcsimex.se.1 = theta.hat.mcsimex.se.1,
                              mcsimex.2 = theta.hat.mcsimex.2, mcsimex.se.2 = theta.hat.mcsimex.se.2,
                              mcsimex.3 = theta.hat.mcsimex.3, mcsimex.se.3 = theta.hat.mcsimex.se.3,
                              mcsimex.4 = theta.hat.mcsimex.4, mcsimex.se.4 = theta.hat.mcsimex.se.4,
                              simex.1 = theta.hat.simex.1, simex.se.1 = theta.hat.simex.se.1,
                              simex.2 = theta.hat.simex.2, simex.se.2 = theta.hat.simex.se.2,
                              simex.3 = theta.hat.simex.3, simex.se.3 = theta.hat.simex.se.3,
                              simex.4 = theta.hat.simex.4, simex.se.4 = theta.hat.simex.se.4,
                              n=n, j=j, seed=seed)
        
        ret <- rbind(ret, ret.tmp)
      } else {
        ret <- data.frame(true = theta.hat.true, true.se = theta.hat.true.se,
                          observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                          crude = theta.hat.crude, crude.se = theta.hat.crude.se, 
                          mcsimex.1 = theta.hat.mcsimex.1, mcsimex.se.1 = theta.hat.mcsimex.se.1,
                          mcsimex.2 = theta.hat.mcsimex.2, mcsimex.se.2 = theta.hat.mcsimex.se.2,
                          mcsimex.3 = theta.hat.mcsimex.3, mcsimex.se.3 = theta.hat.mcsimex.se.3,
                          mcsimex.4 = theta.hat.mcsimex.4, mcsimex.se.4 = theta.hat.mcsimex.se.4,
                          simex.1 = theta.hat.simex.1, simex.se.1 = theta.hat.simex.se.1,
                          simex.2 = theta.hat.simex.2, simex.se.2 = theta.hat.simex.se.2,
                          simex.3 = theta.hat.simex.3, simex.se.3 = theta.hat.simex.se.3,
                          simex.4 = theta.hat.simex.4, simex.se.4 = theta.hat.simex.se.4,
                          n=n, j=j, seed=seed)
      }
    })
    
  }
}

# save data format
# sim.{scenario}.{eta}.{mis_prob}.{n}.{project}
# eta-1: eta_a3=c(0,-0.02,-0.04,-0.06,-0.08,-0.1)
# eta-2: eta_t3=c(0,0.2,0.4,0.6,0.8,1)
# TODO: add a dictionary to save the parameters settings
sim.2.1.1.2500.pa <- ret
save(sim.2.1.1.2500.pa, file="../sim.data/sim.2.1.1.2500.pa.RData")
sim.2.2.1.2500.pa <- ret
save(sim.2.2.1.2500.pa, file="../sim.data/sim.2.2.1.2500.pa.RData")
# sim.{scenario}.{gen}.{eta}.{mis_prob}.{A}.{n}.{project}
# A-1 0
# A-2 0.1
# A-3 0.5
# A-4 1
sim.2.surv.4.4.3.2500.pa <- ret
save(sim.2.surv.4.4.3.2500.pa, file="../sim.data/sim.2.surv.4.4.3.2500.pa.RData")
sim.2.surv.4.4.1.2500.pa <- ret
save(sim.2.surv.4.4.1.2500.pa, file="../sim.data/sim.2.surv.4.4.1.2500.pa.RData")
sim.2.surv.4.4.4.2500.pa <- ret
save(sim.2.surv.4.4.4.2500.pa, file="../sim.data/sim.2.surv.4.4.4.2500.pa.RData")
sim.2.surv.4.4.2.2500.pa <- ret
save(sim.2.surv.4.4.2.2500.pa, file="../sim.data/sim.2.surv.4.4.2.2500.pa.RData")
# summary

sim.2.surv.4.4.3.2500.pa %>%
  head()
load(file = "../sim.data/sim.2.surv.4.4.1.2500.pa.RData")
load(file = "../sim.data/sim.2.surv.4.4.2.2500.pa.RData")
load(file = "../sim.data/sim.2.surv.4.4.3.2500.pa.RData")
load(file = "../sim.data/sim.2.surv.4.4.4.2500.pa.RData")
colnames <- colnames(sim.2.surv.4.4.3.2500.pa)
ret.3.est <- sim.2.surv.4.4.3.2500.pa %>%
              gather(key='type', value='est', colnames[1:22]) %>%
              filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4",
                                 "mcsimex.1", "mcsimex.2", "mcsimex.3", "mcsimex.4")) %>%
              mutate(bias=est-(0.5))
ret.3.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

ret.4.est <- sim.2.surv.4.4.4.2500.pa %>%
  gather(key='type', value='est', colnames[1:22]) %>%
  filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4",
                     "mcsimex.1", "mcsimex.2", "mcsimex.3", "mcsimex.4")) %>%
  mutate(bias=est-(1))

ret.4.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

ret.1.est <- sim.2.surv.4.4.1.2500.pa %>%
  gather(key='type', value='est', colnames[1:22]) %>%
  filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4",
                     "mcsimex.1", "mcsimex.2", "mcsimex.3", "mcsimex.4")) %>%
  mutate(bias=est-(0))

ret.1.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

ret.2.est <- sim.2.surv.4.4.2.2500.pa %>%
  gather(key='type', value='est', colnames[1:22]) %>%
  filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4",
                     "mcsimex.1", "mcsimex.2", "mcsimex.3", "mcsimex.4")) %>%
  mutate(bias=est-(0.1))

ret.2.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))
sim.2.surv.1 <- ret.1.est %>%
  filter(type %in% c("observed","simex.1", "mcsimex.1")) %>%
  mutate(setting="No treatment effect null")
sim.2.surv.2 <- ret.2.est %>%
  filter(type %in% c("observed","simex.1", "mcsimex.1")) %>%
  mutate(setting="Weak alternative")
sim.2.surv.3 <- ret.3.est %>%
  filter(type %in% c("observed","simex.1", "mcsimex.1")) %>%
  mutate(setting="Moderate alternative")
sim.2.surv.4 <- ret.4.est %>%
  filter(type %in% c("observed","simex.1", "mcsimex.1")) %>%
  mutate(setting="Strong alternative")
   



fig.sim.2.surv.bias.1 <- bind_rows(
          sim.2.surv.1,
          sim.2.surv.2,
          sim.2.surv.3,
          sim.2.surv.4)

fig.sim.2.surv.bias.1$setting <- factor(fig.sim.2.surv.bias.1$setting,
                        levels = c("No treatment effect null",
                                   "Weak alternative", "Moderate alternative", "Strong alternative"))
fig.sim.2.surv.bias.1.plot <- fig.sim.2.surv.bias.1%>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme(legend.position="bottom")

ggsave(file="../figures/Fig_me_surv_bias.eps", width = 290,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig.sim.2.surv.bias.1.plot) #saves g

##################
sim.2.surv.1 <- ret.1.est %>%
  filter(type %in% c("observed", "mcsimex.1", "mcsimex.2", "mcsimex.3")) %>%
  mutate(setting="No treatment effect null")
sim.2.surv.2 <- ret.2.est %>%
  filter(type %in% c("observed", "mcsimex.1", "mcsimex.2", "mcsimex.3")) %>%
  mutate(setting="Weak alternative")
sim.2.surv.3 <- ret.3.est %>%
  filter(type %in% c("observed", "mcsimex.1", "mcsimex.2", "mcsimex.3")) %>%
  mutate(setting="Moderate alternative")
sim.2.surv.4 <- ret.4.est %>%
  filter(type %in% c("observed", "mcsimex.1", "mcsimex.2", "mcsimex.3")) %>%
  mutate(setting="Strong alternative")


sim.2.surv.1 <- ret.1.est %>%
  filter(type %in% c("observed", "simex.1", "simex.2", "simex.3")) %>%
  mutate(setting="No treatment effect null")
sim.2.surv.2 <- ret.2.est %>%
  filter(type %in% c("observed", "simex.1", "simex.2", "simex.3")) %>%
  mutate(setting="Weak alternative")
sim.2.surv.3 <- ret.3.est %>%
  filter(type %in% c("observed", "simex.1", "simex.2", "simex.3")) %>%
  mutate(setting="Moderate alternative")
sim.2.surv.4 <- ret.4.est %>%
  filter(type %in% c("observed", "simex.1", "simex.2", "simex.3")) %>%
  mutate(setting="Strong alternative")

fig.sim.2.surv.bias.2 <- bind_rows(
  sim.2.surv.1,
  sim.2.surv.2,
  sim.2.surv.3,
  sim.2.surv.4)

fig.sim.2.surv.bias.2 <- fig.sim.2.surv.bias.2 %>%
  mutate(model=case_when((type == "observed") ~ "observed",
    (type == "simex.1") ~ "simex (correct)", 
                         (type == "simex.2") ~ "simex (mild)",
                         (type == "simex.3") ~ "simex (severe)"
  ))

fig.sim.2.surv.bias.2 <- fig.sim.2.surv.bias.2 %>%
  mutate(model=case_when((type == "observed") ~ "observed", 
                         (type == "mcsimex.1") ~ "mcsimex (correct)", 
                         (type == "mcsimex.2") ~ "mcsimex (mild)",
                         (type == "mcsimex.3") ~ "mcsimex (severe)"
  ))

fig.sim.2.surv.bias.2$setting <- factor(fig.sim.2.surv.bias.2$setting,
                                        levels = c("No treatment effect null",
                                                   "Weak alternative", "Moderate alternative", "Strong alternative"))
fig.sim.2.surv.bias.2.plot.simex <- fig.sim.2.surv.bias.2%>%
  ggplot(aes(x=bias, colour = model)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme_bw() + 
  theme(legend.position="bottom")

ggsave(file="../figures/Fig_me_surv_bias_simex.eps", width = 240,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig.sim.2.surv.bias.2.plot.simex) #saves g


fig.sim.2.surv.bias.2.plot.mcsimex <- fig.sim.2.surv.bias.2%>%
  ggplot(aes(x=bias, colour = model)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme_bw() + 
  theme(legend.position="bottom")

ggsave(file="../figures/Fig_me_surv_bias_mcsimex.eps", width = 240,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig.sim.2.surv.bias.2.plot.mcsimex) #saves g


fig2a.plot <- fig2a %>%
  filter(!setting %in% c("Treatment effect homogeneity null")) %>%
  mutate(ci.type = factor(ci.type, levels=c("Wald", "Bootstrap"))) %>%
  rename(CI.type = ci.type) %>%
  ggplot() +
  geom_line(aes(log(n), coverage.c, colour=CI.type, linetype=CI.type)) +
  scale_linetype_manual(values=c("dotdash", "solid")) +
  scale_color_manual(values=c('#009900','#FF0000')) +
  scale_size_manual(values=c(2, 1)) +
  geom_hline(yintercept=.95, colour="#A0A0A0", linetype="longdash") +
  ylim(0.4, 1) +
  facet_grid(~ setting) +
  labs(y = expression(Empirical~coverage~of~confidence~intervals~"for"~psi[0])) +
  scale_x_continuous(
    name="Log sample size",
    limits = c(log(200), log(5500)),
    breaks=c(log(250), log(500), log(1000), log(2500), log(5000)),
    labels=c("log(250)", "log(500)", "log(1K)", "log(2.5K)", "log(5K)")) +
  theme_bw() +
    theme(legend.position="bottom")

ggsave(file="Figures/Fig20_coverage_Sim1.eps", width = 290,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig2a.plot) #saves g


ret.3.se <- sim.2.surv.4.4.3.2500.pa %>%
  gather(key='se.type', value='se', colnames[1:22]) %>%
  filter(se.type %in% c("true.se", "observed.se","simex.se.1", "simex.se.2", "simex.se.3", "simex.se.4",
                     "mcsimex.se.1", "mcsimex.se.2", "mcsimex.se.3", "mcsimex.se.4"))
ret.1.se <- sim.2.surv.4.4.1.2500.pa %>%
  gather(key='se.type', value='se', colnames[1:22]) %>%
  filter(se.type %in% c("true.se", "observed.se","simex.se.1", "simex.se.2", "simex.se.3", "simex.se.4",
                        "mcsimex.se.1", "mcsimex.se.2", "mcsimex.se.3", "mcsimex.se.4"))
ret.2.se <- sim.2.surv.4.4.2.2500.pa %>%
  gather(key='se.type', value='se', colnames[1:22]) %>%
  filter(se.type %in% c("true.se", "observed.se","simex.se.1", "simex.se.2", "simex.se.3", "simex.se.4",
                        "mcsimex.se.1", "mcsimex.se.2", "mcsimex.se.3", "mcsimex.se.4"))
ret.4.se <- sim.2.surv.4.4.4.2500.pa %>%
  gather(key='se.type', value='se', colnames[1:22]) %>%
  filter(se.type %in% c("true.se", "observed.se","simex.se.1", "simex.se.2", "simex.se.3", "simex.se.4",
                        "mcsimex.se.1", "mcsimex.se.2", "mcsimex.se.3", "mcsimex.se.4"))
ret.3 <- cbind(ret.3.est, ret.3.se)[c(1,2,3,4,5,6,10,11)]
ret.1 <- cbind(ret.1.est, ret.1.se)[c(1,2,3,4,5,6,10,11)]
ret.2 <- cbind(ret.2.est, ret.2.se)[c(1,2,3,4,5,6,10,11)]
ret.4 <- cbind(ret.4.est, ret.4.se)[c(1,2,3,4,5,6,10,11)]
# why we don't need to divided by n?
ret.3 %>%
  mutate(ll=est-qnorm(1-(1-0.95)/2)*se,
         ul=est+qnorm(1-(1-0.95)/2)*se) %>%
  mutate(coverage=(ll <= 0.5 & 0.5 <= ul)) %>%
  group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)),
            coverage = mean(coverage))

ret.1 %>%
  mutate(ll=est-qnorm(1-(1-0.95)/2)*se,
         ul=est+qnorm(1-(1-0.95)/2)*se) %>%
  mutate(coverage=(ll <= 0 & 0 <= ul)) %>%
  group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)),
            coverage = mean(coverage))

ret.2 %>%
  mutate(ll=est-qnorm(1-(1-0.95)/2)*se,
         ul=est+qnorm(1-(1-0.95)/2)*se) %>%
  mutate(coverage=(ll <= 0.1 & 0.1 <= ul)) %>%
  group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)),
            coverage = mean(coverage))

ret.4 %>%
  mutate(ll=est-qnorm(1-(1-0.95)/2)*se,
         ul=est+qnorm(1-(1-0.95)/2)*se) %>%
  mutate(coverage=(ll <= 1 & 1 <= ul)) %>%
  group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)),
            coverage = mean(coverage))


qnorm(1-(1-0.95)/2)
head(ret)
ret$ll = ret$est - qnorm(1-(1-control$conf.level)/2) * ret$se / sqrt(n)
ret$ul = ret$est + qnorm(1-(1-control$conf.level)/2) * ret$se / sqrt(n)

psi.summaries.1 <- ddply(subset(ests.sim.1, (type %in% 'psi.est')), .(n, type), summarize, na = sum(is.na(est)),
                         coverage = mean(ll <= psi0 & psi0 <= ul, na.rm=TRUE),
                         bias = mean(est - psi0, na.rm=TRUE),
                         var = var(est, na.rm=TRUE),
                         mse = mean((est - psi0)^2, na.rm=TRUE))

# plot
# check the bias distribution
ret %>%
  gather(key='type', value='est', colnames(ret)[1:22]) %>%
  filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4",
                     "mcsimex.1", "mcsimex.2", "mcsimex.3", "mcsimex.4")) %>%
  mutate(bias=est-(0.1)) %>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash")

#
# 2)
generate.data.gaussian <- function(n, theta = c(0,1), p_ij=p_ij,
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
  X3.star <- factor(unlist(lapply(X3, gen_chd_given_parents, category=levels(X3), prob_matrix=p_ij)))
  sim.data <- data.frame(X1=X1, X2=X2, X3=X3, X3.star=X3.star)
  
  # treatment 
  glm_a <- with(sim.data, eta_a$eta_a1[X1] + eta_a$eta_a2[X2] 
                + eta_a$eta_a3[X3])
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample Y conditional on A and X
  sim.data$mu <- with(sim.data, 0.5 + theta[A] + eta_t$eta_t1[X1] + eta_t$eta_t2[X2] + eta_t$eta_t3[X3])
  sim.data$Y <- rnorm(n, mean=sim.data$mu, 
                      sd=0.1)
  
  return(sim.data)
}
sim.data <- generate.data.gaussian(n=n, theta = c(0,0.5), p_ij = p_ij,
                          eta_a = list(eta_a1=c(0,-0.1),
                                       eta_a2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_t3=c(0,-0.2,-0.4,-0.6,-0.8,-1)))

ret <- data.frame()
for (n in 2500) {
  for (j in 1:500) {
    
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)
    cat(n, j, seed, '\n')
    
    # generate data
    sim.data <- generate.data.gaussian(n=n, theta = c(0,0.5), p_ij = p_ij,
                                       eta_a = list(eta_a1=c(0,-0.1),
                                                    eta_a2=c(0,-0.01,-0.02,-0.03,-0.04),
                                                    eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                                       eta_t = list(eta_t1=c(0,-0.1),
                                                    eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                                    eta_t3=c(0,-0.2,-0.4,-0.6,-0.8,-1)))
    
    try({
      # observed adjusted causal hazard ratio
      fit.obs <- glm(Y ~ A + X1 + X2 + X3.star, family = "gaussian", data=sim.data)
      theta.hat.obs <- fit.obs$coefficients[2]
      theta.hat.obs.se <- sqrt(vcov(fit.obs)[2,2])
      
      fit.obs.cont <- glm(Y ~ A + X1 + X2 + as.numeric(X3.star), family = "gaussian", data=sim.data)
      theta.hat.obs.cont <- fit.obs.cont$coefficients[2]
      theta.hat.obs.cont.se <- sqrt(vcov(fit.obs.cont)[2,2])
      
      # simex
      sim.data$X3.star.cont <- as.numeric(sim.data$X3.star)-1
      naive.model.continuous <- glm(Y ~ A + X1 + X2 + X3.star.cont, family = "gaussian", data=sim.data, x=T, y=T)
      var.true <- var(as.numeric(sim.data$X3)-as.numeric(sim.data$X3.star))
      fit.simex.1 <- simex(naive.model.continuous, measurement.error = var.true, 
                           SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.1 <- fit.simex.1$coefficients[2]
      theta.hat.simex.se.1 <- sqrt(fit.simex.1$variance.jackknife[2,2])
      
      fit.simex.2 <- simex(naive.model.continuous, measurement.error = var.true*1.1, 
                           SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.2 <- fit.simex.2$coefficients[2]
      theta.hat.simex.se.2 <- sqrt(fit.simex.2$variance.jackknife[2,2])
      
      fit.simex.3 <- simex(naive.model.continuous, measurement.error = var.true*1.6, 
                           SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.3 <- fit.simex.3$coefficients[2]
      theta.hat.simex.se.3 <- sqrt(fit.simex.3$variance.jackknife[2,2])
      
      fit.simex.4 <- simex(naive.model.continuous, measurement.error = var.true*0.6, 
                           SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.4 <- fit.simex.4$coefficients[2]
      theta.hat.simex.se.4 <- sqrt(fit.simex.4$variance.jackknife[2,2])
      
      # mecor
      fit.mecor.1 <- mecor(Y~MeasErrorRandom(as.numeric(X3.star), variance = var.true)+A+X1+X2, data=sim.data)
      theta.hat.mecor.1 <- fit.mecor.1$corfit$coef[3]
      theta.hat.mecor.se.1 <- sqrt(fit.mecor.1$corfit$zerovar_vcov[3,3])
      
      fit.mecor.2 <- mecor(Y~MeasErrorRandom(as.numeric(X3.star), variance = var.true*1.1)+A+X1+X2, data=sim.data)
      theta.hat.mecor.2 <- fit.mecor.2$corfit$coef[3]
      theta.hat.mecor.se.2 <- sqrt(fit.mecor.2$corfit$zerovar_vcov[3,3])
      
      fit.mecor.3 <- mecor(Y~MeasErrorRandom(as.numeric(X3.star), variance = var.true*1.6)+A+X1+X2, data=sim.data)
      theta.hat.mecor.3 <- fit.mecor.3$corfit$coef[3]
      theta.hat.mecor.se.3 <- sqrt(fit.mecor.3$corfit$zerovar_vcov[3,3])
      
      fit.mecor.4 <- mecor(Y~MeasErrorRandom(as.numeric(X3.star), variance = var.true*0.6)+A+X1+X2, data=sim.data)
      theta.hat.mecor.4 <- fit.mecor.4$corfit$coef[3]
      theta.hat.mecor.se.4 <- sqrt(fit.mecor.4$corfit$zerovar_vcov[3,3])
      
      if (length(ret)>0){
        ret.tmp <- data.frame(
                              observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                              observed.cont = theta.hat.obs.cont, observed.cont.se = theta.hat.obs.cont.se,
                              simex.1 = theta.hat.simex.1, simex.se.1 = theta.hat.simex.se.1,
                              simex.2 = theta.hat.simex.2, simex.se.2 = theta.hat.simex.se.2,
                              simex.3 = theta.hat.simex.3, simex.se.3 = theta.hat.simex.se.3,
                              simex.4 = theta.hat.simex.4, simex.se.4 = theta.hat.simex.se.4,
                              mecor.1 = theta.hat.mecor.1, mecor.se.1 = theta.hat.mecor.se.1,
                              mecor.2 = theta.hat.mecor.2, mecor.se.2 = theta.hat.mecor.se.2,
                              mecor.3 = theta.hat.mecor.3, mecor.se.3 = theta.hat.mecor.se.3,
                              mecor.4 = theta.hat.mecor.4, mecor.se.4 = theta.hat.mecor.se.4,
                              n=n, j=j, seed=seed)
        
        ret <- rbind(ret, ret.tmp)
      } else {
        ret <- data.frame(observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                          observed.cont = theta.hat.obs.cont, observed.cont.se = theta.hat.obs.cont.se,
                          simex.1 = theta.hat.simex.1, simex.se.1 = theta.hat.simex.se.1,
                          simex.2 = theta.hat.simex.2, simex.se.2 = theta.hat.simex.se.2,
                          simex.3 = theta.hat.simex.3, simex.se.3 = theta.hat.simex.se.3,
                          simex.4 = theta.hat.simex.4, simex.se.4 = theta.hat.simex.se.4,
                          mecor.1 = theta.hat.mecor.1, mecor.se.1 = theta.hat.mecor.se.1,
                          mecor.2 = theta.hat.mecor.2, mecor.se.2 = theta.hat.mecor.se.2,
                          mecor.3 = theta.hat.mecor.3, mecor.se.3 = theta.hat.mecor.se.3,
                          mecor.4 = theta.hat.mecor.4, mecor.se.4 = theta.hat.mecor.se.4,
                          n=n, j=j, seed=seed)
      }
    })
    
  }
}

sim.2.gaussian.4.4.3.2500.pa <- ret
save(sim.2.gaussian.4.4.3.2500.pa, file="../sim.data/sim.2.gaussian.4.4.3.2500.pa.RData")
sim.2.gaussian.4.4.1.2500.pa <- ret
save(sim.2.gaussian.4.4.1.2500.pa, file="../sim.data/sim.2.gaussian.4.4.1.2500.pa.RData")
sim.2.gaussian.4.4.4.2500.pa <- ret
save(sim.2.gaussian.4.4.4.2500.pa, file="../sim.data/sim.2.gaussian.4.4.4.2500.pa.RData")
sim.2.gaussian.4.4.2.2500.pa <- ret
save(sim.2.gaussian.4.4.2.2500.pa, file="../sim.data/sim.2.gaussian.4.4.2.2500.pa.RData")

# plot
# check the bias distribution
sim.2.gaussian.4.4.1.2500.pa %>%
  gather(key='type', value='est', colnames(ret)[1:20]) %>%
  filter(type %in% c("observed", "observed.cont","simex.1", "simex.2", "simex.3", "simex.4",
                     "mecor.1", "mecor.2", "mecor.3", "mecor.4")) %>%
  mutate(bias=est-(0)) %>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash")

colnames.2 <- colnames(ret)[1:20]
ret.1.3.est <- sim.2.gaussian.4.4.3.2500.pa %>%
  gather(key='type', value='est', colnames.2[1:20]) %>%
  filter(type %in% c("observed", "observed.cont","simex.1", "simex.2", "simex.3", "simex.4",
                     "mecor.1", "mecor.2", "mecor.3", "mecor.4")) %>%
  mutate(bias=est-(0.5))
ret.1.3.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

ret.1.3.est <- sim.2.gaussian.4.4.3.2500.pa %>%
  gather(key='type', value='est', colnames.2[1:20]) %>%
  filter(type %in% c("observed", "observed.cont","simex.1", "simex.2", "simex.3", "simex.4",
                     "mecor.1", "mecor.2", "mecor.3", "mecor.4")) %>%
  mutate(bias=est-(0.5))
ret.1.3.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

ret.1.3.est <- sim.2.gaussian.4.4.3.2500.pa %>%
  gather(key='type', value='est', colnames.2[1:20]) %>%
  filter(type %in% c("observed", "observed.cont","simex.1", "simex.2", "simex.3", "simex.4",
                     "mecor.1", "mecor.2", "mecor.3", "mecor.4")) %>%
  mutate(bias=est-(0.5))
ret.1.3.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

load(file = "../sim.data/sim.2.gaussian.4.4.1.2500.pa.RData")
ret.1.1.est <- sim.2.gaussian.4.4.1.2500.pa %>%
  gather(key='type', value='est', colnames.2[1:20]) %>%
  filter(type %in% c("observed", "observed.cont","simex.1", "simex.2", "simex.3", "simex.4",
                     "mecor.1", "mecor.2", "mecor.3", "mecor.4")) %>%
  mutate(bias=est-(0))
ret.1.1.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))
ret.1.1.se <- sim.2.gaussian.4.4.1.2500.pa %>%
  gather(key='se.type', value='se', colnames.2[1:20]) %>%
  filter(se.type %in% c("observed.se", "observed.cont.se","simex.se.1", "simex.se.2", "simex.se.3", "simex.se.4",
                        "mecor.se.1", "mecor.se.2", "mecor.se.3", "mecor.se.4"))

ret.1.1 <- cbind(ret.1.1.est, ret.1.1.se)[c(1,2,3,4,5,6,10,11)]
ret.1.1 %>%
  mutate(ll=est-qnorm(1-(1-0.95)/2)*se,
         ul=est+qnorm(1-(1-0.95)/2)*se) %>%
  mutate(coverage=(ll <= 0 & 0 <= ul)) %>%
  group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)),
            coverage = mean(coverage))
ret.1.2.est <- sim.2.gaussian.4.4.2.2500.pa %>%
  gather(key='type', value='est', colnames.2[1:20]) %>%
  filter(type %in% c("observed", "observed.cont","simex.1", "simex.2", "simex.3", "simex.4",
                     "mecor.1", "mecor.2", "mecor.3", "mecor.4")) %>%
  mutate(bias=est-(0.1))
ret.1.2.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

ret.1.3.est <- sim.2.gaussian.4.4.3.2500.pa %>%
  gather(key='type', value='est', colnames.2[1:20]) %>%
  filter(type %in% c("observed", "observed.cont","simex.1", "simex.2", "simex.3", "simex.4",
                     "mecor.1", "mecor.2", "mecor.3", "mecor.4")) %>%
  mutate(bias=est-(0.5))
ret.1.3.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

ret.1.4.est <- sim.2.gaussian.4.4.4.2500.pa %>%
  gather(key='type', value='est', colnames.2[1:20]) %>%
  filter(type %in% c("observed", "observed.cont","simex.1", "simex.2", "simex.3", "simex.4",
                     "mecor.1", "mecor.2", "mecor.3", "mecor.4")) %>%
  mutate(bias=est-(1))
ret.1.4.est %>% group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

sim.2.gaussian.1 <- ret.1.1.est %>%
  filter(type %in% c("observed", "simex.1", "simex.2", "simex.3")) %>%
  mutate(setting="No treatment effect null")
sim.2.gaussian.2 <- ret.1.2.est %>%
  filter(type %in% c("observed", "simex.1", "simex.2", "simex.3")) %>%
  mutate(setting="Weak alternative")
sim.2.gaussian.3 <- ret.1.3.est %>%
  filter(type %in% c("observed", "simex.1", "simex.2", "simex.3")) %>%
  mutate(setting="Moderate alternative")
sim.2.gaussian.4 <- ret.1.4.est %>%
  filter(type %in% c("observed", "simex.1", "simex.2", "simex.3")) %>%
  mutate(setting="Strong alternative")

sim.2.gaussian.1 <- ret.1.1.est %>%
  filter(type %in% c("observed", "mecor.1", "mecor.2", "mecor.3")) %>%
  mutate(setting="No treatment effect null")
sim.2.gaussian.2 <- ret.1.2.est %>%
  filter(type %in% c("observed", "mecor.1", "mecor.2", "mecor.3")) %>%
  mutate(setting="Weak alternative")
sim.2.gaussian.3 <- ret.1.3.est %>%
  filter(type %in% c("observed", "mecor.1", "mecor.2", "mecor.3")) %>%
  mutate(setting="Moderate alternative")
sim.2.gaussian.4 <- ret.1.4.est %>%
  filter(type %in% c("observed", "mecor.1", "mecor.2", "mecor.3")) %>%
  mutate(setting="Strong alternative")

fig.sim.2.surv.bias.2 <- bind_rows(
  sim.2.surv.1,
  sim.2.surv.2,
  sim.2.surv.3,
  sim.2.surv.4)

fig.sim.2.gaussian.bias.2 <- fig.sim.2.gaussian.bias.2 %>%
  mutate(model=case_when((type == "observed") ~ "observed",
                         (type == "simex.1") ~ "simex (correct)", 
                         (type == "simex.2") ~ "simex (mild)",
                         (type == "simex.3") ~ "simex (severe)"
  ))

fig.sim.2.gaussian.bias.2 <- fig.sim.2.gaussian.bias.2 %>%
  mutate(model=case_when((type == "observed") ~ "observed",
                         (type == "mecor.1") ~ "mecor (correct)", 
                         (type == "mecor.2") ~ "mecor (mild)",
                         (type == "mecor.3") ~ "mecor (severe)"
  ))

fig.sim.2.surv.bias.2$setting <- factor(fig.sim.2.surv.bias.2$setting,
                                        levels = c("No treatment effect null",
                                                   "Weak alternative", "Moderate alternative", "Strong alternative"))
fig.sim.2.surv.bias.2.plot.simex <- fig.sim.2.surv.bias.2%>%
  ggplot(aes(x=bias, colour = model)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme_bw() + 
  theme(legend.position="bottom")


fig.sim.2.gaussian.bias.2 <- bind_rows(
  sim.2.gaussian.1,
  sim.2.gaussian.2,
  sim.2.gaussian.3,
  sim.2.gaussian.4)

fig.sim.2.gaussian.bias.2$setting <- factor(fig.sim.2.gaussian.bias.2$setting,
                                        levels = c("No treatment effect null",
                                                   "Weak alternative", "Moderate alternative", "Strong alternative"))
fig.sim.2.gaussian.bias.2.plot <- fig.sim.2.gaussian.bias.2%>%
  ggplot(aes(x=bias, colour = model)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  # scale_x_continuous(limits = c(-0.05,0.1)) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme_bw() +  
  theme(legend.position="bottom")

ggsave(file="../figures/Fig_me_gaussian_bias_simex.eps", width = 240,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig.sim.2.gaussian.bias.2.plot) #saves g
