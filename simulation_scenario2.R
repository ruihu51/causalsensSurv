rm(list=ls())

# install packages
library(tidyr)
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

sim.data <- generate.data(n=n, theta = c(0,0.5), p_ij = p_ij.4,
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

ret <- data.frame()
for (n in 2500) {
  for (j in 19:500) {
    
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)
    cat(n, j, seed, '\n')
    
    # generate data
    sim.data <- generate.data(n=n, theta = c(0,0.5), p_ij = p_ij,
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

# plot
# check the bias distribution
ret %>%
  gather(key='type', value='est', colnames(ret)[1:22]) %>%
  filter(type %in% c("true", "observed","simex.1", "simex.2", "simex.3", "simex.4",
                     "mcsimex.1", "mcsimex.2", "mcsimex.3", "mcsimex.4")) %>%
  mutate(bias=est-(0.5)) %>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash")

name <- colnames(ret)[1:22]
name[name %in% c("n", "j")]