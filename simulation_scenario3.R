rm(list=ls())

library(survival)
library(simex)
library(tidyr)
library(dplyr)
library(ggplot2)

generate.data.hybrid <- function(n, alpha=0, beta=0, theta = c(0,1), p_ij=p_ij,
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
  U <- rbinom(n, 1, p=0.5)
  
  sim.data <- data.frame(X1=X1, X2=X2, X3=X3, X3.star=X3.star, U=U)
  
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

p_ij <- matrix(c(0.8, 0.1, 0.01, 0.01, 0.01, 0.01,
                 0.16, 0.8, 0.1, 0.01, 0.01, 0.01,
                 0.01, 0.07, 0.8, 0.1, 0.01, 0.01,
                 0.01, 0.01, 0.07, 0.8, 0.07, 0.01,
                 0.01, 0.01, 0.01, 0.07, 0.8, 0.16,
                 0.01, 0.01, 0.01, 0.01, 0.1, 0.8), nrow = 6, byrow = TRUE)
p_ij.5 <- matrix(c(0.6, 0.14, 0.02, 0.02, 0.02, 0.02,
                   0.3, 0.6, 0.14, 0.02, 0.02, 0.02,
                   0.04, 0.2, 0.6, 0.14, 0.02, 0.02,
                   0.02, 0.02, 0.2, 0.6, 0.2, 0.04,
                   0.02, 0.02, 0.02, 0.2, 0.6, 0.3,
                   0.02, 0.02, 0.02, 0.02, 0.14, 0.6), nrow = 6, byrow = TRUE) 

n <- 2500
sim.data <- generate.data.hybrid(n=n, alpha = 0.5, beta = 0.5, theta = c(0,1), p_ij = p_ij.5,
                          eta_a = list(eta_a1=c(0,-0.1),
                                       eta_a2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_t3=c(0,-0.2,-0.4,-0.6,-0.8,-1)))

dimnames(p_ij) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))
dimnames(p_ij.5) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))

# true
fit.true <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3 + U, data=sim.data)
theta.hat.true <- fit.true$coefficients[1]
theta.hat.true.se <- sqrt(fit.true$var[1,1])

# observed adjusted causal hazard ratio - no unmeasured & residual confounding adjustment
fit.obs <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data)
theta.hat.obs <- fit.obs$coefficients[1]
theta.hat.obs.se <- sqrt(fit.obs$var[1,1])

# only unmeasured confounding adjustment
alpha <- -0.5
beta <- -0.5
X <- cbind(sim.data$X1, sim.data$X2, sim.data$X3.star)
A <- as.numeric(sim.data$A)-1
fit.stoEM_reg <- survSensitivity(sim.data$Y, sim.data$D, A, X, "stoEM_reg",
                                 zetaT = beta, zetaZ = alpha,
                                 B = 5)
theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se

# only residual confounding adjustment
sim.data$X3.star.cont <- as.numeric(sim.data$X3.star)-1
naive.model.continuous <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont, data=sim.data, model = TRUE)
fit.simex.1 <- simex(naive.model.continuous, measurement.error = 0.8, 
                     SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
theta.hat.simex.1 <- fit.simex.1$coefficients[1]
theta.hat.simex.se.1 <- sqrt(fit.simex.1$variance.jackknife[1,1])

alpha <- -0.5
beta <- -0.5
theta <- 1
p_ij.gen <- p_ij.5
ret <- data.frame()
for (n in 2500) {
  for (j in 1:500) {
    
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)
    cat(n, j, seed, '\n')
    
    # generate data
    sim.data <- generate.data.hybrid(n=n, alpha = alpha, beta = beta, theta = c(0,theta), p_ij = p_ij.gen,
                                     eta_a = list(eta_a1=c(0,-0.1),
                                                  eta_a2=c(0,-0.01,-0.02,-0.03,-0.04),
                                                  eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                                     eta_t = list(eta_t1=c(0,-0.1),
                                                  eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                                  eta_t3=c(0,-0.2,-0.4,-0.6,-0.8,-1)))
    
    try({
      # true
      fit.true <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3 + U, data=sim.data)
      theta.hat.true <- fit.true$coefficients[1]
      theta.hat.true.se <- sqrt(fit.true$var[1,1])
      
      # observed adjusted causal hazard ratio - no unmeasured & residual confounding adjustment
      fit.obs <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data)
      theta.hat.obs <- fit.obs$coefficients[1]
      theta.hat.obs.se <- sqrt(fit.obs$var[1,1])
      
      # only unmeasured confounding adjustment
      X <- cbind(sim.data$X1, sim.data$X2, sim.data$X3.star)
      A <- as.numeric(sim.data$A)-1
      fit.stoEM_reg <- survSensitivity(sim.data$Y, sim.data$D, A, X, "stoEM_reg",
                                       zetaT = beta, zetaZ = alpha,
                                       B = 5)
      theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
      theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se
      
      # only residual confounding adjustment
      sim.data$X3.star.cont <- as.numeric(sim.data$X3.star)-1
      naive.model.continuous <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont, data=sim.data, model = TRUE)
      fit.simex.3 <- simex(naive.model.continuous, measurement.error = 0.8, 
                           SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.3 <- fit.simex.3$coefficients[1]
      theta.hat.simex.se.3 <- sqrt(fit.simex.3$variance.jackknife[1,1])
      
      if (length(ret)>0){
        ret.tmp <- data.frame(true = theta.hat.true, true.se = theta.hat.true.se,
                              observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                              stoEM_reg = theta.hat.stoEM_reg, stoEM_reg.se = theta.hat.stoEM_reg.se, 
                              simex.3 = theta.hat.simex.3, simex.se.3 = theta.hat.simex.se.3,
                              n=n, j=j, seed=seed)
        
        ret <- rbind(ret, ret.tmp)
      } else {
        ret <- data.frame(true = theta.hat.true, true.se = theta.hat.true.se,
                          observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                          stoEM_reg = theta.hat.stoEM_reg, stoEM_reg.se = theta.hat.stoEM_reg.se, 
                          simex.3 = theta.hat.simex.3, simex.se.3 = theta.hat.simex.se.3,
                          n=n, j=j, seed=seed)
      }
    })
    
  }
}

# save data format
# sim.{scenario}.{alpha}.{beta}.{theta}.{mis_prob}.{n}.{project}
sim.3.05.05.1.5.2500.pa <- ret
save(sim.3.05.05.1.5.2500.pa, file="../sim.data/sim.3.05.05.1.5.2500.pa.RData")