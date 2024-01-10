# rm(list=ls())

library(survival)
library(simex)

# continuous smoking measure
n <- 1000
S <- rgamma(n, 2, 1)
table(round(S))
eta <- rnorm(n, 0, 1)
S_report <- floor(S-eta)
trans_cpd <- function(x){
  if (x<=0) { 
    return(0)
  } else if (x>=5) {
    return(5)
  } else {
    return(x)
  }
}
S_obs <- unlist(lapply(S_report, trans_cpd))
table(S_obs)

generate.data <- function(n, theta = c(0,-0.2), p_ij=p_ij,
                          alpha = list(a_w1=c(0,-0.1),
                                       a_s=-0.05), 
                          beta = list(b_w1=c(0,0.01),
                                      b_s=0.5,
                                      b_eta=-0.02)){
  
  # covariates
  W1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
  S <- rgamma(n, 2, 1)
  eta <- rnorm(n, 0, 1)
  S_report <- floor(S-eta)
  S_obs <- factor(unlist(lapply(S_report, trans_cpd))+1)
  sim.data <- data.frame(W1=W1, S=S, S_obs=S_obs, eta=eta)
  
  # treatment 
  glm_a <- with(sim.data, alpha$a_w1[W1]
                + alpha$a_s*S)
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  rate <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                   + beta$b_s*S + beta$b_eta*eta))
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

sim.data <- generate.data(n=n)

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S + eta, data=sim.data)

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S, data=sim.data)

fit.obs <- coxph(Surv(Y, D) ~ A + W1 + S_obs, data=sim.data)



table(sim.data$S_obs, floor(sim.data$S))
table(floor(sim.data$S))
p_ij_emp <- prop.table(table(sim.data$S_obs, floor(sim.data$S)), margin = 2)
p_ij_emp <- rbind(p_ij_emp, rep(0,5))
p_ij_emp <- cbind(p_ij_emp, c(0.01,0.01,0.01,0.01,0.01,0.95))
dimnames(p_ij_emp) <- list(seq(0,5), seq(0,5))
p_ij.1.1 <- build.mc.matrix(p_ij_emp, method = "jlt") # do it twice
check.mc.matrix(list(p_ij.1.1))
p_ij.1.1 <- build.mc.matrix(p_ij.1.1, method = "jlt")

sim.data$S_obs <- factor(sim.data$S_obs)
naive.model <- coxph(Surv(Y, D) ~ A + W1 + S_obs, data=sim.data, model = TRUE)
fit.mcsimex.1 <- mcsimex(naive.model, mc.matrix = list(S_obs=p_ij_emp), 
                         SIMEXvariable = c("S_obs"), asymptotic = FALSE)

fit.mcsimex.1 <- mcsimex(naive.model, mc.matrix = list(S_obs=p_ij.1.1), 
                         SIMEXvariable = c("S_obs"), asymptotic = FALSE)

fit.mcsimex.random <- mcsimex(naive.model, mc.matrix = list(S_obs=p_ij.random), 
                         SIMEXvariable = c("S_obs"), asymptotic = FALSE)


load("p_ij.random.RData")
load("../output.data/p_ij.3.1.RData")
load("../output.data/p_ij.3.2.RData")
load("../output.data/p_ij.3.3.RData")

n <- 1000
library(plyr)
res <- ldply(1:500, function(j) {
  seed <- sample(1e3:1e8, 1)
  set.seed(seed)
  cat(n, j, seed, '\n')
  sim.data <- generate.data(n=n)
  
  fit.true <- coxph(Surv(Y, D) ~ A + W1 + S + eta, data=sim.data)
  theta.hat.true <- fit.true$coefficients[1]
  theta.hat.true.se <- sqrt(fit.true$var[1,1])
  
  fit.noeta <- coxph(Surv(Y, D) ~ A + W1 + S, data=sim.data)
  theta.hat.noeta <- fit.noeta$coefficients[1]
  theta.hat.noeta.se <- sqrt(fit.noeta$var[1,1])
  
  fit.obs <- coxph(Surv(Y, D) ~ A + W1 + S_obs, data=sim.data)
  theta.hat.obs <- fit.obs$coefficients[1]
  theta.hat.obs.se <- sqrt(fit.obs$var[1,1])
  
  naive.model <- coxph(Surv(Y, D) ~ A + W1 + S_obs, data=sim.data, model = TRUE)
  
  fit.mcsimex.random <- mcsimex(naive.model, mc.matrix = list(S_obs=p_ij.random), 
                                SIMEXvariable = c("S_obs"), asymptotic = FALSE)
  theta.hat.mcsimex.random <- fit.mcsimex.random$coefficients[1]
  theta.hat.mcsimex.se.random  <- sqrt(fit.mcsimex.random$variance.jackknife[1,1])
  
  fit.mcsimex.3.1 <- mcsimex(naive.model, mc.matrix = list(S_obs=p_ij.3.1), 
                                SIMEXvariable = c("S_obs"), asymptotic = FALSE)
  theta.hat.mcsimex.3.1 <- fit.mcsimex.3.1$coefficients[1]
  theta.hat.mcsimex.se.3.1  <- sqrt(fit.mcsimex.3.1$variance.jackknife[1,1])
  
  fit.mcsimex.3.2 <- mcsimex(naive.model, mc.matrix = list(S_obs=p_ij.3.2), 
                             SIMEXvariable = c("S_obs"), asymptotic = FALSE)
  theta.hat.mcsimex.3.2 <- fit.mcsimex.3.2$coefficients[1]
  theta.hat.mcsimex.se.3.2  <- sqrt(fit.mcsimex.3.2$variance.jackknife[1,1])
  
  fit.mcsimex.3.3 <- mcsimex(naive.model, mc.matrix = list(S_obs=p_ij.3.3), 
                             SIMEXvariable = c("S_obs"), asymptotic = FALSE)
  theta.hat.mcsimex.3.3 <- fit.mcsimex.3.3$coefficients[1]
  theta.hat.mcsimex.se.3.3  <- sqrt(fit.mcsimex.3.3$variance.jackknife[1,1])
  
  data.frame(
    n=n, j=j, seed = seed,
    est.true = theta.hat.true, se.true = theta.hat.true.se,
    est.noeta = theta.hat.noeta, se.noeta = theta.hat.noeta.se,
    est.obs = theta.hat.obs, se.obs = theta.hat.obs.se,
    est.random = theta.hat.mcsimex.random, se.random = theta.hat.mcsimex.se.random,
    est.3.1 = theta.hat.mcsimex.3.1, se.3.3 = theta.hat.mcsimex.se.3.1,
    est.3.2 = theta.hat.mcsimex.3.2, se.3.3 = theta.hat.mcsimex.se.3.2,
    est.3.3 = theta.hat.mcsimex.3.3, se.3.3 = theta.hat.mcsimex.se.3.3)
})

#######################
#

hist(res$est.true, breaks='FD', col = "blue")
hist(res$est.noeta, breaks='FD', add=T)
hist(res$est.obs, breaks='FD', col = "brown", add=T)

hist(res$est.obs, breaks='FD', col = "brown")
hist(res$est.3.3, breaks='FD', col = "yellow")
hist(res$est.random, breaks='FD', col = "green", add=T)
hist(res$est.3.1, breaks='FD', col = "yellow", add=T)
hist(res$est.3.2, breaks='FD', col = "yellow", add=T)
hist(res$est.3.3, breaks='FD', col = "yellow", add=T)

hist(res$est.3.3 - res$est.obs, freq = FALSE)
hist(res$est.3.2 - res$est.obs, add=T, freq = FALSE)
hist(res$est.3.1 - res$est.obs, add=T, freq = FALSE)
hist(res$est.random - res$est.obs, add=T, freq = FALSE)

hist(res$est.true - res$est.obs)
res$est.obs


####################

dimnames(p_ij.random) <- list(levels(sim.data$S_obs), levels(sim.data$S_obs))
dimnames(p_ij.random) <- list(levels(sim.data$S_obs), levels(sim.data$S_obs))
p_ij.5 <- matrix(c(0.6, 0.14, 0.02, 0.02, 0.02, 0.02,
                   0.3, 0.6, 0.14, 0.02, 0.02, 0.02,
                   0.04, 0.2, 0.6, 0.14, 0.02, 0.02,
                   0.02, 0.02, 0.2, 0.6, 0.2, 0.04,
                   0.02, 0.02, 0.02, 0.2, 0.6, 0.3,
                   0.02, 0.02, 0.02, 0.02, 0.14, 0.6), nrow = 6, byrow = TRUE) # why 0.7,0.2 doesn't work

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