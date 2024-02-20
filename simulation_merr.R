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

sim.data <- generate.data(n=n, sigma_eta=2,
                                      beta = list(b_w1=c(0,0.01),
                                                  b_s=0.5,
                                                  b_eta=-0.5))

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S + epsilon, data=sim.data)

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S, data=sim.data)

fit.obs <- coxph(Surv(Y, D) ~ A + W1 + S_obs, data=sim.data)

plot(density(sim.data[sim.data$eta>0, ]$T), col="red")
lines(density(sim.data[sim.data$eta<0, ]$T), col="blue")

################

trans_cpd_q <- function(x){
  if (x<=quantile(S, 0.25)) { 
    return(0)
  } else if (x<=quantile(S, 0.56)) {
    return(1)
  } else if (x<=quantile(S, 0.76)) {
    return(2)
  } else if (x<=quantile(S, 0.89)) {
    return(3)
  } else if (x<=quantile(S, 0.97)) {
    return(4)
  } else {
    return(5)
  }
}

quantile(S, 0.25)
generate.data <- function(n, theta = c(0,-0.2), sigma_eta=1,
                          alpha = list(a_w1=c(0,-0.1),
                                       a_s=-0.05), 
                          beta = list(b_w1=c(0,0.01),
                                      b_s=0.5,
                                      b_eta=-0.1,
                                      b_eta2=0.25,
                                      b_eta3=-0.5)){
  
  # covariates
  W1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
  # S <- rgamma(n, 2, 1)
  S <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, 1)
  # eta <-rnorm(n, epsilon, sigma_eta)
  eta <- epsilon
  # S_report <- floor(S-eta)
  S_report <- S-eta
  S_obs <- factor(unlist(lapply(S_report, trans_cpd_q))+1)
  sim.data <- data.frame(W1=W1, S=S, S_obs=S_obs, eta=eta, epsilon=epsilon)
  
  # treatment 
  glm_a <- with(sim.data, alpha$a_w1[W1]
                + alpha$a_s*S)
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  rate <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                   + beta$b_s*S + beta$b_eta*epsilon))
  
  rate2 <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                   + beta$b_s*S + beta$b_eta2*epsilon))
  
  rate3 <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                   + beta$b_s*S + beta$b_eta3*epsilon))
  # outcome
  T1 <- rexp(n, rate)
  T2 <- rexp(n, rate2)
  T3 <- rexp(n, rate3)
  C <- runif(n, 0.5, 1)
  
  # observed outcome
  Y <- pmin(T, C)
  D <- (T<=C)
  sim.data$T1 <- T1
  sim.data$T2 <- T2
  sim.data$T3 <- T3
  sim.data$C <- C
  
  sim.data$Y1 <- pmin(T1, C)
  sim.data$Y2 <- pmin(T2, C)
  sim.data$Y3 <- pmin(T3, C)
  sim.data$D1 <- (T1<=C)
  sim.data$D2 <- (T2<=C)
  sim.data$D3 <- (T3<=C)
  
  return(sim.data)
}

sim.data <- generate.data(n=n)
plot(survfit(Surv(sim.data[sim.data$eta>0, ]$Y1, sim.data[sim.data$eta>0, ]$D1)~1), 
     col=1, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta>0, ]$Y2, sim.data[sim.data$eta>0, ]$D2)~1), 
      col=2, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta>0, ]$Y3, sim.data[sim.data$eta>0, ]$D3)~1),
      col=3, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta<0, ]$Y1, sim.data[sim.data$eta<0, ]$D1)~1),
      lty=2, col=1, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta<0, ]$Y2, sim.data[sim.data$eta<0, ]$D2)~1),
      lty=2, col=2, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta<0, ]$Y3, sim.data[sim.data$eta<0, ]$D3)~1),
      lty=2, col=3, conf.int = FALSE)

legend(0.6, 0.8, legend=c("b_eps_weak_under", "b_eps_weak_over", 
                       "b_eps_moderate_under", "b_eps_moderate_over", 
                       "b_eps_strong_under", "b_eps_strong_over"),
       col=c(1,1,2,2,3,3), lty=c(1,2,1,2,1,2), cex=0.8)


n <- 1000
library(plyr)
res <- ldply(1:500, function(j) {
  seed <- sample(1e3:1e8, 1)
  set.seed(seed)
  cat(n, j, seed, '\n')
  sim.data <- generate.data(n=n)

  fit.obs <- coxph(Surv(Y1, D1) ~ A + W1 + S_obs, data=sim.data)
  theta.hat.obs.1 <- fit.obs$coefficients[1]
  theta.hat.obs.se.1 <- sqrt(fit.obs$var[1,1])
  
  fit.obs <- coxph(Surv(Y2, D2) ~ A + W1 + S_obs, data=sim.data)
  theta.hat.obs.2 <- fit.obs$coefficients[1]
  theta.hat.obs.se.2 <- sqrt(fit.obs$var[1,1])
  
  fit.obs <- coxph(Surv(Y3, D3) ~ A + W1 + S_obs, data=sim.data)
  theta.hat.obs.3 <- fit.obs$coefficients[1]
  theta.hat.obs.se.3 <- sqrt(fit.obs$var[1,1])
  
  data.frame(
    n=n, j=j, seed = seed,
    est.obs.1 = theta.hat.obs.1, se.obs.1 = theta.hat.obs.se.1,
    est.obs.2 = theta.hat.obs.2, se.obs.2 = theta.hat.obs.se.2,
    est.obs.3 = theta.hat.obs.3, se.obs.3 = theta.hat.obs.se.3)
})


mean(abs(res$est.obs.1 - (-0.2)))
mean(abs(res$est.obs.2 - (-0.2)))
mean(abs(res$est.obs.3 - (-0.2)))

#######################################3
# Jan 16

generate.data <- function(n, theta = c(0,-0.2), sigma_eta=1,
                          alpha = list(a_w1=c(0,-0.1),
                                       a_s=-0.5), 
                          beta = list(b_w1=c(0,0.01),
                                      b_s=0.5,
                                      b_eta=-0.25,
                                      b_eta2=-0.1,
                                      b_eta3=0,
                                      b_eta4=0.1,
                                      b_eta5=0.25)){
  
  # covariates
  W1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
  # S <- rgamma(n, 2, 1)
  S <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, 1)
  eta <-rnorm(n, epsilon, sigma_eta)
  # eta <- epsilon
  # S_report <- floor(S-eta)
  S_report <- S-eta
  S_obs <- factor(unlist(lapply(S_report, trans_cpd_q))+1)
  sim.data <- data.frame(W1=W1, S=S, S_obs=S_obs, eta=eta, epsilon=epsilon)
  
  # treatment 
  glm_a <- with(sim.data, alpha$a_w1[W1]
                + alpha$a_s*S)
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  rate <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                   + beta$b_s*S + beta$b_eta*epsilon))
  
  rate2 <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                    + beta$b_s*S + beta$b_eta2*epsilon))
  
  rate3 <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                    + beta$b_s*S + beta$b_eta3*epsilon))
  
  rate4 <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                    + beta$b_s*S + beta$b_eta4*epsilon))
  
  rate5 <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                    + beta$b_s*S + beta$b_eta5*epsilon))
  
  # outcome
  T1 <- rexp(n, rate)
  T2 <- rexp(n, rate2)
  T3 <- rexp(n, rate3)
  T4 <- rexp(n, rate4)
  T5 <- rexp(n, rate5)
  C <- runif(n, 0.5, 1)
  
  # observed outcome
  Y <- pmin(T, C)
  D <- (T<=C)
  sim.data$T1 <- T1
  sim.data$T2 <- T2
  sim.data$T3 <- T3
  sim.data$T4 <- T4
  sim.data$T5 <- T5
  sim.data$C <- C
  
  sim.data$Y1 <- pmin(T1, C)
  sim.data$Y2 <- pmin(T2, C)
  sim.data$Y3 <- pmin(T3, C)
  sim.data$Y4 <- pmin(T4, C)
  sim.data$Y5 <- pmin(T5, C)
  sim.data$D1 <- (T1<=C)
  sim.data$D2 <- (T2<=C)
  sim.data$D3 <- (T3<=C)
  sim.data$D4 <- (T4<=C)
  sim.data$D5 <- (T5<=C)
  
  return(sim.data)
}

sim.data <- generate.data(n=n)
plot(survfit(Surv(sim.data[sim.data$eta>0, ]$Y1, sim.data[sim.data$eta>0, ]$D1)~1), 
     col=1, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta>0, ]$Y2, sim.data[sim.data$eta>0, ]$D2)~1), 
      col=2, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta>0, ]$Y3, sim.data[sim.data$eta>0, ]$D3)~1),
      col=3, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta<0, ]$Y1, sim.data[sim.data$eta<0, ]$D1)~1),
      lty=2, col=1, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta<0, ]$Y2, sim.data[sim.data$eta<0, ]$D2)~1),
      lty=2, col=2, conf.int = FALSE)
lines(survfit(Surv(sim.data[sim.data$eta<0, ]$Y3, sim.data[sim.data$eta<0, ]$D3)~1),
      lty=2, col=3, conf.int = FALSE)

legend(0.6, 0.8, legend=c("b_eps_weak_under", "b_eps_weak_over", 
                          "b_eps_moderate_under", "b_eps_moderate_over", 
                          "b_eps_strong_under", "b_eps_strong_over"),
       col=c(1,1,2,2,3,3), lty=c(1,2,1,2,1,2), cex=0.8)


n <- 1000
library(plyr)
res <- ldply(1:1000, function(j) {
  seed <- sample(1e3:1e8, 1)
  set.seed(seed)
  cat(n, j, seed, '\n')
  sim.data <- generate.data(n=n)
  
  fit.obs <- coxph(Surv(Y1, D1) ~ A + W1 + S_obs, data=sim.data)
  theta.hat.obs.1 <- fit.obs$coefficients[1]
  theta.hat.obs.se.1 <- sqrt(fit.obs$var[1,1])
  
  fit.obs <- coxph(Surv(Y2, D2) ~ A + W1 + S_obs, data=sim.data)
  theta.hat.obs.2 <- fit.obs$coefficients[1]
  theta.hat.obs.se.2 <- sqrt(fit.obs$var[1,1])
  
  fit.obs <- coxph(Surv(Y3, D3) ~ A + W1 + S_obs, data=sim.data)
  theta.hat.obs.3 <- fit.obs$coefficients[1]
  theta.hat.obs.se.3 <- sqrt(fit.obs$var[1,1])
  
  fit.obs <- coxph(Surv(Y4, D4) ~ A + W1 + S_obs, data=sim.data)
  theta.hat.obs.4 <- fit.obs$coefficients[1]
  theta.hat.obs.se.4 <- sqrt(fit.obs$var[1,1])
  
  fit.obs <- coxph(Surv(Y5, D5) ~ A + W1 + S_obs, data=sim.data)
  theta.hat.obs.5 <- fit.obs$coefficients[1]
  theta.hat.obs.se.5 <- sqrt(fit.obs$var[1,1])
  
  data.frame(
    n=n, j=j, seed = seed,
    est.obs.1 = theta.hat.obs.1, se.obs.1 = theta.hat.obs.se.1,
    est.obs.2 = theta.hat.obs.2, se.obs.2 = theta.hat.obs.se.2,
    est.obs.3 = theta.hat.obs.3, se.obs.3 = theta.hat.obs.se.3,
    est.obs.4 = theta.hat.obs.4, se.obs.4 = theta.hat.obs.se.4,
    est.obs.5 = theta.hat.obs.5, se.obs.5 = theta.hat.obs.se.5)
})


mean(abs(res$est.obs.1 - (-0.2)))
mean(abs(res$est.obs.2 - (-0.2)))
mean(abs(res$est.obs.3 - (-0.2)))
mean(abs(res$est.obs.4 - (-0.2)))
mean(abs(res$est.obs.5 - (-0.2)))
merr.d3.r1.as2 <- res
save(merr.d3.r1.as2, file="merr.d3.r1.as2.RData")


res <- merr.d3.r1.as2

#################
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

## bias distribution of est.obs
plot(density(res$est.obs - (-0.2)))
mean(res$est.obs - (-0.2)) # [1] -0.01060851
mean(res$est.obs - (-0.2)) /0.2

# coverage
mean(res$est.obs - 1.96*res$se.obs <= (-0.2) & (-0.2) <= res$est.obs + 1.96*res$se.obs)
# actually pretty good 0.946

## bias distribution of est.obs
plot(density(res$est.obs - res$est.3.1))
mean(res$est.obs - res$est.3.1) # [1] -0.002467619
plot(density(res$est.obs - res$est.3.2))
mean(res$est.obs - res$est.3.2) # [1] -0.007857603
plot(density(res$est.obs - res$est.3.3))
mean(res$est.obs - res$est.3.3) # [1] -0.01260562
plot(density(res$est.obs - res$est.random))
mean(res$est.obs - res$est.random) # [1] -0.01948067

####################
plot(density(res$est.3.3 - (-0.2)))
mean(res$est.3.3 - (-0.2)) # [1] 0.001997113
plot(density(res$est.random - (-0.2)))
mean(res$est.random - (-0.2)) # [1] 0.008872156
####################
sim.data <- generate.data(n=n)
sim.data$S_obs
cor(as.numeric(sim.data$S_obs), as.numeric(sim.data$A), method = "spearman")
cor(as.numeric(sim.data$S-sim.data$eta), as.numeric(sim.data$A), method = "spearman")
cor(as.numeric(sim.data$S), as.numeric(sim.data$A), method = "spearman")
n <- 1000
library(plyr)
res$seed
res.corr <- ldply(1:500, function(j) {
  seed <- res$seed[j]
  set.seed(seed)
  cat(n, j, seed, '\n')
  sim.data <- generate.data(n=n)
  
  cor.A.S_obs <- cor(as.numeric(sim.data$S_obs), as.numeric(sim.data$A), method = "spearman")
  cor.A.S_report <- cor(as.numeric(sim.data$S-sim.data$eta), as.numeric(sim.data$A), method = "spearman")
  cor.A.S <- cor(as.numeric(sim.data$S), as.numeric(sim.data$A), method = "spearman")

  data.frame(
    n=n, j=j, seed = seed,
    cor.A.S_obs = cor.A.S_obs,
    cor.A.S_report = cor.A.S_report,
    cor.A.S = cor.A.S)
})

density(abs(res.corr$cor.A.S) - abs(res.corr$cor.A.S_obs))
density(res.corr$cor.A.S)
density(res.corr$cor.A.S_obs)
density(res.corr$cor.A.S_report)
density(abs(res.corr$cor.A.S_report) - abs(res.corr$cor.A.S_obs))

var(sim.data$S)
var(sim.data$S-sim.data$eta)
# res$cor.A.S larger than res$cor.A.S_obs a little bit
###################

dimnames(p_ij.random) <- list(levels(sim.data$S_obs), levels(sim.data$S_obs))
dimnames(p_ij.random) <- list(levels(sim.data$S_obs), levels(sim.data$S_obs))

####################

generate.data <- function(n, theta = c(0,-0.2), sigma_eta=2,
                          alpha = list(a_w1=c(0,-0.1),
                                       a_s=-0.05), 
                          beta = list(b_w1=c(0,0.01),
                                      b_s=0.5,
                                      b_eta=-0.02)){
  
  # covariates
  W1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
  S <- rgamma(n, 2, 1)
  epsilon <- rnorm(n, 0, 1)
  eta <-rnorm(n, epsilon, sigma_eta)
  S_report <- floor(S-eta)
  S_obs <- factor(unlist(lapply(S_report, trans_cpd))+1)
  sim.data <- data.frame(W1=W1, S=S, S_obs=S_obs, eta=eta, epsilon=epsilon)
  
  # treatment 
  glm_a <- with(sim.data, alpha$a_w1[W1]
                + alpha$a_s*S)
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  rate <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                   + beta$b_s*S + beta$b_eta*epsilon))
  # outcome
  T <- rexp(n, rate)
  C <- runif(n, 0.5, 1)
  
  # observed outcome
  Y <- pmin(T, C)
  D <- (T<=C)
  sim.data$T <- T
  sim.data$C <- C
  sim.data$Y <- Y
  sim.data$D <- D
  
  return(sim.data)
}

sim.data <- generate.data(n=n)

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S + epsilon, data=sim.data)

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S, data=sim.data)

fit.obs <- coxph(Surv(Y, D) ~ A + W1 + S_obs, data=sim.data)

n <- 1000
library(plyr)
res <- ldply(1:500, function(j) {
  seed <- sample(1e3:1e8, 1)
  set.seed(seed)
  cat(n, j, seed, '\n')
  sim.data <- generate.data(n=n, sigma_eta=4,
                            beta = list(b_w1=c(0,0.01),
                                             b_s=0.5,
                                             b_eta=-0.25))
  
  fit.true <- coxph(Surv(Y, D) ~ A + W1 + S + epsilon, data=sim.data)
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
    n=n, j=j, seed = seed, sigma_eta = 4, b_epsilon = -0.25,
    est.true = theta.hat.true, se.true = theta.hat.true.se,
    est.noeta = theta.hat.noeta, se.noeta = theta.hat.noeta.se,
    est.obs = theta.hat.obs, se.obs = theta.hat.obs.se,
    est.random = theta.hat.mcsimex.random, se.random = theta.hat.mcsimex.se.random,
    est.3.1 = theta.hat.mcsimex.3.1, se.3.3 = theta.hat.mcsimex.se.3.1,
    est.3.2 = theta.hat.mcsimex.3.2, se.3.3 = theta.hat.mcsimex.se.3.2,
    est.3.3 = theta.hat.mcsimex.3.3, se.3.3 = theta.hat.mcsimex.se.3.3)
})

merr.1.1.2.3 <- res
save(merr.1.1.2.3, file="../output.data/merr.1.1.2.3.RData")
#######################
# organize data
#######################
# name the data
# merr.{S_dist}.{g}.{beta_epsilon}.{sigma_eta}
# S_dist - 1: Gamma
# g - 1: floor
# beta_epsilon - 1: -0.1; 2: -0.25; 3: -0.5
# sigma_eta - 1: 1; 2: 2; 3: ? # reliability ratio?

merr.1.1.1.1 <- res
merr.1.1.1.2 <- merr.1.1.1.1
load("res_0114_1.RData")
merr.1.1.3.2 <- res

load("res_0111_2.RData") # merr.1.1.1.2
load("res_0112.RData") # merr.1.1.1.1

# Controlling $\sigma_{\eta}^2$
merr.1.1.1.2 # already have
merr.1.1.2.2 # running
merr.1.1.3.2 # already have

head(merr.1.1.3.2)
res <- merr.1.1.1.2
## bias distribution of est.obs
plot(density(res$est.obs - (-0.2)))
mean(res$est.obs - (-0.2)) # [1] -0.01060851
var(res$est.obs - (-0.2))
mean(res$est.obs - (-0.2)) /0.2

mean(abs(merr.1.1.1.2$est.obs - (-0.2)))
mean(abs(merr.1.1.3.2$est.obs - (-0.2)))
mean(abs(merr.1.1.2.2$est.obs - (-0.2)))

plot(density(merr.1.1.1.2$est.obs - (-0.2)))
lines(density(merr.1.1.3.2$est.obs - (-0.2)), col="red")
hist(merr.1.1.1.2$est.obs - (-0.2))
hist(merr.1.1.3.2$est.obs - (-0.2), breaks='FD', col = "brown", add=T)
# coverage
mean(res$est.obs - 1.96*res$se.obs <= (-0.2) & (-0.2) <= res$est.obs + 1.96*res$se.obs)
# actually pretty good 0.946

## bias distribution of est.obs
plot(density(res$est.obs - res$est.3.1))
mean(abs(res$est.obs - res$est.3.1)) # [1] -0.002467619
plot(density(res$est.obs - res$est.3.2))
mean(res$est.obs - res$est.3.2) # [1] -0.007857603
plot(density(res$est.obs - res$est.3.3))
mean(res$est.obs - res$est.3.3) # [1] -0.01260562
plot(density(res$est.obs - res$est.random))
mean(res$est.obs - res$est.random)


mean(abs(merr.1.1.3.2$est.obs - merr.1.1.3.2$est.3.1))
mean(abs(merr.1.1.3.2$est.obs - merr.1.1.3.2$est.3.2))
mean(abs(merr.1.1.3.2$est.obs - merr.1.1.3.2$est.3.3))
mean(abs(merr.1.1.3.2$est.obs - merr.1.1.3.2$est.random))
# Controlling $\beta_\epsilon$
# beta_epsilon = -0.25
merr.1.1.2.1 # sigma_eta=1
merr.1.1.2.2 # already have
merr.1.1.2.3 # sigma_eta=4











#########################

generate.data <- function(n, theta = c(0,-0.2), sigma_eta=5,
                          alpha = list(a_w1=c(0,-0.1),
                                       a_s=-0.05), 
                          beta = list(b_w1=c(0,0.01),
                                      b_s=0.5,
                                      b_eta=-0.02)){
  
  # covariates
  W1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
  S <- runif(n, 0, 8)
  epsilon <- rnorm(n, 0, 1)
  eta <-rnorm(n, epsilon, sigma_eta)
  S_report <- floor(S-eta)
  S_obs <- factor(unlist(lapply(S_report, trans_cpd))+1)
  sim.data <- data.frame(W1=W1, S=S, S_obs=S_obs, eta=eta, epsilon=epsilon)
  
  # treatment 
  glm_a <- with(sim.data, alpha$a_w1[W1]
                + alpha$a_s*S)
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  rate <- exp(with(sim.data, theta[A] + beta$b_w1[W1] 
                   + beta$b_s*S + beta$b_eta*epsilon))
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

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S + epsilon, data=sim.data)

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S, data=sim.data)

fit.obs <- coxph(Surv(Y, D) ~ A + W1 + S_obs, data=sim.data)

# standardized the variables??

sim.data <- generate.data(n=n)

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S + eta, data=sim.data)

fit.true <- coxph(Surv(Y, D) ~ A + W1 + S, data=sim.data)

fit.obs <- coxph(Surv(Y, D) ~ A + W1 + S_obs, data=sim.data)