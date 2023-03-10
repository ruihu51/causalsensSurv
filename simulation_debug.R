rm(list=ls())

# install packages
devtools::install_github("Rong0707/survSens")
library(survSens)
library(survival)
library(plyr)
library(dplyr)

####################
# debug area
####################
set.seed(34087494) # no idea
# sample sizes
n <- 2500
alpha <- 0
beta <- 0
theta <- c(0,-1) # treatment effect
eta_a <- list(eta_a1=c(0,0.1),
              eta_a2=c(0,0.01,0.02,0.03,0.04),
              eta_a3=c(0,-0.02,-0.04,-0.06,-0.08,-0.1))
eta_t <- list(eta_t1=c(0,-0.1),
              eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
              eta_t3=c(0,0.02,0.04,0.06,0.08,0.1))

# covariates
# X1: binary - eg. gender
X1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
# X2: ordinal - eg. age
X2 <- factor(sample(0:4, size = n, replace = TRUE, 
                    prob = c(0.12, 0.2, 0.29, 0.35, 0.04)))
# X3: possible mismeasured confounders - eg. smoking behavior - CPDs
X3 <- factor(sample(0:5, size = n, replace = TRUE, 
                    prob = c(0.21, 0.30, 0.21, 0.15, 0.10, 0.03)))
# U: unmeasured confounder
U <- rbinom(n, 1, p=0.5)
sim.data <- data.frame(X1=X1, X2=X2, X3=X3, U=U)

# treatment - physical activity
# TODO: for now binary
glm_a <- with(sim.data, eta_a$eta_a1[X1] + eta_a$eta_a2[X2] 
              + eta_a$eta_a3[X3] + alpha*U)
prob_a <- pnorm(glm_a)
A <- rbinom(n, 1, prob_a)
sim.data$A <- factor(A)

# sample T conditional on A and X
rate <- exp(with(sim.data, theta[A] + eta_t$eta_t1[X1] + eta_t$eta_t2[X2] 
                 + eta_t$eta_t3[X3] + beta*U))
# outcome
# survival time T
# censoring time C
T <- rexp(n, rate)
C <- runif(n, 0.5, 1)

# observed outcome
# observed event time
# censoring indicator
Y <- pmin(T, C)
cat("follow-up years: ", mean(Y), sd(Y)) 
D <- (T<=C)
cat("percentage of deaths: ", mean(D)) 
sim.data$Y <- Y
sim.data$D <- D


# naive causal hazard ratio
fit.naive <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3, data=sim.data)
theta.hat.naive <- fit.naive$coefficients[1]
theta.hat.naive.se <- sqrt(fit.naive$var[1,1])

# stochastic EM with regression
# cannot handle factor A but can handle factor X?
X <- cbind(sim.data$X1, sim.data$X2, sim.data$X3)
A <- as.numeric(sim.data$A)-1
Y <- sim.data$Y
D <- as.numeric(sim.data$D)
fit.stoEM_reg <- survSensitivity(Y, D, A, X, "stoEM_reg",
                                 zetaT = beta, zetaZ = alpha,
                                 B = 5)
coxph(Surv(Y, D) ~ A + X + offset(rep(0,n)), data=sim.data)
sim.data
theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se
