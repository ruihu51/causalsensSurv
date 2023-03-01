rm(list=ls())

# install packages
devtools::install_github("Rong0707/survSens")
library(survSens)
library(survival)
library(plyr)
library(dplyr)

# simulation scenarios
# 1) only one unmeasured confounder U
# 2) only one misclassified ordinal variable X3
# 3) TODO:

# Scenario I
# Aim:
# 1) What is the performance of stoEM when U~Ber(0.5)
# 2) What is the performance of stoEM when U~norm(0,1)
# 3) TODO: what about multinomial?

# generate data

# sample sizes
n <- 2500
alpha <- -1
beta <- 1
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
fit.stoEM_reg <- survSensitivity(sim.data$Y, sim.data$D, A, X, "stoEM_reg",
                                 zetaT = beta, zetaZ = alpha,
                                 B = 5)
theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se

# Scenario II
# misclassification probabilities (K * K matrix)
# p_ij=P(X*=i|X=j)
# meet tapered assumption
mis_prob <- matrix(c(0.8, 0.1, 0, 0, 0, 0,
                     0.2, 0.8, 0.1, 0, 0, 0,
                     0, 0.1, 0.8, 0.1, 0, 0,
                     0, 0, 0.1, 0.8, 0.1, 0,
                     0, 0, 0, 0.1, 0.8, 0.2,
                     0, 0, 0, 0, 0.1, 0.8), nrow = 6, byrow = TRUE)
# generate X* given X based on misclassification probabilities
X3 <- factor(sample(0:5, size = n, replace = TRUE, 
                    prob = c(0.21, 0.30, 0.21, 0.15, 0.10, 0.03)))
var_category <- levels(X3)
gen_obs_given_unobs <- function(x, category){
  j <- which(category==x)
  sample_prob <- mis_prob[,j]
  x.star <- sample(0:(length(category)-1), size = 1, replace = TRUE, prob = sample_prob)
  return(x.star)
}
X3.star <- factor(unlist(lapply(X3, gen_obs_given_unobs, category=var_category)))
# small trick: check correlation between X3 and X3.star
data.frame(table(X3,X3.star)) %>%
  ggplot(aes(x=X3, y=Freq, fill=X3.star)) + geom_bar(stat="identity")

# sample sizes
n <- 2500
theta <- c(0,-1) # treatment effect
# eta_a3 and eta_t3 meet monotonicity assumption
eta_a <- list(eta_a1=c(0,0.1),
              eta_a2=c(0,0.01,0.02,0.03,0.04),
              eta_a3=c(0,-0.02,-0.04,-0.06,-0.08,-0.1))
eta_t <- list(eta_t1=c(0,-0.1),
              eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
              eta_t3=c(0,0.02,0.04,0.06,0.08,0.1))

# covariates
X1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
X2 <- factor(sample(0:4, size = n, replace = TRUE, 
                    prob = c(0.12, 0.2, 0.29, 0.35, 0.04)))
X3 <- factor(sample(0:5, size = n, replace = TRUE, 
                    prob = c(0.21, 0.30, 0.21, 0.15, 0.10, 0.03)))
var_category <- levels(X3)
X3.star <- factor(unlist(lapply(X3, gen_obs_given_unobs, category=var_category)))
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
cat("follow-up years: ", mean(Y), sd(Y)) 
D <- (T<=C)
cat("percentage of deaths: ", mean(D)) 
sim.data$Y <- Y
sim.data$D <- D

# observed adjusted causal hazard ratio
fit.obs <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data)
theta.hat.obs <- fit.obs$coefficients[1]
theta.hat.obs.se <- sqrt(fit.obs$var[1,1])

# crude effect
fit.crude <- coxph(Surv(Y, D) ~ A, data=sim.data)
theta.hat.crude <- fit.crude$coefficients[1]
theta.hat.crude.se <- sqrt(fit.crude$var[1,1])
