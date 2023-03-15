# bias correction methods
# methods summary
# MCSIMEX - simex
# Regression Calibration - mecor


# data
# K=2
p_ij <- matrix(c(0.8,0.1,0.2,0.9), nrow = 2, byrow = TRUE)

generate.data <- function(n, theta = c(0,-1), 
                          eta_a = list(eta_a1=c(0,0.1),
                                       eta_a2=c(0,0.01,0.02,0.03,0.04),
                                       eta_a3=c(0,-0.02,-0.04,-0.06,-0.08,-0.1)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_t3=c(0,0.02,0.04,0.06,0.08,0.1))){
  
  # covariates
  X1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
  X2 <- factor(sample(0:1, size = n, replace = TRUE, 
                      prob = c(0.4, 0.6)))
  X2.star <- factor(unlist(lapply(X2, gen_chd_given_parents, category=levels(X2), prob_matrix=p_ij)))
  sim.data <- data.frame(X1=X1, X2=X2, X2.star=X2.star)
  
  # treatment 
  glm_a <- with(sim.data, eta_a$eta_a1[X1] + eta_a$eta_a2[X2])
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  rate <- exp(with(sim.data, theta[A] + eta_t$eta_t1[X1] + eta_t$eta_t2[X2]))
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



# MCSIMEX
# 1) survival outcome Y, binary A, binary X1 and X2, misclassified X2*
sim.data <- generate.data(n=n, theta = c(0,-1), 
                          eta_a = list(eta_a1=c(0,0.1),
                                       eta_a2=c(0,0.8)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-1.2)))
# estimate
library(simex)
naive.model <- coxph(Surv(Y, D) ~ A + X1 + X2.star, data=sim.data, model = TRUE)
dimnames(p_ij) <- list(levels(sim.data$X2.star), levels(sim.data$X2.star))
fit.mcsimex <- mcsimex(naive.model, mc.matrix = list(X2.star=p_ij), 
                       SIMEXvariable = c("X2.star"), asymptotic = FALSE)
summary(fit.mcsimex)

# K>2
p_ij <- matrix(c(0.8, 0.1, 0, 0, 0, 0,
                 0.2, 0.8, 0.1, 0, 0, 0,
                 0, 0.1, 0.8, 0.1, 0, 0,
                 0, 0, 0.1, 0.8, 0.1, 0,
                 0, 0, 0, 0.1, 0.8, 0.2,
                 0, 0, 0, 0, 0.1, 0.8), nrow = 6, byrow = TRUE)

generate.data <- function(n, theta = c(0,-1), p_ij=p_ij,
                          eta_a = list(eta_a1=c(0,0.1),
                                       eta_a2=c(0,0.01,0.02,0.03,0.04),
                                       eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_t3=c(0,0.2,0.4,0.6,0.8,1))){
  
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

sim.data <- generate.data(n=n, theta = c(0,-1), 
                          eta_a = list(eta_a1=c(0,0.1),
                                       eta_a2=c(0,0.01,0.02,0.03,0.04),
                                       eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_t3=c(0,0.2,0.4,0.6,0.8,1)))
p_ij.star <- matrix(c(0.8, 0.1, 0.01, 0.01, 0.01, 0.01,
                      0.16, 0.8, 0.1, 0.01, 0.01, 0.01,
                      0.01, 0.07, 0.8, 0.1, 0.01, 0.01,
                      0.01, 0.01, 0.07, 0.8, 0.07, 0.01,
                      0.01, 0.01, 0.01, 0.07, 0.8, 0.16,
                      0.01, 0.01, 0.01, 0.01, 0.1, 0.8), nrow = 6, byrow = TRUE)
naive.model <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data, model = TRUE)
dimnames(p_ij) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))
dimnames(p_ij.star) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))
fit.mcsimex <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.star), 
                       SIMEXvariable = c("X3.star"), asymptotic = FALSE)
fit.mcsimex$coefficients[1]

ret <- data.frame()
for (n in 2500) {
  for (j in 1:500) {
    
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)
    cat(n, j, seed, '\n')
    
    # generate data
    sim.data <- generate.data(n=n, theta = c(0,-1), p_ij = p_ij,
                              eta_a = list(eta_a1=c(0,0.1),
                                           eta_a2=c(0,0.01,0.02,0.03,0.04),
                                           eta_a3=c(0,-1.2,-1.4,-1.6,-1.8,-2)), 
                              eta_t = list(eta_t1=c(0,-0.1),
                                           eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                           eta_t3=c(0,1.2,1.4,1.6,1.8,2)))
    
    try({
      # observed adjusted causal hazard ratio
      fit.obs <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data)
      theta.hat.obs <- fit.obs$coefficients[1]
      theta.hat.obs.se <- sqrt(fit.obs$var[1,1])
      
      # crude effect
      fit.crude <- coxph(Surv(Y, D) ~ A, data=sim.data)
      theta.hat.crude <- fit.crude$coefficients[1]
      theta.hat.crude.se <- sqrt(fit.crude$var[1,1])
      
      # mcsimex
      naive.model <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data, model = TRUE)
      fit.mcsimex <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.star), 
                             SIMEXvariable = c("X3.star"), asymptotic = FALSE)
      theta.hat.simex <- fit.mcsimex$coefficients[1]
      
      if (length(ret)>0){
        ret.tmp <- data.frame(observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                              crude = theta.hat.crude, crude.se = theta.hat.crude.se, 
                              simex = theta.hat.simex,
                              n=n, j=j, seed=seed)
        
        ret <- rbind(ret, ret.tmp)
      } else {
        ret <- data.frame(observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                          crude = theta.hat.crude, crude.se = theta.hat.crude.se,
                          simex = theta.hat.simex,
                          n=n, j=j, seed=seed)
      }
    })
    
  }
}

ret %>%
  gather(key='type', value='est', c("observed", "crude", "observed.se", "crude.se", "simex")) %>%
  filter(type %in% c("observed","crude", "simex")) %>%
  mutate(bias=est-(-1)) %>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash")

ret %>%
  gather(key='type', value='est', c("observed", "simex", "mecor")) %>%
  # filter(type %in% c("observed","crude", "simex")) %>%
  mutate(bias=est-(-1)) %>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash")

ret %>%
  gather(key='type', value='est', c("observed", "simex", "simex.1", "mecor")) %>%
  # filter(type %in% c("observed","crude", "simex")) %>%
  mutate(bias=est-(-1)) %>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash")

save(ret, file = "../sim.data/tmp.ret.RData")

# save data format
# sim.{scenario}.{eta}.{mis_prob}.{n}.{project}
# eta-1: eta_a3=c(0,-0.02,-0.04,-0.06,-0.08,-0.1)
# eta-2: eta_t3=c(0,0.2,0.4,0.6,0.8,1)
# update
# sim.{scenario}.{gen}.{eta}.{mis_prob}.{n}.{project}
# gen-null: generate.data (survival outcome)
# gen-gaussian: 

sim.2.1.1.2500.pa <- ret
save(sim.2.1.1.2500.pa, file="../sim.data/sim.2.1.1.2500.pa.RData")
sim.2.2.1.2500.pa <- ret
save(sim.2.2.1.2500.pa, file="../sim.data/sim.2.2.1.2500.pa.RData")
sim.2.2.2.2500.pa <- ret
save(sim.2.2.2.2500.pa, file="../sim.data/sim.2.2.2.2500.pa.RData")
sim.2.2.3.2500.pa <- ret
save(sim.2.2.3.2500.pa, file="../sim.data/sim.2.2.3.2500.pa.RData")
sim.2.3.2.2500.pa <- ret
save(sim.2.3.2.2500.pa, file="../sim.data/sim.2.3.2.2500.pa.RData")
sim.2.gaussian.2.3.2500.pa.1 <- ret
save(sim.2.gaussian.2.3.2500.pa.1, file="../sim.data/sim.2.gaussian.2.3.2500.pa.1.RData")

# 2) continuous outcome Y, binary A, binary X1 and X2, misclassified X2*
generate.data.gaussian <- function(n, theta = c(0,-1), 
                          eta_a = list(eta_a1=c(0,0.1),
                                       eta_a2=c(0,0.01,0.02,0.03,0.04),
                                       eta_a3=c(0,-0.02,-0.04,-0.06,-0.08,-0.1)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_t3=c(0,0.02,0.04,0.06,0.08,0.1))){
  
  # covariates
  X1 <- factor(sample(0:1, size = n, replace = TRUE, prob = c(0.3, 0.7)))
  X2 <- factor(sample(0:1, size = n, replace = TRUE, 
                      prob = c(0.4, 0.6)))
  X2.star <- factor(unlist(lapply(X2, gen_chd_given_parents, category=levels(X2), prob_matrix=p_ij)))
  sim.data <- data.frame(X1=X1, X2=X2, X2.star=X2.star)
  
  # treatment 
  glm_a <- with(sim.data, eta_a$eta_a1[X1] + eta_a$eta_a2[X2])
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  sim.data$mu <- with(sim.data, 0.5+theta[A] + eta_t$eta_t1[X1] + eta_t$eta_t2[X2])
  sim.data$Y <- rnorm(n, mean=sim.data$mu, 
                      sd=0.1)

  return(sim.data)
}
sim.data <- generate.data.gaussian(n=n, theta = c(0,-1), 
                          eta_a = list(eta_a1=c(0,0.1),
                                       eta_a2=c(0,0.8)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-1.2)))
naive.model <- glm(Y ~ A + X1 + X2.star, family = "gaussian", data=sim.data, x=T, y=T)
dimnames(p_ij) <- list(levels(sim.data$X2.star), levels(sim.data$X2.star))
fit.mcsimex <- mcsimex(naive.model, mc.matrix = list(X2.star=p_ij), 
                       SIMEXvariable = c("X2.star"))

glm(Y ~ A + X1 + X2, family = "gaussian", data=sim.data)


# mecor
library(mecor)
generate.data.mecor <- function(n, theta = c(0,-1), 
                                   eta_a = list(eta_a1=c(0,0.1),
                                                eta_a2=c(0,0.01,0.02,0.03,0.04),
                                                eta_a3=c(0,-0.02,-0.04,-0.06,-0.08,-0.1)), 
                                   eta_t = list(eta_t1=c(0,-0.1),
                                                eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                                eta_t3=c(0,0.02,0.04,0.06,0.08,0.1))){
  
  # covariates
  X1 <- runif(n,0,1)
  X2 <- rnorm(n,1,2)
  # X2.star <- X2 + rnorm(n,0,1)
  X2.star <- 0.2 + 0.3*X2 - 1.2*X1 + rnorm(n,0,1)
  sim.data <- data.frame(X1=X1, X2=X2, X2.star=X2.star)
  
  # treatment 
  glm_a <- 0.1 + 0.3*X1 - 0.5*X2
  prob_a <- pnorm(glm_a)
  A <- rbinom(n, 1, prob_a)
  sim.data$A <- factor(A)
  
  # sample T conditional on A and X
  sim.data$mu <- with(sim.data, 0.5 + theta[A] + 0.2*X1 - 0.8*X2)
  sim.data$Y <- rnorm(n, mean=sim.data$mu, 
                      sd=0.1)
  
  return(sim.data)
}
sim.data <- generate.data.mecor(n=n)
sim.data <- generate.data.gaussian(n=n, theta = c(0,-1), 
                                   eta_a = list(eta_a1=c(0,0.1),
                                                eta_a2=c(0,0.8)), 
                                   eta_t = list(eta_t1=c(0,-0.1),
                                                eta_t2=c(0,-1.2)))
naive.model <- glm(Y ~ A + X1 + X2.star, family = "gaussian", data=sim.data)
glm(Y ~ A + X1 + X2.star, family = "gaussian", data=sim.data)
glm(Y ~ A + X1 + X2, family = "gaussian", data=sim.data)
fit.mecor <- mecor(Y~MeasErrorRandom(X2.star, variance = 1)+A+X1, data=sim.data)
measerror.fit <- lm(X2~X2.star+A+X1,data=sim.data)
val.data <- generate.data.mecor(n=n)
measerror.fit <- lm(X2~X2.star+A+X1,data=val.data)

fit.mecor <- mecor(Y~MeasErrorExt(X2.star, model = measerror.fit)+A+X1, data=sim.data)
fit.mecor <- mecor(Y~MeasErrorExt(X2.star, list(coef = c(1.3, 0.6, -1.8, 1)))+A+X1, data=sim.data)

glm(Y ~ A + X1 + as.numeric(X2.star), family = "gaussian", data=sim.data)
glm(Y ~ A + X1 + as.numeric(X2), family = "gaussian", data=sim.data)
fit.mecor <- mecor(Y~MeasErrorRandom(as.numeric(X2.star), variance = 0.14)+A+X1, data=sim.data)
var(as.numeric(sim.data$X2)-as.numeric(sim.data$X2.star))

glm(Y ~ A + X1 + as.numeric(X2), family = "gaussian", data=sim.data)

generate.data.gaussian <- function(n, theta = c(0,-1), p_ij=p_ij,
                                   eta_a = list(eta_a1=c(0,0.1),
                                                eta_a2=c(0,0.01,0.02,0.03,0.04),
                                                eta_a3=c(0,-0.02,-0.04,-0.06,-0.08,-0.1)), 
                                   eta_t = list(eta_t1=c(0,-0.1),
                                                eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                                eta_t3=c(0,0.02,0.04,0.06,0.08,0.1))){
  
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
  sim.data$mu <- with(sim.data, 0.5+theta[A] + eta_t$eta_t1[X1] + eta_t$eta_t2[X2] + eta_t$eta_t3[X3])
  sim.data$Y <- rnorm(n, mean=sim.data$mu, 
                      sd=0.1)
  
  return(sim.data)
}

sim.data <- generate.data.gaussian(n=n, theta = c(0,-1), p_ij = p_ij.star,
                                   eta_a = list(eta_a1=c(0,0.1),
                                                     eta_a2=c(0,0.01,0.02,0.03,0.04),
                                                     eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                                   eta_t = list(eta_t1=c(0,-0.1),
                                                eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                                eta_t3=c(0,0.2,0.4,0.6,0.8,1)))
glm(Y ~ A + X1 + X2 + X3.star, family = "gaussian", data=sim.data)
glm(Y ~ A + X1 + X2 + X3, family = "gaussian", data=sim.data)
fit.mecor <- mecor(Y~MeasErrorRandom(as.numeric(X3.star), variance = 0.46)+A+X1+X2, data=sim.data)
var(as.numeric(sim.data$X3)-as.numeric(sim.data$X3.star))
fit.mecor
naive.model <- glm(Y ~ A + X1 + X2 + X3.star, family = "gaussian", data=sim.data, x=T, y=T)
dimnames(p_ij) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))
fit.mcsimex <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.star), 
                       SIMEXvariable = c("X3.star"))

# 1) continuous Y + continuous X + random measurement error
# the method is actually not stable when the variance is large
# var=1, the estimate is close to the true effect
# var=3, the estimate is biased
# 2) continuous Y + continuous X + external measurement error
# it works when the validation data follows totally the same generation process as the simulation data
# when the measurement error model is misspecified, the estimate is biased
# 3) continuous Y + discrete X
# mecor cannot handle the factor X directly
# treat X as continuous variable
# a) K=2 + random measurement error
# the estimate is close to the true effect
# b) K=6
# the estimate is also close to the true effect, so it is worthwhile doing some comparison between mecor and simex
# 4) mecor not available for survival outcome but
# TODO: maybe can extend mecor to survival outcome
# but for now it is good to compare SIMEX and mecor in continuous Y and discrete X
# TODO: it is actually also interesting to see what if SIMEX misclassified

p_ij.star1 <- matrix(c(0.76, 0.1, 0.02, 0.02, 0.02, 0.02,
                      0.16, 0.75, 0.1, 0.02, 0.02, 0.02,
                      0.02, 0.08, 0.75, 0.1, 0.02, 0.02,
                      0.02, 0.02, 0.08, 0.75, 0.08, 0.02,
                      0.02, 0.02, 0.02, 0.08, 0.75, 0.16,
                      0.02, 0.02, 0.02, 0.02, 0.1, 0.76), nrow = 6, byrow = TRUE)

ret <- data.frame()
for (n in 2500) {
  for (j in 317:500) {
    
    seed <- sample(1e3:1e8, 1)
    # seed <- sim.2.gaussian.2.3.2500.pa[(sim.2.gaussian.2.3.2500.pa$n==n) & (sim.2.gaussian.2.3.2500.pa$j==j), "seed"][1]
    set.seed(seed)
    cat(n, j, seed, '\n')
    
    # generate data
    sim.data <- generate.data.gaussian(n=n, theta = c(0,-1), p_ij = p_ij.star,
                                       eta_a = list(eta_a1=c(0,0.1),
                                                    eta_a2=c(0,0.01,0.02,0.03,0.04),
                                                    eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                                       eta_t = list(eta_t1=c(0,-0.1),
                                                    eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                                    eta_t3=c(0,0.2,0.4,0.6,0.8,1)))
    
    try({
      # observed adjusted causal hazard ratio
      fit.obs <- glm(Y ~ A + X1 + X2 + X3.star, family = "gaussian", data=sim.data)
      theta.hat.obs <- fit.obs$coefficients[2]
      # theta.hat.obs.se <- sqrt(fit.obs$var[1,1])
      
      # mcsimex
      naive.model <- glm(Y ~ A + X1 + X2 + X3.star, family = "gaussian", data=sim.data, x=T, y=T)
      dimnames(p_ij.star) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))
      dimnames(p_ij.star1) <- list(levels(sim.data$X3.star), levels(sim.data$X3.star))
      fit.mcsimex <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.star), 
                             SIMEXvariable = c("X3.star"))
      theta.hat.simex <- fit.mcsimex$coefficients[2]
      fit.mcsimex.mis <- mcsimex(naive.model, mc.matrix = list(X3.star=p_ij.star1), 
                             SIMEXvariable = c("X3.star"))
      theta.hat.simex.mis <- fit.mcsimex.mis$coefficients[2]
      
      # mecor
      var.true <- round(var(as.numeric(sim.data$X3)-as.numeric(sim.data$X3.star)),2)
      fit.mecor <- mecor(Y~MeasErrorRandom(as.numeric(X3.star), variance = var.true)+A+X1+X2, data=sim.data)
      theta.hat.mecor <- fit.mecor$corfit$coef[3]
      
      if (length(ret)>0){
        ret.tmp <- data.frame(observed = theta.hat.obs, simex = theta.hat.simex.mis,
                              simex = theta.hat.simex, mecor = theta.hat.mecor,
                              n=n, j=j, seed=seed)
        
        ret <- rbind(ret, ret.tmp)
      } else {
        ret <- data.frame(observed = theta.hat.obs, simex = theta.hat.simex.mis,
                          simex = theta.hat.simex, mecor = theta.hat.mecor,
                          n=n, j=j, seed=seed)
      }
    })
    
  }
}
save(ret, file = "../sim.data/tmp.ret.RData")




