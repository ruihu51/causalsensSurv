rm(list=ls())

# install packages
library(tidyr)
library(ggplot2)

# Scenario II
# Aim:
# 1) validate whether the observed adjusted effect lie between the true and crude effect


# 1)
# specified misclassification probabilities + monotonicity in eta_a3 + eta_t3
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
gen_obs_given_unobs <- function(x, category){
  j <- which(category==x)
  sample_prob <- mis_prob[,j]
  x.star <- sample(0:(length(category)-1), size = 1, replace = TRUE, prob = sample_prob)
  return(x.star)
}


generate.data <- function(n, theta = c(0,-1), 
                          eta_a = list(eta_a1=c(0,0.1),
                                       eta_a2=c(0,0.01,0.02,0.03,0.04),
                                       eta_a3=c(0,-0.02,-0.04,-0.06,-0.08,-0.1)), 
                          eta_t = list(eta_t1=c(0,-0.1),
                                       eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                       eta_t3=c(0,0.02,0.04,0.06,0.08,0.1))){
  
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
  D <- (T<=C)
  sim.data$Y <- Y
  sim.data$D <- D
  
  return(sim.data)
}


ret <- data.frame()
for (n in 2500) {
  for (j in 1:500) {
    
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)
    cat(n, j, seed, '\n')
    
    # generate data
    sim.data <- generate.data(n=n, theta = c(0,-1), 
                              eta_a = list(eta_a1=c(0,0.1),
                                           eta_a2=c(0,0.01,0.02,0.03,0.04),
                                           eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                              eta_t = list(eta_t1=c(0,-0.1),
                                           eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                           eta_t3=c(0,0.2,0.4,0.6,0.8,1)))
    
    try({
      # observed adjusted causal hazard ratio
      fit.obs <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data)
      theta.hat.obs <- fit.obs$coefficients[1]
      theta.hat.obs.se <- sqrt(fit.obs$var[1,1])
      
      # crude effect
      fit.crude <- coxph(Surv(Y, D) ~ A, data=sim.data)
      theta.hat.crude <- fit.crude$coefficients[1]
      theta.hat.crude.se <- sqrt(fit.crude$var[1,1])
      
      if (length(ret)>0){
        ret.tmp <- data.frame(observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                              crude = theta.hat.crude, crude.se = theta.hat.crude.se, 
                              n=n, j=j, seed=seed)
        
        ret <- rbind(ret, ret.tmp)
      } else {
        ret <- data.frame(observed = theta.hat.obs, observed.se = theta.hat.obs.se, 
                          crude = theta.hat.crude, crude.se = theta.hat.crude.se, 
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
  gather(key='type', value='est', c("observed", "crude", "observed.se", "crude.se")) %>%
  filter(type %in% c("observed","crude")) %>%
  mutate(bias=est-(-1)) %>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash")

# change to a more significant eta_a3 and eta_t3
ret %>%
  gather(key='type', value='est', c("observed", "crude", "observed.se", "crude.se")) %>%
  filter(type %in% c("observed","crude")) %>%
  mutate(bias=est-(-1)) %>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") +
  ggtitle(expression(paste("eta_2")))






