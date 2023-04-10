rm(list=ls())

# install packages
library(tidyr)
library(ggplot2)

# Scenario I
# Aim:
# 1) What is the performance of stoEM when U~Ber(0.5)
# 2) What is the performance of stoEM when U~norm(0,1)
# 3) TODO: what about multinomial?

# 1) U~Ber(0.5)
# alpha=beta=0, theta=-1, n=2500, lambda0=1
# est: naive, stoEM

generate.data <- function(n, alpha=0, beta=0, theta = c(0,1), 
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
  # U: unmeasured confounder
  U <- rbinom(n, 1, p=0.5) # dist-1
  # U <- rnorm(n, 0, 1) # dist-2
  # U <- runif(n, 0, 1) # dist-3
  sim.data <- data.frame(X1=X1, X2=X2, X3=X3, U=U)
  
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

alpha <-  2
beta <- 2
ret <- data.frame()
for (n in 2500) {
  for (j in 1:500) {
    
    seed <- sample(1e3:1e8, 1)
    set.seed(seed)
    cat(n, j, seed, '\n')
    
    # generate data
    sim.data <- generate.data(n=n, alpha=alpha, beta=beta, theta = c(0,0))
    
    try({
      fit.naive <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3, data=sim.data)
      theta.hat.naive <- fit.naive$coefficients[1]
      theta.hat.naive.se <- sqrt(fit.naive$var[1,1])
    
      # stochastic EM with regression
      X <- cbind(sim.data$X1, sim.data$X2, sim.data$X3)
      A <- as.numeric(sim.data$A)-1
      fit.stoEM_reg <- survSensitivity(sim.data$Y, sim.data$D, A, X, "stoEM_reg",
                                      zetaT = beta, zetaZ = alpha,
                                      B = 5)
      theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
      theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se
      
      if (length(ret)>0){
        ret.tmp <- data.frame(naive = theta.hat.naive, naive.se = theta.hat.naive.se, 
                              stoEM_reg = theta.hat.stoEM_reg, stoEM_reg.se = theta.hat.stoEM_reg.se, 
                              n=n, j=j, seed=seed)
        
        ret <- rbind(ret, ret.tmp)
      } else {
        ret <- data.frame(naive = theta.hat.naive, naive.se = theta.hat.naive.se, 
                          stoEM_reg = theta.hat.stoEM_reg, stoEM_reg.se = theta.hat.stoEM_reg.se, 
                          n=n, j=j, seed=seed)
      }
    })
    
  }
}

# save data format
# sim.{scenario}.{alpha}.{beta}.{Udist}.{n}.{project}
# sim.1.0.0.1.2500
# Udist: 1-Ber(0.5)
sim.1.0.0.1.2500.pa <- ret
save(sim.1.0.0.1.2500.pa, file="../sim.data/sim.1.0.0.1.2500.pa.RData")
sim.1.n1.1.1.2500.pa <- ret
save(sim.1.n1.1.1.2500.pa, file="../sim.data/sim.1.n1.1.1.2500.pa.RData")
sim.1.n2.2.1.2500.pa <- ret
save(sim.1.n2.2.1.2500.pa, file="../sim.data/sim.1.n2.2.1.2500.pa.RData")
sim.1.n1.1.2.2500.pa <- ret
save(sim.1.n1.1.2.2500.pa, file="../sim.data/sim.1.n1.1.2.2500.pa.RData")
sim.1.n1.1.3.2500.pa <- ret
save(sim.1.n1.1.3.2500.pa, file="../sim.data/sim.1.n1.1.3.2500.pa.RData")

sim.1.1.1.1.2500.pa <- ret
save(sim.1.1.1.1.2500.pa, file="../sim.data/sim.1.1.1.1.2500.pa.RData")
sim.1.2.2.1.2500.pa <- ret
save(sim.1.2.2.1.2500.pa, file="../sim.data/sim.1.2.2.1.2500.pa.RData")
sim.1.2.2.1.2500.pa <- ret
save(sim.1.2.2.1.2500.pa, file="../sim.data/sim.1.2.2.1.2500.pa.RData")
sim.1.0.0.1.2500.pa <- ret
save(sim.1.0.0.1.2500.pa, file="../sim.data/sim.1.0.0.1.2500.pa.RData")
sim.1.1.1.2.2500.pa <- ret
save(sim.1.1.1.2.2500.pa, file="../sim.data/sim.1.1.1.2.2500.pa.RData")
sim.1.1.1.3.2500.pa <- sim.1.1.1.1.2500.pa
save(sim.1.1.1.3.2500.pa, file="../sim.data/sim.1.1.1.3.2500.pa.RData")

load(file="../sim.data/sim.1.1.1.1.2500.pa.RData")

# TODO: add simulation when theta=0
sim.1.0.0.1.1.2500.pa <- ret
save(sim.1.0.0.1.1.2500.pa, file="../sim.data/sim.1.0.0.1.1.2500.pa.RData")


# Fig 1a
sim.1.0 <- sim.1.0.0.1.2500.pa %>%
  gather(key='type', value='est', c("naive", "stoEM_reg", "naive.se", "stoEM_reg.se")) %>%
  filter(type %in% c("naive","stoEM_reg")) %>%
  mutate(bias=est-(1),
         setting="No unmeasured confounding")

sim.1.1 <- sim.1.1.1.1.2500.pa %>%
  gather(key='type', value='est', c("naive", "stoEM_reg", "naive.se", "stoEM_reg.se")) %>%
  filter(type %in% c("naive","stoEM_reg")) %>%
  mutate(bias=est-(1),
         setting="Moderate unmeasured confounder")

sim.1.2 <- sim.1.2.2.1.2500.pa %>%
  gather(key='type', value='est', c("naive", "stoEM_reg", "naive.se", "stoEM_reg.se")) %>%
  filter(type %in% c("naive","stoEM_reg")) %>%
  mutate(bias=est-(1),
         setting="Strong unmeasured confounder")

fig.sim.1.ber.bias <- bind_rows(
  sim.1.0,
  sim.1.1,
  sim.1.2)

fig.sim.1.ber.bias$setting <- factor(fig.sim.1.ber.bias$setting,
                                        levels = c("No unmeasured confounding",
                                                   "Moderate unmeasured confounder", 
                                                   "Strong unmeasured confounder"))
fig.sim.1.ber.bias.plot <- fig.sim.1.ber.bias%>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme_bw() + 
  theme(legend.position="bottom")

ggsave(file="../figures/Fig_um_ber_bias.eps", width = 290,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig.sim.1.ber.bias.plot) #saves g

# Fig 1b
sim.1.ber <- sim.1.1.1.1.2500.pa %>%
  gather(key='type', value='est', c("naive", "stoEM_reg", "naive.se", "stoEM_reg.se")) %>%
  filter(type %in% c("naive","stoEM_reg")) %>%
  mutate(bias=est-(1),
         setting="Bernoulli distribution")

sim.1.normal <- sim.1.1.1.2.2500.pa %>%
  gather(key='type', value='est', c("naive", "stoEM_reg", "naive.se", "stoEM_reg.se")) %>%
  filter(type %in% c("naive","stoEM_reg")) %>%
  mutate(bias=est-(1),
         setting="Normal distribution")

sim.1.uniform <- sim.1.1.1.3.2500.pa %>%
  gather(key='type', value='est', c("naive", "stoEM_reg", "naive.se", "stoEM_reg.se")) %>%
  filter(type %in% c("naive","stoEM_reg")) %>%
  mutate(bias=est-(1),
         setting="Uniform distribution")

fig.sim.1.dist <- bind_rows(
  sim.1.ber,
  sim.1.normal,
  sim.1.uniform)

fig.sim.1.dist$setting <- factor(fig.sim.1.dist$setting,
                                     levels = c("Bernoulli distribution",
                                                "Normal distribution", 
                                                "Uniform distribution"))
fig.sim.1.dist.plot <- fig.sim.1.dist%>%
  ggplot(aes(x=bias, colour = type)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme_bw() + 
  theme(legend.position="bottom")

ggsave(file="../figures/Fig_um_dist_bias.eps", width = 290,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig.sim.1.dist.plot) #saves g
