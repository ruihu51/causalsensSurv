rm(list=ls())

library(survival)
library(survSens)
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
alpha <- 0.5
beta <- 0.5
X <- cbind(sim.data$X1, sim.data$X2, sim.data$X3.star)
A <- as.numeric(sim.data$A)-1
fit.stoEM_reg <- survSensitivity(sim.data$Y, sim.data$D, A, X, "stoEM_reg",
                                 zetaT = beta, zetaZ = alpha,
                                 B = 5)
theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se

#
X <- cbind(sim.data$X1, sim.data$X2, sim.data$X3.star.cont)

# only residual confounding adjustment
sim.data$X3.star.cont <- as.numeric(sim.data$X3.star)-1
naive.model.continuous <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont, data=sim.data, model = TRUE)
fit.simex.1 <- simex(naive.model.continuous, measurement.error = 0.8, 
                     SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
theta.hat.simex.1 <- fit.simex.1$coefficients[1]
theta.hat.simex.se.1 <- sqrt(fit.simex.1$variance.jackknife[1,1])

# adjusted both
# source("https://raw.githubusercontent.com/Rong0707/survSens/master/R/SimulateU_surv.R")
sim.data$Usim <- SimulateU_surv(sim.data$Y, sim.data$D, A, X, zetat = beta, zetaz = alpha, theta = 0.5, offset = TRUE)$U 
# SimulateU_surv include iter=20 but not B=5
naive.model.U <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont + Usim, data=sim.data, model = TRUE) # simex cannot take offset
fit.simex.U.3 <- simex(naive.model.U, measurement.error = 0.8, 
                     SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)

# add B
B <- 5
theta.hat.simex.UB.3 <- c()
for (i in 1:B){
  sim.data$Usim <- SimulateU_surv(sim.data$Y, sim.data$D, A, X, zetat = beta, zetaz = alpha, theta = 0.5, offset = TRUE)$U 
  # SimulateU_surv include iter=20 but not B=5
  naive.model.U <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont + Usim, data=sim.data, model = TRUE) # simex cannot take offset
  fit.simex.U.3 <- simex(naive.model.U, measurement.error = 0.8, 
                         SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
  theta.hat.simex.UB.3[i] <- fit.simex.U.3$coefficients[1]
}
############
alpha <- 0.5
beta <- 0.5
theta <- 1
p_ij.gen <- p_ij.5
seed.data <- sim.3.05.05.1.5.2500.pa
ret <- data.frame()
for (n in 2500) {
  for (j in 1:500) {
    seed <- seed.data[seed.data$j==j,"seed"]
    if (identical(seed, integer(0))){
      seed <- sample(1e3:1e8, 1)
    }
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
      # # true
      # fit.true <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3 + U, data=sim.data)
      # theta.hat.true <- fit.true$coefficients[1]
      # theta.hat.true.se <- sqrt(fit.true$var[1,1])
      # 
      # # observed adjusted causal hazard ratio - no unmeasured & residual confounding adjustment
      # fit.obs <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star, data=sim.data)
      # theta.hat.obs <- fit.obs$coefficients[1]
      # theta.hat.obs.se <- sqrt(fit.obs$var[1,1])
      # 
      # # only unmeasured confounding adjustment
      # X <- cbind(sim.data$X1, sim.data$X2, sim.data$X3.star)
      # A <- as.numeric(sim.data$A)-1
      # fit.stoEM_reg <- survSensitivity(sim.data$Y, sim.data$D, A, X, "stoEM_reg",
      #                                  zetaT = beta, zetaZ = alpha,
      #                                  B = 5)
      # theta.hat.stoEM_reg <- fit.stoEM_reg$tau1
      # theta.hat.stoEM_reg.se <- fit.stoEM_reg$tau1.se
      # 
      # # only residual confounding adjustment
      sim.data$X3.star.cont <- as.numeric(sim.data$X3.star)-1
      # naive.model.continuous <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont, data=sim.data, model = TRUE)
      # fit.simex.3 <- simex(naive.model.continuous, measurement.error = 0.8, 
      #                      SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      # theta.hat.simex.3 <- fit.simex.3$coefficients[1]
      # theta.hat.simex.se.3 <- sqrt(fit.simex.3$variance.jackknife[1,1])
      
      # adjusted both
      sim.data$Usim <- SimulateU_surv(sim.data$Y, sim.data$D, A, X, zetat = beta, zetaz = alpha, theta = 0.5, offset = TRUE)$U 
      # SimulateU_surv include iter=20 but not B=5
      naive.model.U <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3.star.cont + Usim, data=sim.data, model = TRUE) # simex cannot take offset
      fit.simex.U.3 <- simex(naive.model.U, measurement.error = 0.8, 
                             SIMEXvariable = c("X3.star.cont"), asymptotic = FALSE)
      theta.hat.simex.U.3 <- fit.simex.U.3$coefficients[1]
      theta.hat.simex.se.U.3 <- sqrt(fit.simex.U.3$variance.jackknife[1,1])
      
      # 
      
      if (length(ret)>0){
        ret.tmp <- data.frame(
                              simex.U.3 = theta.hat.simex.U.3, simex.se.U.3 = theta.hat.simex.se.U.3,
                              n=n, j=j, seed=seed)
        
        ret <- rbind(ret, ret.tmp)
      } else {
        ret <- data.frame( 
                          simex.U.3 = theta.hat.simex.U.3, simex.se.U.3 = theta.hat.simex.se.U.3,
                          n=n, j=j, seed=seed)
      }
    })
    
  }
}

library(parallel)

tm <- system.time(
  ret <- mclapply(1:500, sim.3.hybrid, n_range=2500,
                                      seed.data=sim.3.05.05.1.5.2500.pa,
                                      alpha=0.1, beta=0.1,
                                      theta=1, p_ij=p_ij.5, B=1,
                                      mc.cores = 4)
)
sim.3.hybrid(n_range=2500, j_range=1:2, seed.data=sim.3.05.05.1.5.2500.pa,
             alpha=0.1, beta=0.1,
             theta=0.1, p_ij=p_ij.5, B=20)
load("../sim.data/sim.3.05.05.1.5.2500.pa.RData")
############
tm <- system.time(
  sim.3.01.01.1.5.200 <- mclapply(1:200, sim.3.hybrid, n_range=2500,
                    seed.data=sim.3.05.05.1.5.2500.pa,
                    alpha=0.1, beta=0.1,
                    theta=0, p_ij=p_ij.5, B=20,
                    mc.cores = 4)
)
save(sim.3.01.01.1.5.200, file="../sim.data/sim.3.01.01.1.5.200.RData")
tm
############
# check data
ret <- do.call(rbind.data.frame, sim.3.02.02.3.0.100)
theta <- 0.5
data.est <- ret %>%
  gather(key='type', value='est', c("true", "observed", 
                                    "stoEM_reg", "simex.3", 
                                    "simex.UB.3")) %>%
  filter(type %in% c("true", "observed", 
                     "stoEM_reg", "simex.3",
                     "simex.UB.3")) %>%
  mutate(bias=est-(theta))

# how about CI coverage rate?
data.se <- ret %>%
  gather(key='type', value='se', c("true.se", "observed.se", 
                                   "stoEM_reg.se", "simex.se.3", 
                                   "simex.se.UB.3")) %>%
  filter(type %in% c("true.se", "observed.se", 
                     "stoEM_reg.se", "simex.se.3", 
                     "simex.se.UB.3")) 
data <- cbind(data.est, data.se)[-c(6,7,8)]
# head(cbind(data.est, data.se))

data %>%
  mutate(ll=est-qnorm(1-(1-0.95)/2)*se,
         ul=est+qnorm(1-(1-0.95)/2)*se) %>%
  mutate(coverage=(ll <= theta & theta <= ul)) %>%
  group_by(type) %>%
  dplyr::summarise(est.mean = mean(est),
                   bias.mean = mean(abs(bias)),
                   coverage = mean(coverage))

###########
# figures
##########
# 0. integrate data
load("../sim.data/sim.3.01.01.3.5.100.RData")
load("../sim.data/sim.3.01.01.2.5.200.RData")

d1 <- do.call(rbind.data.frame, sim.3.01.01.3.5.100)
d2 <- do.call(rbind.data.frame, sim.3.01.01.3.5.200)
sim.3.01.01.3.5.2500.pa <- bind_rows(d1, d2)
save(sim.3.01.01.1.5.2500.pa, file="../sim.data/sim.3.01.01.1.5.2500.pa.RData")

sim.3.01.01.1.5.2500.pa <- do.call(rbind.data.frame, sim.3.01.01.1.5.200)

# 1. b)
for (setting in file_keys[1:3]){
  index <- as.character(as.integer(file_par_list[[setting]][1]))
  theta <- file_par_list[[setting]][2]
  file_name = paste("../sim.data/sim.3.02.02.", index, ".5.2500.pa.RData", sep="")
  data_name = paste("sim.3.02.02.", index, ".5.2500.pa", sep="")
  load(file = file_name)
  ret <- get(data_name)
  data.est <- ret %>%
    gather(key='type', value='est', c("true", "observed", 
                                      "stoEM_reg", "simex.3", 
                                      "simex.UB.3")) %>%
    filter(type %in% c("true", "observed", 
                       "stoEM_reg", "simex.3",
                       "simex.UB.3")) %>%
    mutate(bias=est-(theta)) %>%
    mutate(setting = setting)
  data.se <- ret %>%
    gather(key='type', value='se', c("true.se", "observed.se", 
                                     "stoEM_reg.se", "simex.se.3", 
                                     "simex.se.UB.3")) %>%
    filter(type %in% c("true.se", "observed.se", 
                       "stoEM_reg.se", "simex.se.3", 
                       "simex.se.UB.3"))
  ret_name <- paste("ret.", index, sep="")
  data <- cbind(data.est, data.se)[-c(6,7,8)]
  assign(ret_name, data)
}
summary.data <- bind_rows(
  ret.1,
  ret.2,
  ret.3)

summary.data$setting <- factor(summary.data$setting,
                               levels = c("No treatment effect",
                                          "Weak", "Moderate"))

fig.sim.3.bias.2.plot <- summary.data %>%
  filter(type %in% c("observed", "stoEM_reg", "simex.3", 
                     "simex.UB.3")) %>%
  mutate(model=case_when((type == "observed") ~ "observed",
                         (type == "stoEM_reg") ~ "est_unmeasured", 
                         (type == "simex.3") ~ "est_residual",
                         (type == "simex.UB.3") ~ "est_hybrid")) %>%
  mutate(model=factor(model, levels=c("observed", "est_unmeasured", "est_residual", "est_hybrid"))) %>%
  ggplot(aes(x=bias, colour = model)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  # scale_x_continuous(limits = c(-0.05,0.1)) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme_bw() +  
  theme(legend.position="bottom")

ggsave(file="../figures/Fig_hybrid_bias_2.eps", width = 290,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig.sim.3.bias.2.plot)

ret.3 %>%
  mutate(ll=est-qnorm(1-(1-0.95)/2)*se,
         ul=est+qnorm(1-(1-0.95)/2)*se) %>%
  mutate(coverage=(ll <= 0.5 & 0.5 <= ul)) %>%
  group_by(type) %>%
  summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)),
            coverage = mean(coverage))

############
# appendix
###########
# c)
for (setting in file_keys[1:3]){
  index <- as.character(as.integer(file_par_list[[setting]][1]))
  theta <- file_par_list[[setting]][2]
  file_name = paste("../sim.data/sim.3.02.02.", index, ".0.2500.pa.RData", sep="")
  data_name = paste("sim.3.02.02.", index, ".0.2500.pa", sep="")
  load(file = file_name)
  ret <- get(data_name)
  data.est <- ret %>%
    gather(key='type', value='est', c("true", "observed", 
                                      "stoEM_reg", "simex.3", 
                                      "simex.UB.3")) %>%
    filter(type %in% c("true", "observed", 
                       "stoEM_reg", "simex.3",
                       "simex.UB.3")) %>%
    mutate(bias=est-(theta)) %>%
    mutate(setting = setting)
  data.se <- ret %>%
    gather(key='type', value='se', c("true.se", "observed.se", 
                                     "stoEM_reg.se", "simex.se.3", 
                                     "simex.se.UB.3")) %>%
    filter(type %in% c("true.se", "observed.se", 
                       "stoEM_reg.se", "simex.se.3", 
                       "simex.se.UB.3"))
  ret_name <- paste("ret.", index, sep="")
  data <- cbind(data.est, data.se)[-c(6,7,8)]
  assign(ret_name, data)
}
summary.data <- bind_rows(
  ret.1,
  ret.2,
  ret.3)

summary.data$setting <- factor(summary.data$setting,
                               levels = c("No treatment effect",
                                          "Weak", "Moderate"))

fig.sim.3.bias.3.plot <- summary.data %>%
  filter(type %in% c("observed", "stoEM_reg", "simex.3", 
                     "simex.UB.3")) %>%
  mutate(model=case_when((type == "observed") ~ "observed",
                         (type == "stoEM_reg") ~ "est_unmeasured", 
                         (type == "simex.3") ~ "est_residual",
                         (type == "simex.UB.3") ~ "est_hybrid")) %>%
  mutate(model=factor(model, levels=c("observed", "est_unmeasured", "est_residual", "est_hybrid"))) %>%
  ggplot(aes(x=bias, colour = model)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  # scale_x_continuous(limits = c(-0.05,0.1)) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme_bw() +  
  theme(legend.position="bottom")

ggsave(file="../figures/Fig_hybrid_bias_3.eps", width = 290,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig.sim.3.bias.3.plot)

# a)
for (setting in file_keys[1:3]){
  index <- as.character(as.integer(file_par_list[[setting]][1]))
  theta <- file_par_list[[setting]][2]
  file_name = paste("../sim.data/sim.3.01.01.", index, ".5.2500.pa.RData", sep="")
  data_name = paste("sim.3.01.01.", index, ".5.2500.pa", sep="")
  load(file = file_name)
  ret <- get(data_name)
  data.est <- ret %>%
    gather(key='type', value='est', c("true", "observed", 
                                      "stoEM_reg", "simex.3", 
                                      "simex.UB.3")) %>%
    filter(type %in% c("true", "observed", 
                       "stoEM_reg", "simex.3",
                       "simex.UB.3")) %>%
    mutate(bias=est-(theta)) %>%
    mutate(setting = setting)
  data.se <- ret %>%
    gather(key='type', value='se', c("true.se", "observed.se", 
                                     "stoEM_reg.se", "simex.se.3", 
                                     "simex.se.UB.3")) %>%
    filter(type %in% c("true.se", "observed.se", 
                       "stoEM_reg.se", "simex.se.3", 
                       "simex.se.UB.3"))
  ret_name <- paste("ret.", index, sep="")
  data <- cbind(data.est, data.se)[-c(6,7,8)]
  assign(ret_name, data)
}
summary.data <- bind_rows(
  ret.1,
  ret.2,
  ret.3)

summary.data$setting <- factor(summary.data$setting,
                               levels = c("No treatment effect",
                                          "Weak", "Moderate"))

fig.sim.3.bias.1.plot <- summary.data %>%
  filter(type %in% c("observed", "stoEM_reg", "simex.3", 
                     "simex.UB.3")) %>%
  mutate(model=case_when((type == "observed") ~ "observed",
                         (type == "stoEM_reg") ~ "est_unmeasured", 
                         (type == "simex.3") ~ "est_residual",
                         (type == "simex.UB.3") ~ "est_hybrid")) %>%
  mutate(model=factor(model, levels=c("observed", "est_unmeasured", "est_residual", "est_hybrid"))) %>%
  ggplot(aes(x=bias, colour = model)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="blue", linetype = "longdash") + 
  facet_grid(~ setting) +
  # scale_x_continuous(limits = c(-0.05,0.1)) +
  labs(x="Bias", y = "Bias distribution density") + 
  theme_bw() +  
  theme(legend.position="bottom")

ggsave(file="../figures/Fig_hybrid_bias_1.eps", width = 290,
       height = 100, units="mm", device=cairo_ps, limitsize = FALSE, fig.sim.3.bias.1.plot)




############
tm <- system.time(
  ret.2 <- mclapply(1:500, sim.3.hybrid, n_range=2500,
                  seed.data=sim.3.05.05.1.5.2500.pa,
                  alpha=0.2, beta=0.2,
                  theta=1, p_ij=p_ij.5, B=1,
                  mc.cores = 4)
)
save(ret.2, file="../sim.data/ret.2.RData")
sim.3.hybrid.B(n_range=2500, j_range=1:2, seed.data=sim.3.02.02.4.5.2500.pa,
               alpha=0.2, beta=0.2,
               theta=1, p_ij=p_ij.5, B=40)
tm <- system.time(
  ret.2.B1 <- mclapply(5:100, sim.3.hybrid.B, n_range=2500,
                    seed.data=sim.3.02.02.4.5.2500.pa,
                    alpha=0.2, beta=0.2,
                    theta=1, p_ij=p_ij.5, B=20,
                    mc.cores = 4)
)
save(ret.2.B1, file="../sim.data/ret.2.B1.RData")

d1 <- do.call(rbind.data.frame, ret.2.B)
d2 <- do.call(rbind.data.frame, ret.2.B1)
sim.3.02.02.4.5.2500.pa.B <- bind_rows(d1, d2)
sim.3.02.02.4.5.2500.pa.B <- join(sim.3.02.02.4.5.2500.pa, sim.3.02.02.4.5.2500.pa.B, by=c("n", "j", "seed"), type="inner")

ret.2.B <- ret.2
dim(sim.3.02.02.4.5.2500.pa.B)

tm <- system.time(
  ret.3 <- mclapply(1:500, sim.3.hybrid, n_range=2500,
                    seed.data=sim.3.05.05.1.5.2500.pa,
                    alpha=0.1, beta=0.1,
                    theta=0, p_ij=p_ij.5, B=1,
                    mc.cores = 4)
)
save(ret.3, file="../sim.data/ret.3.RData")


# save data format
# sim.{scenario}.{alpha}.{beta}.{theta}.{mis_prob}.{n}.{project}
sim.3.01.01.4.5.2500.pa <- ret
save(sim.3.01.01.1.5.2500.pa, file="../sim.data/sim.3.01.01.1.5.2500.pa.RData")

# result
sim.3.01.01.1.5.2500.pa <- do.call(rbind.data.frame, ret.3)
head(sim.3.01.01.1.5.2500.pa)
data.est <- sim.3.02.02.4.5.2500.pa.B %>%
  gather(key='type', value='est', c("true", "observed", 
                                    "stoEM_reg", "simex.3", "simex.U.3",
                                    "stoEM_reg.B", "simex.UB.3")) %>%
  filter(type %in% c("true", "observed", 
                     "stoEM_reg", "simex.3", "simex.U.3",
                     "stoEM_reg.B", "simex.UB.3")) %>%
  mutate(bias=est-(1))

data.est %>% group_by(type) %>%
  dplyr::summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)))

# how about CI coverage rate?
data.se <- sim.3.02.02.4.5.2500.pa.B %>%
  gather(key='type', value='se', c("true.se", "observed.se", 
                                    "stoEM_reg.se", "simex.se.3", "simex.se.U.3",
                                   "stoEM_reg.se.B", "simex.se.UB.3")) %>%
  filter(type %in% c("true.se", "observed.se", 
                     "stoEM_reg.se", "simex.se.3", "simex.se.U.3",
                     "stoEM_reg.se.B", "simex.se.UB.3")) 
data <- cbind(data.est, data.se)[-c(6,7,8)]
head(cbind(ret.4.est, ret.4.se))

data %>%
  mutate(ll=est-qnorm(1-(1-0.95)/2)*se,
         ul=est+qnorm(1-(1-0.95)/2)*se) %>%
  mutate(coverage=(ll <= 1 & 1 <= ul)) %>%
  group_by(type) %>%
  dplyr::summarise(est.mean = mean(est),
            bias.mean = mean(abs(bias)),
            coverage = mean(coverage))

ret.4.est$type
library(plyr)

detach_package("plyr", TRUE)
sim.3.05.05.4.5.2500.pa <- join(sim.3.05.05.4.5.2500.pa, ret1, by=c("n", "j", "seed"), type="inner")
head(sim.3.05.05.4.5.2500.pa)

# the range of alpha and beta
sim.data$X1.cont <- as.numeric(sim.data$X1)-1
sim.data$X2.cont <- as.numeric(sim.data$X2)-1
sim.data$X3.cont <- as.numeric(sim.data$X3)-1
coxph(Surv(Y, D) ~ A + X1.cont + X2.cont + X3.cont, data=sim.data, model = TRUE)

# alpha, beta = 0.1, 0.1

