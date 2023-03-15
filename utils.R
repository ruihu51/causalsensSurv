# a dictionary save the parameters settings
par_list <- list(sim1=list(),
                 sim2=list(
                   eta=list("2"=list(
                     eta_a = list(eta_a1=c(0,0.1),
                                  eta_a2=c(0,0.01,0.02,0.03,0.04),
                                  eta_a3=c(0,-0.2,-0.4,-0.6,-0.8,-1)), 
                     eta_t = list(eta_t1=c(0,-0.1),
                                  eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                  eta_t3=c(0,0.2,0.4,0.6,0.8,1))),
                     "3"=list(eta_a = list(eta_a1=c(0,0.1),
                                           eta_a2=c(0,0.01,0.02,0.03,0.04),
                                           eta_a3=c(0,-1.2,-1.4,-1.6,-1.8,-2)), 
                              eta_t = list(eta_t1=c(0,-0.1),
                                           eta_t2=c(0,-0.01,-0.02,-0.03,-0.04),
                                           eta_t3=c(0,1.2,1.4,1.6,1.8,2)))
                     ),
                   mis_prob=list("2"=list(
                     p_ij=matrix(c(0.8, 0.1, 0, 0, 0, 0,
                                   0.2, 0.8, 0.1, 0, 0, 0,
                                   0, 0.1, 0.8, 0.1, 0, 0,
                                   0, 0, 0.1, 0.8, 0.1, 0,
                                   0, 0, 0, 0.1, 0.8, 0.2,
                                   0, 0, 0, 0, 0.1, 0.8), nrow = 6, byrow = TRUE),
                     p_ij.star=matrix(c(0.8, 0.1, 0.01, 0.01, 0.01, 0.01,
                                       0.16, 0.8, 0.1, 0.01, 0.01, 0.01,
                                       0.01, 0.07, 0.8, 0.1, 0.01, 0.01,
                                       0.01, 0.01, 0.07, 0.8, 0.07, 0.01,
                                       0.01, 0.01, 0.01, 0.07, 0.8, 0.16,
                                       0.01, 0.01, 0.01, 0.01, 0.1, 0.8), nrow = 6, byrow = TRUE)),
                     "3"=list(
                       p_ij=matrix(c(0.8, 0.1, 0.01, 0.01, 0.01, 0.01,
                                     0.16, 0.8, 0.1, 0.01, 0.01, 0.01,
                                     0.01, 0.07, 0.8, 0.1, 0.01, 0.01,
                                     0.01, 0.01, 0.07, 0.8, 0.07, 0.01,
                                     0.01, 0.01, 0.01, 0.07, 0.8, 0.16,
                                     0.01, 0.01, 0.01, 0.01, 0.1, 0.8), nrow = 6, byrow = TRUE),
                       p_ij.star=matrix(c(0.8, 0.1, 0.01, 0.01, 0.01, 0.01,
                                          0.16, 0.8, 0.1, 0.01, 0.01, 0.01,
                                          0.01, 0.07, 0.8, 0.1, 0.01, 0.01,
                                          0.01, 0.01, 0.07, 0.8, 0.07, 0.01,
                                          0.01, 0.01, 0.01, 0.07, 0.8, 0.16,
                                          0.01, 0.01, 0.01, 0.01, 0.1, 0.8), nrow = 6, byrow = TRUE))
                   )
                 ))

# simulation data generation coding
# gaussian
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