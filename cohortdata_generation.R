rm(list=ls())

# install packages
library(survival)

# generate data following the data structure of the cohort data
# sample sizes
n <- 160000

# treatment - physical activity
A <- sample(0:2, size = n, replace = TRUE, prob = c(0.15, 0.35, 0.5))
cat("distribution of the treatment: ", prop.table(table(A)))
data <- data.frame(sort(A))
colnames(data) <- c("A")

# covariates
# binary - eg. gender
# 0-male # 1-female
sex  <- c()
prob_list <- matrix(c(0.72, 0.28,
                      0.71, 0.29,
                      0.59, 0.41), nrow = 3, byrow = TRUE)
size_list <- table(A)
for (i in 1:length(size_list)){
  sex  <- append(sex, 
                   sample(0:1, size = size_list[i], replace = TRUE, 
                          prob = prob_list[i,]))
}
data$sex <- factor(sex)

# ordinal - eg. age
# sample age based on A (but here assume age is independent of A)
age <- c()
for (size in table(A)){
  age <- append(age, 
               sample(0:4, size = size, replace = TRUE, 
                      prob = c(0.12, 0.2, 0.29, 0.35, 0.04)))
}
data$age <- factor(age)

# categorical - eg. health status
# sample health status based on A
# better health status indicates higher physical activity score
health <- c()
prob_list <- matrix(c(0.24, 0.42, 0.33, 0.01,
                      0.15, 0.38, 0.46, 0.01,
                      0.06, 0.25, 0.68, 0.01), nrow = 3, byrow = TRUE)
for (i in 1:length(size_list)){
  health <- append(health, 
                sample(0:3, size = size_list[i], replace = TRUE, 
                       prob = prob_list[i,])) # - as.factor
}
data$health <- factor(health)

# smoking behavior
# years since cessation
years_quit <- c()
prob_list <- matrix(c(0.66, 0.18, 0.11, 0.05,
                      0.71, 0.15, 0.10, 0.4,
                      0.82, 0.11, 0.05, 0.02), nrow = 3, byrow = TRUE)
for (i in 1:length(size_list)){
  years_quit <- append(years_quit, 
                   sample(0:3, size = size_list[i], replace = TRUE, 
                          prob = prob_list[i,])) # - as.factor
}
data$years_quit <- factor(years_quit)

# CPDs
cpd<- c()
prob_list <- matrix(c(0.18, 0.26, 0.21, 0.16, 0.14, 0.05,
                      0.21, 0.30, 0.21, 0.15, 0.10, 0.03,
                      0.33, 0.33, 0.17, 0.10, 0.05, 0.02), nrow = 3, byrow = TRUE)
for (i in 1:length(size_list)){
  cpd <- append(cpd, 
             sample(0:5, size = size_list[i], replace = TRUE, 
                    prob = prob_list[i,])) # - as.factor
}
data$cpd <- factor(cpd)

data$A <- factor(data$A)

# try to mimic Table 2
prop.table(table(data$A, data$sex), margin=1)
prop.table(table(data$A, data$age), margin=1)
prop.table(table(data$A, data$health), margin=1)
prop.table(table(data$A, data$years_quit), margin=1)
prop.table(table(data$A, data$cpd), margin=1)

# sample T conditional on A and X
# small trick: multiply factor variables
rate <- (1/30)*exp(with(data, c(0, -0.12, -0.19)[A] + c(0, -0.05)[sex] 
                        + c(0, 0.01, 0.01, 0.01, 0.01)[age] + c(0, -0.03, -0.01, 0)[health] 
                        + c(0, 0.03, 0.05, 0.07)[years_quit] + c(0, 0.02, 0.03, 0.05, 0.07, 0.09)[cpd]))
# outcome
# survival time T
# censoring time C
T <- rexp(n, rate)
C <- runif(n, 24, 26)

# observed outcome
# observed event time
# censoring indicator
Y <- pmin(T, C)
cat("follow-up years: ", mean(Y), sd(Y)) # mean - 18.9, sd - 6.3
D <- (T<=C)
cat("percentage of deaths: ", mean(D)) # around 0.5
data$Y <- Y
data$D <- D

# naive causal hazard ratio
# include all possible treatment values
fit.naive <- coxph(Surv(Y, D) ~ A + sex + age + health
                   + years_quit + cpd, data=data)
exp(fit.naive$coefficients[1:2])
# 1 vs 0 - 0.88
# 2 vs 0 - 0.83
fit.naive.1 <- coxph(Surv(Y, D) ~ A + sex + age + health
                   + years_quit + cpd, data = data[data$A!=2,])
exp(fit.naive.1$coefficients[1:2])
fit.naive.2 <- coxph(Surv(Y, D) ~ A + sex + age + health
                     + years_quit + cpd, data = data[data$A!=1,])
exp(fit.naive.2$coefficients[1:2])
# coefficient slightly different
# TODO: can we do separate models?


# possible measurement error (misclassification in this case because the variables are categorical)
# TODO: measurement error model for ordinal variables

# possible unmeasurement variable
U <- rbinom(n, 1, 0.5) # the simplest setup
# generate continuous variable and transform it into categorical variables
U <- runif(n, 0, 1)
U <- rnorm(n, 0, 1)


