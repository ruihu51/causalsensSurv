# generate data following the data structure of the cohort data
# sample sizes
n <- 160000

# outcome
# survival time T
# censoring time C
T <- rexp(n, 1/35)
C <- runif(n, 24, 26)

# observed outcome
# observed event time
# censoring indicator
Y <- pmin(T, C)
cat("follow-up years: ", mean(Y), sd(Y)) # mean - 18.9, sd - 6.3
D <- (T<=C)
cat("percentage of deaths: ", mean(D)) # around 0.5

# treatment - physical activity
A <- sample(0:2, size = n, replace = TRUE, prob = c(0.15, 0.35, 0.5))
cat("distribution of the treatment: ", prop.table(table(A)))

# covariates
# binary - eg. gender
X1 <- rbinom(n, 1, 0.33) # 0-male # 1-female
# ordinal - eg. age
X2 <- sample(0:4, size = n, replace = TRUE, prob = c(0.12, 0.2, 0.29, 0.35, 0.04)) # - as.factor
# TODO: sample X2 based on A
# categorical - eg. health status
X3 <- sample(0:4, size = n, replace = TRUE, prob = c(0.12, 0.2, 0.29, 0.35, 0.04)) # - as.factor

# try to mimic Table 2
# TODO: generate A conditional on X or not
prop.table(table(A, X2))

# smoking behavior
# years since cessation
years_quit <- sample(0:3, size = n, replace = TRUE, prob = c(0.12, 0.2, 0.29, 0.35, 0.04))
# CPDs
cpd <- sample(0:5, size = n, replace = TRUE, prob = c(0.12, 0.2, 0.29, 0.35, 0.04))
# possible measurement error (misclassification in this case because the variables are categorical)
# TODO: measurement error model for ordinal variables
  
# possible unmeasurement variable
U <- rbinom(n, 1, 0.5) # the simplest setup
# generate continuous variable and transform it into categorical variables
U <- runif(n, 0, 1)
U <- rnorm(n, 0, 1)

# naive causal hazard ratio
# include all possible treatment values
fit.naive <- coxph(Surv(Y, D) ~ A + X1 + X2 + X3 + years_quit + cpd)
tau.naive  <- fit.naive$coefficients[1]
# 1 vs 0 - 0.88
# 2 vs 0 - 0.83
# TODO: separate the data
