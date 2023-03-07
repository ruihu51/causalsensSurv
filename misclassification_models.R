# Notations
# The misclassification model of a categorical variable is defined as a K by K matrix, mis_prob.
# Let p_ij = P(C*=i|C=j) denote the elements of the matrix. i=1,...K and j=1,...,K.
# To generate C from C*, (K-1)*K parameters needed to be specified.
# Let n_i denote the number of C=i which is unknown and m_i denote the number of C*=i which is known. 
# Let n_ij denote the number of (C*=i, C=j). Let N denote the total number of observations.
# Define q_ij = P(C=i|C*=j)

# utils
gen_chd_given_parents <- function(parents, category, prob_matrix){
  j <- which(category==parents)
  sample_prob <- prob_matrix[,j]
  chd <- sample(0:(length(category)-1), size = 1, replace = TRUE, prob = sample_prob)
  return(chd)
}

# parameters
n <- 1000
p_ij <- matrix(c(0.8,0.1,0.2,0.9), nrow = 2, byrow = TRUE)
C <- factor(sample(0:1, size = n, replace = TRUE, 
                   prob = c(0.4, 0.6)))

# 1)
# generate C* given C based on misclassification probabilities p_ij (applicable to K>=2)
C.star <- factor(unlist(lapply(C, gen_chd_given_parents, category=levels(C), prob_matrix=p_ij)))
# small trick: check correlation between X3 and X3.star
data.frame(table(C,C.star)) %>%
  ggplot(aes(x=C, y=Freq, fill=C.star)) + geom_bar(stat="identity")

# 2)
# generate C given C*
# CASE I: K=2
gen_unobs <- function(N, C.star, p_11, p_22){
  m1 <- sum(C.star==levels(C.star)[1]) # number of C.star==0
  m2 <- N - m1
  # calculate p_ij
  p_ij <- matrix(c(p_11, 1-p_22, 1-p_11, p_22), nrow = 2, byrow = TRUE)
  # calculate n_i using m_i and mis_prob
  n2 <- (m2 - N*(1-p_ij[1,1]))/(p_ij[2,2]-(1-p_ij[1,1]))
  n1 <- N - n2
  # calculate n_ij using n_i and p_ij
  n_col <- c(n1, n2)
  n_ij <- t(n_col*t(p_ij))
  cat(n_col,"\n")
  # calculate q_ij
  m_row <- c(m1, m2)
  q_ij <- t((1/m_row)*(n_ij))
  cat(q_ij)
  cat(apply(q_ij,2,sum))
  
  C <- numeric(N)
  for (i in 1:N){
    if(C.star[i]=="1"){
      u <- runif(1,0,1)
      C[i] <- ifelse(u<q_ij[2,2],"1","0")
    } else {
      u <- runif(1,0,1)
      C[i] <- ifelse(u<q_ij[1,1],"0","1")
    }
  }
  # generate C from C* based on q_ij
  # C <- factor(unlist(lapply(C.star, gen_chd_given_parents,
  #                                 category=levels(C.star), prob_matrix=q_ij)))
  return(C)
}
C.check <- gen_unobs(N = 1000, C.star = C.star, p_11 = 0.8, p_22=0.9)

# validation for generate C based on C* - simulation study
# 1) sample true C, then generate C* given C
# 2) use the above function to generate C given C.star, to check whether the function is correct
n <- 1000
p_ij <- matrix(c(0.9,0.1,0.1,0.9), nrow = 2, byrow = TRUE)
C <- factor(sample(0:1, size = n, replace = TRUE, 
                   prob = c(0.4, 0.6)))
C.star <- factor(unlist(lapply(C, gen_chd_given_parents, category=levels(C), prob_matrix=p_ij)))
C.check <- gen_unobs(N = n, C.star = C.star, p_11 = p_ij[1,1], p_22=p_ij[2,2])

table(C.star, C.check)
table(C.check,C)
table(C)
table(C.star)

chisq.test(C, C.check) # H0: independent
chisq.test(C, C.star)

mean(C==C.check)
mean(C==C.star)

# results we want to see:
# C is close to C.check
# but here we see C==C.star larger than C==C.check
# so I think this might just due to the randomness
# TODO: not sure

# CASE II: K=3
# TODO:

# CASE III: K>=2 with additional assumptions
# TODO: derivation


############################
# scripts for generating C given C*
# 0. Known parameters
N <- 1664
m1 <- 1449
m2 <- N - m1

# Specify the parameters: p_ij (2-1)*2=2 parameters needed to be specified
p_22 <- 0.78
p_11 <- 0.99

# 1. Calculate the mis_prob using known p_ij

p_ij <- matrix(c(p_11, 1-p_22, 1-p_11, p_22), nrow = 2, byrow = TRUE)

# 2. Calculate n_i using m_i and mis_prob
n2 <- (m2 - N*(1-p_ij[1,1]))/(p_ij[2,2]-(1-p_ij[1,1]))
n1 <- N - n2

# 3. Calculate n_ij using n_i and p_ij
n_col <- c(n1, n2)
n_ij <- t(n_col*t(p_ij))
# t(c(1,2)*t(matrix(c(3,4,5,6), nrow = 2, byrow = TRUE)))

# 4. Calculate q_ij
m_row <- c(m1, m2)
q_ij <- t((1/m_row)*(n_ij)) # n_ij transpose twice

# 5. Generate C from C* based on q_ij
gen_chd_given_parents <- function(parents, category, prob_matrix){
  j <- which(category==parents)
  sample_prob <- prob_matrix[,j]
  chd <- sample(0:(length(category)-1), size = 1, replace = TRUE, prob = sample_prob)
  return(chd)
}
var_category <- levels(C.star)
C.check <- factor(unlist(lapply(C.star, gen_chd_given_parents, category=var_category, prob_matrix=q_ij)))
# small trick: check correlation between X3 and X3.star
data.frame(table(C,C.check)) %>%
  ggplot(aes(x=C, y=Freq, fill=C.check)) + geom_bar(stat="identity")
table(C.check)
table(C.star)
table(C)
