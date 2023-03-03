# Quantitative bias analysis pipeline

# Notations
# The misclassification model of a categorical variable is defined as a K by K matrix, mis_prob.
# Let p_ij = P(C*=i|C=j) denote the elements of the matrix. i=1,...K and j=1,...,K.
# To generate C from C*, (K-1)*K parameters needed to be specified.
# Let n_i denote the number of C=i which is unknown and m_i denote the number of C*=i which is known. 
# Let n_ij denote the number of (C*=i, C=j). Let N denote the total number of observations.
# Define q_ij = P(C=i|C*=j)

# CASE I: K=2
# 0. Known parameters
N <- 1000
m1 <- 800
m2 <- N - m1

# Specify the parameters: p_ij (2-1)*2=2 parameters needed to be specified
p_22 <- 0.8
p_11 <- 0.9

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
q_ij <- t((1/m_row)*t(n_ij))

# 5. Generate C from C* based on q_ij
gen_chd_given_parents <- function(parents, category, prob_matrix){
  j <- which(category==parents)
  sample_prob <- prob_matrix[,j]
  chd <- sample(0:(length(category)-1), size = 1, replace = TRUE, prob = sample_prob)
  return(chd)
}
var_category <- levels(C.star)
C <- factor(unlist(lapply(C.star, gen_chd_given_parents, category=var_category, prob_matrix=q_ij)))

# CASE II: K=3
# TODO:

# CASE III: K>=2 with additional assumptions
# TODO: derivation


# validation for generate C based on C* - simulation study
