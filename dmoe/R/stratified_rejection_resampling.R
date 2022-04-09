#----------------------------------------------------------------------
# Function for resampling particles using stratified sampling with rejection scheme
#----------------------------------------------------------------------
#  The following function implements the stratified sampling algorithm with rejection
#  The resampling scheme is based on the paper " Online inference for hidden Markov models via particle filters"
#----------------------------------------------------------------------
# Function for resampling particles using stratified sampling with rejection scheme
#----------------------------------------------------------------------

# The partition function implements the partition algorithm used to search particles that should be rejected/ discarded
# It is based on the cumulative distribution of the importance weights.
# Input:
# w: a vector of importance weights (of length M).
# n < M: the deisred size of output vector.
# output:
# k: the threshold of the largest value that is less than all the importance weights of the n best particles.
partition <- function(w, n) {
  w <- sort(w)
  i <- length(w)
  k <- w[i]
  bk <- sum(w[which(w < k)])
  ak <- length(which(w >= k))
  f <- (bk / k) + ak - n
  if (f > 0) {
    return(k = 1 / n)
  }
  else {
    i <- i - 1
    while (f < 0) {
      k <- w[i]
      bk <- sum(w[which(w < k)])
      ak <- length(which(w >= k))
      f <- (bk / k) + ak - n
      i <- i - 1
    }
  }
  return(k)
}

# implementing the stratified sampling algorithm. Same input as above.
# Input:
# w: a vector of importance weights (of length M).
# n < M: the deisred size of output vector.
# output:
# index: a vector of selected indexes of length n.
strat_resample <- function(w, n) {
  m <- length(w)
  k <- sum(w / n)
  u <- runif(1, 0, k)
  i <- 1
  index <- c()
  while (i <= m) {
    u <- u - w[i]
    if (u < 0) {
      index <- c(index, i)
      u <- u + k
    }
    i <- i + 1
  }
  return(index)
}
#  The following function implements the stratified sampling algorithm with rejection
#  The resampling scheme is based on the paper " Online inference for hidden Markov models via particle filters"
#  by Paul Fearnhead and Peter Clifford.
# The aim of this function is to shuffle the particles.
# Input:
# w: a vector of importance weights (of length M).
# size < M: the desired size of output vector.
# output:
# index: vector of particles indexes that were resampled
# weights: new normalized weights

reject_strat_resample <- function(w, size) {
  probtreshold <- partition(w, size)
  index1 <- which(w >= probtreshold)
  w1 <- w[index1]
  w2 <- w[which(w < probtreshold)]
  n2 <- size - length(w1)
  if(n2 > 0){
  index2 <- strat_resample(w = w2, n = n2)
  new_weights <- c(w1, rep(probtreshold, length(index2)))
  new_weights <- new_weights / sum(new_weights)
  index <- c(index1, index2)
  }else{ # no resampling performed
    index <- seq_len(length(w))
    new_weights <- w
  }
  return(list(index = index , weights = new_weights))
}

# this function normalizes any set empirical densities.
# input:
# values: a vector of unnormalized densities
# output:
# exp_val:a vector of normalized densities
normalize <- function(values) {
  n <- length(values)
  precision <- log(1e-16) - log(n) # all values less than this precision will be set to precision
  m <- max(values)
  values <- values - m
  exp_val <- ifelse(values >= precision, exp(values), exp(precision))
  return(exp_val / sum(exp_val))
}
