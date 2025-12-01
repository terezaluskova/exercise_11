
library(Biostrings)

viterbi <- function(obs, HMM1) {
  N  <- HMM1$N          # set of hidden states
  M  <- HMM1$M          # set of emitted characters
  A  <- HMM1$A          # transition probability matrix
  B  <- HMM1$B          # emission probability matrix
  pi <- HMM1$pi         # initial probability distribution vector
  
  obs <- unlist(strsplit(as.character(obs), ""))
  
  obslen <- length(obs)
  nlen <- length(N)
  
  obs_idx <- match(obs, M)
  
  V <- matrix(-Inf, nlen, obslen)
  back <- matrix(NA, nrow = nlen, ncol = obslen)
 
  for (i in 1:nlen){
    V[i,1] <- pi[i] + B[i, obs_idx[1]]
    back[i, 1] <- 0
  }
  for (t in 2:obslen){
    for (j in 1:nlen){
      scores <- V[, t-1] + A[, j] + B[j, obs_idx[t]]
      V[j, t] <- max(scores)
      back[j, t] <- which.max(scores)
    }
  }
  best_last <- which.max(V[, obslen])
  path_idx <- numeric(obslen)
  path_idx[obslen] <- best_last
  
  for (t in (obslen-1):1) {
    path_idx[t] <- back[path_idx[t+1], t+1]
  }
  path <- N[path_idx]
  
  return(list(V = V, path = path))
}

load("V:/MPA-PRG/exercise_11/HMM1.RData")

seq <- AAString("TGA")
result <- viterbi(seq, HMM1)
result$V      
result$path
