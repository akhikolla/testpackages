# Script to simulate the simdiss object included in the package
# (the number of clusters is fixed and set to 3)
make_simdata <- function(prob, n, S, seed, groups = c(0.3, 0.5, 0.2)) {
  set.seed(seed)
  
  sample1 <- sample(1:3, n, replace = TRUE, prob = c(0.4, 0.3, 0.3))
  sample2 <- sample(1:3, n, replace = TRUE, prob = c(0.5, 0.2, 0.3))
  sample3 <- sample(1:3, n, replace = TRUE, prob = c(0.3, 0.2, 0.5))
  
  diss1 <- as.matrix(dist(sample1))
  diss1[diss1 > 0] <- 1
  diss2 <- as.matrix(dist(sample2))
  diss2[diss2 > 0] <- 1
  diss3 <- as.matrix(dist(sample3))
  diss3[diss3 > 0] <- 1
  diss <- list()
  diss[[1]] <- diss1
  diss[[2]] <- diss2
  diss[[3]] <- diss3
  
  set.seed(seed)
  counts <- round(S*groups, digits = 0)
  counts[length(counts)] <- S - sum(counts[-length(counts)])
  diss_ran <- list()
  for (i in 1:3) {
    diss_ran[[i]] <- list()
    for (j in 1:counts[i]) {
      diss_ran[[i]][[j]] <- diss[[i]]
    }
  }
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      for (g in 1:3) {
        if (diss[[g]][i, j] == 1) {
          for (k in 1:counts[g]) {
            diss_ran[[g]][[k]][i, j] <- diss_ran[[g]][[k]][j, i] <-
              rbinom(1, 1, prob)
          }
        } else {
          for (k in 1:counts[g]) {
            diss_ran[[g]][[k]][i, j] <- diss_ran[[g]][[k]][j, i] <-
              rbinom(1, 1, 0.1)
          }
        }
      }
    }
  }
  
  subjects <- list()
  for (i in 1:3) {
    subjects[[i]] <- list()
    for (j in 1:counts[i]) {
      subjects[[i]][[j]] <- paste0("diss", i, "_", j)
    }
  }
  subjects <- unlist(subjects)
  
  all_diss <- list()
  s <- 1
  for (i in 1:3) {
    for (j in 1:counts[i]) {
      all_diss[[s]] <- as.dist(diss_ran[[i]][[j]])
      attr(all_diss[[s]], "call") <- NULL
      s <- s + 1
    }
  }
  names(all_diss) <- subjects
  
  return(all_diss)
}

library(dmbc)

n <- 16
S <- 10

all_diss_9 <- make_simdata(prob = 0.9, n = n, S = S, seed = 120)
all_diss_8 <- make_simdata(prob = 0.8, n = n, S = S, seed = 120)
all_diss_7 <- make_simdata(prob = 0.7, n = n, S = S, seed = 120)
all_diss_6 <- make_simdata(prob = 0.6, n = n, S = S, seed = 120)
all_diss_5 <- make_simdata(prob = 0.5, n = n, S = S, seed = 120)

simdata <- new("dmbc_data", diss = all_diss_9, n = n, S = S, 
  family = "binomial")
