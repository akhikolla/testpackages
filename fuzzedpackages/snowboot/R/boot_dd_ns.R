#Two internal functions for bootstrapping non-seeds in boot_dd:
#boot_dd_ns -- usual bootstrap;
#boot_dd_ns_cl -- bootstrap using a cluster.

boot_dd_ns <- function(dns, k, B) {
  if (length(dns) < 2) {
    tmp <- matrix(0, length(k), B)
    if (length(dns) == 0){
      tmp
    } else {
      tmp[dns,] <- 1
      tmp
    }
  } else {
    sapply(1:B, FUN = function(b)
      table(c(sample(dns, replace = TRUE, prob = dns^{-1}), k)) - 1
    )
  }
}


boot_dd_ns_cl <- function(dns, k, B, cl) {
  if (length(dns) < 2) {
    tmp <- matrix(0, length(k), B)
    if (length(dns) == 0){
      tmp
    } else {
      tmp[dns,] <- 1
      tmp
    }
  } else {
    parallel::parSapply(cl, X = 1:B, FUN = function(b)
      table(c(sample(dns, replace = TRUE, prob = dns^{-1}), k)) - 1
    )
  }
}
