#' Generate synthetic functional ANOVA dataset
#'
#' This function generates a synthetic 3-factor functional ANOVA dataset.
#'
#' @param st.dev The standard deviation of the error.
#' @param n.replicates The number of replicates for each factor combination.
#' 
#' @return A list containing the data without noise, the data with noise, and the 
#' design matrix.
#' @export
#' @examples
#' data <- GenerateSyntheticAnova(st.dev = 5, n.replicates = 10)
#' ix <- 1
#' plot(data$clean.Y[ix, ], type = "l", col = "red", ylab = "")
#' lines(data$noisy.Y[ix, ], col = "blue")
GenerateSyntheticAnova <- function(st.dev = 10, n.replicates = 5) {
  
  m <- 512
  # First factor
  factor.A.1 <- DJ.EX(n = m, noisy = FALSE)$doppler
  factor.A.2 <- factor.A.1
  # Second factor
  factor.B.1 <- DJ.EX(n = m, noisy = FALSE)$bumps
  factor.B.2  <- factor.B.1 
  factor.B.2[310 : 360] = 0.5 * factor.B.1[310 : 360]
  factor.B.3  <- factor.B.1 
  factor.B.3[310 : 360] = 0.2 * factor.B.1[310 : 360]
  # Third factor
  factor.C.1 <- DJ.EX(n = m, noisy = FALSE)$blocks
  factor.C.2 <- factor.C.1
  factor.C.2[45 : 60]  <- 1.5 * factor.C.1[45 : 60] 
 
  sex <- as.factor(c("M", "F"))
  snp <- as.factor(c("aa", "ab", "bb"))
  ttt <- as.factor(c("T", "U"))
  X <- expand.grid(factorA = sex,
                   factorB = snp,
                   factorC = ttt)
  X <- X[rep(1 : nrow(X), each = n.replicates), ]
  n_tot <- nrow(X)
  
  small.signal <- data.frame(factor.A.1 = factor.A.1, 
                             factor.A.2 = factor.A.2, 
                             factor.B.1 = factor.B.1, 
                             factor.B.2 = factor.B.2, 
                             factor.B.3 = factor.B.3,
                             factor.C.1 = factor.C.1,
                             factor.C.2 = factor.C.2)

  clean.Y <- matrix(NA, nrow = n_tot, ncol= m)
  noisy.Y <- matrix(NA, nrow = n_tot, ncol= m)
  
  frml <- formula(~ 1 + factorA + factorB + factorC)
  Z <- model.matrix(frml, X)
  
  for (i in 1 : nrow(X)) { # for each factor combination
    clean.Y[i, ] <- small.signal$factor.A.1 + 
      (small.signal$factor.A.2 - small.signal$factor.A.1) * Z[i, "factorAM"] + 
      small.signal$factor.B.1 + 
      (small.signal$factor.B.2 - small.signal$factor.B.1 ) * Z[i, "factorBab"] + 
      (small.signal$factor.B.3 - small.signal$factor.B.1 ) * Z[i, "factorBbb"] +
      small.signal$factor.C.1 + 
      (small.signal$factor.C.2 - small.signal$factor.C.1 ) * Z[i, "factorCU"]
    noisy.Y[i, ] <- clean.Y[i, ] + rnorm(m, 0, st.dev)
  }
  
  return(list(noisy.Y = noisy.Y,
              X = X,
              clean.Y = clean.Y))
}

