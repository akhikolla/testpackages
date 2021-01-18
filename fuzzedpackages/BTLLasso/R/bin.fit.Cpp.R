


bin.fit.Cpp <- function(resp, design, kat, epsilon = 1e-05, penalty, 
  lambda, max.iter = 200, start = NULL, adaptive = NULL, norm = "L1", 
  control = list(c = 1e-06, gama = 20, index = NULL), m, hat.matrix = FALSE, 
  lambda2 = 1e-04) {

  N <- length(resp)
  q <- kat - 1
  n <- N/q
  
  acoefs <- penalty$acoefs
  
  if (is.null(start)) {
    start <- rep(0,ncol(design))
    if(any(which(rowSums(abs(acoefs)) == 0))){
      start[which(rowSums(abs(acoefs)) == 0)] <- coef(glm.fit(y = resp, x = design[,which(rowSums(abs(acoefs)) == 0)], family = binomial()))
    }
    if (any(is.na(start))) {
      start[which(is.na(start))] <- 0
    }
  }
  
  if (is.null(adaptive)) {
    weight <- as.vector(rep(1, ncol(acoefs)))
  } else {
    weight <- abs(t(acoefs) %*% adaptive)
    if (any(weight == 0)) 
      weight[which(weight == 0)] <- epsilon
    weight <- as.vector(weight^(-1))
  }
  
  pen.nums <- c(penalty$numpen.order, penalty$numpen.intercepts, 
    penalty$numpen.X, penalty$numpen.Z1, penalty$numpen.Z2)
  
  if (sum(pen.nums) > 0) {
    if (penalty$weight.penalties) {
      pen.nums.scaled <- c(penalty$numpen.order/penalty$n.order, 
        penalty$numpen.intercepts/(m - 1), penalty$numpen.X/penalty$p.X/(m - 
          1), penalty$numpen.Z1/penalty$p.Z1/m, penalty$numpen.Z2/penalty$p.Z2)
      weight <- weight/rep(pen.nums.scaled, pen.nums)
    }
  }
  
  beta.old <- beta.new <- start
  diff <- 1
  delta.new <- delta.old <- 1
  
  
  rcpp.out <- binfit(matrix(beta.new, ncol = 1), epsilon, max.iter, 
    acoefs, lambda, matrix(weight, ncol = 1), control, design, 
    N, n, q, matrix(resp, ncol = 1), control$index, control$c, 
    control$gama, norm, as.numeric(hat.matrix), lambda2)
  
  
  beta.new <- rcpp.out$beta.new
  start <- rcpp.out$start
  df <- rcpp.out$df
  df2 <- rcpp.out$df2
  
  rownames(beta.new) <- names(start)
  
  
  if (norm == "grouped") {
    beta.pen <- matrix(beta.new[rowSums(acoefs) != 0], nrow = m - 
      1)
    norm.col <- sqrt(colSums((beta.pen)^2))/(m - 1)
    beta.pen[, norm.col < epsilon] <- 0
    beta.new[rowSums(acoefs) != 0] <- c(beta.pen)
  }
  
  
  return(list(coefficients = beta.new, start = start, df = df, weight = weight, df2 = df2))
  
}
