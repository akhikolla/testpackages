kNormFun <- function(X1, x2, par, fun) { 

    ## X1, x2 : matrices with same number of columns 'd' (dimension)
    n1 <- nrow(X1)
    n2 <- nrow(x2)
    d <- ncol(X1)
    Grad <- array(NA, dim = c(n1, n2, d + 1L),
                  dimnames = list(rownames(X1), rownames(x2), names(par)))
    SS2 <- 0  
    for (j in 1L:d){
        Hj <- outer(X1[ , j], x2[ , j], "-")
        Hj2 <- (Hj / par[j])^2
        SS2 <- SS2 + Hj2
        Grad[ , , j] <- Hj2 / par[j]
    }
    
    H <- sqrt(SS2)
    K <- fun(x1 = H, x2 = 0, par = c(1.0, 1.0))
    Grad[ , , d + 1L] <- K
    dK <- attr(K, "gradient")

    ## now put in matrix/array format. Keep only the 'range' part
    ## of the gradient position #1
    K <- array(K, dim = c(n1, n2))
    dK <- array(dK[ , 1L, 1L],  dim = c(n1, n2))
    K <- par[d + 1L] * K 

    fac <- par[d + 1L] * dK  / SS2
    fac[SS2 < 1e-8] <- 0.0
    
    Grad[ , , 1L:d] <- sweep(Grad[ , , 1L:d, drop = FALSE],
                             MARGIN = c(1L, 2L), STATS = fac, FUN = "*")
    
    attr(K, "gradient") <- Grad
    return(K)
}


kMatern <- function(d = 1, nu = "5/2"){
  
  if (!is.element(nu, c("1/2", "3/2", "5/2"))) {
      stop("the possible values of nu are \"1/2\", \"3/2\", \"5/2\"")
  }
  if (nu == "1/2") {
    fun <- k1FunExp
    label <- "exponential kernel"}
  if (nu == "3/2") {
    fun <- k1FunMatern3_2
    label <- "Matern kernel with nu = 3/2"}
  if (nu == "5/2") {
    fun <- k1FunMatern5_2
    label <- "Matern kernel with nu = 5/2"}
  
  kernFun <- function(x1, x2, par){
    kNormFun(x1, x2, par, fun)
  }
  k <- covMan( 
    kernel = kernFun,
    hasGrad = TRUE,
    acceptMatrix = TRUE,
    label = label,
    d = d,
    par = c(rep(1, d), 1),
    parLower = rep(1e-8, d + 1L),
    parUpper = rep(Inf, d + 1L),
    parNames = c(paste("theta", 1L:d, sep = "_"), "sigma2")
  )
#  k@kernel <- get("fun", envir = environment(k@kernel))
  return(k)
}

kGauss <- function(d = 1){
  kernFun <- function(x1, x2, par){
    kNormFun(x1, x2, par, k1FunGauss)
  }
  k <- covMan( 
    kernel = kernFun,
    hasGrad = TRUE,
    acceptMatrix = TRUE,
    d = d,
    par = c(rep(1, d), 1),
    parLower = rep(1e-8, d + 1L),
    parUpper = rep(Inf, d + 1L),
    parNames = c(paste("theta", 1L:d, sep = "_"), "sigma2"),
    label = "Gauss kernel"
  )
  return(k)
}
