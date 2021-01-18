#### In order not to call geoR, just copy this simple function here to be used
"materncopy" <-
  function (u, phi, kappa)
  {
    if(is.vector(u)) names(u) <- NULL
    if(is.matrix(u)) dimnames(u) <- list(NULL, NULL)
    uphi <- u/phi
    uphi <- ifelse(u > 0,
                   (((2^(-(kappa-1)))/ifelse(0, Inf,gamma(kappa))) *
                      (uphi^kappa) *
                      besselK(x=uphi, nu=kappa)), 1)
    uphi[u > 600*phi] <- 0
    return(uphi)
  }


#### Poisson marginal
poisson.gc <- function(link = "log", lambda = NULL){
  ans <- list()
  ans$discrete <- TRUE
  if(is.null(lambda)){
  fm <- poisson( substitute( link ) )
  ans$start <- function(y, x, effort) {
    mfit <- suppressWarnings(glm.fit(x, y/effort, family = fm))
    mu0 <- fitted(mfit)
    reg0 <- coef(mfit)
    glb <- qnorm(ppois(y-1, lambda = mu0))
    gub <- qnorm(ppois(y, lambda = mu0))
    names(reg0)[1] <- "Intercept"
    return(reg0)
  }
  ans$nod <- 0
  ans$bounds <- function(y, x, pars, effort) {
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    a <- qnorm(ppois( y-1, lambda = M));  b <- qnorm(ppois( y, lambda = M))
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    pdf <- dpois(y, lambda = M, log = FALSE)
    return(pdf)
  }
  ans$cdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    cdf <- ppois( y, lambda = M)
    return(cdf)
  }
  ans$fm <- fm
  }
  if(is.numeric(lambda)){
  q <- function(p) qpois(p = p, lambda = lambda)
  ans$margvar <- lambda
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rpois(n = nrep, lambda = lambda)
  ans$q <- q
  }
  class(ans) <- c("marginal.gc")
  return(ans)
}




#### NB marginal
negbin.gc <- function(link = "log", mu = NULL, od = NULL){
  ans <- list()
  ans$discrete <- TRUE
  if(is.null(mu) & is.null(od)){
  fm <- poisson( substitute( link ) )
  ans$start <- function(y, x, effort) {
    mfit <- suppressWarnings(glm.fit(x, y/effort, family = fm))
    reg0 <- coef(mfit); mu0 <- fitted(mfit)
    od0 <- max(10*.Machine$double.eps, mean(((y-mu0)^2-mu0)/mu0^2))
    glb <- qnorm(pnbinom(y-1, size = 1/od0, mu = mu0))
    gub <- qnorm(pnbinom(y, size = 1/od0, mu = mu0))
    names( od0 ) <- "overdispersion"
    names(reg0)[1] <- "Intercept"
    return(c(reg0, od0))
  }
  ans$nod <- 1
  ans$bounds <- function(y, x, pars, effort) {
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    S <- 1/pars[ncol(x)+1]
    a <- qnorm(pnbinom( y-1, size = S, mu = M))
    b <- qnorm(pnbinom( y, size = S, mu = M))
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    S <- 1/pars[ncol(x)+1]
    pdf <- dnbinom(y, size = S, mu = M, log = FALSE)
    return(pdf)
  }
  ans$cdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    S <- 1/pars[ncol(x)+1]
    cdf <- pnbinom(y, size = S, mu = M)
    return(cdf)
  }
  ans$fm <- fm
  }
  if(is.numeric(mu) & is.numeric(od)){
  q <- function(p) qnbinom(p = p, size = 1/od, mu = mu)
  ans$margvar <- mu*(1+mu*od)
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rnbinom(n = nrep, size = 1/od, mu = mu)
  ans$q <- q
  }
  class(ans) <- c("marginal.gc")
  return(ans)
}


#### Binomial marginal

binomial.gc <- function(link = "logit", size = NULL, prob = NULL){
  ans <- list()
  ans$discrete <- TRUE
  if(is.null(size) & is.null(prob)){
  fm <- binomial( substitute( link ) )
  ans$start <- function(y, x, effort) {
    mfit <- suppressWarnings(glm.fit(x, y/(effort), family = fm))
    reg0 <- coef(mfit)
    glb <- qnorm(pbinom(y-1, size = effort, prob = fitted(mfit)))
    gub <- qnorm(pbinom(y, size = effort, prob = fitted(mfit)))
    names(reg0)[1] <- "Intercept"
    return(reg0)
  }
  ans$nod <- 0
  ans$bounds <- function(y, x, pars, effort) {
    p <- fm$linkinv(pars[1:ncol(x)]%*%t(x))
    a <- qnorm(pbinom( y-1, size = effort, prob = p))
    b <- qnorm(pbinom( y, size = effort, prob = p))
    return(list(lower = a, upper = b))
  }

  ans$pdf <- function(y, x, pars, effort){
    p <- fm$linkinv(pars[1:ncol(x)]%*%t(x))
    pdf <- dbinom(y, size = effort, prob = p, log = FALSE)
    return(pdf)
  }

  ans$cdf <- function(y, x, pars, effort){
    p <- fm$linkinv(pars[1:ncol(x)]%*%t(x))
    cdf <- pbinom(y, size = effort, prob = p)
    return(cdf)
  }
  ans$fm <- fm
  }
  if(is.numeric(size) & is.numeric(prob))
    {
    q <- function(p) qbinom(p = p, size = size, prob = prob)
    ans$margvar <- size*prob*(1-prob)
    ans$int.marg <- function (order) {
      if (requireNamespace("EQL", quietly = TRUE)) {
        integrate(function(x, order)
          ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                  q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
          order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
      }else{
        stop("Please install {EQL} first!")
      }
    }
    ans$rsample <- function(nrep) rbinom(n = nrep, size = size, prob = prob)
    ans$q <- q
  }
  class(ans) <- c("marginal.gc")
  return(ans)
}


#### Zero-inflated Poisson


zip.gc <- function(link = "log", mu = NULL, od = NULL){
  ans <- list()
  ans$discrete <- TRUE
  if(is.null(mu) & is.null(od)){
  fm <- poisson( substitute( link ) )
  ans$start <- function(y, x, effort) {
    mfit <- suppressWarnings(glm.fit(x, y/effort, family = fm))
    reg0 <- coef(mfit); mu0 <- fitted(mfit)
    od0 <- max(10*.Machine$double.eps, mean(((y-mu0)^2-mu0)/mu0^2))
    glb <- qnorm((y >= 1)*od0/(1+od0)+ppois(y-1, (1+od0)*mu0)/(1+od0))
    gub <- qnorm((y >= 0)*od0/(1+od0)+ppois(y, (1+od0)*mu0)/(1+od0))
    names(reg0)[1] <- "Intercept"
    names(od0) <- "overdispersion"
    return(c(reg0, od0))
  }
  ans$nod <- 1
  ans$bounds <- function(y, x, pars, effort) {
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    od <- pars[ncol(x)+1]
    a <- qnorm( (y>=1)*od/(1+od) + ppois(y-1, lambda = (1+od)*M)/(1+od) )
    b <- qnorm( (y>=0)*od/(1+od) + ppois(y, lambda = (1+od)*M)/(1+od) )
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    od <- pars[ncol(x)+1]
    pdf <-  (y==0)*od/(1+od) + dpois(y, (1+od)*M)/(1+od)
    return(pdf)
  }
  ans$cdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    M <- ifelse(M == Inf, .Machine$double.xmax, M)
    od <- pars[ncol(x)+1]
    cdf <- (y >= 0)*od/(1+od)+ppois(y, (1+od)*M)/(1+od)
    return(cdf)
  }
  ans$fm <- fm
 }
  if(is.numeric(mu) & is.numeric(od))
  {
  q <- function(p) qpois(pmax( 0, (p-od/(1+od))/(1-od/(1+od)) ), (1+od)*mu)
  ans$margvar <- mu*(1+mu*od)
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) ifelse(rbinom(n=nrep, 1, od/(1+od)), 0, rpois(n=nrep, (1+od)*mu))
  ans$q <- q
 }
  class(ans) <- c("marginal.gc")
  return(ans)
}



#### Correlation Functions


matern.gc <- function(range = NULL, kappa = 0.5, nugget = TRUE){
  ans <- list()
  if(is.null(range)){
  if(nugget == TRUE){
    ans$nug <- 1
    ans$npar.cor <- 2
    ans$start <- function(D) {
      corstart <- c(median(D)/2,0.2)
      names(corstart) <- c("range","nugget")
      return(corstart)
    }
    ans$corr <- function(corrpar,D){
      S <- (1 - corrpar[2]) * materncopy(D, corrpar[1], kappa) + corrpar[2]*diag(NROW(D))
      return(S)
    }
  }else{
    ans$nug <- 0
    ans$npar.cor <- 1
    ans$start <- function(D)  {
      corstart <- median(D)/2
      names(corstart) <- "range"
      corstart
    }
    ans$corr <- function(corrpar,D){
      S <- materncopy(D, corrpar, kappa)
      return(S)
    }
   }
  }
  if(is.numeric(range) & is.numeric(nugget)){
  ans$S <- function(D) (1-nugget)*materncopy(D, range, kappa) + nugget*diag(NROW(D))
  }
  class(ans) <- c("corr.gc")
  return(ans)
}



powerexp.gc <- function(range = NULL, kappa = 1, nugget = TRUE){
  ans <- list()
  if(kappa > 2)
    stop("'Kappa' must be between 0 and 2")
  if(is.null(range)){
  if(nugget == TRUE){
    ans$nug <- 1
    ans$npar.cor <- 2
    ans$start <- function(D) {
      corstart <- c(median(D)/2, 0.2)
      names(corstart) <- c("range","nugget")
      return(corstart)
    }
    ans$corr <- function(corrpar, D ){
      S <- (1 - corrpar[2])*exp(-abs( (D/corrpar[1])^(kappa))) + corrpar[2]*diag(NROW(D))
      return(S)
    }
  }else{
    ans$nug <- 0
    ans$npar.cor <- 1
    ans$start <- function(D) {
      corstart <- median(D)/2
      names(corstart) <- "range"
      return(corstart)
    }
    ans$corr <- function(corrpar, D ){
      S <- exp(-abs((D/corrpar[1])^(kappa)))
      return(S)
    }
   }
  }
  if(is.numeric(range) & is.numeric(nugget)){
    ans$S <- function(D) (1-nugget)*exp(-abs( (D/range)^(kappa) )) + nugget*diag(NROW(D))
  }
  class(ans) <- c("corr.gc")
  return(ans)
}


spherical.gc<- function(range = NULL, nugget = TRUE){
  ans <- list()
  if(is.null(range)){
  if(nugget == TRUE){
    ans$nug <- 1
    ans$npar.cor <- 2
    ans$start <- function(D) {
      corstart <- c(median(D)/2, 0.2)
      names(corstart) <- c("range","nugget")
      return(corstart)
    }
    ans$corr <- function(corrpar, D ){
      S <- (1 - corrpar[2])*(1 - 1.5*D/corrpar[1] +
                               0.5*(D/corrpar[1])^3)*(D<=corrpar[1]) + corrpar[2]*diag(NROW(D))
      return(S)
    }
   }else{
     ans$nug <- 0
     ans$npar.cor <- 1
     ans$start <- function(D) {
       corstart <- median(D)/2
       names(corstart) <- "range"
       return(corstart)
    }
    ans$corr <- function(corrpar, D ){
       S <- (1 - 1.5*D/corrpar[1] + 0.5*(D/corrpar[1])^3)*(D<=corrpar[1])
       return(S)
    }
   }
  }
  if(is.numeric(range) & is.numeric(nugget)){
    ans$S <- function(D) (1-nugget)*(1 - 1.5*D/range + 0.5*(D/range)^3)*(D<=range)+nugget*diag(NROW(D))
  }
  class(ans) <- c("corr.gc")
  return(ans)
}



#### Continuous Marginals for the purpose of correlation matching and data simulation

gaussian.gc <- function(mean = 0, sd = 1){
  q <- function(p) qnorm(p = p, mean = mean, sd = sd)
  ans <- list()
  ans$discrete <- FALSE
  ans$margvar <- sd^2
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rnorm(n = nrep, mean = mean, sd = sd)
  ans$q <- q
  class(ans) <- c("marginal.gc")
  return(ans)
}




gm.gc <- function(shape = 1, rate = 1){
  q <- function(p) qgamma(p = p, shape = shape, rate = rate)
  ans <- list()
  ans$discrete <- FALSE
  ans$margvar <- shape/(rate^2)
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rgamma(n = nrep, shape = shape, rate = rate)
  ans$q <- q
  class(ans) <- c("marginal.gc")
  return(ans)
}




beta.gc <- function(shape1 = 1, shape2 = 1){
  q <- function(p) qbeta(p = p, shape1 = shape1, shape2 = shape2)
  ans <- list()
  ans$discrete <- FALSE
  ans$margvar <- shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1))
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rbeta(n = nrep, shape1 = shape1, shape2 = shape2)
  ans$q <- q
  class(ans) <- c("marginal.gc")
  return(ans)
}





weibull.gc <- function(shape = 1, scale = 1){
  q <- function(p) qweibull(p = p, shape = shape, scale = scale)
  ans <- list()
  ans$discrete <- FALSE
  ans$margvar <- (scale^2)*(base::gamma(1+2/shape) - base::gamma(1+1/shape)^2)
  ans$int.marg <- function (order) {
    if (requireNamespace("EQL", quietly = TRUE)) {
      integrate(function(x, order)
        ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
        order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
    }else{
      stop("Please install {EQL} first!")
    }
  }
  ans$rsample <- function(nrep) rweibull(n = nrep, shape = shape, scale = scale)
  ans$q <- q
  class(ans) <- c("marginal.gc")
  return(ans)
}
