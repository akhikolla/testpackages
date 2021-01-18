#### create cumulative response vector
cumul.response <- function(Y) {
  get.resp <- function(x) {
    as.numeric(as.numeric(x) <= 1:(Y$q))
  }
  
  cum.resp <- c()
  for (i in 1:length(Y$response)) {
    cum.resp <- c(cum.resp, get.resp(Y$response[i]))
  }
  
  cum.resp
}



#### check functions for matrices Z1 and Z2

fun.check <- function(x, object.names) {
  identical(sort(x), sort(object.names))
}

check <- function(Z, object.names, subject) {
  
  ## which Z matrix are we talking about
  Z.name <- deparse(substitute(Z))
  
  ## names of the objects
  m <- length(object.names)
  
  ## check if ncol(Z) is a multiple of m
  if (!is.matrix(Z)) 
    stop(paste(Z.name, "has to be a matrix"))
  if (ncol(Z)%%m != 0) 
    stop(paste("Number of columns of", Z.name, "has to be a multiple of the number of objects"))
  
  ## names of objects in Z matrix
  objects.Z <- matrix(word(colnames(Z), 2, sep = fixed(".")), 
    ncol = m, byrow = TRUE)
  
  ## number of Z-covariates
  p.Z <- nrow(objects.Z)
  
  ## order of the objects compared to general order from
  ## response
  order.Z <- apply(objects.Z, 1, order)
  
  ## names of the Z-covariates
  var.Z <- matrix(word(colnames(Z), 1, sep = fixed(".")), ncol = m, 
    byrow = TRUE)
  
  
  ## 1. names of the objects have to be identical 2. within Z,
  ## objects have to be equal for all covariates, in the same
  ## order
  if (!all(apply(objects.Z, 1, fun.check, object.names)) | 
    !all(apply(var.Z[, -1, drop = FALSE], 2, function(x) {
      identical(unique(c(x)), unique(c(var.Z[, 1, drop = FALSE])))
    }))) {
    text <- paste("The matrix", Z.name, "has to be arranged so that all columns 
        corresponding to one covariate are put next to each other.\n
        The colnames of", 
      Z.name, "have to be named according to the scheme 
        'firstvar.object1',...,'firstvar.objectm',...,'lastvar.objectm'.\n The object names
        'object1',...,'objectm' have to be identical to the object names used in the response object.\n
        The variable names and the object names have to be separated by '.'.")
    
    stop(text)
  }
  
  ## every subject needs to have its individual line in the Z
  ## matrix, identified with rownames
  if (!all(subject %in% rownames(Z))) {
    text <- paste("The rownames of the matrix", Z.name, "have to be equal to the subjects specified in the response object.")
    stop(text)
  }
  
  
  list(vars.Z = unique(c(var.Z)), order.Z = order.Z)
}

#### log likelihood function
loglik <- function(coef, y, design, kat) {
  eta <- matrix(design %*% coef, ncol = kat - 1, byrow = TRUE)
  pi.help <- matrix(exp(eta)/(1 + exp(eta)), ncol = kat - 1)
  # print(pi.help)
  
  pi <- pi.help
  if (kat > 2) {
    for (i in 2:(kat - 1)) {
      pi[, i] <- pi.help[, i] - pi.help[, i - 1]
    }
  }
  pi <- cbind(pi, 1 - pi.help[, kat - 1])
  
  
  yhelp <- rep(y, each = kat)
  yhelp <- matrix(as.numeric(yhelp == rep(1:kat, length(y))), 
    byrow = T, ncol = kat)
  
  
  loglik <- sum(yhelp * log(pi))
  loglik
}

#### predict function
predBTLLasso <- function(coef, q, design) {
  k <- q + 1
  
  eta <- matrix(design %*% coef, ncol = q, byrow = TRUE)
  pi.help <- matrix(exp(eta)/(1 + exp(eta)), ncol = q)
  # print(pi.help)
  
  pi <- pi.help
  if (k > 2) {
    for (i in 2:(q)) {
      pi[, i] <- pi.help[, i] - pi.help[, i - 1]
    }
  }
  pi <- cbind(pi, 1 - pi.help[, q])
  
  c(t(pi))
}

#### function to reparameterize from reference category to
#### symmetric side constraints
reparam <- function(x) {

  z <- ncol(x) + 1
  K <- matrix((-1/z), ncol = z - 1, nrow = z - 1)
  diag(K) <- (z - 1)/z
  x.sym <- x %*% K
  
  
  x2 <- x.sym %*% matrix(rep(-1, z - 1), nrow = z - 1)
  
  x.sym <- cbind(x.sym, x2)
  
  x.sym
}

#### function to create complete coefficient matrix either
#### together with zero-columns for reference categories or with
#### symmetric side constraints
expand.coefs <- function(coef, D, Y, symmetric = TRUE, name.order = "Order") {
  
  n.theta <- D$n.theta
  n.order <- D$n.order
  n.intercepts <- D$n.intercepts
  p.X <- D$p.X
  p.Z1 <- D$p.Z1
  p.Z2 <- D$p.Z2
  m <- Y$m
  object.names <- Y$object.names
  vars.X <- D$vars.X
  
  ## initialize new coefficient vector
  coef.new <- c()
  
  ## if threshold parameters exist, leave them as they are
  if (n.theta > 0) {
    coef.new <- coef[, 1:n.theta, drop = FALSE]
    colnames(coef.new) <- paste("theta", 1:n.theta, sep = ".")
  }
  
  ## if order parameters exist, leave them as they are
  if (n.order > 0) {
    order.effects <- coef[, (n.theta + 1):(n.theta + n.order), 
      drop = FALSE]
    if(n.order==m){
      colnames(order.effects) <- paste(name.order, object.names, sep=".")
    }else{
      colnames(order.effects) <- name.order
    }
    
    coef.new <- cbind(coef.new, order.effects)
  }
  
  ## if intercepts parameters exist, reparameterize them
  if (n.intercepts > 0) {
    intercepts <- coef[, (n.theta + n.order + 1):(n.theta + 
      n.order + n.intercepts), drop = FALSE]
    if (symmetric) {
      intercepts.new <- reparam(intercepts)
    } else {
      intercepts.new <- cbind(intercepts, 0)
    }
    colnames(intercepts.new) <- object.names
    coef.new <- cbind(coef.new, intercepts.new)
  }
  
  ## if parameters for subject-specific covariates exist,
  ## reparameterize them
  if (p.X > 0) {
    gamma <- coef[, (n.theta + n.intercepts + n.order + 1):(n.theta + 
      n.intercepts + n.order + p.X * (m - 1)), drop = FALSE]
    index <- 1
    for (i in 1:p.X) {
      if (symmetric) {
        coef.X <- reparam(gamma[, index:(index + m - 
          2), drop = FALSE])
      } else {
        coef.X <- cbind(0, gamma[, index:(index + m - 
          2), drop = FALSE])
      }
      colnames(coef.X) <- paste(rep(vars.X[i], each = m), 
        object.names, sep = ".")
      coef.new <- cbind(coef.new, coef.X)
      index <- index + m - 1
    }
  }
  
  ## if parameters for object-specific covariates exist, leave
  ## them as they are
  if (p.Z1 + p.Z2 > 0) {
    rest <- coef[, (n.theta + n.intercepts + n.order + p.X * 
      (m - 1) + 1):ncol(coef), drop = FALSE]
    if (p.Z1 > 0) {
      colnames(rest) <- c(paste(rep(D$vars.Z1, each = m), 
        object.names, sep = "."), D$vars.Z2)
    } else {
      colnames(rest) <- D$vars.Z2
    }
    
    coef.new <- cbind(coef.new, rest)
  }
  
  ## return reparameterized parameter vector
  coef.new
}


#### create subset of a response object for cross validation
subsetY <- function(Y, id.ex) {
  Y$response <- Y$response[-id.ex]
  Y$first.object <- Y$first.object[-id.ex]
  Y$second.object <- Y$second.object[-id.ex]
  Y$subject <- Y$subject[-id.ex]
  Y$subject.names <- levels(as.factor(Y$subject))
  # Y$subject.names <- Y$subject.names[-id.ex]
  Y$n <- length(Y$subject.names)
  
  Y
}

bootY <- function(Y, id.vec) {
  Y$response <- Y$response[id.vec]
  Y$first.object <- Y$first.object[id.vec]
  Y$second.object <- Y$second.object[id.vec]
  Y$subject <- Y$subject[id.vec]
  Y$subject.names <- levels(as.factor(Y$subject))
  Y$n <- length(Y$subject.names)
  
  Y
}