
#### function to create complete design matrix

create.design <- function(X, Z1, Z2, first.object, second.object, 
  m, subject, control, order.Z1, order.Z2, with.order) {

  design.X <- design.X.repar <- design.Z1 <- design.Z2 <- design.order <- c()
  
  I <- m * (m - 1)/2
  n <- length(first.object)
  
  design.help <- matrix(0, ncol = m, nrow = n)
  
  for (i in 1:n) {
    design.help[i, first.object[i]] <- 1
    design.help[i, second.object[i]] <- -1
  }
  
  if (!is.null(X)) {
    design.X <- create.design.X(design.help[, -m], X[subject,,drop = FALSE])
    design.X.repar <- create.design.X(design.help, X[subject,,drop = FALSE])
  }
  
  if (!is.null(Z1)) {
    p.Z1 <- ncol(Z1)/m
    index.Z1 <- (0:(m - 1)) * (m) + (1:(m))
    for (pp in 1:p.Z1) {
      Z1.help <- Z1[, ((pp - 1) * m) + (1:m)]
      design.help.Z1 <- create.design.X(design.help, Z1.help[subject, 
        order.Z1[, pp]])
      design.Z1 <- cbind(design.Z1, design.help.Z1[, index.Z1])
    }
  }
  
  if (!is.null(Z2)) {
    p.Z2 <- ncol(Z2)/m
    for (pp in 1:p.Z2) {
      Z2.help <- rowSums(design.help * (Z2[subject, ((pp - 
        1) * m) + (1:m)])[, order.Z2[, pp]])
      design.Z2 <- cbind(design.Z2, Z2.help)
    }
  }
  
  ### inclusion of order (home) effects
  order.effect <- control$order.effect
  object.order.effect <- control$object.order.effect
  
  if (order.effect | object.order.effect) {
    if (object.order.effect) {
      design.order <- c(design.help)
      if (control$order.center) {
        design.order[design.order == -1] <- 1
        design.order <- design.order * 0.5
      } else {
        design.order[design.order == -1] <- 0
      }
      design.order <- matrix(design.order, ncol = m)
    } else {
      design.order <- matrix(1, nrow = n)
    }
    design.order[!with.order,] <- 0
  }



  
  ## inclusion of intercepts
  if (control$include.intercepts) {
    design.help.repar <- design.help
    design.help <- design.help[, -m]
  } else {
    design.help <- design.help.repar <- c()
  }
  

  design <- cbind(design.order, design.help, design.X, design.Z1, 
    design.Z2)
  
  design.repar <- cbind(design.order, design.help.repar, design.X.repar, design.Z1, 
                  design.Z2)
  
  return(list(design = design, design.repar = design.repar))
}


create.design.X <- function(design.help, X) {

  design <- matrix(c(apply(X, 2, function(xx) {
    xx * c(design.help)
  })), nrow = nrow(design.help))
  
  return(design)
}
