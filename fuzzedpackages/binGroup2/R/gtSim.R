##################################################################
# sim.gt() function                                              #
##################################################################
sim.gt <- function (x = NULL, gshape = 20, gscale = 2, par,
                    linkf = c("logit", "probit", "cloglog"),
                    sample.size, group.size, sens = 1, spec = 1, 
                    sens.ind = NULL, spec.ind = NULL){
  if (is.null(sens.ind))
    sens.ind <- sens
  if (is.null(spec.ind))
    spec.ind <- spec
  if (is.null(x)) {
    x <- rgamma(n = sample.size, shape = gshape, scale = gscale)
    X <- cbind(1, x)
  }
  else {
    X <- cbind(1, x)
    sample.size <- nrow(X)
  }
  linkf <- match.arg(linkf)
  pijk <- switch(linkf, logit = plogis(X %*% par),
                 probit = pnorm(X %*% par),
                 cloglog = 1 - exp(-exp(X %*% par)))
  ind <- rbinom(n = sample.size, size = 1, prob = pijk)
  num.g <- ceiling(sample.size/group.size)
  vec <- 1:sample.size
  groupn <- rep(1:num.g, each = group.size)[vec]
  save.sum <- tapply(ind, groupn, sum)
  save.group <- as.vector(ifelse(save.sum > 0, 1, 0))
  save.obs <- rep(NA, num.g)
  ret <- rep(NA, sample.size)
  for (i in 1:num.g)
    save.obs[i] <- ifelse(save.group[i] == 1, rbinom(1, 1, sens),
                          1 - rbinom(1, 1, spec))
  gres <- rep(save.obs, each = group.size)[vec]
  for (i in vec) {
    if (gres[i] == 1)
      ret[i] <- ifelse(ind[i] == 1, rbinom(1,
                                           1, sens.ind), 1 - rbinom(1, 1, spec.ind))
  }   
  grd <- data.frame(gres = gres, x = x, groupn = groupn, ind = ind, retest = ret)
  if (ncol(X) > 2)
    for (i in 2:ncol(X))
      colnames(grd)[i] <- paste("x", i - 1, sep="")
  grd
}



##################################################################
# sim.halving() function                                         #
##################################################################
sim.halving <- function (x = NULL, gshape = 20, gscale = 2, par,
                         linkf = c("logit", "probit", "cloglog"),
                         sample.size, group.size, sens = 1, spec = 1,
                         sens.ind = NULL, spec.ind = NULL){
  if (is.null(sens.ind))
    sens.ind <- sens
  if (is.null(spec.ind))
    spec.ind <- spec
  if (is.null(x)) {
    x <- rgamma(n = sample.size, shape = gshape, scale = gscale)
    X <- cbind(1, x)
  }
  else {
    X <- cbind(1, x)
    sample.size <- nrow(X)
  }
  linkf <- match.arg(linkf)
  pijk <- switch(linkf, logit = plogis(X %*% par),
                 probit = pnorm(X %*% par),
                 cloglog = 1 - exp(-exp(X %*% par)))
  ind <- rbinom(n = sample.size, size = 1, prob = pijk)
  num.g <- ceiling(sample.size/group.size)
  vec <- 1:sample.size
  groupn <- rep(1:num.g, each = group.size)[vec]
  save.sum <- tapply(ind, groupn, sum)
  save.group <- as.vector(ifelse(save.sum > 0, 1, 0))
  save.obs <- rep(NA, num.g)
  subgroup <- ret <- rep(NA, sample.size)
  for (grn in 1:num.g) {
    vec1 <- vec[groupn == grn]
    gs <- length(vec1)
    save.obs[grn] <- ifelse(save.group[grn] == 1, rbinom(1, 1, sens),
                            1 - rbinom(1, 1, spec))
    if (save.obs[grn] == 1) {
      sub1 <- vec1[1:ceiling(gs/2)]
      sub2 <- vec1[(ceiling(gs/2) + 1):gs]
      tZ1 <- sum(ind[sub1])
      tZ2 <- sum(ind[sub2])
      Z1 <- ifelse(tZ1 == 1, rbinom(1,
                                    1, sens), 1 - rbinom(1, 1, spec))
      Z2 <- ifelse(tZ2 == 1, rbinom(1,
                                    1, sens), 1 - rbinom(1, 1, spec))
      if (Z1 == 1) {
        for (i1 in sub1) {
          ret[i1] <- ifelse(ind[i1] == 1, rbinom(1,
                                                 1, sens.ind), 1 - rbinom(1, 1, spec.ind))
        }
      }
      if (Z2 == 1) {
        for (i1 in sub2) {
          ret[i1] <- ifelse(ind[i1] == 1, rbinom(1,
                                                 1, sens.ind), 1 - rbinom(1, 1, spec.ind))
        }
      }
      subgroup[sub1] <- Z1
      subgroup[sub2] <- Z2
    }
  }
  gres <- rep(save.obs, each = group.size)[vec]
  grd <- data.frame(gres = gres, x = x, groupn = groupn, ind = ind,
                    retest = ret, subgroup = subgroup)
  if (ncol(X) > 2)
    for (i in 2:ncol(X))
      colnames(grd)[i] <- paste("x", i - 1, sep="")
  grd
}


##################################################################
# sim.mp() function                                              #
##################################################################
# "mp" refers to matrix pooling, another name for array testing

sim.mp <- function (x = NULL, gshape = 20, gscale = 2, par,
                    linkf = c("logit", "probit", "cloglog"),
                    n.row, n.col, sens = 1, spec = 1, 
                    sens.ind = NULL, spec.ind = NULL){
  if (is.null(sens.ind))
    sens.ind <- sens
  if (is.null(spec.ind))
    spec.ind <- spec
  if (length(n.row) != length(n.col))
    stop("vector n.row and n.col must have the same length")
  linkf <- match.arg(linkf)
  if (is.null(x)) {
    sample.size <- sum(n.col * n.row)
    x <- rgamma(n = sample.size, shape = gshape, scale = gscale)
    X <- cbind(1, x)
  }
  else {
    X <- cbind(1, x)
    sample.size <- nrow(X)
    if (sum(n.col * n.row) != sample.size)
      stop("n.row and n.col not consistent with the sample size")
  }
  len <- length(n.row)
  pijk <- switch(linkf, logit = plogis(X %*% par),
                 probit = pnorm(X %*% par),
                 cloglog = 1 - exp(-exp(X %*% par)))
  ind <- rbinom(n = sample.size, size = 1, prob = pijk)
  individual <- col.groupn <- row.groupn <- numeric(0)
  rowr <- colr <- numeric(0)
  ret <- rep(NA, sample.size)
  for (i in 1:len) {
    if (i > 1)
      index <- seq(max(index) + 1, length = (n.row * n.col)[i])
    else index <- 1:(n.row * n.col)[1]
    indm <- matrix(ind[index], nrow = n.row[i])
    col.resp <- apply(indm, MARGIN = 2, FUN = sum)
    col.resp <- ifelse(col.resp > 0, 1, 0)
    col.err <- rep(NA, n.col[i])
    for (j in 1:n.col[i])
      col.err[j] <- ifelse(col.resp[j] == 1, rbinom(1, 1, sens),
                           1 - rbinom(1, 1, spec))
    row.resp <- apply(indm, MARGIN = 1, FUN = sum)
    row.resp <- ifelse(row.resp > 0, 1, 0)
    row.err <- rep(NA, n.row[i])
    for (j in 1:n.row[i])
      row.err[j] <- ifelse(row.resp[j] == 1, rbinom(1, 1, sens),
                           1 - rbinom(1, 1, spec))
    temp.c <- rep(1:n.col[i], each = n.row[i])
    col.groupn <- c(col.groupn, temp.c)
    temp.r <- rep(1:n.row[i], n.col[i])
    row.groupn <- c(row.groupn, temp.r)
    temp2.c <- rep(col.err, each = n.row[i])
    colr <- c(colr, temp2.c)
    temp2.r <- rep(row.err, n.col[i])
    rowr <- c(rowr, temp2.r)
    if (all(row.err == 0)) {
      for (j in index) {
        if (colr[j] == 1)
          ret[j] <- ifelse(ind[j] == 1, rbinom(1,
                                               1, sens.ind), 1 - rbinom(1, 1, spec.ind))
      }
    }
    else {
      if (all(col.err == 0)) {
        for (j in index) {
          if (rowr[j] == 1)
            ret[j] <- ifelse(ind[j] == 1, rbinom(1,
                                                 1, sens.ind), 1 - rbinom(1, 1, spec.ind))
        }
      }
      else {
        for (j in index) {
          if (rowr[j] == 1 && colr[j] == 1)
            ret[j] <- ifelse(ind[j] == 1, rbinom(1, 1,
                                                 sens.ind), 1 - rbinom(1, 1, spec.ind))
        }
      }
    }
    individual <- c(individual, list(indm))
  }
  sq <- rep(1:len, n.col * n.row)
  if (all(colr == 0) && all(rowr == 0))
    return(NULL)
  grd <- data.frame(x = x, col.resp = colr,
                    row.resp = rowr, coln = col.groupn, rown = row.groupn,
                    arrayn = sq, retest = ret)
  if (ncol(X) > 2)
    for (i in 1:(ncol(X) - 1))
      colnames(grd)[i] <- paste("x", i, sep="")
  list(dframe = grd, ind = individual, prob = as.vector(pijk))
}





##################################################################
# gtSim() function                                               #
##################################################################

#' @title Simulation function for group testing data
#' 
#' @description Simulates data in group testing form ready to be 
#' fit by \code{\link{gtReg}}. 
#' 
#' @param type \kbd{"sp"} for simple pooling, \kbd{"halving"} for 
#' halving protocol, and \kbd{"array"} for array testing 
#' (also known as matrix pooling).
#' @param x a matrix of user-submitted covariates with which to 
#' simulate the data. Default is NULL, in which case a gamma 
#' distribution is used to generate the covariates automatically.
#' @param gshape shape parameter for the gamma distribution. The 
#' value must be non-negative. Default value is set to 20.
#' @param gscale scale parameter for the gamma distribution. The 
#' value must be strictly positive. Default value is set to 2.
#' @param par the true coefficients in the linear predictor.
#' @param size1 sample size of the simulated data (for use with 
#' \kbd{"sp"} and \kbd{"halving"} methods) or a vector that specifies the 
#' number of rows in each matrix (for use with \kbd{"array"} method). If 
#' only one matrix is simulated, this value is a scalar.
#' @param size2 group size in pooling individual samples (for use 
#' with \kbd{"sp"} and \kbd{"halving"} methods) or a vector that specifies 
#' the number of columns in each matrix (for use with \kbd{"array"} method). 
#' If only one matrix is simulated, this value is a scalar. 
#' @param linkf a character string specifying one of the three 
#' link functions to be used: \kbd{"logit"} (default), \kbd{"probit"}, or 
#' \kbd{"cloglog"}.
#' @param sens sensitivity of the group tests. Default value is 
#' set to 1.
#' @param spec specificity of the group tests. Default value is 
#' set to 1.
#' @param sens.ind sensitivity of the individual retests. If NULL, 
#' set to be equal to sens.
#' @param spec.ind specificity of the individual retests. If NULL, 
#' set to be equal to spec. 
#' 
#' @details Generates group testing data in simple pooling form 
#' (\kbd{type="sp"}), for the halving protocol (\kbd{type="halving"}), 
#' or in array testing form (\kbd{type="array"}). 
#' The covariates are either specified by the x argument or they 
#' are generated from a gamma distribution with the given gshape 
#' and gscale parameters. The individual probabilities are 
#' calculated from the covariates, the coefficients given in par, 
#' and the link function specified through linkf. The true binary
#' individual responses are then simulated from the individual 
#' probabilities. 
#' 
#' The true group responses are found from the 
#' individual responses within the groups (i.e., if at least one 
#' response is positive, the group is positive; otherwise, the group 
#' response is negative). Finally, the observed group responses
#' are simulated using the ginve sens and spec. Individual retests
#' are simulated from sens.ind and spec.ind for samples in 
#' observed positive groups. Note that with a given group size 
#' (specified by size2 with method="sp" or method="halving"), 
#' the last group may have fewer individuals.
#' 
#' The true binary individual responses are then 
#' simulated from the individual probabilities. The group, subgroup, 
#' and individual retests are simulated using the given sens 
#' and spec under the halving protocol.
#' 
#' The true binary individual responses are then simulated from 
#' the individual probabilites. The individuals are organized into
#' (by column) one or more matrices specified by n.row and n.col, 
#' and the true group responses are found (i.e., if at least one 
#' response is positive, the group is positive; otherwise, the 
#' group response is negative). The observed row and column group 
#' responses are then simulated using the given sens and spec 
#' values. Individual retests are simulated from sens.ind and 
#' spec.ind for individuals that lie on the intersection of an 
#' observed positive row and and observed positive column. In 
#' the case where no column (row) tests positive in a matrix, 
#' all the individuals in any observed positive rows (columns) 
#' will be assigned a simulated retest result. If no column or 
#' row is observed positive, NULL is returned.
#' 
#' @return For simple pooling (\kbd{type="sp"}) and the halving 
#' protocol (\kbd{type="array"}), a data frame or for array 
#' testing (\kbd{type="array"}), a list, which may include the following:
#' \item{gres}{the group response, for simple pooling and the halving 
#' protocol only.}
#' \item{col.resp}{the column group response, for array testing only.}
#' \item{row.resp}{the row group response, for array testing only.}
#' \item{x}{the covariate.}
#' \item{groupn}{the group number, for simple pooling and the halving 
#' protocol only.}
#' \item{arrayn}{the array number, for array testing only.}
#' \item{coln}{the column group number, for array testing only.}
#' \item{rown}{the row group number, for array testing only.}
#' \item{ind}{the true individual responses. For simple pooling and the 
#' halving protocol, these are included in the data frame of results. 
#' For array testing, these are included in the list of results, with 
#' individual responses presented in matrices.}
#' \item{retest}{the results of individual retests.}
#' \item{subgroup}{the subgroup number, for the halving protocol.}
#' \item{prob}{the individual probabilities, for array testing only.}
#' 
#' @author This function is a combination of \kbd{sim.gt}, \kbd{sim.halving}, 
#' and \kbd{sim.mp} written by Boan Zhang for the \code{binGroup} package. 
#' Minor modifications have been made for inclusion of the functions in the 
#' \code{binGroup2} package.
#' 
#' @seealso \code{\link{gtReg}} to fit simulated group testing data.
#' 
#' @examples 
#' set.seed(46)
#' gt.data <- gtSim(type="sp", par=c(-12, 0.2), 
#'                  size1=700, size2=5)
#' 
#' x1 <- sort(runif(100, 0, 30))
#' x2 <- rgamma(100, shape=17, scale=1.5)
#' gt.data <- gtSim(type="sp", x=cbind(x1, x2), 
#'                  par=c(-14, 0.2, 0.3), size2=4, 
#'                  sens=0.98, spec=0.98)
#'                    
#' set.seed(46)
#' gt.data <- gtSim(type="halving", par=c(-6, 0.1), 
#'                  gshape=17, gscale=1.4, size1=5000, 
#'                  size2=5, sens=0.95, spec=0.95)
#'                    
#' # 5x6 and 4x5 matrix
#' set.seed(9128)
#' sa1a <- gtSim(type="array", par=c(-7, 0.1), 
#'               size1=c(5, 4), size2=c(6, 5), 
#'               sens=0.95, spec=0.95)
#' sa1a$dframe

# Brianna Hitt - 01-06-2020

gtSim <- function(type="sp", x=NULL, gshape=20, gscale=2, par, 
                    linkf=c("logit", "probit", "cloglog"),
                    size1, size2, sens=1, spec=1, 
                    sens.ind=NULL, spec.ind=NULL){
  
  if(type=="sp"){
    # for use with gtreg 
    results <- sim.gt(x=x, gshape=gshape, gscale=gscale, par=par,
                      linkf=linkf, sample.size=size1, group.size=size2, 
                      sens=sens, spec=spec, 
                      sens.ind=sens.ind, spec.ind=spec.ind)
  } else if(type=="halving"){
    # for use with gtreg.halving
    results <- sim.halving(x=x, gshape=gshape, gscale=gscale, par=par,
                           linkf=linkf, sample.size=size1, group.size=size2, 
                           sens=sens, spec=spec,
                           sens.ind=sens.ind, spec.ind=spec.ind)
  } else if(type=="array"){
    # for use with gtreg.mp
    results <- sim.mp(x=x, gshape=gshape, gscale=gscale, par=par,
                      linkf=linkf, n.row=size1, n.col=size2, 
                      sens=sens, spec=spec, 
                      sens.ind=sens.ind, spec.ind=spec.ind)
  }
  results
}

#