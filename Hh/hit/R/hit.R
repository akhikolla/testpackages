#' @title Hierarchical Inference Testing
#'
#' @description Hierarchical inference testing for linear models with high-dimensional 
#' and/or correlated covariates by repeated sample splitting.
#'
#' @param x Design matrix of dimension \code{n * p}, without intercept. Variables not part 
#' of the dendrogram are added to the HO-model, see Details below.
#' @param y Quantitative response variable dimension \code{n}.
#' @param hierarchy Object of class \code{\link{as.hierarchy}}. Must include all variables 
#' of \code{x} which should be tested.
#' @param family Family of response variable distribution. Ether \code{y} is "gaussian" or 
#' "poisson" in which case \code{y} must be a vector or it is "binomial" distibuded and is 
#' either a vecror of zeros and ones, factor with two levels,  or a two-column matrix of 
#' counts or proportions. The second column is treated as the target class. For a factor, 
#' the last level in alphabetical order is the target class. For "binomial" if \code{y} is 
#' presented as a vector, it will be coerced into a factor.
#' @param B Number of sample-splits.
#' @param p.samp1 Fraction of data used for the LASSO. The hierachical ANOVA testing uses 
#' \code{1 - p.samp1}.
#' @param nfolds Number of folds (default is 10). See \code{\link[glmnet]{cv.glmnet}} for 
#' more details.
#' @param overall.lambda Logical, if true, lambda is estimated once, if false, lambda is 
#' estimated for each sample split.
#' @param lambda.opt Criterion for optimum selection of cross-validated lasso. Either 
#' "lambda.1se" (default) or "lambda.min". See \code{\link[glmnet]{cv.glmnet}} for more 
#' details.
#' @param alpha A single value in the range of 0 to 1 for the elastic net mixing parameter.
#' @param gamma Vector of gamma-values.
#' @param max.p.esti Maximum alpha level. All p-values above this value are set to one. 
#' Small \code{max.p.esti} values reduce computing time.
#' @param mc.cores Number of cores for parallelising. Theoretical maximum is 'B'. For 
#' details see \code{\link[parallel]{mclapply}}.
#' @param trace If TRUE it prints current status of the program.
#' @param ... Additional arguments for \code{\link[glmnet]{cv.glmnet}}.
#'
#' @details The H0-model contains variables, with are not tested, like experimental-design 
#' variables. These variables are not penalised in the LASSO model selection and are always 
#' include in the reduced ANOVA model.
#'
#' @references Mandozzi, J. and Buehlmann, P. (2013). \emph{Hierarchical testing in the 
#' high-dimensional setting with correlated variables}. To appear in the Journal of the 
#' American Statistical Association. Preprint arXiv:1312.5556
#'
#' @examples
#'
#' # Simulation:
#' set.seed(123)
#' n <- 80
#' p <- 82
#' ## x with correlated columns
#' corMat <- toeplitz((p:1/p)^2)
#' corMatQ <- chol(corMat)
#' x <- matrix(rnorm(n * p), nrow = n) %*% corMatQ
#' colnames(x) <- paste0("x", 1:p)
#' ## y
#' mu <- x[, c(5, 10, 72)] %*% c(2, -2, 2)
#' y <-  rnorm(n, mu)
#' ## clustering of the clumns of x
#' hc <- hclust(dist(t(x)))
#'
#' # HIT with AF
#' out <- hit(x, y, hc)
#' summary(out)
#'
#' @importFrom parallel mclapply
#' @importFrom glmnet cv.glmnet
#' @importFrom stats reorder
#' @export
hit <- function(x, y, hierarchy, family = "gaussian", B = 50, p.samp1 = 0.5, 
                nfolds = 10, overall.lambda = FALSE, lambda.opt = "lambda.1se", alpha = 1,
                gamma = seq(0.05, 0.99, length.out = 100), max.p.esti = 1,
                mc.cores = 1L, trace = FALSE, ...) {
  #   Mandozzi and Buehlmann (2015), 2 Description of method
  ## Checks args
  if (is.null(colnames(x)))
    stop("column names of 'x' are missing")
  n <- nrow(x)
  p <- ncol(x)
  if (class(hierarchy) == "dendrogram" || class(hierarchy) == "hclust")
    hierarchy <- as.hierarchy(hierarchy)
  stopifnot(class(hierarchy) == "hierarchy")
  lambda.opt <- match.arg(lambda.opt, c("lambda.1se", "lambda.min"))
  ##### Check family and assign test
  family <- match.arg(family, c("gaussian", "binomial", "poisson"))
  test <- "LRT"
  if (family == "gaussian")
    test <- "F"
  family2 <- eval(call(family))
  ##### Checks order of variats
  x.names <- colnames(x)
  hier.names <- names(hierarchy)
  if (!all(hier.names %in% x.names))
    stop("'hierarchy' includs variabels not in 'x'")
  if (identical(hier.names, x.names)) {
    x.notest <- integer(0L)
  } else if (setequal(hier.names, x.names)) {
    hierarchy <- reorder(hierarchy, x.names)
    x.notest <- integer(0L)
  } else {
    hierarchy <- reorder(hierarchy, x.names)
    x.notest <- match(setdiff(x.names, names(hierarchy)), x.names)
  }
  penalty.factor <- rep(1L, p)
  penalty.factor[x.notest] <- 0L
  ##### Sub-sample splits
  if (p.samp1 < .1 | p.samp1 > .9)
    stop("'p.samp1' must be between .1 and 0.9")
  n.samp1 <- as.integer(n * p.samp1)
  n.samp2 <- n - n.samp1
  allSamp1.ids <- replicate(B, sample.int(n, n.samp1), simplify = FALSE)
  ##  2.2 Screening
  if (trace)
    cat("LASSO has started at:\n\t", as.character(Sys.time()), "\n")
  if (overall.lambda) {
    lambdas <- samp1.lambda.overall(x, y, family, alpha, lambda.opt, nfolds,
                                    penalty.factor, n.samp2, ...)
    allActSet.ids <- mclapply(allSamp1.ids, samp1.lasso.overall,
                              x, y, family, alpha, lambdas,
                              penalty.factor, n.samp2, ...,
                              mc.cores = mc.cores, mc.cleanup = TRUE)
  } else {
    allActSet.ids <- mclapply(allSamp1.ids, samp1.lasso.split,
                              x, y, family, alpha, lambda.opt, nfolds,
                              penalty.factor, n.samp2, ...,
                              mc.cores = mc.cores, mc.cleanup = TRUE)
  }

  ##  2.3 Testing and multiplicity adjustmen; and
  ##  2.4 Aggregating and Hierarchical adjustment
  if (trace)
    cat("Significance testing has started at:\n\t",
        as.character(Sys.time()), "\n")
  pValues <- samp2.sigHierarchy(1L, x, y, hierarchy, family2, test, x.notest,
                                B, allSamp1.ids, allActSet.ids, 0,
                                max.p.esti, sort(gamma), 0L, mc.cores)
  ##### Results
  if (trace)
    cat("HIT has finished at:\n\t", as.character(Sys.time()), "\n")
  asi <- sort(unlist(allActSet.ids))
  sel.tab <- rep(0, p)
  sel.tab[unique(asi)] <- table(asi) / B
  tested <- rep(TRUE, p)
  tested[x.notest] <- FALSE
  names(tested) <- x.names
  out <- list(pValues = pValues,
              selectFreq = sel.tab,
              hierarchy = hierarchy,
              tested = tested,
              max.p.esti = max.p.esti)
  class(out) <- "hit"
  out
}


#' @title Cross-validation of LASSO penenlty lambda
#'
#' @description Cross-validation of LASSO penenlty lambda.
#'
#' @param x Design matrix, of dimension n x p.
#' @param y Vector of quantitative response variable.
#' @param family Distribution family of \code{y}.
#' @param alpha A single value n the range of 0 to 1 for the elastic net mixing parameter. 
#' @param lambda.opt Criterion for optimum selection of cross-validated lasso. Either 
#' "lambda.1se" (default) or "lambda.min". See \code{\link[glmnet]{cv.glmnet}} for more 
#' details.
#' @param nfolds Number of folds (default is 10). See
#' @param penalty.factor See glmnet.
#' @param n.samp2 Number of individuals in samp2 which is the max. for non zero 
#' coefficients.
#' @param ... Additional agruments.
#'
#' @importFrom parallel mclapply
#' @importFrom glmnet cv.glmnet
#' @keywords internal
samp1.lambda.overall <- function(x, y, family, alpha, lambda.opt, nfolds, 
                                 penalty.factor, n.samp2, ...) {
  foldid <- sample(rep(1L:nfolds, length.out = nrow(x)))
  suppressWarnings(
    cvfit <- cv.glmnet(x, y, family = family, foldid = foldid, alpha = alpha,
                       penalty.factor = penalty.factor, pmax = n.samp2 - 2L, ...)
  )
  if (lambda.opt == "lambda.min")
    optlambda <- cvfit$lambda[cvfit$lambda >= cvfit$lambda.min]
  else
    optlambda <- cvfit$lambda[cvfit$lambda >= cvfit$lambda.1se]
  optlambda
}



#' @title Variabel Screening
#'
#' @description LASSO function of the HIT algorithem.
#'
#' @param samp1 List of index for subsample (mclapply index).
#' @param x Design matrix, of dimension n x p.
#' @param y Vector of quantitative response variable.
#' @param family Distribution family of \code{y}.
#' @param alpha Mixing value for elnet.
#' @param lambda A vector of lambda values sorted from large to small where the smallest is 
#' the optimal value or NULL.
#' @param penalty.factor See glmnet.
#' @param n.samp2 Number of individuals in samp2 which is the max. for non zero 
#' coefficients.
#' @param ... Additional agruments.
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats coef
#' @keywords internal
samp1.lasso.overall <- function(samp1, x, y, family, alpha, lambda,
                                penalty.factor, n.samp2, ...) {
  suppressWarnings(
    fit <- glmnet(x[samp1, ], y[samp1], family = family,
                  alpha = alpha, lambda = lambda,
                  penalty.factor = penalty.factor, pmax = n.samp2 - 2L, ...)
  )
  beta <- coef(fit, s = lambda[length(lambda)])[-1L]
  actSet <- which(beta != 0 & penalty.factor == 1L)
  actSet
}

#' @title Variabel Screening
#'
#' @description LASSO function of the HIT algorithem.
#'
#' @param samp1 List of index for subsample (mclapply index).
#' @param x Design matrix, of dimension n x p.
#' @param y Vector of quantitative response variable.
#' @param family Distribution family of \code{y}.
#' @param alpha Mixing value for elnet.
#' @param lambda.opt Criterion for optimum selection of cross-validated lasso. Either 
#' "lambda.1se" (default) or "lambda.min". See \code{\link[glmnet]{cv.glmnet}} for more 
#' details.
#' @param nfolds Number of folds (default is 10). See
#' @param penalty.factor See glmnet.
#' @param n.samp2 Number of individuals in samp2 which is the max. for non zero 
#' coefficients.
#' @param ... Additional agruments.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#' @keywords internal
samp1.lasso.split <- function(samp1, x, y, family, alpha, lambda.opt, nfolds, 
                              penalty.factor, n.samp2, ...) {
  suppressWarnings(
    fit <- cv.glmnet(x[samp1, ], y[samp1], family = family,  nfolds = nfolds, 
                     alpha = alpha, penalty.factor = penalty.factor, 
                     pmax = n.samp2 - 2L, ...)
  )
  if (lambda.opt == "lambda.min")
    beta <- coef(fit, s = fit$lambda.min)[-1L]
  else
    beta <- coef(fit, s = fit$lambda.1se)[-1L]
  actSet <- which(beta != 0 & penalty.factor == 1L)
  actSet
}


#' @title Variabel Testing along the hierarchy
#'
#' @description ANOVA Testing, Multiplicity Adjustment, Aggregating and Hierarchical 
#' Adjustment
#'
#' @param j Index for cluster (mclapply index).
#' @param x Design matrix, of dimension n x p.
#' @param y Vector of quantitative response variable.
#' @param hierarchy A hierarchy object.
#' @param family Distribution family of \code{y}.
#' @param test name of test.
#' @param x.notest  Vector of indeces of non tested variabels.
#' @param B Number of sample-splits.
#' @param allSamp1.ids  List of subsampels.
#' @param allActSet.ids List of active sets.
#' @param upper.p P value upper huerarchy level of the clusters variables.
#' @param max.p.esti Maximum alpha level. All p-values above this value are set to one. 
#' Small max.p.esti values reduce computing time.
#' @param gamma Vector of gamma-values.
#' @param hl.count Hierarchy level counter for parallelism.
#' @param cores Number of cores for parallelising.
#'
#' @importFrom stats quantile
#' @keywords internal
samp2.sigHierarchy <- function(j, x, y, hierarchy, family, test, x.notest,
                               B, allSamp1.ids, allActSet.ids, upper.p,
                               max.p.esti, gamma, hl.count, cores) {
  ## 2.3 Testing and multiplicity adjustment
  cluster <- hierarchy[[j]]
  p.cluster <- sapply(1L:B, samp2.sigNode,
                      x, y, cluster, family, test, x.notest,
                      allSamp1.ids, allActSet.ids)
  ##  2.4 Aggregating and Hierarchical adjustment
  ### 2.4-1 Aggregating
  q.aggre <- sapply(gamma,
                    function(i, x) { min(1, quantile(x / i, i)) },
                    x = p.cluster)
  p.aggre <- min(1, (1 - log(min(gamma))) * min(q.aggre))
  ### 2.4-2 Hierarchical adjustment
  p.value <- max(p.aggre, upper.p)
  ##### Estimation at next lower level and find a way to parallelize
  if (!is.null(js <- attr(cluster, "subset"))) {
    if (p.value < max.p.esti) {
      if (hl.count == 0L && length(js) >= cores) {
        hl.count <- as.integer(sqrt(cores))
        pValues <- mclapply(js, samp2.sigHierarchy,
                            x, y, hierarchy, family, test, x.notest,
                            B, allSamp1.ids, allActSet.ids, p.value,
                            max.p.esti, gamma, hl.count + 1L, cores,
                            mc.cores = cores)
      } else if (hl.count <= as.integer(sqrt(cores))) {
        mc.cores <- ifelse(as.integer(sqrt(cores)) == 1L, 1L, 2L)
        pValues <- mclapply(js, samp2.sigHierarchy,
                            x, y, hierarchy, family, test, x.notest,
                            B, allSamp1.ids, allActSet.ids, p.value,
                            max.p.esti, gamma, hl.count + 1L, cores,
                            mc.cores = mc.cores, mc.allow.recursive = TRUE)
      } else {
        pValues <- lapply(js, samp2.sigHierarchy,
                          x, y, hierarchy, family, test, x.notest,
                          B, allSamp1.ids, allActSet.ids, p.value,
                          max.p.esti, gamma, hl.count + 1L, cores)
      }
    } else {
      pOne <- function(j) {
        if (!is.null(js <- attr(hierarchy[[j]], "subset")))
          return(c(1, sapply(js, pOne)))
        return(1)
      }
      pValues <-  sapply(js, pOne)
    }
  } else {
    pValues <- c()
  }
  out <- c(p.value, unlist(pValues))
  out
}


#' @title ANOVA testing and multiplicity adjustment
#'
#' @description ANOVA test (at single node) of the HIT algorithem.
#'
#' @param k Index for subsample (mclapply index).
#' @param x Design matrix, of dimension n x p.
#' @param y Vector of quantitative response variable.
#' @param cluster Clusters to be tested.
#' @param family Distribution family of \code{y}.
#' @param test name of test.
#' @param x.notest Vector of indeces of non tested variabels.
#' @param allSamp1.ids  List of subsampels.
#' @param allActSet.ids List of active sets.
#'
#' @keywords internal
samp2.sigNode <- function(k, x, y, cluster, family, test, x.notest,
                          allSamp1.ids, allActSet.ids) {
  ## 2.3 Testing and multiplicity adjustment
  n <- nrow(x)
  actClust <- intersect(allActSet.ids[[k]], cluster)
  nonActClust <- setdiff(allActSet.ids[[k]], cluster)
  ### 2.3-1 Testing
  p.cluster <- 1
  if (l.actClust <- length(actClust)) {
    ##### ANOVA between active set and active set minus cluster
    samp2 <- (1L:n)[-allSamp1.ids[[k]]]
    if (is.matrix(y))
      y <- y[samp2, ]
    else
      y <- y[samp2]
    if (l.nonActClust <- length(nonActClust)) {
      if (l.nonTested <- length(x.notest)) {
        x  <- cbind(1L, x[samp2, c(x.notest, nonActClust, actClust)])
        assign <- rep(0L:3L, c(1L, l.nonTested, l.nonActClust, l.actClust))
        get.p <- 3L
      } else {
        x  <- cbind(1L, x[samp2, c(nonActClust, actClust)])
        assign <- rep(0L:2L, c(1L, l.nonActClust, l.actClust))
        get.p <- 2L
      }
    } else {
      if (l.nonTested <- length(x.notest)) {
        x  <- cbind(1L, x[samp2, c(x.notest, actClust)])
        assign <- rep(0L:2L, c(1L, l.nonTested, l.actClust))
        get.p <- 2L
      } else {
        x  <- cbind(1L, x[samp2, actClust])
        assign <- rep(0L:1L, c(1L, l.actClust))
        get.p <- 1L
      }
    }
    p.test <- fast.anova(x, y, assign, family, test)[get.p]
    ### 2.3-2 Multiplicity adjustment
    p.cluster <- min(1, p.test * (l.actClust + l.nonActClust) / l.actClust)
  }
  p.cluster
}


#' @title Summary of HIT
#'
#' @description Significant clusters at alpha threshold.
#'
#' @param object A \code{\link{hit}} object.
#' @param alpha A alpha significance threshold.
#' @param max.height max. Height to consider.
#' @param ... Further arguments passed to or from other methods (not used).
#'
#' @method summary hit
#' @export
summary.hit <- function(object, alpha = 0.05, max.height, ...) {
  make.pVal <- function(i) {
    p.value <- object$pValues[i]
    inx <- object$hierarchy[[i]]
    height <- attr(object$hierarchy[[i]], "height")
    if (p.value <= alpha && 
        (all(is.na(P.CLUSTER[inx])) || 
         (all(!is.na(P.CLUSTER[inx])) && 
          all(P.CLUSTER[inx] <= alpha)))) {
      P.CLUSTER[inx] <<-  p.value
      ID.CLUSTER[inx] <<- COUNTER
      H.CLUSTER[inx] <<- height
      COUNTER <<- COUNTER + 1L
    }
    return(NULL)
  }
  stopifnot(inherits(object, "hit"))
  stopifnot(alpha <= 0.999)
  if (missing(max.height))
    max.height <- attr(object$hierarchy[[1L]], "height")
  P.CLUSTER <- rep(NA_real_, length(object$hierarchy[[1L]]))
  ID.CLUSTER <- rep(NA_integer_, length(object$hierarchy[[1L]]))
  H.CLUSTER <- rep(NA_real_, length(object$hierarchy[[1L]]))
  COUNTER <- 1L
  lapply(length(object$hierarchy):1L, make.pVal)
  non.na <- which(!is.na(ID.CLUSTER))
  out <- data.frame(clusters = ID.CLUSTER[non.na],
                    heights = H.CLUSTER[non.na],
                    pValues = P.CLUSTER[non.na])
  rownames(out) <- names(object$hierarchy[[1L]])[non.na]
  out <- out[out$heights <= max.height, ]
  if (ll <- length(unique(out[, 1L])))
    out[, 1L] <- as.integer(factor(out[, 1L], labels = 1L:ll))
  out
}
