#' Fitting vennLasso models
#'
#' @param x input matrix of dimension nobs by nvars. Each row is an observation,
#' each column corresponds to a covariate
#' @param y numeric response vector of length nobs
#' @param groups A list of length equal to the number of groups containing vectors of integers
#' indicating the variable IDs for each group. For example, \code{groups = list(c(1,2), c(2,3), c(3,4,5))} specifies
#' that Group 1 contains variables 1 and 2, Group 2 contains variables 2 and 3, and Group 3 contains
#' variables 3, 4, and 5. Can also be a matrix of 0s and 1s with the number of columns equal to the
#' number of groups and the number of rows equal to the number of variables. A value of 1 in row i and
#' column j indicates that variable i is in group j and 0 indicates that variable i is not in group j.
#' @param family \code{"gaussian"} for least squares problems, \code{"binomial"} for binary response, 
#' and \code{"coxph"} for time-to-event outcomes (not yet available)
#' @param nlambda The number of lambda values. Default is 100.
#' @param lambda A user-specified sequence of lambda values. Left unspecified, the a sequence of lambda values is
#' automatically computed, ranging uniformly on the log scale over the relevant range of lambda values.
#' @param lambda.min.ratio Smallest value for lambda, as a fraction of \code{lambda.max}, the (data derived) entry
#' value (i.e. the smallest value for which all parameter estimates are zero). The default
#' depends on the sample size \code{nobs} relative to the number of variables \code{nvars}. If
#' \code{nobs > nvars}, the default is 0.0001, close to zero. If \code{nobs < nvars}, the default
#' is 0.01. A very small value of \code{lambda.min.ratio} can lead to a saturated fit
#' when \code{nobs < nvars}.
#' @param lambda.fused tuning parameter for fused lasso penalty
#' @param penalty.factor vector of weights to be multiplied to the tuning parameter for the
#' group lasso penalty. A vector of length equal to the number of groups
#' @param group.weights A vector of values representing multiplicative factors by which each group's penalty is to
#' be multiplied. Often, this is a function (such as the square root) of the number of predictors in each group.
#' The default is to use the square root of group size for the group selection methods.
#' @param adaptive.lasso Flag indicating whether or not to use adaptive lasso weights. If set to \code{TRUE} and
#' \code{group.weights} is unspecified, then this will override \code{group.weights}. If a vector is supplied to group.weights,
#' then the \code{adaptive.lasso} weights will be multiplied by the \code{group.weights} vector
#' @param adaptive.fused Flag indicating whether or not to use adaptive fused lasso weights. 
#' @param standardize Should the data be standardized? Defaults to \code{FALSE}.
#' @param intercept Should an intercept be fit? Defaults to \code{TRUE}
#' @param one.intercept Should a single intercept be fit for all subpopulations instead of one
#' for each subpopulation? Defaults to \code{FALSE}.
#' @param compute.se Should standard errors be computed? If \code{TRUE}, then models are re-fit with no penalization and the standard
#' errors are computed from the refit models. These standard errors are only theoretically valid for the
#' adaptive lasso (when \code{adaptive.lasso} is set to \code{TRUE})
#' @param conf.int level for confidence intervals. Defaults to \code{NULL} (no confidence intervals). Should be a value between 0 and 1. If confidence
#' intervals are to be computed, compute.se will be automatically set to \code{TRUE}
#' @param gamma power to raise the MLE estimated weights by for the adaptive lasso. defaults to 1
#' @param rho ADMM parameter. must be a strictly positive value. By default, an appropriate value is automatically chosen
#' @param dynamic.rho \code{TRUE}/\code{FALSE} indicating whether or not the rho value should be updated throughout the course of the ADMM iterations
#' @param maxit integer. Maximum number of ADMM iterations. Default is 500.
#' @param abs.tol absolute convergence tolerance for ADMM iterations for the relative dual and primal residuals.
#' Default is \code{10^{-5}}, which is typically adequate.
#' @param rel.tol relative convergence tolerance for ADMM iterations for the relative dual and primal residuals.
#' Default is \code{10^{-5}}, which is typically adequate.
#' @param irls.maxit integer. Maximum number of IRLS iterations. Only used if \code{family != "gaussian"}. Default is 100.
#' @param irls.tol convergence tolerance for IRLS iterations. Only used if \code{family != "gaussian"}. Default is 10^{-5}.
#' @param model.matrix logical flag. Should the design matrix used be returned?
#' @param ... not used
#' @return An object with S3 class "vennLasso"
#'
#' @useDynLib vennLasso, .registration=TRUE
#'
#' @import Rcpp
#'
#' @export
#' @examples
#' library(Matrix)
#' 
#' # first simulate heterogeneous data using
#' # genHierSparseData
#' set.seed(123)
#' dat.sim <- genHierSparseData(ncats = 2, nvars = 25,
#'                              nobs = 200, 
#'                              hier.sparsity.param = 0.5,
#'                              prop.zero.vars = 0.5,
#'                              family = "gaussian")
#'
#' x          <- dat.sim$x
#' conditions <- dat.sim$group.ind
#' y          <- dat.sim$y
#' 
#' true.beta.mat <- dat.sim$beta.mat
#' 
#' fit <- vennLasso(x = x, y = y, groups = conditions)
#'
#' (true.coef <- true.beta.mat[match(dimnames(fit$beta)[[1]], rownames(true.beta.mat)),])
#' round(fit$beta[,,21], 2)
#'
#' ## fit adaptive version and compute confidence intervals
#' afit <- vennLasso(x = x, y = y, groups = conditions, conf.int = 0.95, adaptive.lasso = TRUE)
#'
#' (true.coef <- true.beta.mat[match(dimnames(fit$beta)[[1]], rownames(true.beta.mat)),])[,1:10]
#' round(afit$beta[,1:10,28], 2)
#' round(afit$lower.ci[,1:10,28], 2)
#' round(afit$upper.ci[,1:10,28], 2)
#'
#' aic.idx <- which.min(afit$aic)
#' bic.idx <- which.min(afit$bic)
#'
#' # actual coverage
#' # actual coverage
#' mean(true.coef[afit$beta[,-1,aic.idx] != 0] >= 
#'              afit$lower.ci[,-1,aic.idx][afit$beta[,-1,aic.idx] != 0] &
#'          true.coef[afit$beta[,-1,aic.idx] != 0] <= 
#'              afit$upper.ci[,-1,aic.idx][afit$beta[,-1,aic.idx] != 0])
#'
#' (covered <- true.coef >= afit$lower.ci[,-1,aic.idx] & true.coef <= afit$upper.ci[,-1,aic.idx])
#' mean(covered)
#'
#'
#' # logistic regression example
#' \dontrun{
#' set.seed(123)
#' dat.sim <- genHierSparseData(ncats = 2, nvars = 25,
#'                              nobs = 200, 
#'                              hier.sparsity.param = 0.5,
#'                              prop.zero.vars = 0.5,
#'                              family = "binomial",
#'                              effect.size.max = 0.5) # don't make any 
#'                                                     # coefficients too big
#'
#' x           <- dat.sim$x
#' conditions  <- dat.sim$group.ind
#' y           <- dat.sim$y
#' true.beta.b <- dat.sim$beta.mat
#'
#' bfit <- vennLasso(x = x, y = y, groups = conditions, family = "binomial")
#'
#' (true.coef.b <- -true.beta.b[match(dimnames(fit$beta)[[1]], rownames(true.beta.b)),])
#' round(bfit$beta[,,20], 2)
#' }
#'
vennLasso <- function(x, y,
                      groups,
                      family           = c("gaussian", "binomial"),
                      nlambda          = 100L,
                      lambda           = NULL,
                      lambda.min.ratio = NULL,
                      lambda.fused     = NULL,
                      penalty.factor   = NULL,
                      group.weights    = NULL,
                      adaptive.lasso   = FALSE,
                      adaptive.fused   = FALSE,
                      gamma            = 1,
                      standardize      = FALSE,
                      intercept        = TRUE,
                      one.intercept    = FALSE,
                      compute.se       = FALSE,
                      conf.int         = NULL,
                      rho              = NULL,
                      dynamic.rho      = TRUE,
                      maxit            = 500L,
                      abs.tol          = 1e-5,
                      rel.tol          = 1e-5,
                      irls.tol         = 1e-5,
                      irls.maxit       = 100L,
                      model.matrix     = FALSE,
                      ...)
{


    family <- match.arg(family)
    this.call = match.call()


    nvars <- ncol(x)
    nobs <- nrow(x)
    y <- drop(y)
    dimy <- dim(y)
    leny <- ifelse(is.null(dimy), length(y), dimy[1])
    stopifnot(leny == nobs)

    stopifnot(nobs == nrow(groups))
    
    if (!is.null(lambda.fused)) stop("fused not available yet")

    ## apply surv
    if (family == "coxph") 
    {
        ## taken from glmnet
        if(!is.matrix(y)||!all(match(c("time","status"),dimnames(y)[[2]],0))) stop("Cox model requires a matrix with columns 'time' (>0) and 'status'  (binary) as a response; a 'Surv' object suffices",call.=FALSE)
        y.surv <- y
        delta = as.double(y[,"status"])
        if(sum(delta) == 0)stop("All data are censored; not permitted for Cox family")
        y = as.double(y[,"time"])
        if(any(y <= 0))stop("negative event times encountered;  not permitted for Cox family")
        maxit=as.integer(maxit)
    }
    else
    {
        delta <- rep(0L, nrow(x))
    }

    if (!is.null(conf.int))
    {
        conf.int <- conf.int[1]
        if (conf.int <= 0 | conf.int >= 1)
        {
            stop("conf.int must be within (0,1)")
        }
        # need to compute SEs for confidence intervals
        compute.se <- TRUE
    }

    ## if coxph, order the data and check the
    ## first one cannot be censored
    if(family == "coxph")
    {
        yorder <- order(y)
        y <- y[yorder]
        delta <- delta[yorder]
        x <- x[yorder,]
        groups <- groups[yorder,]
        if (delta[1] ==0)
        {
            x <- x[-(1:(which.max(delta)-1)),]
            y <- y[-(1:(which.max(delta)-1))]
            groups <- groups[-(1:(which.max(delta)-1)),]
            delta <- delta[-(1:(which.max(delta)-1))]

        }
        intercept <- FALSE
    }

    ## we only care about one.intercept if intercept = TRUE
    one.intercept <- !(!one.intercept & intercept)

    storage.mode(groups) <- "numeric"

    ## starting and ending index of each dataset
    n.conditions <- ncol(groups)
    M            <- 2 ^ n.conditions
    p            <- ncol(x)              # number of covariates
    N            <- nrow(x)              # total number of observations
    vnames       <- colnames(x)
    if(is.null(vnames)) vnames <- paste0("V", seq(p))

    if (!one.intercept) vnames <- c("(Intercept)", vnames)


    if (family == "binomial" ) 
    {
        nc=dim(y)
        maxit=as.integer(maxit)
        if(is.null(nc))
        {
            ## Need to construct a y matrix, and include the weights
            yy=as.factor(y)
            ntab=table(yy)
            classnames=names(ntab)
            nc=as.integer(length(ntab))
            #y=diag(nc)[as.numeric(y),]
            rm(yy)
        } else 
        {
            noo=nc[1]
            if(noo!=N)stop("x and y have different number of rows in call to glmnet",call.=FALSE)
            nc=as.integer(nc[2])
            classnames=colnames(y)
        }
    } else 
    {
        classnames <- NULL
    }




    ######################################################################
    ##
    ##               the following code sets up the data
    ##                according to the specified strata
    ##
    ######################################################################

    # all possible combinations of conditions
    combin.mat <- array(NA, dim = c(M, n.conditions))
    data.indices <- vector(mode = "list", length = M)
    group.vec <- apply(groups, 1, FUN = function(x) {paste(x, collapse = ",")})
    ind.2.remove <- NULL
    for (c in 1:M) {
        combin.mat[c,] <- as.integer(intToBits(c)[1:n.conditions])
        data.indices[[c]] <- which(group.vec == paste(combin.mat[c,], collapse = ","))
        if (length(data.indices[[c]]) == 0) {
            ind.2.remove <- c(c, ind.2.remove)
        }
    }
    # remove combinations that have no observations
    if (!is.null(ind.2.remove)) {
        data.indices[ind.2.remove] <- NULL
        M <- length(data.indices)
        combin.mat <- combin.mat[-ind.2.remove,]
    }
    # sort in order that allows fitting to work
    #ord.idx <- order(rowSums(combin.mat), decreasing = FALSE)
    #order combin.mat by each column
    ord.idx <- do.call(order, split(combin.mat, rep(1:ncol(combin.mat), each = nrow(combin.mat))))
    combin.mat <- combin.mat[ord.idx,]
    data.indices <- data.indices[ord.idx]

    direct.above.idx.list <- above.idx.list <- vector(mode = "list", length = M)


    for (c in 1:M) 
    {
        indi <- which(combin.mat[c,] == 1)

        # get all indices of g terms which are directly above the current var
        # in the hierarchy, ie. for three categories A, B, C, for the A group,
        # this will return the indices for AB and AC. For AB, it will return
        # the index for ABC. for none, it will return the indices for
        # A, B, and C, etc
        inner.loop <- (1:(M))[-c]
        for (j in inner.loop) 
        {
            diffs.tmp <- combin.mat[j,] - combin.mat[c,]
            if (all( diffs.tmp >= 0 )) 
            {
                if (sum( diffs.tmp == 1 ) == 1) 
                {
                    direct.above.idx.list[[c]] <- c(direct.above.idx.list[[c]], j)
                }
                above.idx.list[[c]] <- c(above.idx.list[[c]], j)
            }
        }
        above.idx.list[[c]] <- c(above.idx.list[[c]], c)
    }
    rsc <- rowSums(combin.mat)
    if (any(rsc == 0)) 
    {
        above.idx.list[[which(rsc == 0)]] <- which(rsc == 0)
    }


    ## this is ugly, but it works.
    ## add an intercept term. this will effectively
    ## ensure that there is a separate intercept for
    ## each subpopulation
    if (!one.intercept & (family != "coxph")) x <- cbind(1, x)

    ## find variables that are missing in an entire condition
    which.missing.idx   <- checkWhichVarsMissingForWhichConditions(x, groups)

    ## find which variables are missing in all conditions but one
    which.missing.idx.2 <- checkWhichVarsMissingForAllButOneCondition(x, groups)

    ## find variables with no variation within a strata
    which.zero.sd.idx   <- checkWhichVarsZeroSDForWhichConditions(x, groups, combin.mat)
    if (!one.intercept & (family != "coxph")) which.zero.sd.idx[1,] <- FALSE

    # matrix with dimension (num. vars)x(num. group combinations) with a value
    # 1 in position (i,j) if variable i is not missing in condition combination (strata) j
    which.not.missing.idx.combin <- 1 - ( (t(tcrossprod(combin.mat[rev(1:nrow(combin.mat)),], which.missing.idx)) == 1) |
                                              (t(tcrossprod(combin.mat[rev(1:nrow(combin.mat)),], which.missing.idx.2)) == 1) |
                                              which.zero.sd.idx )

    which.not.missing.idx.combin.vec <- as.vector(which.not.missing.idx.combin)


    ## number of variables gets larger if we have many intercepts
    if (!one.intercept & (family != "coxph")) p <- p + 1

    if (is.null(penalty.factor)) penalty.factor <- rep(1, p * M )

    ## unpenalize intercepts
    if (!one.intercept & (family != "coxph"))
    {
        ## indices of intercepts
        intercept.idx <- seq.int(1, (M - 1) * p + 1, by = p)
        intercept <- FALSE
    } else
    {
        intercept.idx <- numeric(0)
    }

    grp.membership.block <- t(sapply(above.idx.list, function(x){tmp <- numeric(M); tmp[x] <- 1; tmp}))
    #group.mat <- array(0, dim = rep(p * M, 2))
    group.mat <- Matrix(0, nrow = (p * M), ncol = (p * M))
    #for (i in 1:p) group.mat[((i-1)*M+1):(i*M),((1:M)*p - p + i )] <- t(t(grp.membership.block) * which.not.missing.idx.combin[i,])
    for (i in 1:p) group.mat[((1:M) * p - p + i ), ((i - 1) * M + 1):(i * M)] <- t(t(grp.membership.block) * which.not.missing.idx.combin[i,])

    ## unpenalize intercepts
    if (!one.intercept & (family != "coxph"))
    {
        intercept.grp.idx <- which(Matrix::colSums(group.mat[intercept.idx,]) > 0)
        penalty.factor[intercept.grp.idx] <- 0
    }

    cs.nz.idx <- which(Matrix::colSums(group.mat) > 0)
    group.mat <- group.mat[which.not.missing.idx.combin.vec == 1, cs.nz.idx]


    penalty.factor <- penalty.factor[which.not.missing.idx.combin.vec == 1]

    above.not.incl.idx.list <- lapply(1:length(above.idx.list), function(i) above.idx.list[[i]][above.idx.list[[i]] != i])

    fused.D.block <- array(0, dim = c(length(unlist(above.not.incl.idx.list)), nrow(grp.membership.block)))

    k <- 0
    for (i in 1:length(above.not.incl.idx.list))
    {
        ab.idx <- above.not.incl.idx.list[[i]]
        if (length(ab.idx) > 0)
        {
            for (j in 1:length(ab.idx))
            {
                k <- k + 1
                fused.D.block[k, c(i, ab.idx[j])] <- c(1, -1)
            }
        }
    }
    if (!is.null(lambda.fused))
    {
        N.fused   <- nrow(fused.D.block)
        fused.mat <- Matrix(0, nrow = (p * N.fused), ncol = (p * M))

        for (i in 1:p) fused.mat[((i - 1) * N.fused + 1):(i * N.fused), ((1:M) * p - p + i )] <- fused.D.block * which.not.missing.idx.combin[i,]

        rs.nz.idx <- which(Matrix::rowSums(abs(fused.mat) ) > 0)
        fused.mat <- fused.mat[rs.nz.idx, which.not.missing.idx.combin.vec == 1]
    } else
    {
        fused.mat <- NULL
        lambda.fused <- 0
    }


    ## create expanded, block diagonal matrix
    ## with each strata's submatrix as a block
    big.x <- as.matrix(array(0, dim = c(N, ncol(group.mat))))
    col.idx.vec <- c(0,cumsum(colSums(which.not.missing.idx.combin)))
    for (c in 1:M) {
        cur.var.idx <- which(which.not.missing.idx.combin[,c]==1)
        big.x[data.indices[[c]], (col.idx.vec[c]+1):col.idx.vec[c+1]] <- x[data.indices[[c]], cur.var.idx]
    }

    ######################################################################
    ##
    ##                          compute model
    ##
    ######################################################################



    fit <- vennLasso::oglasso(x                    = big.x,
                              y                    = y,
                              delta                = delta,
                              group                = group.mat,
                              fused                = fused.mat,
                              family               = family,
                              nlambda              = nlambda,
                              lambda               = lambda,
                              lambda.min.ratio     = lambda.min.ratio,
                              lambda.fused         = lambda.fused,
                              penalty.factor       = penalty.factor,
                              penalty.factor.fused = penalty.factor,
                              group.weights        = group.weights,
                              adaptive.lasso       = adaptive.lasso,
                              adaptive.fused       = adaptive.fused,
                              gamma                = gamma,
                              standardize          = standardize,
                              intercept            = intercept,
                              compute.se           = compute.se,
                              dynamic.rho          = dynamic.rho,
                              maxit                = maxit,
                              abs.tol              = abs.tol,
                              rel.tol              = rel.tol,
                              irls.tol             = irls.tol,
                              irls.maxit           = irls.maxit,
                              ...)

    nlambda <- fit$nlambda

    class2 <- switch(family,
                     "gaussian" = "venngaussian",
                     "binomial" = "vennbinomial", "coxph" = "venncoxph")

    beta.tmp   <- fit$beta

    ######################################################################
    ##
    ##                    compute baseline harzard
    ##                           if needed
    ##
    ######################################################################

    if (family == "coxph")
    {
        yorder <- order(y)
        yy     <- y[yorder]
        deltaorder <- delta[yorder]
        w          <- array(0,c(N,nlambda))
        for (l in 1:nlambda)
        {
           w[,l] <- exp(big.x[yorder,] %*% fit$beta[-1,l])
        }
    }

    ######################################################################
    ##
    ##                    compute standard errors
    ##                           if needed
    ##
    ######################################################################

    n.free.param <- rep(NA, length(fit$lambda))

    xtx <- fit$XtX
    xty <- fit$XtY

    ## for logistic an extra grp is added at the end with
    ## zero weight for intercept
    if (length(fit$group.weights) > ncol(big.x))
    {
        if (fit$group.weights[length(fit$group.weights)] == 0)
        {
            fit$group.weights <- fit$group.weights[-length(fit$group.weights)]
        }
    }

    if (compute.se)
    {
        nz.idx <- (abs(beta.tmp[-1,]) >= sqrt(.Machine$double.eps))
        se <- beta.refit <- array(0, dim = dim(beta.tmp[-1,]))
        intercept.refit <- numeric(length(beta.tmp[1,]))
        if (!is.null(conf.int))
        {
            z.val <- qnorm((1 + conf.int)/2, 0, 1)
            lower.ci <- upper.ci <- lower.ci.refit <- upper.ci.refit <- array(0, dim = dim(beta.tmp[-1,]))
        }

        do.not.compute.se <- FALSE
        for (l in 1:nlambda)
        {
            fit.it <- FALSE
            if (l == 1)
            {
                if (sum(nz.idx[,l]) > 0)
                {
                    fit.it <- TRUE
                }
            } else
            {
                # refit again if any nonzero are different from last lambda
                if (any(nz.idx[,l] != nz.idx[,l - 1]))
                {
                    fit.it <- TRUE
                }
            }

            if (sum(nz.idx[,l]) >= N)
            {
                ## do not re-fit model if it is saturated
                fit.it <- FALSE
                do.not.compute.se <- TRUE
            }

            if (fit.it)
            {

                if (family == "gaussian")
                {

                    if (intercept)
                    {
                        #xtx <- crossprod(cbind(1, big.x))

                        idx.nz <- (c(TRUE, nz.idx[,l]))
                    } else
                    {
                        #xtx <- crosprod(big.x)

                        idx.nz <- (nz.idx[,l])
                    }

                    compute.mod.cov <- sum(idx.nz) < 0.95 * ncol(xtx)

                    if(compute.mod.cov)
                    {
                        G11 <- as(xtx[ idx.nz,  idx.nz, drop = FALSE], "dpoMatrix")
                        G12 <- xtx[ idx.nz, !idx.nz, drop = FALSE]
                        G21 <- xtx[!idx.nz,  idx.nz, drop = FALSE]
                        G22 <- as(xtx[!idx.nz, !idx.nz, drop = FALSE], "dpoMatrix")
                    } else
                    {
                        G11 <- as(xtx[ idx.nz,  idx.nz, drop = FALSE], "dpoMatrix")
                    }


                    if (intercept)
                    {

                        #beta.refit.tmp <- drop(solve(G11, crossprod(cbind(1, big.x[,nz.idx[,l]]), y)))
                        beta.refit.tmp <- drop(solve(G11, xty[c(TRUE, nz.idx[,l])]))

                        if (sum(nz.idx[,l]) == 1)
                        {
                            resids <- y - drop(big.x[,nz.idx[,l]] * beta.refit.tmp[-1]) - beta.refit.tmp[1]
                        } else
                        {
                            resids <- y - big.x[,nz.idx[,l]] %*% beta.refit.tmp[-1] - beta.refit.tmp[1]
                        }

                    } else
                    {

                        #beta.refit.tmp <- drop(solve(G11, crossprod(big.x[,nz.idx[,l]], y)))
                        beta.refit.tmp <- drop(solve(G11, xty[nz.idx[,l]]))

                        resids <- y - big.x[,nz.idx[,l]] %*% drop(beta.refit.tmp)
                    }

                    ## sigma ^ 2
                    sigma.sq <- sum(resids ^ 2) / (length(resids) - ncol(G11))

                    ## place the coefficient results
                    ## in an array instead of a long vector
                    beta.cur.tmp <- NULL
                    if (intercept)
                    {
                        beta.cur.tmp <- beta.tmp[1, l]
                        #ada.wts.denom <- as.vector(Matrix::rowSums( Matrix::t(Matrix::t(group.mat) * fit$group.weights) ))
                        grp.norms <- as.vector(sqrt( group.mat %*% beta.tmp[-1, l] ^ 2 ))

                        grp.mult <- 1 / (grp.norms * fit$group.weights)
                        grp.mult[is.infinite(grp.mult)] <- 1e15
                        grp.mult[is.nan(grp.mult)] <- 1e15

                        ada.wts.denom <- as.vector(Matrix::rowSums( Matrix::t(Matrix::t(group.mat) * (grp.mult)  ) ))
                        #denom <- ada.wts.denom * other.wts.denom
                        #denom[denom == 0] <- 1e-10
                        ada.wts.denom[ada.wts.denom == 0] <- 1e15
                        A.diag <- c(0, ada.wts.denom)
                    } else
                    {
                        #ada.wts.denom <- as.vector(Matrix::rowSums( Matrix::t(Matrix::t(group.mat) * fit$group.weights) ))
                        grp.norms <- as.vector(sqrt( group.mat %*% beta.tmp[-1, l] ^ 2 ))

                        grp.mult <- 1 / (grp.norms * fit$group.weights)
                        grp.mult[is.infinite(grp.mult)] <- 1e15
                        grp.mult[is.nan(grp.mult)] <- 1e15

                        ada.wts.denom <- as.vector(Matrix::rowSums( Matrix::t(Matrix::t(group.mat) * (grp.mult)  ) ))
                        #denom <- ada.wts.denom * other.wts.denom
                        #denom[denom == 0] <- 1e-10
                        ada.wts.denom[ada.wts.denom == 0] <- 1e15
                        A.diag <- ada.wts.denom
                    }

                    if (intercept)
                    {
                        A.diag.nz <- A.diag[c(TRUE, nz.idx[,l])]
                    } else
                    {
                        A.diag.nz <- A.diag[nz.idx[,l]]
                    }

                    beta.cur.tmp <- c(beta.cur.tmp, beta.tmp[-1,l][nz.idx[,l]])

                    G11.inv <- solve(G11)

                    if(compute.mod.cov)
                    {
                        E <- G22 - G21 %*% solve(G11, G12)

                        G11.tilde <- G11
                        Matrix::diag(G11.tilde) <- Matrix::diag(G11.tilde) + (nobs * fit$lambda[l] * A.diag.nz)
                        G11.tilde.inv <- solve(G11.tilde)

                        cov.mat <- G11.inv + (G11.inv - G11.tilde.inv) %*% G12 %*% solve(E, G21 %*% (G11.inv - G11.tilde.inv))
                    } else
                    {
                        cov.mat <- G11.inv
                    }

                    se[nz.idx[,l], l] <- sqrt(sigma.sq * diag(cov.mat))

                    #n.free.param[l] <- sum(diag( solve(xtx + nobs * fit$lambda[l] * diag(A.diag), xtx) ))

                    n.free.param[l] <- ncol(xtx) - sum(nobs * fit$lambda[l] * diag(A.diag) / (N + nobs * fit$lambda[l] * diag(A.diag)))

                    #} else
                    #{
                    #    df.sub <- data.frame(y = y, big.x[,nz.idx[,l] ])
                    #    re.fit <- lm(y ~ ., data = df.sub)

                    #    # take out intercept
                    #    se[nz.idx[,l], l] <- sqrt(diag(vcov(re.fit)))[-1]
                    #    coef.tmp <- coef(re.fit)
                    #    beta.refit[nz.idx[,l], l] <- coef.tmp[-1]
                    #    intercept.refit[l] <- coef.tmp[1]
                    #}

                } else if (family == "binomial")
                {
                    df.sub <- data.frame(y = y, big.x[,nz.idx[,l] ])
                    re.fit <- glm(y ~ ., data = df.sub, family = binomial())

                    # take out intercept
                    se[nz.idx[,l], l] <- sqrt(diag(vcov(re.fit)))[-1]
                    coef.tmp <- coef(re.fit)
                    beta.refit[nz.idx[,l], l] <- coef.tmp[-1]
                    intercept.refit[l] <- coef.tmp[1]
                } else if (family == "coxph")
                {
                    df.sub <- data.frame(y = y.surv, big.x[,nz.idx[,l] ])
                    re.fit <- coxph(y ~ ., data = df.sub)

                    # no need to take out intercept
                    if (anyNA(re.fit$coefficients))
                    {
                        se[nz.idx[,l-1], l] <- se[nz.idx[,l-1], l-1]
                        beta.refit[nz.idx[,l-1], l] <- beta.refit[nz.idx[,l-1], l-1]
                        print("Results for a larger lambda is returned")
                    } else {
                        se[nz.idx[,l], l] <- sqrt(diag(vcov(re.fit)))
                        beta.refit[nz.idx[,l], l] <- coef(re.fit)
                    }

 #                   w.refit[,l] <- exp(big.x[yorder,] %*% coef(re.fit))
                }
            } else
            {
                ## if the selected vars are the same as
                ## for the last lambda, then just take
                ## the last computed standard errors
                if (l > 1)
                {
                    se[nz.idx[,l-1], l] <- se[nz.idx[,l-1], l-1]
                    n.free.param[l]   <- n.free.param[l-1]
                } else
                {
                    n.free.param[1] <- 1 * intercept
                }
            }
            if (do.not.compute.se)
            {
                se[, l] <- NA
            }
            if (!is.null(conf.int))
            {
                lower.ci[nz.idx[,l], l] <- beta.tmp[-1, ][nz.idx[,l], l] - z.val * se[nz.idx[,l], l]
                upper.ci[nz.idx[,l], l] <- beta.tmp[-1, ][nz.idx[,l], l] + z.val * se[nz.idx[,l], l]
                lower.ci.refit[nz.idx[,l], l] <- beta.refit[nz.idx[,l], l] - z.val * se[nz.idx[,l], l]
                upper.ci.refit[nz.idx[,l], l] <- beta.refit[nz.idx[,l], l] + z.val * se[nz.idx[,l], l]
            }
        }
    } else ## if not compute.se
    {

        for (l in 1:nlambda)
        {
            if (intercept)
            {
                beta.cur.tmp <- beta.tmp[1, l]

                grp.norms <- as.vector(sqrt( group.mat %*% beta.tmp[-1, l] ^ 2 ))

                grp.mult <- 1 / (grp.norms * fit$group.weights)
                grp.mult[is.infinite(grp.mult)] <- 1e15
                grp.mult[is.nan(grp.mult)] <- 1e15

                ada.wts.denom <- as.vector(Matrix::rowSums( Matrix::t(Matrix::t(group.mat) * (grp.mult)  ) ))
                ada.wts.denom[ada.wts.denom == 0] <- 1e15
                A.diag <- c(0, ada.wts.denom)
            } else
            {
                grp.norms <- as.vector(sqrt( group.mat %*% beta.tmp[-1, l] ^ 2 ))

                grp.mult <- 1 / (grp.norms * fit$group.weights)
                grp.mult[is.infinite(grp.mult)] <- 1e15
                grp.mult[is.nan(grp.mult)] <- 1e15

                ada.wts.denom <- as.vector(Matrix::rowSums( Matrix::t(Matrix::t(group.mat) * (grp.mult)  ) ))
                ada.wts.denom[ada.wts.denom == 0] <- 1e15
                A.diag <- ada.wts.denom
            }

            ## estimate number of free params
            n.free.param[l] <- ncol(xtx) - sum(nobs * fit$lambda[l] * diag(A.diag) / (N + nobs * fit$lambda[l] * diag(A.diag)))

        }

    }

    ######################################################################
    ##
    ##                    place betas in a matrix
    ##                   dim = (n groups) x (n vars)
    ##
    ######################################################################
    ##                place the SE and conf.int results
    ##              in an array instead of a long vector
    ##
    beta.array <- array(NA, dim = c(M, p, nlambda))
    if (compute.se)
    {
        se.array <- beta.refit.array <- beta.array
        for (c in 1:M) {
            cur.var.idx <- which(which.not.missing.idx.combin[,c]==1)
            for (l in 1:nlambda) {
                se.array[c, cur.var.idx, l] <- se[(col.idx.vec[c]+1):col.idx.vec[c+1],l]
                beta.refit.array[c, cur.var.idx, l] <- beta.refit[(col.idx.vec[c]+1):col.idx.vec[c+1],l]
            }
        }
        if (!is.null(conf.int))
        {
            lower.array <- upper.array <- lower.refit.array <- upper.refit.array <- se.array
            for (c in 1:M) {
                cur.var.idx <- which(which.not.missing.idx.combin[,c]==1)
                for (l in 1:nlambda) {
                    lower.array[c, cur.var.idx, l] <- lower.ci[(col.idx.vec[c]+1):col.idx.vec[c+1],l]
                    upper.array[c, cur.var.idx, l] <- upper.ci[(col.idx.vec[c]+1):col.idx.vec[c+1],l]
                    lower.refit.array[c, cur.var.idx, l] <- lower.ci.refit[(col.idx.vec[c]+1):col.idx.vec[c+1],l]
                    upper.refit.array[c, cur.var.idx, l] <- upper.ci.refit[(col.idx.vec[c]+1):col.idx.vec[c+1],l]
                }
            }
        } else
        {
            lower.array <- upper.array <- lower.refit.array <- upper.refit.array <- NA
        }
    } else
    {
        lower.array <- upper.array <- lower.refit.array <- upper.refit.array <- NA
        se.array <- NA
        beta.refit.array <- NA
        intercept.refit <- NA
    }

    ## place the coefficient results
    ## in an array instead of a long vector
    for (c in 1:M) 
    {
        cur.var.idx <- which(which.not.missing.idx.combin[,c]==1)
        for (l in 1:nlambda) {
            beta.array[c, cur.var.idx, l] <- beta.tmp[1+(col.idx.vec[c]+1):col.idx.vec[c+1],l]
        }
    }
    combin.names <- apply(combin.mat, 1, function(x) paste(x, collapse = ","))
    dimnames(beta.array) <- list(combin.names, vnames, paste(seq(along=fit$lambda)))


    if (compute.se)
    {
        dimnames(se.array) <- dimnames(beta.refit.array) <- dimnames(beta.array)
        if (!is.null(conf.int))
        {
            dimnames(lower.array) <- dimnames(upper.array) <- dimnames(beta.array)
            dimnames(lower.refit.array) <- dimnames(upper.refit.array) <- dimnames(beta.array)
        }
    }

    ######################################################################
    ##
    ##                  arrange results to be returned
    ##
    ######################################################################

    fit$beta                   <- beta.array
    fit$intercept              <- beta.tmp[1,]
    fit$se                     <- se.array
    fit$lower.ci               <- lower.array
    fit$upper.ci               <- upper.array
    fit$beta.refit             <- beta.refit.array
    fit$lower.ci.refit         <- lower.refit.array
    fit$upper.ci.refit         <- upper.refit.array
    fit$intercept.refit        <- intercept.refit
    fit$n.conditions           <- n.conditions
    fit$n.combinations         <- M
    fit$nobs                   <- N
    fit$nvars                  <- p
    fit$intercept.fit          <- intercept
    fit$standardized           <- standardize
    fit$df.strata              <- predict.vennLasso(fit, type="nvars")

    free.params                <- colSums(fit$df.strata) + 1 * intercept
    free.params                <- n.free.param
    fit$df                     <- N - free.params

    fit$logLik                 <- logLik.vennLasso(fit)
    fit$aic                    <-      2 * free.params - 2 * fit$logLik
    fit$bic                    <- log(N) * free.params - 2 * fit$logLik
    c <- 4
    fit$mbic                   <- fit$bic + 2 * free.params * log( (ncol(big.x) + 1 * intercept)/c - 1 )
    fit$aicc                   <- fit$aic + (2 * free.params) * (free.params + 1) / (N - free.params - 1)
    fit$condition.combinations <- combin.mat
    fit$combin.names           <- combin.names
    fit$var.names              <- vnames
    fit$one.intercept          <- one.intercept
    fit$call                   <- this.call
    if (model.matrix)
    {
        fit$model.matrix <- big.x
    }
    if (family == "coxph") {
        fit$status <- deltaorder
        fit$W <- w
        fit$time <- yy
    }
    beta.idx <- match("beta", names(fit))
    int.idx <- match("intercept", names(fit))
    se.idx <- match("se", names(fit))
    fit <- fit[c(beta.idx, int.idx, se.idx, (1:length(fit))[-c(beta.idx, int.idx, se.idx)])]
    fit$data.indices <- data.indices
    class(fit) <- c("vennLasso", class2)
    fit
}





