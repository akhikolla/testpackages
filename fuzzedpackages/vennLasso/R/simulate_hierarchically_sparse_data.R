#' function to generate data with hierarchical sparsity
#'
#' @param ncats number of categories to stratify on
#' @param nvars number of variables
#' @param nobs number of observations per strata to simulate
#' @param nobs.test number of independent test observations per strata to simulate
#' @param hier.sparsity.param parameter between 0 and 1 which determines how much hierarchical sparsity there is. To achieve
#' a desired total level of sparsity among the variables with hierarchical sparsity, this parameter can be estimated using the
#' function 'estimate.hier.sparsity.param'
#' @param avg.hier.zeros desired percent of zero variables among the variables with hierarchical zero patterns. If this is specified, it will
#' override the given hier.sparsity.param value and estimate it. This takes a while
#' @param prop.zero.vars proportion of all variables that will be zero across all strata
#' @param effect.size.max maximum magnitude of the true effect sizes
#' @param misspecification.prop proportion of variables with hierarchical missingness misspecified
#' @param family family for the response variable
#' @param sd standard devation for gaussian simulations
#' @param snr signal-to-noise ratio (only used for \code{family = "gaussian"})
#' @param beta a matrix of true beta values. If given, then no beta will be created and data will be simulated from the given beta
#' @param tau rate parameter for \code{rexp()} for generating time-to-event outcomes
#' @param covar scalar, pairwise covariance term for covariates
#' @importFrom stats rnorm
#' @importFrom stats rexp
#' @importFrom stats var
#' @importFrom stats approx
#' @importFrom stats runif
#' @import MASS
#' @export
#' @examples
#' set.seed(123)
#'
#' dat.sim <- genHierSparseData(ncats = 3, nvars = 100, nobs = 200)
#'
#' # estimate hier.sparsity.param for 0.15 total proportion of nonzero variables
#' # among vars with hierarchical zero patterns
#' \dontrun{
#' hsp <- estimate.hier.sparsity.param(ncats = 3, nvars = 50, avg.hier.zeros = 0.15, nsims = 100)
#' }
#' # the above results in the following value
#' hsp <- 0.6270698
#'
#' # check that this does indeed achieve the desired level of sparsity
#' mean(replicate(50, mean(genHierSparseBeta(ncats = 3, 
#'                        nvars = 50, hier.sparsity.param = hsp) != 0)  ))
#'                            
#' dat.sim2 <- genHierSparseData(ncats = 3, nvars = 100, nobs = 200, hier.sparsity.param = hsp)
#'
#' sparseBeta <- genHierSparseBeta(ncats = 3, nvars = 100, hier.sparsity.param = hsp)
#'
#' ## generate data with already generated beta
#' dat.sim3 <- genHierSparseData(ncats = 3, nvars = 100, nobs = 200, beta = sparseBeta)
#'
#'
#'
#'
#' ## complete example:
#' ## 50% sparsity:
#' hsp <- 0.2626451
#'
#' dat.sim <- genHierSparseData(ncats = 3, nvars = 25,
#'                              nobs = 150, nobs.test = 1000,
#'                              hier.sparsity.param = hsp,
#'                              prop.zero.vars = 0.5,
#'                              effect.size.max = 0.25,
#'                              family = "gaussian")
#'
#' x        <- dat.sim$x
#' x.test   <- dat.sim$x.test
#' y        <- dat.sim$y
#' y.test   <- dat.sim$y.test
#' grp      <- dat.sim$group.ind
#' grp.test <- dat.sim$group.ind.test
#'
#' fit.adapt <- cv.vennLasso(x, y,
#'                           grp,
#'                           adaptive.lasso = TRUE,
#'                           nlambda        = 25,
#'                           family         = "gaussian",
#'                           abs.tol        = 1e-5,
#'                           rel.tol        = 1e-5,
#'                           maxit          = 1000,
#'                           irls.maxit     = 15L,
#'                           gamma          = 0.2,
#'                           standardize    = FALSE,
#'                           intercept      = TRUE,
#'                           nfolds         = 3,
#'                           model.matrix   = TRUE)
#'
#' preds.a <- predict(fit.adapt$vennLasso.fit, x.test, grp.test, s = fit.adapt$lambda.min,
#'                    type = 'response')
#'
#'
genHierSparseData <- function(ncats,
                              nvars,
                              nobs,
                              nobs.test = 100,
                              hier.sparsity.param = 0.5,
                              avg.hier.zeros = NULL,
                              prop.zero.vars = 0.5,
                              effect.size.max = 0.5,
                              misspecification.prop = 0,
                              family = c("gaussian", "binomial", "coxph"),
                              sd  = 1,
                              snr = NULL,
                              beta = NULL,
                              tau = 10,
                              covar = 0
                              )
{

    family <- match.arg(family)
    x.list <- km.list <- y.list <- vector(length = 2 ^ ncats, mode = "list")

    x.test.list <- y.test.list <- y.test.pred.list  <- z.test.list <- var.idx.list <- x.list

    stopifnot(covar[1] >= 0)

    ncmb <- 2 ^ ncats
    for (i in 1:ncmb) 
    {
        km.list[[i]] <- as.integer(intToBits(i)[1:ncats])
    }
    km.mat <- Reduce(rbind, km.list)
    attr(km.mat, "dimnames") <- NULL
    ord.idx <- do.call(order, split(km.mat, rep(1:ncol(km.mat), each = nrow(km.mat))))
    km.mat <- km.mat[ord.idx,]
    km.list <- km.list[ord.idx]

    n.zero.cols <- floor(prop.zero.vars * nvars)
    n.nzvars    <- nvars - n.zero.cols
    zero.cols   <- (n.nzvars + 1):nvars
    
    if (is.null(beta)) 
    {
        rs <- rowSums(km.mat)
        id1 <- which(rs == 1)
        id2 <- which(rs > 1)
        id0 <- which(rs == 0)

        beta.mat <- genHierSparseBeta(ncats, n.nzvars, hier.sparsity.param,
                                      avg.hier.zeros, effect.size.max, misspecification.prop)
        beta.mat.sparse <- cbind(beta.mat, array(0, dim = c(nrow(beta.mat), n.zero.cols)))
        rownames(beta.mat.sparse) <- rownames(beta.mat)
    } else
    {
        beta.mat.sparse <- beta
    }

    xb.all <- NULL
    if (family == "coxph")
    {
        t.list       <- NULL
        ct.list      <- NULL
        c.list       <- NULL
        t.test.list  <- NULL
        ct.test.list <- NULL
        c.test.list  <- NULL
    }

    for (i in 1:(2 ^ ncats)) 
    {
        if (covar[1] == 0)
        {
            x.list[[i]] <- matrix(rnorm(nobs * nvars), ncol = nvars)
            x.test.list[[i]] <- matrix(rnorm(nobs.test * nvars), ncol = nvars)
        } else
        {
            sigma <- covar[1] ^ abs(outer(1:nvars, 1:nvars, FUN = "-"))

            x.list[[i]] <- mvrnorm(n = nobs, mu = rep(0, nvars), Sigma = sigma)
            x.test.list[[i]] <- mvrnorm(n = nobs.test, mu = rep(0, nvars), Sigma = sigma)
        }

        if (family == "binomial")
        {
            y.list[[i]] <- rbinom( nrow(x.list[[i]]), 1, 1 / (1 + exp(x.list[[i]] %*% (beta.mat.sparse[i,])  )) )
            y.test.list[[i]] <- rbinom(nrow(x.test.list[[i]]), 1, 1 / (1 + exp(x.test.list[[i]] %*% (beta.mat.sparse[i,]) )) )
        } else if (family == "gaussian")
        {
            xb.tmp <- x.list[[i]] %*% beta.mat.sparse[i,]
            if (!is.null(snr))
            {
                xb.all <- c(xb.all, xb.tmp)
            }
            y.list[[i]] <- xb.tmp # + rnorm(nobs, sd = sd)
        } else if (family == "coxph")
        {

            t.list[[i]] <- rexp(nrow(x.list[[i]]), exp(x.list[[i]] %*% (beta.mat.sparse[i,])  ))
            ct.list[[i]] <- rexp(nrow(x.list[[i]]), tau)
            c.list[[i]] <- (t.list[[i]]>ct.list[[i]])
            y.list[[i]] <- c.list[[i]] * t.list[[i]] + (1-c.list[[i]]) * ct.list[[i]]

            t.test.list[[i]] <- rexp(nrow(x.test.list[[i]]), exp(x.test.list[[i]] %*% (beta.mat.sparse[i,])  ))
            ct.test.list[[i]] <- rexp(nrow(x.test.list[[i]]), tau)
            c.test.list[[i]] <- (t.test.list[[i]]>ct.test.list[[i]])
            y.test.list[[i]] <- c.test.list[[i]] * t.test.list[[i]] + (1-c.test.list[[i]]) * ct.test.list[[i]]
            z.test.list[[i]] <- x.test.list[[i]] %*% (beta.mat.sparse[i,])
        }
    }
    if (family == "gaussian")
    {
        if (!is.null(snr))
        {
            var.xb <- var(xb.all)
            sd <- sqrt(var.xb / snr)
        }

        for (i in 1:(2 ^ ncats))
        {
            y.list[[i]] <- y.list[[i]] + rnorm(nobs, sd = sd)
            y.test.list[[i]] <- x.test.list[[i]] %*% beta.mat.sparse[i,] + rnorm(nobs.test, sd = sd)
        }
    }




    for (c in 1:(2 ^ ncats)) 
    {
        indi <- which(km.mat[c,] == 1)
        # this returns indices of all g terms to be included in the
        # product for the beta coefficient
        var.idx.list[[c]] <- which(apply(km.mat, 1, function(x) all(x[indi] == 1)))
    }

    x.all <- Reduce(rbind, x.list)
    y.all <- Reduce(c, y.list)
    if (family == "coxph")
    {
        c.all <- Reduce(c,c.list)
    }
    group.ind <- Reduce(rbind, lapply(km.list, function(x) matrix(rep(x, nobs), ncol = length(x), byrow = TRUE)) )

    #df <- data.frame(y = y.all, x.all, group.ind)


    x.test.all <- Reduce(rbind, x.test.list)
    group.test.ind <- Reduce(rbind, lapply(km.list, function(x) matrix(rep(x, nobs.test), ncol = length(x), byrow = TRUE)) )
    y.test.all <- Reduce(c, y.test.list)
    if (family == "coxph")
    {
        c.test.all <- Reduce(c,c.test.list)
        z.test.all <- Reduce(c,z.test.list)
    }
    #df.test <- data.frame(y = y.test.all, x.test.all, group.test.ind)


    #fm <- as.formula(paste("y ~ (", paste(paste("X", 1:nvars, sep = ""), collapse = "+"), ")*(",
    #                       paste(paste("X", 1:ncats, ".1", sep = ""), collapse = "+"), ") - 1", sep = "" ))
    #x.mm <- model.matrix(fm, data = df)
    #x.test.mm <- model.matrix(fm, data = df.test)




    rownames(beta.mat.sparse) <- apply(km.mat, 1, function(x) paste(x, collapse = ",")) # was rev(x) in the paste before
    beta.mat.sparse2 <- beta.mat.sparse[order(rownames(beta.mat.sparse)),]
    beta.mat.sparse2 <- beta.mat.sparse2[c(2:nrow(beta.mat.sparse), 1),]


    var.idx.list2 <- var.idx.list
    for (c in 1:(2 ^ ncats)) 
    {
        indi <- which(km.mat[c,] == 1)
        # this returns indices of all g terms to be included in the
        # product for the beta coefficient
        var.idx.list2[[c]] <- which(apply(km.mat, 1, function(x) all(x[indi] == 1)))
    }


    gg.mat2 <- beta.mat.sparse2
    completed <- NULL
    for (i in 0:ncats) 
    {
        cur.idx <- which(sapply(var.idx.list2, function(x) length(x)==2^i))
        for (j in 1:length(cur.idx)) {
            cur.vars <- var.idx.list2[[cur.idx[j]]]
            which.compl <- cur.vars %in% completed
            cur.var <- cur.vars[which(!which.compl)]
            div.vars <- cur.vars[which(which.compl)]

            if (sum(1L * which.compl) == 0) {
                gg.mat2[cur.var,] <- beta.mat.sparse2[cur.var,]
            } else {
                multipl <- apply(gg.mat2[div.vars,,drop=FALSE], 2, prod)
                gg.mat2[cur.var,] <- beta.mat.sparse2[cur.var,] / ifelse(multipl != 0,multipl,1)
            }
            completed <- c(completed, cur.var)
        }
    }


    beta.mat2 <- beta.mat.sparse2

    for (cc in 1:(2 ^ ncats)) 
    {
        #indi <- which(km.mat[c,] == 1)
        # this returns indices of all g terms to be included in the
        # product for the beta coefficient
        #var.idx <- which(apply(km.mat, 1, function(x) all(x[indi] == 1)))

        # compute beta
        if (length(var.idx.list2[[cc]]) > 1) 
        {
            beta.mat2[cc, ] <- apply(gg.mat2[var.idx.list2[[cc]],], 2, function(x) prod(x))
        } else 
        {
            beta.mat2[cc, ] <- gg.mat2[var.idx.list2[[cc]],]
        }

    }

    if (family == "coxph")
    {

        y.all <- Surv(y.all,c.all)
        y.test.all <- Surv(y.test.all,c.test.all)
    }

    if (family == "coxph")
    {
      list(beta.mat = beta.mat.sparse, x = x.all, y = y.all, group.ind = group.ind,
           x.test = x.test.all, y.test = y.test.all, z.test = z.test.all, group.ind.test = group.test.ind,
           beta.mat.2 = beta.mat.sparse2)
    } else {
      list(beta.mat = beta.mat.sparse, x = x.all, y = y.all, group.ind = group.ind,
           x.test = x.test.all, y.test = y.test.all, group.ind.test = group.test.ind,
           beta.mat.2 = beta.mat.sparse2)
    }
}


#' function to generate coefficient matrix with hierarchical sparsity
#'
#' @param ncats number of categories to stratify on
#' @param nvars number of variables
#' @param hier.sparsity.param parameter between 0 and 1 which determines how much hierarchical sparsity there is. To achieve
#' a desired total level of sparsity among the variables with hierarchical sparsity, this parameter can be estimated using the
#' function 'estimate.hier.sparsity.param'
#' @param avg.hier.zeros desired percent of zero variables among the variables with hierarchical zero patterns. If this is specified, it will
#' override the given hier.sparsity.param value and estimate it. This takes a while
#' @param effect.size.max maximum magnitude of the true effect sizes
#' @param misspecification.prop proportion of variables with hierarchical missingness misspecified
#'
#' @export
#' @examples
#' set.seed(123)
#'
#' # estimate hier.sparsity.param for 0.15 total proportion of nonzero variables
#' # among vars with hierarchical zero patterns
#' # NOT RUN: Takes a long time
#' # hsp <- estimate.hier.sparsity.param(ncats = 3, nvars = 25, avg.hier.zeros = 0.15, nsims = 100)
#' # the above results in the following value
#' hsp <- 0.6341772
#'
#' # check that this does indeed achieve the desired level of sparsity
#' mean(replicate(100, mean(genHierSparseBeta(ncats = 3, 
#'                            nvars = 25, hier.sparsity.param = hsp) != 0)  ))
#'
#' sparseBeta <- genHierSparseBeta(ncats = 3, nvars = 25, hier.sparsity.param = hsp)
#'
genHierSparseBeta <- function(ncats,
                              nvars,
                              hier.sparsity.param = 0.5,
                              avg.hier.zeros = NULL,
                              effect.size.max = 0.5,
                              misspecification.prop = 0)
{

    km.list <- vector(length = 2 ^ ncats, mode = "list")


    ncmb <- 2 ^ ncats
    for (i in 1:ncmb) 
    {
        km.list[[i]] <- as.integer(intToBits(i)[1:ncats])
    }
    km.mat <- Reduce(rbind, km.list)
    attr(km.mat, "dimnames") <- NULL
    ord.idx <- do.call(order, split(km.mat, rep(1:ncol(km.mat), each = nrow(km.mat))))
    km.mat <- km.mat[ord.idx,]
    km.list <- km.list[ord.idx]

    n.nzvars <- nvars

    beta.mat <- array(NA, dim = c(ncmb, n.nzvars))
    rownames(beta.mat) <- apply(km.mat, 1, function(x) paste(x, collapse = ","))
    rs <- rowSums(km.mat)
    id1 <- which(rs == 1)
    id2 <- which(rs > 1)
    id0 <- which(rs == 0)
    beta.mat[id1,] <- matrix(runif(length(id1) * n.nzvars, min=effect.size.max/2, max=effect.size.max) *
                                 (2 * rbinom(length(id1) * n.nzvars, 1,0.5)-1), ncol = n.nzvars)
    beta.mat[id1[1],] <- runif(n.nzvars, min=effect.size.max/2, max=effect.size.max) *
        (2 * rbinom(n.nzvars, 1,0.5)-1)
    beta.mat[id1[-1],] <- t(t(matrix(runif(length(id1[-1]) * n.nzvars, min=effect.size.max/2, max=effect.size.max) *
                                         (2 * rbinom(length(id1[-1]) * n.nzvars, 1,0.5)-1), ncol = n.nzvars)) )

    beta.mat[which(rs == 0),] <- runif(n.nzvars, min = effect.size.max/2, max = effect.size.max) *
        (2 * rbinom(n.nzvars, 1,0.5)-1)

    for (s in 2:ncats) 
    {
        id2 <- which(rs == s)
        for (i in id2) 
        {
            rmeans <- numeric(n.nzvars)
            #for (k in 1:ncats) {
            #    if (km.mat[i, k] == 1) {
            #        rmeans <- rmeans  + beta.mat[id1[k],] * 0.75
            #    }
            #}
            beta.mat[i,] <- runif(n.nzvars, min=effect.size.max * 0.5, max=effect.size.max) *
                (2 * rbinom(n.nzvars, 1,0.5)-1) #runif(n.nzvars, min=-0.5+rmeans, max=0.5+rmeans)#rnorm(n.nzvars, mean = rmeans, sd = 0.05)
        }
    }


    func.to.solve <- function(pct.var, nsims = 150)
    {
        avg.nonzero <- numeric(nsims)
        for (i in 1:nsims)
        {
            ## compute average number nonzero coefs
            avg.nonzero[i] <- mean(induce.hier.sparsity(beta.mat, pct.var,
                                                        km.mat, ncmb, id0,
                                                        effect.size.max,
                                                        misspecification.prop) != 0)
        }
        avg.hier.zeros - mean(avg.nonzero)
    }

    if (!is.null(avg.hier.zeros))
    {
        est.pct <- uniroot(func.to.solve, interval = c(1e-10, 1-1e-10), tol = 1e-3)
        hier.sparsity.param <- est.pct$root
    }

    beta.mat.sparse <- induce.hier.sparsity(beta.mat, hier.sparsity.param,
                                            km.mat, ncmb, id0,
                                            effect.size.max,
                                            misspecification.prop)
    beta.mat.sparse
}

#' function to estimate the hierarchical sparsity parameter for a desired level of sparsity for simulated hierarchical coefficients
#'
#' @param ncats number of categories to stratify on
#' @param nvars number of variables
#' @param avg.hier.zeros desired percent of zero variables among the variables with hierarchical zero patterns.
#' @param nsims number of simulations to estimate the average sparsity. A larger number will be more accurate but take much longer.
#' @param effect.size.max maximum magnitude of the true effect sizes
#' @param misspecification.prop proportion of variables with hierarchical missingness misspecified
#' @importFrom stats runif
#' @importFrom stats rbinom
#' @importFrom stats uniroot
#' @export
#' @examples
#' set.seed(123)
#'
#'
#' # estimate hier.sparsity.param for 0.15 total proportion of nonzero variables
#' # among vars with hierarchical zero patterns
#' \dontrun{
#' hsp <- estimate.hier.sparsity.param(ncats = 3, nvars = 25, avg.hier.zeros = 0.15, nsims = 100)
#' }
#' # the above results in the following value
#' hsp <- 0.6341772
#' 
#' # check that this does indeed achieve the desired level of sparsity
#' mean(replicate(100, mean(genHierSparseBeta(ncats = 3, 
#'                            nvars = 25, hier.sparsity.param = hsp) != 0)  ))
#'
#' sparseBeta <- genHierSparseBeta(ncats = 3, nvars = 25, hier.sparsity.param = hsp)
#'
#'
#' \dontrun{
#' hsp2 <- estimate.hier.sparsity.param(ncats = 2, nvars = 100, 
#'                         avg.hier.zeros = 0.30, nsims = 50) # 0.5778425
#' hsp3 <- estimate.hier.sparsity.param(ncats = 3, nvars = 100, 
#'                         avg.hier.zeros = 0.30, nsims = 50) # 0.4336312
#' hsp4 <- estimate.hier.sparsity.param(ncats = 4, nvars = 100, 
#'                         avg.hier.zeros = 0.30, nsims = 50) # 0.2670061
#' hsp5 <- estimate.hier.sparsity.param(ncats = 5, nvars = 100, 
#'                         avg.hier.zeros = 0.30, nsims = 50) # 0.146682
#' }
#' # 0.07551241 for hsp6
#'
estimate.hier.sparsity.param <- function(ncats,
                                         nvars,
                                         avg.hier.zeros = 0.30,
                                         nsims = 150,
                                         effect.size.max = 0.5,
                                         misspecification.prop = 0)
{

    km.list <- vector(length = 2 ^ ncats, mode = "list")

    ncmb <- 2 ^ ncats
    for (i in 1:ncmb) 
    {
        km.list[[i]] <- as.integer(intToBits(i)[1:ncats])
    }
    km.mat <- Reduce(rbind, km.list)
    attr(km.mat, "dimnames") <- NULL
    ord.idx <- do.call(order, split(km.mat, rep(1:ncol(km.mat), each = nrow(km.mat))))
    km.mat <- km.mat[ord.idx,]
    km.list <- km.list[ord.idx]

    #n.zero.cols <- floor(prop.zero.vars * nvars)
    n.nzvars    <- nvars

    beta.mat <- array(NA, dim = c(ncmb, n.nzvars))
    rownames(beta.mat) <- apply(km.mat, 1, function(x) paste(x, collapse = ","))
    rs <- rowSums(km.mat)
    id1 <- which(rs == 1)
    id2 <- which(rs > 1)
    id0 <- which(rs == 0)
    beta.mat[id1,] <- matrix(runif(length(id1) * n.nzvars, min=effect.size.max/2, max=effect.size.max) *
                                 (2 * rbinom(length(id1) * n.nzvars, 1,0.5)-1), ncol = n.nzvars)
    beta.mat[id1[1],] <- runif(n.nzvars, min=effect.size.max/2, max=effect.size.max) *
        (2 * rbinom(n.nzvars, 1,0.5)-1)
    beta.mat[id1[-1],] <- t(t(matrix(runif(length(id1[-1]) * n.nzvars, min=effect.size.max/2, max=effect.size.max) *
                                         (2 * rbinom(length(id1[-1]) * n.nzvars, 1,0.5)-1), ncol = n.nzvars)) )

    beta.mat[which(rs == 0),] <- runif(n.nzvars, min = effect.size.max/2, max = effect.size.max) *
        (2 * rbinom(n.nzvars, 1,0.5)-1)

    for (s in 2:ncats) 
    {
        id2 <- which(rs == s)
        for (i in id2) 
        {
            rmeans <- numeric(n.nzvars)
            #for (k in 1:ncats) {
            #    if (km.mat[i, k] == 1) 
            #    {
            #        rmeans <- rmeans  + beta.mat[id1[k],] * 0.75
            #    }
            #}
            beta.mat[i,] <- runif(n.nzvars, min=effect.size.max/2, max=effect.size.max) *
                (2 * rbinom(n.nzvars, 1,0.5)-1) #runif(n.nzvars, min=-0.5+rmeans, max=0.5+rmeans)#rnorm(n.nzvars, mean = rmeans, sd = 0.05)
        }
    }

    #pct.zero.cols <- 0
    #beta.mat <- t(t(beta.mat) * rbinom(n.nzvars, 1, 1 - pct.zero.cols))



    func.to.solve <- function(pct.var)
    {
        avg.nonzero <- numeric(nsims)
        for (i in 1:nsims)
        {
            ## compute average number nonzero coefs
            avg.nonzero[i] <- mean(induce.hier.sparsity(beta.mat, pct.var,
                                                        km.mat, ncmb, id0,
                                                        effect.size.max,
                                                        misspecification.prop) != 0)
        }
        avg.hier.zeros - mean(avg.nonzero)
    }
    est.pct <- uniroot(func.to.solve, interval = c(1e-10, 1-1e-10), tol = 1e-3)
    hier.sparsity.param <- est.pct$root
    hier.sparsity.param
}



induce.hier.sparsity <- function(beta.mat, pct.zero.general,
                                 km.mat, ncmb, id.no.conditions,
                                 effect.size.max, misspecification.prop)
{
    sparsity.inducer <- matrix(rbinom(prod(dim(beta.mat)), 1, 1-pct.zero.general), ncol = ncol(beta.mat))

    for (o in 1:floor(ncol(sparsity.inducer) * (1 - misspecification.prop)))
    {
        x <- sparsity.inducer[,o]

        res <- rep(1, length(x))
        for (i in 1:ncmb) 
        {
            if (x[i] == 0) 
            {
                res[i] <- 0
                cidxx <- (1:ncmb)[!(1:ncmb %in% i)]
                for (j in cidxx) 
                {
                    if (sum(km.mat[j,]) < sum(km.mat[i,])) 
                    {
                        if (any(km.mat[i,] == km.mat[j,] & (km.mat[j,] == 1 | all(km.mat[j,] == 0))   )) 
                        {
                            res[j] <- 0
                        }
                    }
                }
            }

        }
        sparsity.inducer[,o] <- res
    }


    for (o in 1:floor(ncol(sparsity.inducer) * (1 - misspecification.prop)))
    {
        x <- sparsity.inducer[,o]

        res <- rep(1, length(x))
        for (i in 1:ncmb) 
        {
            if (x[i] == 0) 
            {
                res[i] <- 0
                cidxx <- (1:ncmb)[!(1:ncmb %in% i)]
                for (j in cidxx) 
                {
                    if (sum(km.mat[j,]) < sum(km.mat[i,])) 
                    {
                        if (any(km.mat[i,] == km.mat[j,] & (km.mat[j,] == 1 | all(km.mat[j,] == 0))   )) 
                        {
                            res[j] <- 0
                        }
                    }
                }
            }

        }
        sparsity.inducer[,o] <- res
    }



    rownames(sparsity.inducer) <- rownames(beta.mat)
    beta.mat <- beta.mat * sparsity.inducer
    ## missingness for the strata with NO conditions should
    ## not have hierarchical missingness
    beta.mat[id.no.conditions,] <- runif(ncol(beta.mat), min=effect.size.max/2,
                                         max=effect.size.max) * (2 * rbinom(ncol(beta.mat), 1,0.5)-1) *
        rbinom(ncol(beta.mat), 1, 1 - pct.zero.general)
    beta.mat
}


