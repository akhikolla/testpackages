## Matt Galloway


#' @title Parallel CV (uses CV_ADMMc)
#' @description Parallel implementation of cross validation.
#'
#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-2, 2, 0.2)}.
#' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{seq(-1, 1, 0.2)}.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param rho initial step size for ADMM algorithm.
#' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
#' @param tau.inc factor in which to increase step size \code{rho}
#' @param tau.dec factor in which to decrease step size \code{rho}
#' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit = loglik} then iterations will stop when the relative change in log-likelihood is less than \code{tol.abs}. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
#' @param tol.abs absolute convergence tolerance. Defaults to 1e-4.
#' @param tol.rel relative convergence tolerance. Defaults to 1e-4.
#' @param maxit maximum number of iterations. Defaults to 1e3.
#' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param crit.cv cross validation criterion (\code{loglik}, \code{penloglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
#' 
#' @return returns list of returns which includes:
#' \item{lam}{optimal tuning parameter.}
#' \item{alpha}{optimal tuning parameter.}
#' \item{min.error}{minimum average cross validation error (cv.crit) for optimal parameters.}
#' \item{avg.error}{average cross validation error (cv.crit) across all folds.}
#' \item{cv.error}{cross validation errors (cv.crit).}
#' 
#' @keywords internal

# we define the CV_ADMMc function
CVP_ADMM = function(X = NULL, lam = 10^seq(-2, 2, 0.2), alpha = seq(0, 
    1, 0.2), diagonal = FALSE, rho = 2, mu = 10, tau.inc = 2, tau.dec = 2, 
    crit = c("ADMM", "loglik"), tol.abs = 1e-04, tol.rel = 1e-04, 
    maxit = 1000, adjmaxit = NULL, K = 5, crit.cv = c("loglik", 
        "penloglik", "AIC", "BIC"), start = c("warm", "cold"), 
    cores = 1, trace = c("progress", "print", "none")) {
    
    # match values
    crit = match.arg(crit)
    crit.cv = match.arg(crit.cv)
    start = match.arg(start)
    trace = match.arg(trace)
    lam = sort(lam)
    alpha = sort(alpha)
    
    # make cluster and register cluster
    num_cores = detectCores()
    if (cores > num_cores) {
        cat("\nOnly detected", paste(num_cores, "cores...", sep = " "))
    }
    if (cores > K) {
        cat("\nNumber of cores exceeds K... setting cores = K")
        cores = K
    }
    
    cluster = makeCluster(cores)
    registerDoParallel(cluster)
    
    # use cluster for each fold in CV
    n = nrow(X)
    ind = sample(n)
    k = NULL
    CV = foreach(k = 1:K, .packages = "ADMMsigma", .inorder = FALSE) %dopar% 
        {
            
            leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * 
                n/K)]
            
            # training set
            X.train = X[-leave.out, , drop = FALSE]
            X_bar = apply(X.train, 2, mean)
            X.train = scale(X.train, center = X_bar, scale = FALSE)
            
            # validation set
            X.valid = X[leave.out, , drop = FALSE]
            X.valid = scale(X.valid, center = X_bar, scale = FALSE)
            
            # sample covariances
            S.train = crossprod(X.train)/(dim(X.train)[1])
            S.valid = crossprod(X.valid)/(dim(X.valid)[1])
            
            # run foreach loop on CVP_ADMMc
            CVP_ADMMc(n = nrow(X.valid), S_train = S.train, S_valid = S.valid, 
                lam = lam, alpha = alpha, diagonal = diagonal, 
                rho = rho, mu = mu, tau_inc = tau.inc, tau_dec = tau.dec, 
                crit = crit, tol_abs = tol.abs, tol_rel = tol.rel, 
                maxit = maxit, adjmaxit = adjmaxit, crit_cv = crit.cv, 
                start = start, trace = trace)
            
        }
    
    # determine optimal tuning parameters
    CV = array(as.numeric(unlist(CV)), dim = c(length(lam), length(alpha), 
        K))
    AVG = apply(CV, c(1, 2), mean)
    best = which(AVG == min(AVG), arr.ind = TRUE)
    error = min(AVG)
    best_lam = lam[best[1]]
    best_alpha = alpha[best[2]]
    
    # stop cluster
    stopCluster(cluster)
    
    # return best lam and alpha values
    return(list(lam = best_lam, alpha = best_alpha, min.error = error, 
        avg.error = AVG, cv.error = CV))
    
}




##-----------------------------------------------------------------------------------




#' @title Parallel Ridge CV (uses CVP_RIDGEc)
#' @description Parallel implementation of cross validation for RIDGEsigma.
#'
#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param lam positive tuning parameters for ridge penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-2, 2, 0.1)}.
#' @param K specify the number of folds for cross validation.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
#' 
#' @return returns list of returns which includes:
#' \item{lam}{optimal tuning parameter.}
#' \item{min.error}{minimum average cross validation error for optimal parameters.}
#' \item{avg.error}{average cross validation error across all folds.}
#' \item{cv.error}{cross validation errors (negative validation likelihood).}
#' 
#' @keywords internal

# we define the CVP_RIDGE function
CVP_RIDGE = function(X = NULL, lam = 10^seq(-2, 2, 0.1), K = 5, 
    cores = 1, trace = c("none", "progress", "print")) {
    
    # make cluster and register cluster
    num_cores = detectCores()
    if (cores > num_cores) {
        cat("\nOnly detected", paste(num_cores, "cores...", sep = " "))
    }
    if (cores > K) {
        cat("\nNumber of cores exceeds K... setting cores = K")
        cores = K
    }
    
    cluster = makeCluster(cores)
    registerDoParallel(cluster)
    
    # use cluster for each fold in CV
    n = dim(X)[1]
    ind = sample(n)
    lam = sort(lam)
    k = NULL
    CV = foreach(k = 1:K, .packages = "ADMMsigma", .combine = "cbind", 
        .inorder = FALSE) %dopar% {
        
        leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * n/K)]
        
        # training set
        X.train = X[-leave.out, , drop = FALSE]
        X_bar = apply(X.train, 2, mean)
        X.train = scale(X.train, center = X_bar, scale = FALSE)
        
        # validation set
        X.valid = X[leave.out, , drop = FALSE]
        X.valid = scale(X.valid, center = X_bar, scale = FALSE)
        
        # sample covariances
        S.train = crossprod(X.train)/(dim(X.train)[1])
        S.valid = crossprod(X.valid)/(dim(X.valid)[1])
        
        # run foreach loop on CVP_RIDGEc
        CVP_RIDGEc(n = nrow(X.valid), S_train = S.train, S_valid = S.valid, 
            lam = lam, trace = trace)
        
    }
    
    # determine optimal tuning parameters
    AVG = as.matrix(apply(CV, 1, mean))
    best = which(AVG == min(AVG), arr.ind = TRUE)
    error = min(AVG)
    best_lam = lam[best[1]]
    
    # stop cluster
    stopCluster(cluster)
    
    # return best lam and alpha values
    return(list(lam = best_lam, min.error = error, avg.error = AVG, 
        cv.error = CV))
    
}
