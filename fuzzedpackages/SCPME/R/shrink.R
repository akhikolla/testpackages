## Matt Galloway


#' @title Shrinking characteristics of precision matrix estimators
#' 
#' @description Shrinking characteristics of precision matrix estimators. Penalized precision matrix estimation using the ADMM algorithm.
#' Consider the case where \eqn{X_{1}, ..., X_{n}} are iid \eqn{N_{p}(\mu,
#' \Sigma)} and we are tasked with estimating the precision matrix,
#' denoted \eqn{\Omega \equiv \Sigma^{-1}}. This function solves the
#' following optimization problem:
#' \describe{
#' \item{Objective:}{
#' \eqn{\hat{\Omega}_{\lambda} = \arg\min_{\Omega \in S_{+}^{p}}
#' \left\{ Tr\left(S\Omega\right) - \log \det\left(\Omega \right) +
#' \lambda\left\| A \Omega B - C \right\|_{1} \right\}}}
#' }
#' where \eqn{\lambda > 0} and we define
#' \eqn{\left\|A \right\|_{1} = \sum_{i, j} \left| A_{ij} \right|}.
#' 
#' @details For details on the implementation of 'shrink', see the vignette
#' \url{https://mgallow.github.io/SCPME/}.
#' 
#' @param X option to provide a nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param Y option to provide nxr response matrix. Each row corresponds to a single response and each column contains n response of a single feature/response.
#' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
#' @param A option to provide user-specified matrix for penalty term. This matrix must have p columns. Defaults to identity matrix.
#' @param B option to provide user-specified matrix for penalty term. This matrix must have p rows. Defaults to identity matrix.
#' @param C option to provide user-specified matrix for penalty term. This matrix must have nrow(A) rows and ncol(B) columns. Defaults to zero matrix.
#' @param nlam number of \code{lam} tuning parameters for penalty term generated from \code{lam.min.ratio} and \code{lam.max} (automatically generated). Defaults to 10.
#' @param lam.max option to specify the maximum \code{lam} tuning parameter. Defaults to NULL.
#' @param lam.min.ratio smallest \code{lam} value provided as a fraction of \code{lam.max}. The function will automatically generate \code{nlam} tuning parameters from \code{lam.min.ratio*lam.max} to \code{lam.max} in log10 scale. If \code{lam.max = NULL}, \code{lam.max} is calculated to be the smallest \code{lam} such that all off-diagonal entries in \code{Omega} are equal to zero. Defaults to 1e-3.
#' @param lam option to provide positive tuning parameters for penalty term. This will cause \code{nlam} and \code{lam.min.ratio} to be disregarded. If a vector of parameters is provided, they should be in increasing order. Defaults to NULL.
#' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. Alpha must be a single value (cross validation across alpha not supported).
#' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores must be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
#' @param rho initial step size for ADMM algorithm.
#' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
#' @param tau.rho factor in which to increase/decrease step size \code{rho}
#' @param iter.rho step size \code{rho} will be updated every \code{iter.rho} steps
#' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit = loglik} then iterations will stop when the relative change in log-likelihood is less than \code{tol.abs}. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
#' @param tol.abs absolute convergence tolerance. Defaults to 1e-4.
#' @param tol.rel relative convergence tolerance. Defaults to 1e-4.
#' @param maxit maximum number of iterations. Defaults to 1e4.
#' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param crit.cv cross validation criterion (\code{MSE}, \code{loglik}, \code{penloglik}, \code{AIC}, or \code{BIC}). Defaults to \code{MSE}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.

#' 
#' @return returns class object \code{ADMMsigma} which includes:
#' \item{Call}{function call.}
#' \item{Iterations}{number of iterations.}
#' \item{Tuning}{optimal tuning parameter.}
#' \item{Lambdas}{grid of lambda values for CV.}
#' \item{maxit}{maximum number of iterations.}
#' \item{Omega}{estimated penalized precision matrix.}
#' \item{Sigma}{estimated covariance matrix from the penalized precision matrix (inverse of Omega).}
#' \item{Path}{array containing the solution path. Solutions will be ordered in ascending alpha values for each lambda.}
#' \item{Z}{final sparse update of estimated penalized precision matrix.}
#' \item{Y}{final dual update.}
#' \item{rho}{final step size.}
#' \item{Loglik}{penalized log-likelihood for Omega}
#' \item{MIN.error}{minimum average cross validation error (cv.crit) for optimal parameters.}
#' \item{AVG.error}{average cross validation error (cv.crit) across all folds.}
#' \item{CV.error}{cross validation errors (cv.crit).}
#' 
#' @references
#' \itemize{
#' \item Boyd, Stephen, Neal Parikh, Eric Chu, Borja Peleato, Jonathan Eckstein, and others. 2011. 'Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.' \emph{Foundations and Trends in Machine Learning} 3 (1). Now Publishers, Inc.: 1-122. \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}
#' \item Hu, Yue, Chi, Eric C, amd Allen, Genevera I. 2016. 'ADMM Algorithmic Regularization Paths for Sparse Statistical Machine Learning.' \emph{Splitting Methods in Communication, Imaging, Science, and Engineering}. Springer: 433-459.
#' \item Molstad, Aaron J., and Adam J. Rothman. (2017). 'Shrinking Characteristics of Precision Matrix Estimators. \emph{Biometrika.}. \url{https://doi.org/10.1093/biomet/asy023}
#' \item Rothman, Adam. 2017. 'STAT 8931 notes on an algorithm to compute the Lasso-penalized Gaussian likelihood precision matrix estimator.'
#' }
#' 
#' @author Matt Galloway \email{gall0441@@umn.edu}
#' 
#' @seealso \code{\link{plot.shrink}}
#' @export
#' 
#' @examples
#' # generate some data
#' data = data_gen(n = 100, p = 5, r = 1)
#' 
#' # lasso penalized omega (glasso)
#' shrink(X = data$X, Y = data$Y)
#' 
#' # lasso penalized beta (print estimated omega)
#' lam.max = max(abs(t(data$X) %*% data$Y))
#' (shrink = shrink(X = data$X, Y = data$Y, B = cov(data$X, data$Y), lam.max = lam.max))
#' 
#' # print estimated beta
#' shrink$Z

# we define the ADMM covariance estimation function
shrink = function(X = NULL, Y = NULL, S = NULL, A = diag(ncol(S)), 
    B = diag(ncol(S)), C = matrix(0, ncol = ncol(B), nrow = ncol(A)), 
    nlam = 10, lam.max = NULL, lam.min.ratio = 0.001, lam = NULL, 
    alpha = 1, path = FALSE, rho = 2, mu = 10, tau.rho = 2, 
    iter.rho = 10, crit = c("ADMM", "loglik"), tol.abs = 1e-04, 
    tol.rel = 1e-04, maxit = 10000, adjmaxit = NULL, K = 5, 
    crit.cv = c("MSE", "loglik", "penloglik", "AIC", "BIC"), 
    start = c("warm", "cold"), cores = 1, trace = c("progress", 
        "print", "none")) {
    
    
    # checks
    if (is.null(X) && is.null(S)) {
        stop("Must provide entry for X or S!")
    }
    if (!all(lam > 0)) {
        stop("lam must be positive!")
    }
    if (!(all(c(rho, mu, tau.rho, iter.rho, tol.abs, tol.rel, 
        maxit, adjmaxit, K, cores) > 0))) {
        stop("Entry must be positive!")
    }
    if ((alpha < 0) || (alpha > 1)) {
        stop("Alpha must be between 0 and 1!")
    }
    if (!(all(sapply(c(rho, mu, tau.rho, iter.rho, tol.abs, 
        tol.rel, maxit, adjmaxit, K, cores, nlam, lam.min.ratio, 
        alpha), length) <= 1))) {
        stop("Entry must be single value!")
    }
    if (all(c(maxit, adjmaxit, K, cores)%%1 != 0)) {
        stop("Entry must be an integer!")
    }
    if (cores < 1) {
        stop("Number of cores must be positive!")
    }
    if (((length(lam) > 1) & (!path || (crit.cv == "MSE"))) & 
        (is.null(X) || is.null(Y))) {
        stop("Must provide entry for X and Y!")
    }
    if (is.null(Y) && (crit.cv == "MSE")) {
        cat("Matrix Y not detected... will use loglik for crit.cv instead!\n\n")
        crit.cv = "loglik"
    }
    if (cores > 1 && path) {
        cat("Parallelization not possible when producing solution path. Setting cores = 1...\n\n")
        cores = 1
    }
    if (path) {
        if (match.arg(crit.cv) == "MSE") {
            cat("MSE crit.cv not available when path == TRUE... setting crit.cv == loglik\n\n")
            crit.cv = "loglik"
        }
        K = 1
    }
    if (cores > K) {
        cat("Number of cores exceeds K... setting cores = K\n\n")
        cores = K
    }
    if (is.null(adjmaxit)) {
        adjmaxit = maxit
    }
    
    
    # match values
    crit = match.arg(crit)
    crit.cv = match.arg(crit.cv)
    start = match.arg(start)
    trace = match.arg(trace)
    call = match.call()
    MIN.error = AVG.error = CV.error = NULL
    n = ifelse(is.null(X), nrow(S), nrow(X))
    if (is.null(Y)) {
        Y = matrix(0)
    }
    
    # compute sample covariance matrix, if necessary
    if (is.null(S)) {
        S = (nrow(X) - 1)/nrow(X) * cov(X)
        if (nrow(A) == 0) {
            A = diag(ncol(S))
        }
        if (nrow(B) == 0) {
            B = diag(ncol(S))
        }
        if (nrow(C) == 0) {
            C = matrix(0, ncol = ncol(B), nrow = ncol(A))
        }
    }
    
    # calculate tau used in algorithm
    tau = max(eigen(crossprod(A))$values) * max(eigen(tcrossprod(B))$values) + 
        1e-08
    
    # more checks
    if (ncol(A) != ncol(S)) {
        stop("Matrix A has incompatible number of columns!")
    }
    if (nrow(B) != ncol(S)) {
        stop("Matrix B has incompatible number of rows!")
    }
    if ((nrow(C) != nrow(A)) || (ncol(C) != ncol(B))) {
        stop("Matrix C has incompatible dimensions!")
    }
    
    # compute grid of lam values, if necessary
    if (is.null(lam)) {
        if (!((lam.min.ratio <= 1) && (lam.min.ratio > 
            0))) {
            cat("lam.min.ratio must be in (0, 1]... setting to 1e-2!\n\n")
            lam.min.ratio = 0.01
        }
        if (!((nlam > 0) && (nlam%%1 == 0))) {
            cat("nlam must be a positive integer... setting to 10!\n\n")
            nlam = 10
        }
        
        # calculate lam.max and lam.min
        if (is.null(lam.max)) {
            lam.max = max(abs(S - diag(S)))
        }
        lam.min = lam.min.ratio * lam.max
        
        # calculate grid of lambda values
        lam = 10^seq(log10(lam.min), log10(lam.max), length = nlam)
        
    } else {
        
        # sort lambda values
        lam = sort(lam)
        
    }
    
    # perform cross validation, if necessary
    initOmega = diag(diag(S)^(-1))
    init = A %*% initOmega %*% B - C
    zeros = matrix(0, nrow = nrow(C), ncol = ncol(C))
    if ((length(lam) > 1) && (!is.null(X) || path)) {
        
        # run CV in parallel?
        if (cores > 1) {
            
            # execute CVP_ADMM
            ADMM = CVP_ADMM(X = X, Y = Y, A = A, B = B, 
                C = C, lam = lam, alpha = alpha, tau = tau, 
                rho = rho, mu = mu, tau.rho = tau.rho, 
                iter.rho = iter.rho, crit = crit, tol.abs = tol.abs, 
                tol.rel = tol.rel, maxit = maxit, adjmaxit = adjmaxit, 
                K = K, crit.cv = crit.cv, start = start, 
                cores = cores, trace = trace)
            MIN.error = ADMM$min.error
            AVG.error = ADMM$avg.error
            CV.error = ADMM$cv.error
            
        } else {
            
            # execute CV_ADMMc
            if (is.null(X)) {
                X = matrix(0)
            }
            ADMM = CV_ADMMc(X = X, S = S, Y = Y, A = A, 
                B = B, C = C, lam = lam, alpha = alpha, 
                path = path, tau = tau, rho = rho, mu = mu, 
                tau_rho = tau.rho, iter_rho = iter.rho, 
                crit = crit, tol_abs = tol.abs, tol_rel = tol.rel, 
                maxit = maxit, adjmaxit = adjmaxit, K = K, 
                crit_cv = crit.cv, start = start, trace = trace)
            MIN.error = ADMM$min.error
            AVG.error = ADMM$avg.error
            CV.error = ADMM$cv.error
            Path = ADMM$path
            
        }
        
        # print warning if lam on boundary
        if (((ADMM$lam == lam[1]) || ADMM$lam == lam[length(lam)]) && 
            ((length(lam) != 1) && (!path))) {
            cat("\nOptimal tuning parameter on boundary...!\n")
        }
        
        # compute final estimate at best tuning parameters
        ADMM = ADMMc(S = S, A = A, B = B, C = C, initOmega = initOmega, 
            initZ = init, initY = zeros, lam = ADMM$lam, 
            alpha = alpha, tau = tau, rho = rho, mu = mu, 
            tau_rho = tau.rho, iter_rho = iter.rho, crit = crit, 
            tol_abs = tol.abs, tol_rel = tol.rel, maxit = maxit)
        
        
    } else {
        
        # execute ADMM_sigmac
        if (length(lam) > 1) {
            stop("Must set specify X, set path = TRUE, or provide single value for lam.")
        }
        
        ADMM = ADMMc(S = S, A = A, B = B, C = C, initOmega = initOmega, 
            initZ = init, initY = zeros, lam = lam, alpha = alpha, 
            tau = tau, rho = rho, mu = mu, tau_rho = tau.rho, 
            iter_rho = iter.rho, crit = crit, tol_abs = tol.abs, 
            tol_rel = tol.rel, maxit = maxit)
        
    }
    
    # compute penalized loglik
    loglik = (-n/2) * (sum(ADMM$Omega * S) - determinant(ADMM$Omega, 
        logarithm = TRUE)$modulus[1] + ADMM$lam * ((1 - 
        alpha)/2 * sum((A %*% ADMM$Omega %*% B - C)^2) + 
        alpha * sum(abs(A %*% ADMM$Omega %*% B - C))))
    
    
    # return values
    tuning = matrix(c(log10(ADMM$lam), ADMM$lam), ncol = 2)
    colnames(tuning) = c("log10(lam)", "lam")
    if (!path) {
        Path = NULL
    }
    
    returns = list(Call = call, Iterations = ADMM$Iterations, 
        Tuning = tuning, Lambdas = lam, maxit = maxit, 
        Omega = ADMM$Omega, Sigma = qr.solve(ADMM$Omega), 
        Path = Path, Z = ADMM$Z, Y = ADMM$Y, rho = ADMM$rho, 
        Loglik = loglik, MIN.error = MIN.error, AVG.error = AVG.error, 
        CV.error = CV.error)
    
    class(returns) = "shrink"
    return(returns)
    
}






##-----------------------------------------------------------------------------------



#' @title Print shrink object
#' @description Prints shrink object and suppresses output if needed.
#' @param x class object shrink
#' @param ... additional arguments.
#' @keywords internal
#' @export
print.shrink = function(x, ...) {
    
    # print warning if maxit reached
    if (x$maxit <= x$Iterations) {
        cat("\nMaximum iterations reached...!")
    }
    
    # print call
    cat("\nCall: ", paste(deparse(x$Call), sep = "\n", 
        collapse = "\n"), "\n", sep = "")
    
    # print iterations
    cat("\nIterations: ", paste(x$Iterations, sep = "\n", 
        collapse = "\n"), "\n", sep = "")
    
    # print optimal tuning parameters
    cat("\nTuning parameters:\n")
    print.default(round(x$Tuning, 3), print.gap = 2L, quote = FALSE)
    
    # print loglik
    cat("\nLog-likelihood: ", paste(round(x$Loglik, 5), 
        sep = "\n", collapse = "\n"), "\n", sep = "")
    
    # print Omega if dim <= 10
    if (nrow(x$Omega) <= 10) {
        cat("\nOmega:\n")
        print.default(round(x$Omega, 5))
    } else {
        cat("\n(...output suppressed due to large dimension!)\n")
    }
    
}



##-----------------------------------------------------------------------------------




#' @title Plot shrink object
#' @description Produces a plot for the cross validation errors, if available.
#' @param x class object shrink.
#' @param type produce either 'heatmap' or 'line' graph
#' @param footnote option to print footnote of optimal values. Defaults to TRUE.
#' @param ... additional arguments.
#' @export
#' @examples
#' # generate some data
#' data = data_gen(n = 100, p = 5, r = 1)
#' 
#' # lasso penalized beta (print estimated omega)
#' lam.max = max(abs(t(data$X) %*% data$Y))
#' (shrink = shrink(X = data$X, Y = data$Y, B = cov(data$X, data$Y), lam.max = lam.max))
#' 
#' # print estimated beta
#' shrink$Z
#' 
#' # create heatmap of CV errors
#' plot(shrink, type = 'heatmap')
#' 
#' # create line graph of CV errors
#' plot(shrink)

plot.shrink = function(x, type = c("line", "heatmap"), 
    footnote = TRUE, ...) {
    
    # check
    type = match.arg(type)
    Means = NULL
    if (is.null(x$CV.error)) {
        stop("No cross validation errors to plot!")
    }
    
    if (type == "line") {
        
        # gather values to plot
        cv = cbind(expand.grid(lam = x$Lambdas, alpha = 0), 
            Errors = as.data.frame.table(x$CV.error)$Freq)
        
        # produce line graph
        graph = ggplot(summarise(group_by(cv, lam), Means = mean(Errors)), 
            aes(log10(lam), Means)) + geom_jitter(width = 0.2, 
            color = "navy blue") + theme_minimal() + geom_line(color = "red") + 
            labs(title = "Cross-Validation Errors", y = "Error") + 
            geom_vline(xintercept = x$Tuning[1], linetype = "dotted")
        
    } else {
        
        # augment values for heat map (helps visually)
        lam = x$Lambdas
        cv = expand.grid(lam = lam, alpha = 0)
        Errors = 1/(c(x$AVG.error) + abs(min(x$AVG.error)) + 
            1)
        cv = cbind(cv, Errors)
        
        # design color palette
        bluetowhite <- c("#000E29", "white")
        
        # produce ggplot heat map
        graph = ggplot(cv, aes(alpha, log10(lam))) + geom_raster(aes(fill = Errors)) + 
            scale_fill_gradientn(colours = colorRampPalette(bluetowhite)(2), 
                guide = "none") + theme_minimal() + labs(title = "Heatmap of Cross-Validation Errors") + 
            theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
                axis.ticks.x = element_blank())
        
    }
    
    if (footnote) {
        
        # produce with footnote
        graph + labs(caption = paste("**Optimal: log10(lam) = ", 
            round(x$Tuning[1], 3), sep = ""))
        
    } else {
        
        # produce without footnote
        graph
        
    }
}
