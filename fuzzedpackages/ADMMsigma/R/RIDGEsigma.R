## Matt Galloway


#' @title Ridge penalized precision matrix estimation
#' 
#' @description Ridge penalized matrix estimation via closed-form solution. If one is only interested in the ridge penalty, this function will be faster and provide a more precise estimate than using \code{ADMMsigma}. \cr
#' Consider the case where
#' \eqn{X_{1}, ..., X_{n}} are iid \eqn{N_{p}(\mu, \Sigma)}
#' and we are tasked with estimating the precision matrix,
#' denoted \eqn{\Omega \equiv \Sigma^{-1}}. This function solves the
#' following optimization problem:
#' \describe{
#' \item{Objective:}{
#' \eqn{\hat{\Omega}_{\lambda} = \arg\min_{\Omega \in S_{+}^{p}}
#' \left\{ Tr\left(S\Omega\right) - \log \det\left(\Omega \right) +
#' \frac{\lambda}{2}\left\| \Omega \right\|_{F}^{2} \right\}}}
#' }
#' where \eqn{\lambda > 0} and \eqn{\left\|\cdot \right\|_{F}^{2}} is the Frobenius
#' norm.
#'
#' @param X option to provide a nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
#' @param lam positive tuning parameters for ridge penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-2, 2, 0.1)}.
#' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores will be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
#' @param K specify the number of folds for cross validation.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.

#' 
#' @return returns class object \code{RIDGEsigma} which includes:
#' \item{Lambda}{optimal tuning parameter.}
#' \item{Lambdas}{grid of lambda values for CV.}
#' \item{Omega}{estimated penalized precision matrix.}
#' \item{Sigma}{estimated covariance matrix from the penalized precision matrix (inverse of Omega).}
#' \item{Path}{array containing the solution path. Solutions are ordered dense to sparse.}
#' \item{Gradient}{gradient of optimization function (penalized gaussian likelihood).}
#' \item{MIN.error}{minimum average cross validation error (cv.crit) for optimal parameters.}
#' \item{AVG.error}{average cross validation error (cv.crit) across all folds.}
#' \item{CV.error}{cross validation errors (cv.crit).}
#' 
#' @references
#' \itemize{
#' \item Rothman, Adam. 2017. 'STAT 8931 notes on an algorithm to compute the Lasso-penalized Gaussian likelihood precision matrix estimator.'
#' }
#' 
#' @author Matt Galloway \email{gall0441@@umn.edu}
#' 
#' @seealso \code{\link{plot.RIDGE}}, \code{\link{ADMMsigma}}
#' 
#' @export
#' 
#' @examples
#' # generate data from a sparse matrix
#' # first compute covariance matrix
#' S = matrix(0.7, nrow = 5, ncol = 5)
#' for (i in 1:5){
#'  for (j in 1:5){
#'    S[i, j] = S[i, j]^abs(i - j)
#'  }
#'  }
#'
#' # generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
#' set.seed(123)
#' Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' out = eigen(S, symmetric = TRUE)
#' S.sqrt = out$vectors %*% diag(out$values^0.5)
#' S.sqrt = S.sqrt %*% t(out$vectors)
#' X = Z %*% S.sqrt
#'
#' # ridge penalty no ADMM
#' RIDGEsigma(X, lam = 10^seq(-5, 5, 0.5))

# we define the ADMM covariance estimation function
RIDGEsigma = function(X = NULL, S = NULL, lam = 10^seq(-2, 2, 0.1), 
    path = FALSE, K = 5, cores = 1, trace = c("none", "progress", 
        "print")) {
    
    # checks
    if (is.null(X) && is.null(S)) {
        stop("Must provide entry for X or S!")
    }
    if (!all(lam > 0)) {
        stop("lam must be positive!")
    }
    if (all(c(K, cores)%%1 != 0)) {
        stop("Entry must be an integer!")
    }
    if (cores < 1) {
        stop("Number of cores must be positive!")
    }
    if (cores > 1 && path) {
        cat("Parallelization not possible when producing solution path. Setting cores = 1...\n\n")
        cores = 1
    }
    if (path) {
        K = 1
    }
    
    # match values
    call = match.call()
    trace = match.arg(trace)
    MIN.error = AVG.error = CV.error = NULL
    Lambdas = lam = sort(lam)
    n = ifelse(is.null(X), nrow(S), nrow(X))
    
    
    # compute sample covariance matrix, if necessary
    if (is.null(S)) {
        S = (nrow(X) - 1)/nrow(X) * cov(X)
    }
    
    # perform cross validation, if necessary
    if ((length(lam) > 1) & (!is.null(X) || path)) {
        
        # run CV in parallel?
        if (cores > 1) {
            
            # execute ParallelCV
            RIDGE = CVP_RIDGE(X = X, lam = lam, K = K, cores = cores, 
                trace = trace)
            MIN.error = RIDGE$min.error
            AVG.error = RIDGE$avg.error
            CV.error = RIDGE$cv.error
            lam = RIDGE$lam
            
        } else {
            
            # execute CV_RIDGEsigma
            if (is.null(X)) {
                X = matrix(0)
            }
            RIDGE = CV_RIDGEc(X = X, S = S, lam = lam, path = path, 
                K = K, trace = trace)
            MIN.error = RIDGE$min.error
            AVG.error = RIDGE$avg.error
            CV.error = RIDGE$cv.error
            lam = RIDGE$lam
            Path = RIDGE$path
            
        }
        
        # print warning if lam on boundary
        if ((lam == Lambdas[1]) || (lam == Lambdas[length(Lambdas)]) && 
            (length(Lambdas) != 1)) {
            cat("\nOptimal tuning parameter on boundary...!")
        }
        
        # compute final estimate at best tuning parameters
        Omega = RIDGEc(S = S, lam = lam)
        
        
    } else {
        
        # execute RIDGEsigmac
        if (length(lam) > 1) {
            stop("Must specify X, set path = TRUE, or provide single value for lam.")
        }
        Omega = RIDGEc(S = S, lam = lam)
        
    }
    
    # compute gradient
    grad = S - qr.solve(Omega) + lam * Omega
    
    # compute penalized loglik
    loglik = (-n/2) * (sum(Omega * S) - determinant(Omega, logarithm = TRUE)$modulus[1] + 
        lam * sum(Omega^2))
    
    # return values
    tuning = matrix(c(log10(lam), lam), ncol = 2)
    colnames(tuning) = c("log10(lam)", "lam")
    if (!path) {
        Path = NULL
    }
    
    returns = list(Call = call, Lambda = tuning, Lambdas = Lambdas, 
        Omega = Omega, Sigma = qr.solve(Omega), Path = Path, Gradient = grad, 
        Loglik = loglik, MIN.error = MIN.error, AVG.error = AVG.error, 
        CV.error = CV.error)
    
    class(returns) = "RIDGE"
    return(returns)
    
}





##-----------------------------------------------------------------------------------



#' @title Print RIDGE object
#' @description Prints RIDGE object and suppresses output if needed.
#' @param x class object RIDGE.
#' @param ... additional arguments.
#' @keywords internal
#' @export
print.RIDGE = function(x, ...) {
    
    # print call
    cat("\nCall: ", paste(deparse(x$Call), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
    
    # print optimal tuning parameter
    cat("\nTuning parameter:\n")
    print.default(round(x$Lambda, 3), print.gap = 2L, quote = FALSE)
    
    # print loglik
    cat("\nLog-likelihood: ", paste(round(x$Loglik, 5), sep = "\n", 
        collapse = "\n"), "\n", sep = "")
    
    # print Omega if dim <= 10
    if (nrow(x$Omega) <= 10) {
        cat("\nOmega:\n")
        print.default(round(x$Omega, 5))
    } else {
        cat("\n(...output suppressed due to large dimension!)\n")
    }
    
}



##-----------------------------------------------------------------------------------




#' @title Plot RIDGE object
#' @description Produces a heat plot for the cross validation errors, if available.
#' @param x class object RIDGE
#' @param type produce either 'heatmap' or 'line' graph
#' @param footnote option to print footnote of optimal values. Defaults to TRUE.
#' @param ... additional arguments.
#' @export
#' @examples
#' # generate data from a sparse matrix
#' # first compute covariance matrix
#' S = matrix(0.7, nrow = 5, ncol = 5)
#' for (i in 1:5){
#'  for (j in 1:5){
#'    S[i, j] = S[i, j]^abs(i - j)
#'  }
#'  }
#'
#' # generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
#' set.seed(123)
#' Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' out = eigen(S, symmetric = TRUE)
#' S.sqrt = out$vectors %*% diag(out$values^0.5)
#' S.sqrt = S.sqrt %*% t(out$vectors)
#' X = Z %*% S.sqrt
#'
#' # produce CV heat map for RIDGEsigma
#' plot(RIDGEsigma(X, lam = 10^seq(-5, 5, 0.5)))
#' 
#' # produce line graph for RIDGEsigma
#' plot(RIDGEsigma(X), type = 'line')

plot.RIDGE = function(x, type = c("heatmap", "line"), footnote = TRUE, 
    ...) {
    
    # check
    type = match.arg(type)
    Means = NULL
    if (is.null(x$CV.error)) {
        stop("No cross validation errors to plot!")
    }
    
    if (type == "line") {
        
        # gather values to plot
        cv = cbind(expand.grid(lam = x$Lambdas, alpha = 0), Errors = as.data.frame.table(x$CV.error)$Freq)
        
        # produce line graph
        graph = ggplot(summarise(group_by(cv, lam), Means = mean(Errors)), 
            aes(log10(lam), Means)) + geom_jitter(width = 0.2, 
            color = "navy blue") + theme_minimal() + geom_line(color = "red") + 
            labs(title = "Cross-Validation Errors", y = "Error")
        
    } else {
        
        # augment values for heat map (helps visually)
        lam = x$Lambdas
        cv = expand.grid(lam = lam, alpha = 0)
        Errors = 1/(c(x$AVG.error) + abs(min(x$AVG.error)) + 1)
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
            x$Lambda[1], sep = ""))
        
    } else {
        
        # produce without footnote
        graph
        
    }
}
