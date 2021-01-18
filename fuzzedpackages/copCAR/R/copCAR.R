
#' Return an adjacency matrix for a square lattice.
#'
#' @details This function builds the adjacency matrix for the \code{m} by \code{n} square lattice.
#'
#' @param m the number of rows in the lattice.
#' @param n the number of columns in the lattice. Defaults to \code{NULL}. If missing, the lattice is assumed to be \code{m} by \code{m}. 
#
#' @return A matrix \eqn{A} of 0s and 1s, where \eqn{A_{ij}} is equal to 1 iff vertices \eqn{i} and \eqn{j} are adjacent.
#'
#' @export

adjacency.matrix = function(m, n = NULL)
{
    if (missing(n))
    {
        A = matrix(0, m^2, m^2)
        for (i in 1:m^2)
        {
            up = i - m
            down = i + m
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= m^2)
                A[i, down] = 1
            if (left %% m != 0)
                A[i, left] = 1
            if (i %% m != 0)
                A[i, right] = 1
        }
    }
    else
    {
        A = matrix(0, m * n, m * n)
        for (i in 1:(m * n))
        {
            up = i - n
            down = i + n
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= (m * n))
                A[i, down] = 1
            if (left %% n != 0)
                A[i, left] = 1
            if (i %% n != 0)
                A[i, right] = 1
        }
    }
    A
}

build.M = function(A, d, rho.max, epsilon)
{
    W = A / d
    cG = max(diag(W %*% W) / d)
    k =  ceiling((log(1 - rho.max) + log(epsilon) - log(cG)) / log(rho.max) - 1)
    eig = eigen(W)
    E = eig$vectors
    lam = eig$values
    Einv = solve(E)
    B = E * t(Einv)
    M = suppressWarnings(Re(buildM_(B, k, lam)))
    M
}

neighbor.list = function(A) 
{
	n = nrow(A)
	N = vector("list", n)
    for (i in 1:n)
    {
        temp = which(as.logical(A[i, ]))
        N[[i]] = temp[temp > i]
    }
    N
}

#' Family function for negative binomial GLMs.
#'
#' @description Provides the information required to apply copCAR with negative binomial marginal distributions.
#'
#' @usage negbinomial(theta = stop("'theta' must be specified."), link = "log")
#'
#' @param theta the dispersion parameter (must be positive).
#' @param link the link function, as a character string, name, or one-element character vector, specifying one of \code{log}, \code{sqrt}, or \code{identity}, or an object of class \dQuote{\code{\link[=family]{link-glm}}}
#'
#' @return An object of class \dQuote{\code{family}}, a list of functions and expressions needed to fit a negative binomial GLM.
#'
#' @export

negbinomial = function(theta = stop("'theta' must be specified."), link = "log") 
{
    if (! is.vector(theta, mode = "numeric") || length(theta) > 1)
        stop("You must supply a scalar value for 'theta'.")
    if (theta <= 0)
        stop("'theta' must be positive.")
    linktemp = substitute(link)
    if (! is.character(linktemp)) 
        linktemp = deparse(linktemp)
    if (linktemp %in% c("log", "identity", "sqrt")) 
        stats = make.link(linktemp)
    else if (is.character(link))
    {
        stats = make.link(link)
        linktemp = link
    }
    else
    {
        if (inherits(link, "link-glm"))
        {
            stats = link
            if (! is.null(stats$name)) 
                linktemp = stats$name
        }
        else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"", linktemp))
    }
    variance = function(mu) mu + mu^2 / theta
    validmu = function(mu) all(mu > 0)
    dev.resids = function(y, mu, wt) 2 * wt * (y * log(pmax(1, y) / mu) - (y + theta) * log((y + theta) / (mu + theta)))
    aic = function(y, n, mu, wt, dev)
    {
        term = (y + theta) * log(mu + theta) - y * log(mu) + lgamma(y + 1) - theta * log(theta) + lgamma(theta) - lgamma(theta + y)
        2 * sum(term * wt)
    }
    initialize = expression({
        if (any(y < 0))
            stop("negative values not allowed for the negative binomial family")
        n = rep(1, nobs)
        mustart = y + (y == 0) / 6
    })
    famname = "negbinomial"
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize, theta = theta,
        validmu = validmu, valideta = stats$valideta), 
        class = "family")
}

#' Simulate areal data.
#'
#' @description \code{rcopCAR} simulates areal data from the copCAR model.
#'
#' @details This function simulates data from the copCAR model with the given spatial dependence parameter \eqn{\rho}, regression coefficients \eqn{\beta}, design matrix \eqn{X}, and adjacency structure \eqn{A}. For negative binomial marginal distributions, a value for the dispersion parameter \eqn{\theta>0} is also required; this value must be passed to the \code{\link{negbinomial}} family function. For more details on the copCAR model, see \code{\link{copCAR}}.
#'
#' @param rho the spatial dependence parameter \eqn{\rho} such that \eqn{\rho \in [0, 1)}.
#' @param beta the vector of regression coefficients \eqn{\beta = (\beta_1, \dots, \beta_p)'}.
#' @param X the \eqn{n} by \eqn{p} design matrix \eqn{X}.
#' @param A the symmetric binary adjacency matrix for the underlying graph.
#' @param family the marginal distribution of the observations and link function to be used in the model. This can be a character string naming a family function, a family function, or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) Supported familes are \code{binomial}, \code{poisson}, and \code{negbinomial}.
#'
#' @return A vector of length \eqn{n} distributed according to the specified copCAR model.
#'
#' @export
#'
#' @examples
#'
#' # Use the 20 x 20 square lattice as the underlying graph.
#'
#' m = 20
#' A = adjacency.matrix(m)
#'
#' # Create a design matrix by assigning coordinates to each vertex
#' # such that the coordinates are restricted to the unit square.
#'
#' x = rep(0:(m - 1) / (m - 1), times = m)
#' y = rep(0:(m - 1) / (m - 1), each = m)
#' X = cbind(x, y)
#'
#' # Set the dependence parameter and regression coefficients.
#'
#' rho = 0.995      # strong dependence
#' beta = c(1, 1)   # the mean surface increases in the direction of (1, 1)
#'
#' # Simulate Poisson data from the corresponding copCAR model.
#'
#' z = rcopCAR(rho, beta, X, A, family = poisson(link = "log"))
#'
#' # Simulate Bernoulli outcomes.
#'
#' z = rcopCAR(rho, beta, X, A, family = binomial(link = "logit"))
#'
#' # Set the dispersion parameter.
#'
#' theta = 10
#'
#' # Simulate negative binomial outcomes.
#'
#' z = rcopCAR(rho, beta, X, A, family = negbinomial(theta))

rcopCAR = function(rho, beta, X, A, family)
{
    if (missing(family))
        stop("You must supply a family.")
    if (is.character(family))
        family = get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family = family()
    if (is.null(family$family))
        stop("'family' not recognized.")
    if (! family$family %in% c("binomial", "poisson", "negbinomial"))
        stop("'family' must be binomial, poisson, or negbinomial.")
    if (missing(A) || ! is.matrix(A) || ! isSymmetric(A) || ! all(A == 0 | A == 1))
        stop("You must supply a symmetric binary adjacency matrix.")
    diag(A) = 0
    if (missing(X) || ! is.matrix(X) || ! is.numeric(X))
        stop("You must supply a suitable design matrix.")
    if (nrow(X) != nrow(A))
        stop("The supplied design matrix and adjacency matrix are not conformable.")
    if (missing(rho) || ! is.vector(rho, mode = "numeric") || length(rho) > 1)
        stop("You must supply a scalar value for 'rho'.")
    if (rho <= 0 || rho >= 1)
        stop("'rho' must be between 0 and 1.")
    if (missing(beta) || ! is.vector(beta, mode = "numeric"))
        stop("You must supply a numeric vector for 'beta'.")
    if (ncol(X) != length(beta))
        stop("The given design matrix and vector of regression coefficients are not conformable.")
    Z = simulate.copCAR(rho = rho, beta = beta, X = X, A = A, family = family)
    as.numeric(Z)
}

simulate.copCAR = function(rho, beta, X, A, family)
{
    d = rowSums(A)
    D = diag(d, length(d))
    Q = D - rho * A
	C = chol(as.spam(Q))
    R = chol2inv.spam(C)
    rootR = t(chol(R))
    Y = rootR %*% rnorm(length(d))
    v = diag(R)
    U = pnorm(Y, 0, sqrt(v))
    linkinv = family$linkinv
    eta = X %*% beta
    mu = linkinv(eta)
    if (family$family == "binomial")
        Z = qbinom(U, 1, mu)
    else if (family$family == "poisson")
        Z = qpois(U, mu)
    else
        Z = qnbinom(U, size = family$theta, mu = mu)
    Z
}

copCAR.DT = function(Z, X, A, rho.max, epsilon, offset, family, maxit)
{
    d = rowSums(A)
    M = build.M(A, d, rho.max, epsilon)
    fit.glm = glm(Z ~ X - 1, offset = offset, family = poisson(link = family$link))
    D = diag(d, length(d))
    C = chol(as.spam(D - 0.5 * A))
    init = c(qnorm(0.5), fit.glm$coef)
    lower = c(qnorm(1e-6), rep(-Inf, ncol(X)))
    upper = c(qnorm(0.999), rep(Inf, ncol(X)))
    if (family$family == "negbinomial")
    {
        init = c(init, 1)
        lower = c(lower, 1e-6)
        upper = c(upper, Inf)
    }
    fit = optim(init, lik.DT, Z = Z, X = X, A = A, D = D, d = d, C = C, M = M, offset = offset, family = family,
		        method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE, control = list(maxit = maxit))
    fit
}

copCAR.CE = function(Z, X, A, rho.max, epsilon, offset, family, m, maxit)
{
    d = rowSums(A)
    M = build.M(A, d, rho.max, epsilon)
    U = matrix(runif(m * length(d)), length(d), m)
    fit.glm = glm(Z ~ X - 1, offset = offset, family = poisson(link = family$link))
    D = diag(d, length(d))
    C = chol(as.spam(D - 0.5 * A))
    init = c(qnorm(0.5), fit.glm$coef)
    lower = c(qnorm(1e-6), rep(-Inf, ncol(X)))
    upper = c(qnorm(0.999), rep(Inf, ncol(X)))
    if (family$family == "negbinomial")
    {
        init = c(init, 1)
        lower = c(lower, 1e-6)
        upper = c(upper, Inf)
    }
    fit = optim(init, lik.CE, Z = Z, X = X, A = A, D = D, d = d, C = C, M = M, U = U, offset = offset, family = family,
		        method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE, control = list(maxit = maxit))
    fit
}

copCAR.CML = function(Z, X, A, offset, family, maxit)
{
    d = rowSums(A)
    D = diag(d, length(d))
    C = chol(as.spam(D - 0.5 * A))
    N = neighbor.list(A)
    if (family$family == "binomial")
        fit.glm = glm(Z ~ X - 1, offset = offset, family = binomial(link = family$link))
    else
        fit.glm = glm(Z ~ X - 1, offset = offset, family = poisson(link = family$link))
    init = c(qnorm(0.5), fit.glm$coef)
    lower = c(qnorm(1e-6), rep(-Inf, ncol(X)))
    upper = c(qnorm(0.999), rep(Inf, ncol(X)))
    if (family$family == "negbinomial")
    {
        init = c(init, 1)
        lower = c(lower, 1e-6)
        upper = c(upper, Inf)
    }
	K = matrix(NA, 4, 2)
    fit = optim(init, lik.CML, Z = Z, X = X, A = A, D = D, C = C, offset = offset, family = family, N = N, K = K,
		        method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE, control = list(maxit = maxit))
    fit
}

lik.DT = function(params, Z, X, A, D, d, C, M, offset, family)
{
    rho = pnorm(params[1])
    beta = params[-1]
    if (family$family == "negbinomial")
    {
        theta = beta[length(beta)]
        beta = beta[-length(beta)]
    }
    k = ncol(M) - 1
    v = as.vector(M %*% rho^(0:k) / d)
    Q = as.spam(D - rho * A)
    C = try(update.spam.chol.NgPeyton(C, Q), silent = TRUE)
    if (any(class(C) == "try-error"))
        return(1e6)
    linkinv = family$linkinv
    eta = offset + X %*% beta
    mu = linkinv(eta)
    if (family$family == "poisson")
        U = (ppois(Z, mu) + ppois(Z - 1, mu)) / 2
    else
        U = (pnbinom(Z, mu = mu, size = theta) + pnbinom(Z - 1, mu = mu, size = theta)) / 2
    Y = qnorm(U, 0, sqrt(v))
    qform =  try(t(Y) %*% Q  %*% Y, silent = TRUE)
    if (any(class(qform) == "try-error"))
        return(1e6)
    if (family$family == "poisson")
        loglik = -sum(log(diag(C))) + 0.5 * (qform - sum(Y^2 / v) - sum(log(v))) - sum(dpois(Z, mu, log = TRUE))
    else
        loglik = -sum(log(diag(C))) + 0.5 * (qform - sum(Y^2 / v) - sum(log(v))) - sum(dnbinom(Z, mu = mu, size = theta, log = TRUE))
    loglik
}

lik.CE = function(params, Z, X, A, D, d, C, M, U, offset, family)
{
    rho = pnorm(params[1])
    beta = params[-1]
    if (family$family == "negbinomial")
    {
        theta = beta[length(beta)]
        beta = beta[-length(beta)]
    }
    k = ncol(M) - 1
    v = as.vector(M %*% rho^(0:k) / d)
    Q = as.spam(D - rho * A)
    C = try(update.spam.chol.NgPeyton(C, Q), silent = TRUE)
    if (any(class(C) == "try-error"))
        return(1e6)
    linkinv = family$linkinv
    eta = offset + X %*% beta
    mu = linkinv(eta)
    m = ncol(U)
    w = numeric(m)
    for (j in 1:m)
    {
        z_ = Z - U[, j]
        if (family$family == "poisson")
            u_ = ppois.CE(z_, mu)
        else
            u_ = pnbinom.CE(z_, mu, theta)
        y_ = qnorm(u_, 0, sqrt(v))
        qform = try(t(y_) %*% Q %*% y_, silent = TRUE)
        if (any(class(qform) == "try-error"))
            return(1e6)
        w[j] = -0.5 * (qform - sum(y_^2 / v))
    }
    s = max(w)
    if (family$family == "poisson")
        loglik = -log(sum(exp(w - s))) - s - sum(dpois(Z, mu, log = TRUE)) - sum(log(diag(C))) - 0.5 * sum(log(v))
    else
        loglik = -log(sum(exp(w - s))) - s - sum(dnbinom(Z, mu = mu, size = theta, log = TRUE)) - sum(log(diag(C))) - 0.5 * sum(log(v))
    if (loglik == Inf)
       return(1e6)
    loglik
}

lik.CML = function(params, Z, X, A, D, C, offset, family, N, K)
{
    rho = pnorm(params[1])
    beta = params[-1]
    if (family$family == "negbinomial")
    {
        theta = beta[length(beta)]
        beta = beta[-length(beta)]
    }
    linkinv = family$linkinv
    eta = offset + X %*% beta
    mu = linkinv(eta)
    Q = as.spam(D - rho * A)
    C = try(update.spam.chol.NgPeyton(C, Q), silent = TRUE)
    if (any(class(C) == "try-error"))
        return(1e6)
    R = chol2inv.spam(C)
    v = diag(R)
    if (family$family == "binomial")
    {
        Y.0 = qnorm(pbinom(Z, 1, mu), 0, sqrt(v)) / sqrt(v)
        Y.1 = qnorm(pbinom(Z - 1, 1, mu), 0, sqrt(v)) / sqrt(v)
    }
    else if (family$family == "poisson")
    {
        Y.0 = qnorm(ppois(Z, mu), 0, sqrt(v)) / sqrt(v)
        Y.1 = qnorm(ppois(Z - 1, mu), 0, sqrt(v)) / sqrt(v)
    }
    else
    {
        Y.0 = qnorm(pnbinom(Z, mu = mu, size = theta), 0, sqrt(v)) / sqrt(v)
        Y.1 = qnorm(pnbinom(Z - 1, mu = mu, size = theta), 0, sqrt(v)) / sqrt(v)
    }
    Y.0 = ifelse(Y.0 == -Inf, -1e6, Y.0)
    Y.1 = ifelse(Y.1 == -Inf, -1e6, Y.1)
    Y.0 = ifelse(Y.0 == Inf, 1e6, Y.0)
    Y.1 = ifelse(Y.1 == Inf, 1e6, Y.1)
    loglik = 0
	u = c(1, -1, -1, 1)
    for (i in 1:nrow(A))
    {
		for (j in N[[i]])
		{
            corr = R[i, j] / sqrt(v[i] * v[j])
	    	K[, 1] = c(Y.0[i], Y.0[i], Y.1[i], Y.1[i])
	    	K[, 2] = c(Y.0[j], Y.1[j], Y.0[j], Y.1[j])
            s = pbivnormf(K, rho = corr)
    	    s = t(s) %*% u
    	    if (s != 0)
    	        loglik = loglik + log(s)
	    }
    }
    -loglik
}

dpois.CE = function(x, mu)
{
    dpois(floor(x + 1), mu)
}

ppois.CE = function(q, mu)
{
    ppois(floor(q), mu) + (q - floor(q)) * dpois.CE(q, mu)
}

dnbinom.CE = function(x, mu, theta)
{
    dnbinom(floor(x + 1), mu = mu, size = theta)
}

pnbinom.CE = function(q, mu, theta)
{
    pnbinom(floor(q), mu = mu, size = theta) + (q - floor(q)) * dnbinom.CE(q, mu, theta)
}

pbivnormf = function (K, rho = 0)
{
    correl = rep(rho, nrow(K))
    lower = as.double(c(0, 0))
    infin = as.integer(c(0, 0))
    uppera = as.double(K[, 1])
    upperb = as.double(K[, 2])
    lt = as.integer(nrow(K))
    prob = double(lt)
    correl = as.double(correl)
    ans = .Fortran("PBIVNORM", prob, lower, uppera, upperb, infin, correl, lt, PACKAGE = "copCAR")[[1]]
    ans
}

copCAR.bootstrap.helper = function(dummy, mu, theta, v, family, rootR, method, control, X, A, offset)
{
	n = length(mu)
	Y = rootR %*% rnorm(n)
	U = pnorm(Y, sd = sqrt(v))
    if (family$family == "binomial")
        Z = qbinom(U, 1, mu)
    else if (family$family == "poisson")
        Z = qpois(U, mu)
    else
        Z = qnbinom(U, size = theta, mu = mu)
    if (method == "CE")
        fit = copCAR.CE(Z = Z, X = X, A = A, rho.max = control$rho.max, epsilon = control$epsilon, offset = offset, family = family, m = control$m, maxit = control$maxit)
    else if (method == "DT")
        fit = copCAR.DT(Z = Z, X = X, A = A, rho.max = control$rho.max, epsilon = control$epsilon, offset = offset, family = family, maxit = control$maxit)
    else
        fit = copCAR.CML(Z = Z, X = X, A = A, offset = offset, family = family, maxit = control$maxit)
	fit
}

copCAR.bootstrap = function(params, X, A, offset, family, method, control, verbose)
{
    rho = pnorm(params[1])
    beta = params[-1]
    if (family$family == "negbinomial")
    {
        theta = beta[length(beta)]
        beta = beta[-length(beta)]
    }
	else
		theta = NULL
    linkinv = family$linkinv
    eta = offset + X %*% beta
    mu = linkinv(eta)
    d = rowSums(A)
    D = diag(d, length(d))
    Q = D - rho * A
	C = chol(as.spam(Q))
	R = chol2inv.spam(C)
    rootR = t(chol(R))
    v = diag(R)
    boot.sample = data.frame(matrix(, control$bootit, length(params)))
    if (! control$parallel)
    {
		if (verbose && requireNamespace("pbapply", quietly = TRUE))
		{
			cat("\n")
			gathered = pbapply::pblapply(1:control$bootit, copCAR.bootstrap.helper, mu, theta, v, family, rootR, method, control, X, A, offset)
		}
		else
		{
			gathered = vector("list", control$bootit)
            for (j in 1:control$bootit)
                gathered[[j]] = copCAR.bootstrap.helper(NULL, mu, theta, v, family, rootR, method, control, X, A, offset)
        }
    }
    else
    {
        cl = parallel::makeCluster(control$nodes, control$type)
        parallel::clusterEvalQ(cl, library(copCAR))
        if (verbose && requireNamespace("pbapply", quietly = TRUE))
		{
		    cat("\n")
			gathered = pbapply::pblapply(1:control$bootit, copCAR.bootstrap.helper, mu, theta, v, family, rootR, method, control, X, A, offset, cl = cl)
		}
		else
            gathered = parallel::clusterApplyLB(cl, 1:control$bootit, copCAR.bootstrap.helper, mu, theta, v, family, rootR, method, control, X, A, offset)
        parallel::stopCluster(cl)      
    }
    for (j in 1:control$bootit)
    {
        fit = gathered[[j]]
        temp = rep(NA, length(params))
        if (fit$convergence != 0)
            warning(fit$message)
        else
            temp = fit$par
        boot.sample[j, ] = temp
	}
    boot.sample
}

copCAR.sandwich.helper = function(dummy, params, mu, theta, M, v, family, C, D, d, N, rootR, method, control, X, A, offset)
{
	n = length(mu)
	Y = rootR %*% rnorm(n)
	U = pnorm(Y, sd = sqrt(v))
    if (family$family == "binomial")
        Z = qbinom(U, 1, mu)
    else if (family$family == "poisson")
        Z = qpois(U, mu)
    else
        Z = qnbinom(U, size = theta, mu = mu)
    if (method == "CML")
	{
		K = matrix(NA, 4, 2)
        gr = - grad(lik.CML, params, Z = Z, X = X, A = A, D = D, C = C, offset = offset, family = family, N = N, K = K)
	}
    else
        gr = - grad(lik.DT, params, Z = Z, X = X, A = A, D = D, d = d, C = C, M = M, offset = offset, family = family)
	gr
}

copCAR.asymptotic = function(params, X, A, offset, family, method, bread, control, verbose)
{
    rho = pnorm(params[1])
    beta = params[-1]
    if (family$family == "negbinomial")
    {
        theta = beta[length(beta)]
        beta = beta[-length(beta)]
    }
	else
		theta = NULL
    linkinv = family$linkinv
    eta = offset + X %*% beta
    mu = linkinv(eta)
    d = rowSums(A)
    D = diag(d, length(d))
    Q = D - rho * A
	C = chol(as.spam(Q))
	R = chol2inv.spam(C)
    rootR = t(chol(R))
    v = diag(R)
	M = NULL
	if (method == "DT")
		M = build.M(A, d, control$rho.max, control$epsilon)
	N = NULL
	if (method == "CML")
		N = neighbor.list(A)
    if (! control$parallel)
    {
		if (verbose && requireNamespace("pbapply", quietly = TRUE))
		{
			cat("\n")
			gathered = pbapply::pblapply(1:control$bootit, copCAR.sandwich.helper, params, mu, theta, M, v, family, C, D, d, N, rootR, method, control, X, A, offset)
		}
		else
		{
			gathered = vector("list", control$bootit)
            for (j in 1:control$bootit)
                gathered[[j]] = copCAR.sandwich.helper(NULL, params, mu, theta, M, v, family, C, D, d, N, rootR, method, control, X, A, offset)
        }
    }
    else
    {
        cl = parallel::makeCluster(control$nodes, control$type)
        parallel::clusterEvalQ(cl, library(copCAR))
        if (verbose && requireNamespace("pbapply", quietly = TRUE))
		{
		    cat("\n")
			gathered = pbapply::pblapply(1:control$bootit, copCAR.sandwich.helper, params, mu, theta, M, v, family, C, D, d, N, rootR, method, control, X, A, offset, cl = cl)
		}
		else
            gathered = parallel::clusterApplyLB(cl, 1:control$bootit, copCAR.sandwich.helper, params, mu, theta, M, v, family, C, D, d, N, rootR, method, control, X, A, offset)
        parallel::stopCluster(cl)      
    }
	meat = matrix(0, length(params), length(params))
    for (j in 1:control$bootit)
    {
		gr = gathered[[j]]
        meat = meat + gr %o% gr / control$bootit
	}   
    cov.hat = bread %*% meat %*% bread
	cov.hat
}

copCAR.control = function(method, confint, control, verbose)
{
	nms = match.arg(names(control), c("maxit", "epsilon", "rho.max", "m", "bootit", "parallel", "type", "nodes"), several.ok = TRUE)
	control = control[nms]
	maxit = control$maxit
    if (is.null(maxit) || ! is.numeric(maxit) || length(maxit) > 1 || maxit != as.integer(maxit) || maxit < 100)
    {
        if (verbose)
            cat("\nControl parameter 'maxit' must be a positive integer >= 100. Setting it to the default value of 1,000.\n")
        control$maxit = 1000
    }
	if (method != "CML")
	{
		epsilon = control$epsilon
        if (is.null(epsilon) || ! is.numeric(epsilon) || length(epsilon) > 1 || epsilon <= 0 || epsilon > 0.1)
        {
            if (verbose)
                cat("\nControl parameter 'epsilon' must be a small positive number. Setting it to the default value of 0.01.\n")
            control$epsilon = 0.01
        }
		rho.max = control$rho.max
        if (is.null(rho.max) || ! is.numeric(rho.max) || length(rho.max) > 1 || rho.max <= 0 || rho.max >= 1)
        {
            if (verbose)
                cat("\nControl parameter 'rho.max' must be from (0, 1), and should be close to 1. Setting it to the default value of 0.999.\n")
            control$rho.max = 0.999
        }
	}
	else
		control$epsilon = control$rho.max = NULL
	if (method == "CE")
	{
		m = control$m
        if (is.null(m) || ! is.numeric(m) || length(m) > 1 || m != as.integer(m) || m < 10)
        {
            if (verbose)
                cat("\nControl parameter 'm' must be a positive integer >= 10. Setting it to the default value of 1,000.\n")
            control$m = 1000
        }
	}
	else
		control$m = NULL
	if (confint != "none")
	{
		if (confint == "bootstrap" || (confint == "asymptotic" && method != "CE"))
		{
            bootit = control$bootit
            if (is.null(bootit) || ! is.numeric(bootit) || length(bootit) > 1 || bootit != as.integer(bootit) || bootit < 100)
            {
                if (verbose)
                    cat("\nControl parameter 'bootit' must be a positive integer >= 100. Setting it to the default value of 1,000.\n")
                control$bootit = 1000
            }
            if (is.null(control$parallel) || ! is.logical(control$parallel) || length(control$parallel) > 1)
            {
                if (verbose)
                    cat("\nControl parameter 'parallel' must be a logical value. Setting it to the default value of TRUE.\n")
                control$parallel = TRUE
            }
            if (control$parallel)
            {
                if (requireNamespace("parallel", quietly = TRUE))
                {
                    if (is.null(control$type) || length(control$type) > 1 || ! control$type %in% c("SOCK", "PVM", "MPI", "NWS"))
                    {
                        if (verbose)
                            cat("\nControl parameter 'type' must be \"SOCK\", \"PVM\", \"MPI\", or \"NWS\". Setting it to \"SOCK\".\n")
                        control$type = "SOCK"
                    }
                    nodes = control$nodes
                    if (is.null(control$nodes) || ! is.numeric(nodes) || length(nodes) > 1 || nodes != as.integer(nodes) || nodes < 2)
                        stop("Control parameter 'nodes' must be a whole number greater than 1.")
                }
                else
                {
                    if (verbose)
                        cat("\nParallel computation requires package parallel. Setting control parameter 'parallel' to FALSE.\n")
                    control$parallel = FALSE
				    control$type = control$nodes = NULL 
                }
            }
	    }
		else
			control$bootit = control$parallel = control$type = control$nodes = NULL
    }
	else
		control$bootit = control$parallel = control$type = control$nodes = NULL
    control
}

#' Fit copCAR model to discrete areal data.
#'
#' @description Fit the copCAR model to Bernoulli, negative binomial, or Poisson observations.
#'
#' @details This function performs frequentist inference for the copCAR model proposed by Hughes (2015). copCAR is a copula-based areal regression model that employs the proper conditional autoregression (CAR) introduced by Besag, York, and Mollié (1991). Specifically, copCAR uses the CAR copula, a Caussian copula based on the proper CAR.
#'
#'
#' The spatial dependence parameter \eqn{\rho \in [0, 1)}, regression coefficients \eqn{\beta = (\beta_1, \dots, \beta_p)' \in R^p}, and, for negative binomial margins, dispersion parameter \eqn{\theta>0} can be estimated using the continous extension (CE) (Madsen, 2009), distributional transform (DT) (Kazianka and Pilz, 2010), or composite marginal likelihood (CML) (Varin, 2008) approaches.
#'
#'
#' The CE approach transforms the discrete observations to continous outcomes by convolving them with independent standard uniforms (Denuit and Lambert, 2005). The true likelihood for the discrete outcomes is the expected likelihood for the transformed outcomes. An estimate (sample mean) of the expected likelihood is optimized to estimate the copCAR parameters. The number of standard uniform vectors, \eqn{m}, can be chosen by the user. The default value is 1,000. The CE approach is exact up to Monte Carlo standard error but is computationally intensive (the computational burden grows rapidly with increasing \eqn{m}). The CE approach tends to perform poorly when applied to Bernoulli outcomes, and so that option is not permitted. 
#'
#'
#' The distributional transform stochastically "smoothes" the jumps of a discrete distribution function (Ferguson, 1967). The DT-based approximation (Kazianka and Pilz, 2010) for copCAR performs well for Poisson and negative binomial marginals but, like the CE approach, tends to perform poorly for Bernoulli outcomes.
#'
#'
#' The CML approach optimizes a composite marginal likelihood formed as the product of pairwise likelihoods of adjacent observations. This approach performs well for Bernoulli, negative binomial, and Poisson outcomes.
#'
#'
#' In the CE and DT approaches, the CAR variances are approximated. The quality of the approximation is determined by the values of control parameters \eqn{\epsilon > 0} and \eqn{\rho^{\max} \in [0, 1)}. The default values are 0.01 and 0.999, respectively.
#'
#'
#' When \code{confint = "bootstrap"}, a parametric bootstrap is carried out, and confidence intervals are computed using the quantile method. Monte Carlo standard errors (Flegal et al., 2008) of the quantile estimators are also provided.
#'
#'
#' When \code{confint = "asymptotic"}, confidence intervals are computed using an estimate of the asymptotic covariance matrix of the estimator. For the CE method, the inverse of the observed Fisher information matrix is used. For the CML and DT methods, the objective function is misspecified, and so the asymptotic covariance matrix is the inverse of the Godambe information matrix (Godambe, 1960), which has a sandwich form. The "bread" is the inverse of the Fisher information matrix, and the "meat" is the covariance matrix of the score function. The former is estimated using the inverse of the observed Fisher information matrix. The latter is estimated using a parametric bootstrap.
#'
#'
#' @param formula an object of class \dQuote{\code{\link{formula}}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of the model specification are given under "Details".
#' @param family the marginal distribution of the observations at the areal units and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) Supported families are \code{binomial}, \code{negbinomial}, and \code{poisson}. When the negative binomial family is used, an initial value for \eqn{\theta} must be passed to the \code{\link{negbinomial}} family function.
#' @param data an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{copCAR} is called.
#' @param offset this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of observations. One or more \code{\link{offset}} terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param A the symmetric binary adjacency matrix for the underlying graph.
#' @param method the method for inference. \code{copCAR} supports the continous extension (\dQuote{\code{CE}}), distributional transform (\dQuote{\code{DT}}), and composite marginal likelihood (\dQuote{\code{CML}}).
#' @param confint the method for computing confidence intervals. This defaults to \dQuote{\code{none}}. The other options are \dQuote{\code{bootstrap}} (for parametric bootstrap intervals using the quantile method) and \dQuote{\code{asymptotic}} (for intervals computed using an estimate of the asymptotic covariance matrix).
#' @param model a logical value indicating whether the model frame should be included as a component of the returned value.
#' @param x a logical value indicating whether the model matrix used in the fitting process should be returned as a component of the returned value.
#' @param y a logical value indicating whether the response vector used in the fitting process should be returned as a component of the returned value.
#' @param verbose a logical value indicating whether to print various messages to the screen, including progress updates. Defaults to \code{FALSE}.
#' @param control a list of parameters for controlling the fitting process.
#' \describe{
#'        \item{\code{bootit}}{the size of the (parametric) bootstrap sample. This applies when \code{confint = "bootstrap"}, or when \code{confint = "asymptotic"} and \code{method = "CML"} or \code{method = "DT"}. Defaults to \code{500}.}
#'        \item{\code{m}}{the number of vectors of standard uniforms used to approximate the expected likelhood when \code{method = "CE"}. Defaults to \code{1000}.}
#'        \item{\code{rho.max}}{the value \eqn{\rho^{\max}}, which is the maximum value of \eqn{\rho} used to approximate the CAR variances when \code{method = "CE"} or \code{method = "DT"}. If missing, assumed to be \code{0.999}.}
#'        \item{\code{epsilon}}{the tolerance \eqn{\epsilon > 0} used to approximate the CAR variances when \code{method = "CE"} or \code{method = "DT"}. If missing, assumed to be \code{0.01}.}
#'        item{\code{maxit}}{the maximum number of iterations to be used by \code{\link{optim}} when optimizing the objective function. Defaults to \code{1000}.}
#'        \item{\code{parallel}}{a logical value indicating whether to parallelize the bootstrap. This defaults to \code{TRUE} if the \pkg{parallel} package can be loaded.}
#'        \item{\code{type}}{the cluster type, one of \dQuote{\code{SOCK}} (default), \dQuote{\code{PVM}}, \dQuote{\code{MPI}}, or \dQuote{\code{NWS}}.}
#'        \item{\code{nodes}}{the number of slave nodes to create.}}
#'
#' @return \code{copCAR} returns an object of class \code{"copCAR"}, which is a list containing the following components:
#'         \item{boot.sample}{(if \code{confint = "bootstrap"}) the bootstrap sample.}
#'         \item{call}{the matched call.}
#'         \item{coefficients}{a named vector of parameter estimates.}
#'         \item{confint}{the value of \code{confint} supplied in the function call.}
#'         \item{control}{a list containing the names and values of the control parameters.}
#'         \item{convergence}{the integer code returned by \code{\link{optim}}.}
#'         \item{cov.hat}{(if \code{confint = "asymptotic"}) the estimate of the asymptotic covariance matrix of the parameter estimator.}
#'         \item{data}{the \code{data} argument.}
#'         \item{family}{the \code{\link{family}} object used.}
#'         \item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
#'         \item{formula}{the formula supplied.}
#'         \item{linear.predictors}{the linear fit on link scale.}
#'         \item{message}{A character string giving any additional information returned by the optimizer, or \code{NULL}.}
#'         \item{method}{the method (CE, CML, or DT) used for inference.}
#'         \item{model}{if requested (the default), the model frame.}
#'         \item{npar}{the number of model parameters.}
#'         \item{offset}{the offset vector used.}
#'         \item{residuals}{the response residuals, i.e., the outcomes minus the fitted values.}
#'         \item{terms}{the \code{\link{terms}} object used.}
#'         \item{value}{the value of the objective function at its minimum.}
#'         \item{x}{if requested, the model matrix.}
#'         \item{xlevels}{(where relevant) a record of the levels of the factors used in fitting.}
#'         \item{y}{if requested (the default), the response vector used.}
#'
#' @references
#' Besag, J., York, J., and Mollié, A. (1991) Bayesian image restoration, with two applications in spatial statistics. \emph{Annals of the Institute of Statistical Mathematics}, \bold{43}(1), 1--20.
#' @references
#' Denuit, M. and Lambert, P. (2005) Constraints on concordance measures in bivariate discrete data. \emph{Journal of Multivariate Analysis}, \bold{93}, 40--57.
#' @references
#' Ferguson, T. (1967) \emph{Mathematical statistics: a decision theoretic approach}, New York: Academic Press.
#' @references
#' Flegal, J., Haran, M., and Jones, G. (2008) Markov Chain Monte Carlo: can we trust the third significant figure? \emph{Statistical Science}, 23(2), 250--260.
#' @references
#' Godambe, V. (1960) An optimum property of regular maximum likelihood estimation. \emph{The Annals of Mathmatical Statistics}, \bold{31}(4), 1208--1211.
#' @references
#' Hughes, J. (2015) copCAR: A flexible regression model for areal data.  \emph{Journal of Computational and Graphical Statistics}, \bold{24}(3), 733--755.
#' @references
#' Kazianka, H. and Pilz, J. (2010) Copula-based geostatistical modeling of continuous and discrete data including covariates. \emph{Stochastic Environmental Research and Risk Assessment}, \bold{24}(5), 661--673.
#' @references
#' Madsen, L. (2009) Maximum likelihood estimation of regression parameters with spatially dependent discrete data. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, \bold{14}(4), 375--391.
#' @references
#' Varin, C. (2008) On composite marginal likelihoods. \emph{Advances in Statistical Analysis}, \bold{92}(1), 1--28.
#'
#' @export
#'
#' @examples \dontrun{
#' # Simulate data and fit copCAR to them.
#'
#' # Use the 20 x 20 square lattice as the underlying graph.
#'
#' m = 20
#' A = adjacency.matrix(m)
#'
#' # Create a design matrix by assigning coordinates to each vertex
#' # such that the coordinates are restricted to the unit square.
#'
#' x = rep(0:(m - 1) / (m - 1), times = m)
#' y = rep(0:(m - 1) / (m - 1), each = m)
#' X = cbind(x, y)
#'
#' # Set the dependence parameter, regression coefficients, and dispersion parameter.
#'
#' rho = 0.995      # strong dependence
#' beta = c(1, 1)   # the mean surface increases in the direction of (1, 1)
#' theta = 2        # dispersion parameter
#'
#' # Simulate negative binomial data from the model.
#'
#' z = rcopCAR(rho, beta, X, A, family = negbinomial(theta))
#'
#' # Fit the copCAR model using the continous extension, and compute 95% (default)
#' # asymptotic confidence intervals. Give theta the initial value of 1. Use m equal to 100.
#'
#' fit.ce = copCAR(z ~ X - 1, A = A, family = negbinomial(1), method = "CE", confint = "asymptotic",
#'                 control = list(m = 100))
#' summary(fit.ce)
#'
#' # Fit the copCAR model using the DT approximation, and compute 90% confidence
#' # intervals. Bootstrap the intervals, based on a bootstrap sample of size 100.
#' # Do the bootstrap in parallel, using ten nodes.
#'
#' fit.dt = copCAR(z ~ X - 1, A = A, family = negbinomial(1), method = "DT", confint = "bootstrap",
#'                 control = list(bootit = 100, nodes = 10))
#' summary(fit.dt, alpha = 0.9)
#'
#' # Fit the copCAR model using the composite marginal likelihood approach.
#' # Do not compute confidence intervals.
#'
#' fit.cml = copCAR(z ~ X - 1, A = A, family = negbinomial(1), method = "CML", confint = "none")
#' summary(fit.cml)
#' }

copCAR = function(formula, family, data, offset, A, method = c("CML", "DT", "CE"), confint = c("none", "bootstrap", "asymptotic"),
                  model = TRUE, x = FALSE, y = TRUE, verbose = FALSE, control = list())
{
    call = match.call()
    if (missing(formula))
        stop("You must supply a formula.")
    if (missing(A) || ! is.matrix(A) || ! isSymmetric(A) || ! all(A == 0 | A == 1))
        stop("You must supply a symmetric binary adjacency matrix.")
	diag(A) = 0
    if (is.character(family))
        family = get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family = family()
    if (is.null(family$family))
        stop("'family' not recognized.")
    if (! family$family %in% c("binomial", "poisson", "negbinomial"))
        stop("'family' must be binomial, poisson, or negbinomial.")
    if (missing(data))
        data = environment(formula)
    mf = match.call(expand.dots = FALSE)
    m = match(c("formula", "data", "offset"), names(mf), 0L)
    mf = mf[c(1L, m)]
    mf[[1L]] = quote(stats::model.frame)
    mf = eval(mf, parent.frame())
    mt = attr(mf, "terms")
    Z = model.response(mf, "numeric")
	if (! is.vector(Z))
		stop("The response must be a vector.")
    if (family$family == "binomial" && ! (Z == 0 || Z == 1))
        stop("The response must be binary if 'family' is binomial.")
	if (family$family == "binomial" && method != "CML")
	    stop("Only the CML method is appropriate for binary outcomes.")
	if (family$family %in% c("poisson", "negbinomial") && ! all(Z >= 0 & (Z == as.integer(Z))))
		stop("The outcomes must be non-negative integers if 'family' is poisson or negbinomial.")
    X = model.matrix(mt, mf)
    offset = as.vector(model.offset(mf))
    if (! is.null(offset) && length(offset) != length(Z))
        stop("The number of offsets should equal the number of observations.")
    if (is.null(offset))
        offset = rep(0, length(Z))
    if (sum(c(nrow(X), nrow(A)) != length(Z)) > 0)
        stop("The supplied response vector/design matrix/adjacency matrix are not conformable.")
    method = match.arg(method)
	confint = match.arg(confint)
    if (! is.logical(verbose) || length(verbose) > 1)
        stop("'verbose' must be a logical value.")
    if (! is.list(control))
        stop("'control' must be a list.")
    control = copCAR.control(method, confint, control, verbose)
    if (method == "CE")
        fit = copCAR.CE(Z = Z, X = X, A = A, rho.max = control$rho.max, epsilon = control$epsilon, offset = offset, family = family, m = control$m, maxit = control$maxit)
    else if (method == "DT")
        fit = copCAR.DT(Z = Z, X = X, A = A, rho.max = control$rho.max, epsilon = control$epsilon, offset = offset, family = family, maxit = control$maxit)
    else
        fit = copCAR.CML(Z = Z, X = X, A = A, offset = offset, family = family, maxit = control$maxit)
    result = list()
    class(result) = "copCAR"
    if (fit$convergence != 0)
    {
        warning("Optimization failed.")
        result$convergence = fit$convergence
        result$message = fit$message
		result$call = call
        return(result)
    }
    npar = length(fit$par)
    linkinv = family$linkinv
    beta.hat = fit$par[-1]
    if (family$family == "negbinomial")
        beta.hat = beta.hat[-length(beta.hat)]
    eta = X %*% beta.hat
    mu = linkinv(eta)
	if (confint != "none")
	{
		if (confint == "bootstrap")
		{
			boot.sample = copCAR.bootstrap(fit$par, X, A, offset, family, method, control, verbose)
			colnames(boot.sample) = c("gamma", names(fit$par[-1]))
		    result$boot.sample = boot.sample
		}
		else
		{
			cov.hat = try(solve(fit$hessian), silent = TRUE)
			if (any(class(cov.hat) == "try-error"))
				result$cov.hat = NULL
			else if (method == "CE")
			    result$cov.hat = cov.hat
			else
			{
				cov.hat = copCAR.asymptotic(fit$par, X, A, offset, family, method, cov.hat, control, verbose)
				result$cov.hat = cov.hat
			}
		}
	}
	result$npar = npar
    result$linear.predictors = eta
    result$fitted.values = mu
	result$residuals = Z - mu
    result$coefficients = c(pnorm(fit$par[1]), fit$par[-1])
    result$confint = confint
	result$family = family
    if (model)
        result$model = mf
    if (x)
        result$x = X
    if (y)
        result$y = Z
    result$convergence = fit$convergence
    result$message = fit$message
    result$value = fit$value
    result$data = data
    result$method = method
    result$xlevels = .getXlevels(mt, mf)
    result$call = call
    result$terms = mt
    result$formula = formula
	result$offset = offset
	result$control = control
    temp = c("rho", colnames(X))
    if (family$family == "negbinomial")
        temp = c(temp, "theta")
    names(result$coefficients) = temp
    result
}

#' Print a summary of a copCAR model fit.
#'
#' @details This function displays (1) the call to \code{\link{copCAR}}, (2) the values of the control parameters, (3) a table of estimates, and (when applicable) (4) confidence intervals and (5) Monte Carlo standard errors.
#'
#' Each row of the table of estimates shows a parameter estimate and (when applicable) the confidence interval for the parameter. If \code{\link{copCAR}} was called with \code{confint = "bootstrap"}, Monte Carlo standard errors are provided.
#'
#' @param object an object of class \code{copCAR}, the result of a call to \code{\link{copCAR}}.
#' @param alpha the significance level for the confidence intervals. The default is 0.05.
#' @param digits the number of significant digits to display. The default is 4.
#' @param \dots additional arguments.
#'
#' @seealso \code{\link{copCAR}}
#'
#' @method summary copCAR
#'
#' @references
#' Flegal, J., Haran, M., and Jones, G. (2008) Markov Chain Monte Carlo: can we trust the third significant figure? \emph{Statistical Science}, \bold{23}(2), 250--260.
#'
#' @export

summary.copCAR = function(object, alpha = 0.05, digits = 4, ...)
{
    cat("\nCall:\n\n")
	print(object$call)
    cat("\nConvergence:\n")
    if (object$convergence != 0)
    {
    	cat("\nOptimization failed =>", object$message, "\n")
		return(invisible())
    }
	else
		cat("\nOptimization converged at", signif(-object$value, digits = 4), "\n")
    cat("\nControl parameters:\n")
    if (length(object$control) > 0)
    {
        control.table = cbind(unlist(object$control))
        colnames(control.table) = ""
        print(control.table, quote = FALSE)
    }
    else
        print("None specified." )
    npar = object$npar
    if (object$confint == "none")
        confint = matrix(rep(NA, 2 * npar), ncol = 2)
	else
	{
	    boot.sample = object$boot.sample
	    if (! is.null(boot.sample))
	    {
	        boot.sample[, 1] = pnorm(boot.sample[, 1])
			boot.sample = boot.sample[complete.cases(boot.sample), ]
	        lo = mcse.q.mat(boot.sample, alpha / 2)
			hi = mcse.q.mat(boot.sample, 1 - alpha / 2)
			confint = cbind(lo, hi)
	    }
		else
		{
			cov.hat = object$cov.hat
			if (! is.null(cov.hat))
			{
			    se = sqrt(diag(cov.hat))
			    scale = qnorm(1 - alpha / 2)
			    coef = object$coef
			    coef[1] = qnorm(coef[1])
			    confint = cbind(coef - scale * se, coef + scale * se)
			    confint[1, ] = pnorm(confint[1, ])
		    }
			else
				confint = matrix(rep(NA, 2 * npar), ncol = 2)
		}
	}
    coef.table = cbind(object$coef, confint)
	if (is.null(object$boot.sample))
        colnames(coef.table) = c("Estimate", "Lower", "Upper")
	else
		colnames(coef.table) = c("Estimate", "Lower", "MCSE (Lower)", "Upper", "MCSE (Upper)")
    rownames(coef.table) = names(object$coef)
    cat("\nCoefficients:\n\n")
	print(signif(coef.table, digits = 4))
	cat("\n")
}

#' Return the estimated covariance matrix for a \code{copCAR} model object.
#'
#' @details Unless \code{\link{copCAR}} was called with \code{confint = "none"}, this function returns an estimate of the covariance matrix of the CE/CML/DT estimator of the parameters. If \code{confint = "bootstrap"}, \code{\link{cov}} is applied to the bootstrap sample to compute the estimate. If \code{confint = "asymptotic"}, an estimate of the asymptotic covariance matrix is returned; this is an estimate of the inverse Fisher information matrix if \code{method = "CE"}, or an estimate of the inverse of the Godambe information matrix if \code{method = "CML"} or \code{method = "DT"}. Note that the entries involving the spatial dependence parameter are for \eqn{\gamma=\Phi^{-1}(\rho)} rather than for \eqn{\rho} (Hughes, 2015).
#'
#' @param object a fitted \code{copCAR} model object.
#' @param \dots additional arguments.
#'
#' @return An estimate of the covariance matrix of the CE/CML/DT estimator of the parameters.
#'
#' @method vcov copCAR
#'
#' @references
#' Hughes, J. (2015) copCAR: A flexible regression model for areal data.  \emph{Journal of Computational and Graphical Statistics}, \bold{24}(3), 733--755.
#'
#' @export

vcov.copCAR = function(object, ...)
{
    if (object$confint == "none")
        stop("Function copCAR was called with 'confint' equal to \"none\".")
    if (! is.null(object$cov.hat))
        cov.hat = object$cov.hat
    else
        cov.hat = cov(object$boot.sample, use = "complete.obs")
	rownames(cov.hat) = colnames(cov.hat) = c("gamma", names(object$coef[-1]))
    cov.hat
}

#' Extract model residuals.
#'
#' @param object an object of class \code{copCAR}, typically the result of a call to \code{\link{copCAR}}.
#' @param type the type of residuals that should be returned. The alternatives are \dQuote{\code{deviance}} (default), \dQuote{\code{pearson}}, and \dQuote{\code{response}}.
#' @param \dots additional arguments.
#'
#' @return A vector of residuals.
#'
#' @seealso \code{\link{copCAR}}, \code{\link{residuals.glm}}
#'
#' @method residuals copCAR
#'
#' @export

residuals.copCAR = function(object, type = c("deviance", "pearson", "response"), ...)
{
    type = match.arg(type)
    y = object$y
    r = object$residuals
    mu = object$fitted.values
	if (type == "response")
		return(r)
	else if (type == "pearson")
		return(r / sqrt(object$family$variance(mu)))
	else
	{
        if (is.null(y))
            y = mu + r
        return(sqrt(pmax((object$family$dev.resids)(y, mu, rep(1, length(y))), 0)))
    }
}

