## Matt Galloway


#' @title Normal regression data generator
#'
#' @description True beta values are generated from p*r independent draws from N(0, 1/p) distribution. X are n independent draws from p multivariate normal N(0, SigmaX). Y is then generated using X and true beta values with an iid error term that follows r multivariate normal distribution N(0, Sigma).
#'
#' @param n desired sample size
#' @param p desired dimension
#' @param r number of responses
#' @param sparsity desired sparsity for beta
#' @param Sigma covariance matrix structure used to generate Y | X
#' @param s option to specify diagonal elements in Sigma
#' @param SigmaX covariance matrix structure used to generate data X
#' @param sx option to specify diagonal elements in SigmaX
#' @param ... additional arguments to pass to data generating functions
#' @return Y, X, betas, Sigma, SigmaX
#'
#' @author Matt Galloway \email{gall0441@@umn.edu}
#'
#' @export
#'
#' @examples
#' # generate 100 observations with predictor dimension equal to 10 and response dimension equal to 5
#' data = data_gen(n = 100, p = 10, r = 5)

# we define the data generation function
data_gen = function(n, p, r = 1, sparsity = 0.5, Sigma = c("tridiag", 
    "dense", "denseQR", "compound"), s = NULL, SigmaX = c("tridiag", 
    "dense", "denseQR", "compound"), sx = NULL, ...) {
    
    # checks
    SigmaX = match.arg(SigmaX)
    Sigma = match.arg(Sigma)
    
    # randomly generate betas
    betas = matrix(rnorm(p * r, 0, sqrt(1/p)), nrow = p, 
        ncol = r)
    betas = betas * matrix(rbinom(p * r, 1, prob = sparsity), 
        nrow = p, ncol = r)
    
    # generate data X
    SigmaX = switch(SigmaX, tridiag = tridiag(p = p, n = n, 
        ...), dense = dense(p = p, n = n, ...), denseQR = denseQR(p = p, 
        n = n, ...), compound = compound(p = p, n = n))
    X = SigmaX$X
    
    SigmaX = SigmaX$S
    if (!is.null(sx)) {
        diag(SigmaX) = sx
    }
    
    # generate sigma matrix
    Sigma = switch(Sigma, tridiag = tridiag(p = r, n = n, 
        ...), dense = dense(p = r, n = n, ...), denseQR = denseQR(p = r, 
        n = n, ...), compound = compound(p = r, n = n))
    
    Sigma = Sigma$S
    if (!is.null(s)) {
        diag(Sigma) = s
    }
    
    # create data
    out = eigen(Sigma, symmetric = TRUE)
    if (length(Sigma > 1)) {
        Sigma.sqrt = out$vectors %*% diag(out$values^0.5) %*% 
            t(out$vectors)
    } else {
        Sigma.sqrt = sqrt(Sigma)
    }
    
    Z = matrix(rnorm(n * r), nrow = n, ncol = r)
    Y = X %*% betas + Z %*% Sigma.sqrt
    
    
    returns = list(Y = Y, X = X, betas = betas, Sigma = Sigma, 
        SigmaX = SigmaX)
    return(returns)
    
}



##-----------------------------------------------------



#' @title Generate tri-diagonal matrices
#' @description Generate p-dimensional matrices so that its inverse is tri-diagonal.
#' @param p desired dimension
#' @param base base multiplier
#' @param n option to generate n observations from covariance matrix S
#' @return Omega, S
#'
#' @keywords internal

# we define the tridiag function
tridiag = function(p = 8, base = 0.7, n = NULL) {
    
    # generate tapered matrices
    S = matrix(0, nrow = p, ncol = p)
    
    for (i in 1:p) {
        for (j in 1:p) {
            S[i, j] = base^abs(i - j)
        }
    }
    
    # oracle
    Omega = qr.solve(S)
    
    # create data, if specified
    if (!is.null(n)) {
        
        # generate n by p matrix X with rows drawn iid N_p(0,
        # sigma)
        out = eigen(S, symmetric = TRUE)
        if (length(S > 1)) {
            S.sqrt = out$vectors %*% diag(out$values^0.5) %*% 
                t(out$vectors)
        } else {
            S.sqrt = sqrt(S)
        }
        
        Z = matrix(rnorm(n * p), nrow = n, ncol = p)
        X = Z %*% S.sqrt
        
        return(list(Omega = Omega, S = S, X = X))
        
    } else {
        
        return(list(Omega = Omega, S = S))
        
    }
    
}




##-----------------------------------------------------



#' @title Generate dense matrices
#' @description Generate p-dimensional matrices so that its inverse is dense.
#' @param p desired dimension
#' @param base base multiplier
#' @param n option to generate n observations from covariance matrix S
#' @return Omega, S
#'
#' @keywords internal

# we define the dense function
dense = function(p = 8, base = 0.9, n = NULL) {
    
    # generate matrix
    S = matrix(base, nrow = p, ncol = p)
    diag(S) = 1
    
    # oracle
    Omega = qr.solve(S)
    
    # create data, if specified
    if (!is.null(n)) {
        
        # generate n by p matrix X with rows drawn iid N_p(0,
        # sigma)
        out = eigen(S, symmetric = TRUE)
        if (length(S > 1)) {
            S.sqrt = out$vectors %*% diag(out$values^0.5) %*% 
                t(out$vectors)
        } else {
            S.sqrt = sqrt(S)
        }
        
        Z = matrix(rnorm(n * p), nrow = n, ncol = p)
        X = Z %*% S.sqrt
        
        return(list(Omega = Omega, S = S, X = X))
        
    } else {
        
        return(list(Omega = Omega, S = S))
        
    }
    
}




##-----------------------------------------------------



#' @title Generate dense matrices (via spectral decomposition)
#' @description Generate p-dimensional matrices so that its inverse is dense. The matrix will be generated so its first 'num' eigen values are 1000 and the remaining are 1. The orthogonal basis is generated via QR decomposition of
#' @param p desired dimension
#' @param num number of 'large' eigen values. Note num must be less than p
#' @param n option to generate n observations from covariance matrix S
#' @return Omega, S
#'
#' @keywords internal

# we define the denseQR function
denseQR = function(p = 8, num = 5, n = NULL) {
    
    # generate eigen values
    eigen = c(rep(1000, num), rep(1, p - num))
    
    # randomly generate orthogonal basis (via QR)
    Q = qr.Q(qr(matrix(rnorm(p * p), nrow = p, ncol = p)))
    
    # generate matrix
    S = Q %*% diag(eigen) %*% t(Q)
    
    # oracle
    Omega = qr.solve(S)
    
    # create data, if specified
    if (!is.null(n)) {
        
        # generate n by p matrix X with rows drawn iid N_p(0,
        # sigma)
        out = eigen(S, symmetric = TRUE)
        if (length(S > 1)) {
            S.sqrt = out$vectors %*% diag(out$values^0.5) %*% 
                t(out$vectors)
        } else {
            S.sqrt = sqrt(S)
        }
        
        Z = matrix(rnorm(n * p), nrow = n, ncol = p)
        X = Z %*% S.sqrt
        
        return(list(Omega = Omega, S = S, X = X))
        
    } else {
        
        return(list(Omega = Omega, S = S))
        
    }
    
}





##-----------------------------------------------------



#' @title Generate compound symmetric matrices
#' @description Generate a p-dimensional compound symmetric matrix.
#' @param p desired dimension
#' @param n option to generate n observations from covariance matrix S
#' @return Omega, S
#'
#' @keywords internal

# we define the compound function
compound = function(p = 8, n = NULL) {
    
    # generate precision matrix
    Omega = matrix(0.9, nrow = p, ncol = p)
    diag(Omega) = 1
    
    # generate covariance matrix
    S = qr.solve(Omega)
    
    # create data, if specified
    if (!is.null(n)) {
        
        # generate n by p matrix X with rows drawn iid N_p(0,
        # sigma)
        out = eigen(S, symmetric = TRUE)
        if (length(S > 1)) {
            S.sqrt = out$vectors %*% diag(out$values^0.5) %*% 
                t(out$vectors)
        } else {
            S.sqrt = sqrt(S)
        }
        
        Z = matrix(rnorm(n * p), nrow = n, ncol = p)
        X = Z %*% S.sqrt
        
        return(list(Omega = Omega, S = S, X = X))
        
    } else {
        
        return(list(Omega = Omega, S = S))
        
    }
    
}
