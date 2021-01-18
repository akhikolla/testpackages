#### Generate autoregressive correlation structure
CorrAR <- function(dim, rho){
  Corr <- abs(matrix(rep(1:dim, length = dim) - rep(1:dim, each = dim), nrow = dim))
  Corr <- rho^Corr
  return(Corr)
}

#### Generate compound sysmetry correlation structure
CorrCS <- function(dim, rho){
  Corr <- matrix(rho, nrow = dim, ncol = dim)
  diag(Corr) <- 1
  return(Corr)
}


#' @title
#' Simulation for functional composition data.
#'
#' @description
#' simulate functional compositional data.
#'
#' @usage
#' Fcomp_Model(n, p, m = 0, intercept = TRUE,
#'             interval = c(0, 1), n_T = 100, obs_spar = 0.6, discrete = FALSE,
#'             SNR = 1, sigma = 2, Nzero_group = 4,
#'             rho_X, Corr_X = c("CorrCS", "CorrAR"),
#'             rho_T, Corr_T = c("CorrAR", "CorrCS"),
#'             range_beta = c(0.5, 1), beta_c = 1, beta_C ,
#'             theta.add = c(1, 2, 5, 6), gamma = 0.5,
#'             basis_beta = c("bs", "OBasis", "fourier"), df_beta = 5, degree_beta = 3,
#'             insert = c("FALSE", "X", "basis"), method = c("trapezoidal", "step"))
#'
#'
#' @param n sample size.
#'
#' @param p number of the components in the functional compositional data.
#'
#' @param intercept whether to include an intercept.
#'                  Default is \code{TRUE}.
#'
#' @param m size of unpenalized variables.
#'          The first \code{ceiling(m/2)} ones are generated with independent \code{bin(1,0.5)} entries;
#'          while the last \code{(m - ceiling(m/2))} ones are generated with independent \code{norm(0, 1)} entries.
#'          Default is 0.
#'
#' @param n_T an integer specifying length of the equally spaced time sequence on domian \code{interval}.
#'
#' @param discrete logical (default is \code{FALSE}) specifying whether the functional compositional data
#'                 \eqn{X} is generated at different time points.
#'                 If \code{distrete = TRUE}, generate \eqn{X} on dense sequence created by
#'                 \code{max(ns_dense = 200 * diff(interval), 5 * n_T)} and then for
#'                 each subject, randomly sample \code{n_T} points.
#'
#' @param interval a vector of length 2 indicating the time domain. Default is \code{c(0, 1)}.
#'
#' @param obs_spar a percentage used to get sparse ovbservation. Each time point is with
#'                 probability \code{obs_spar} to be observed. It allows different subject to be observed on
#'                 different time points.
#'                 \code{obs_spar * n_T > 5} is required.
#'
#' @param SNR signal to noise ratio.
#'
#' @param Corr_X,Corr_T character string specifying correlation structure
#'                      bewteen components and between time points, respectively.
#'                      \itemize{
#'                      \item \code{"CorrCS"}(Default for \code{Corr_X}) compound symmetry.
#'                      \item \code{"CorrAR"}(Default for \code{Corr_T}) autoregressive.
#'                      }
#'
#' @param rho_X,rho_T parameters used to generate correlation matrices.
#'
#' @param sigma variance used to generate the covariance matrix
#'              \code{CovMIX = sigma^2 * kronecker(T.Sigma, X.Sigma)}.
#'              The "non-normalized" data \eqn{w_{i}}{w_i} for each subject is
#'              genearted from multivariate normal distribution with covariance \code{CovMIX}.
#'              \code{T.Sigma} and \code{X.Sigma} are correlation matrices for
#'              time points and components, respectively.
#'
#' @param Nzero_group an even integer specifying that the first \code{Nzero_group} compositional predictors
#'                    are with non-zero effects. Default is 4.
#'
#' @param range_beta  a sorted vector of length 2, specifying the range of coefficient
#'                    matrix \code{B} of demension \eqn{p \times k}{p*k}.
#'                    Specifically, each column of \code{B} is filled with \code{Nzero_group/2} values
#'                    from the unifom distribution over \code{range_beta} and their negative counterparts.
#'                    Default is \code{c(0.5, 1)}.
#'
#' @param beta_c value of coefficients for beta0 and beta_c (coefficients for intercept and time-invariant predictors).
#'               Default is 1.
#'
#' @param beta_C vectorized coefficient matrix.
#'               If missing, the program will generate \code{beta_C} according to \code{range_beta}
#'               and \code{Nzero_group}.
#'
#' @param theta.add logical or integer(s).
#'                  \itemize{
#'                  \item If integer(s), a vector with value(s) in \code{[1,p]},
#'                        indicating which component(s) of compostions is of high
#'                        level mean curve.
#'                  \item If \code{TRUE},
#'                        the components \code{c(1:ceiling(Nzero_group/2)} and
#'                        \code{Nzero_group + (1:ceiling(Nzero_group/2)))}
#'                        are set to with high level mean.
#'                  \item if \code{FALSE}, all mean curves are set to 0's.
#'                  }
#'
#' @param gamma for the high-level mean groups, log(p * gamma) is added on the "non-normalized"
#'              data \eqn{w_{i}}{w_i} before the data are converted to be compositional.
#'
#' @param basis_beta,df_beta,degree_beta \code{basis_fun}, \code{k} and \code{degree}
#'                                       in \code{\link{FuncompCGL}} respectively.
#'
#' @inheritParams FuncompCGL
#'
#' @details
#' The setup of this simulation follows
#' Sun, Z., Xu, W., Cong, X., Li G. and Chen K. (2020) \emph{Log-contrast regression with
#' functional compositional predictors: linking preterm infant's gut microbiome trajectories
#' to neurobehavioral outcome}, \href{https://arxiv.org/abs/1808.02403}{https://arxiv.org/abs/1808.02403}
#' \emph{Annals of Applied Statistics}.
#' \cr
#' Specifically, we first generate correlation matrix \code{X.sigma} for components of a composition
#' based on \code{rho_X} and \code{Corr_X}, and correlation matrix \code{T.sigma}
#' for time points based on \code{rho_T} and \code{Corr_T}. Then, the "non-normalized"
#' data \eqn{w_i=[w_i(t_1)^T,...,w_i(t_{n_T})^T]}
#' for each subject are generated from multivariate normal
#' distrubtion with covariance \code{CovMIX = sigma^2 * kronecker(T.Sigma, X.Sigma)}, and
#' the mean vector is determined by \code{theta.add} and \code{gamma}.
#' Each \eqn{w_i(t_v)} is a \code{p}-vector for each time point \eqn{v =1,...,T_n}.
#' Finally, the compositional data are obtained as
#' \deqn{
#' x_{ij}(t_v) = exp(w_{ij}(t_v))/sum_{k=1}^{p} exp(w_{ik}(t_v)),
#' }
#' for each subject \eqn{i=1,...,n}, component of a composition \eqn{j=1,...,p}
#' and time point \eqn{v=1,...,n_T}.
#'
#'
#' @return a list including
#' \item{data}{a list of observed data,
#'             \itemize{
#'             \item \code{y} a vector of response variable,
#'             \item \code{Comp} a data frame of observed functional compositional
#'                               data, a column of \code{Subject_ID}, and a
#'                               column of \code{TIME},
#'             \item \code{Zc} a matrix of unpenalized variables with dimension
#'                             \eqn{n \times m}{n*m},
#'             \item \code{intercept} whether an \code{intercept}
#'                                    is included.
#'             }}
#'
#' \item{beta}{a length \code{p*df_beta + m + 1} vector of coefficients}
#' \item{basis.info}{matrix of the basis function to generate the coefficient curves}
#'
#' \item{data.raw}{ a list consisting of
#'                 \itemize{
#'                 \item \code{Z_t.full} the functional compositional data.
#'                 \item \code{Z_ITG} integrated functional compositional data.
#'                 \item \code{Y.tru} true response vector without noise.
#'                 \item \code{X} functional "non-normalized" data \code{W}.
#'                 }}
#'
#' \item{parameter}{a list of parameters used in the simulation.}
#'
#'
#' @examples
#' Data <- Fcomp_Model(n = 50, p = 30, m = 0, intercept = TRUE, Nzero_group = 4,
#'                     n_T = 20, SNR = 3, rho_X = 0, rho_T = 0.6,
#'                     df_beta = 5, obs_spar = 1, theta.add = FALSE)
#'
#' @export
#'
#' @inherit coef.FuncompCGL author references
#'
#' @importFrom MASS mvrnorm
#' @importFrom plyr alply ldply
#' @importFrom utils tail head

Fcomp_Model <- function(n, p, m = 0, intercept = TRUE,
                        interval = c(0, 1), n_T = 100, obs_spar = 0.6, discrete = FALSE,
                        SNR = 1, sigma = 2, Nzero_group = 4,
                        rho_X, Corr_X = c("CorrCS", "CorrAR"),
                        rho_T, Corr_T= c("CorrAR", "CorrCS"),
                        range_beta = c(0.5, 1), beta_c = 1, beta_C ,
                        theta.add = c(1, 2, 5, 6), gamma = 0.5,
                        basis_beta = c("bs", "OBasis", "fourier"), df_beta = 5, degree_beta = 3,
                        insert = c("FALSE", "X", "basis"), method = c("trapezoidal", "step")) {



  # >>> Warnings <<<
  Corr_X <- match.arg(Corr_X)
  Corr_T <- match.arg(Corr_T)
  insert <- match.arg(insert)
  method <- match.arg(method)
  basis_beta <- match.arg(basis_beta)
  this.call <- match.call()

  if(obs_spar > 1 || obs_spar < 0) stop("obs_spar is a percentage")
  if(basis_beta %in% c("bs", "OBasis") && df_beta < 4) {
    warning("With B-spline basis including intercept DF is at least 4. Set df_beta to be 4.")
    df_beta <- 4 }
  if(basis_beta == "fourier" && df_beta %% 2 == 0) {
    warning("With Fourier basis, DF needs to be odd. Set df_beta = df_beta - 1.")
    df_beta <- df_beta - 1 }
  # <<< Warnings >>>

  # >>> Coefficients:  beta_C generation <<<
  ## beta_C with dimension p by df_beta, stored as a long-vector
  if(missing(beta_C)) {
    C <- matrix(0, nrow = p, ncol = df_beta)
    for(j in 1:df_beta){
      col_input <- round(runif((Nzero_group/2), min = range_beta[1], max = range_beta[2]),
                         digits = 2)
      col_input <- c(col_input, -col_input)
      col_permuation <- sample(Nzero_group)
      C[1:Nzero_group, j] <- col_input[col_permuation]
    }
    beta_C <- as.vector(t(C))
  }
  # <<< Coefficients: beta_C generation >>>


  # >>> Mean pattern <<<
  add.on <- 0
  if(is.numeric(theta.add)) {
    if(length(theta.add) > p) stop("None zero theta must be smaller than p")
    add.on <- theta.add
  } else if (theta.add) {
    ## If true and non-numeric,
    ## first ceiling(Nzero_group/2) of the None-zero groups are set mean with log(p*gamma)
    ## and first ceiling(Nzero_group/2) of the zero groups are set mean with log(p*gamma)
    add.on <- c(1:ceiling(Nzero_group/2), Nzero_group + (1:ceiling(Nzero_group/2)))
  }
  # <<< Mean pattern  >>>


  # >>> time sequence and basis functions <<<
  ## time sequence
  ns_dense <- n_T
  if(discrete) {
    ns_dense <- max(200 * diff(interval), 5 * n_T) # max(1000, 5*n_T)
  }
  sseq <- round(seq(interval[1], interval[2], length.out = ns_dense), 10) # sseq <- round(seq(0, 1, length.out = ns_dense), 10)
  if(n_T* obs_spar < 5) stop("Too Sparse")

  ## basis function
  beta.basis <- switch(basis_beta,
                       "bs" = bs(x = sseq, df = df_beta, degree = degree_beta,
                                 Boundary.knots = interval, intercept = TRUE),
                       "fourier" = eval.basis(sseq,
                                              basisobj = create.fourier.basis(rangeval = interval, nbasis = df_beta )),
                       "OBasis" = {
                         beta_nknots <- df_beta - (degree_beta + 1)
                         if(beta_nknots > 0) {
                           beta_knots <- ((1:beta_nknots) / (beta_nknots + 1)) * diff(interval) + interval[1]
                         } else {
                           beta_knots <- NULL
                         }
                         evaluate(OBasis(expand.knots(c(interval[1], beta_knots, interval[2])),
                                         order = degree_beta + 1),
                                  sseq)
                       })
  # <<< time sequence and basis functions >>>


  # >>> Covariance matrix: CovMIX <<<
  ## X.sigma dependence between compositions
  ## T.sigma dependence between time points
  X.Sigma <- do.call(Corr_X, list(dim = p, rho = rho_X))
  T.Sigma <- do.call(Corr_T, list(dim = ns_dense, rho = rho_T))
  CovMIX <- sigma^2 * kronecker(T.Sigma, X.Sigma)
  # <<<  Covariance matrix: CovMIX >>>



  # >>> multinormal vector: W_T generation. Mean variation included <<<
  mean_curve <- rep(0, times = p * ns_dense)
  group <- matrix(1:(n_T * p), nrow = n_T)
  mean_curve[as.vector(group[, add.on])] <- log(p*gamma)
  W <- mvrnorm(n, mean_curve, CovMIX)
  #### multinormal vector: W_T generation ####


  # >>> Functional compositional Data: X.comp_full generation <<<
  X <- array(dim = c(n, p, ns_dense), NA)
  dimnames(X)[[2]] <- paste0("X", 1:p)
  X.comp <- X
  I <- diag(p)
  for(s in 1:n) X[s, ,] <- matrix(W[s, ], byrow = FALSE, nrow = p, ncol = ns_dense)

  for(s in 1:n_T) {
    ## X.comp follows a logistic normal distribution for each time point after exp(x(s))/sum(exp(x(s)))
    X.comp[, , s] <- exp(X[, , s])
    X.comp[, , s] <- X.comp[, , s] / rowSums(X.comp[, , s])

    # ## Threshold feature
    # ## Under threshold, measurement could not be detected and recorded as 0.
    # ## If threshold=0, all measurement detected.
    # X.comp[, , s][X.comp[, , s] <= threshold] <- 0
  }
  # <<< longitinal composition: X.comp_full >>>


  # >>> Integration:Z_ITG  <<<
  D <- alply(X.comp, .margins = 1,
             function(x, sseq) data.frame(t(rbind(TIME = sseq, x))),
             sseq = sseq )

  if(discrete) {
    D <- lapply(D, function(x, n_T, ns_dense) {
      T.obs <- sample(ns_dense, n_T)
      T.obs <- sort(T.obs)
      x <- x[T.obs, ]
      return(x)
    }, n_T = n_T, ns_dense = ns_dense)
  }
  Z_t.full <- ldply(D[1:n], data.frame, .id = "Subject_ID")


  Z_ITG <- sapply(D, function(x, sseq, beta.basis, insert, method){
    x[, -1] <- log(x[, -1])
    ITG(X = x, basis = beta.basis, sseq = sseq, insert = insert, method = method, interval = interval)$integral
  } ,sseq = sseq, beta.basis = beta.basis, insert = insert, method = method)

  Z_ITG <- t(Z_ITG)
  # <<< Integration:Z_ITG >>>

  # >>> Control varaible combining Intercept: Z_c, generation <<<
  ## beta_c and beta0
  beta0 <- beta_c

  if(m > 0) {
    Z_control <- matrix(NA, nrow = n, ncol = m)
    ## ceiling(m/2) is generated by norm;
    ## (m - ceiling(m/2)) is generated by binary with prob = 0.5
    Z_control[, (floor(m/2) + 1):m] <- matrix(rnorm(n * length((floor(m/2) + 1):m)),
                                              nrow = n)
    Z_control[, -((floor(m/2) + 1):m)] <- matrix(sample(c(0,1), n * floor(m/2), replace = TRUE, prob = c(0.5, 0.5)),
                                                 nrow = n)
    ## Combine Intercept with control variables
    beta_c <- rep(beta_c, times = m)
  } else {
    Z_control <- matrix(0, nrow = n, ncol = 1)
    beta_c <- 0
  }
  # <<< Control varaible combining Intercept >>>

  # >>>  Response: y, generation <<<
  Y.tru <- Z_control %*% beta_c + Z_ITG %*% beta_C + ifelse(intercept, beta0, 0)
  sd_ratio <- sd( Z_control * beta_c) / sd( Z_ITG %*% beta_C)
  error <- rnorm(n, 0, 1)
  sigma_e <- sd(Y.tru) / (sd(error) * SNR)
  error <- sigma_e*error
  Y.obs <- Y.tru + error
  y_sd <- sd(Y.obs)
  # <<< Response: y >>>


  # >>> Compositioanal Data output: Z_t.obs <<<
  ## Each time point is with prob obs_spar to be observed,
  ## Different Subject with varying #observations and observed time points
  if(obs_spar == 1) {
    Z_t.obs <- Z_t.full
  } else {
    Z_t.obs <- lapply(D, function(x, obs_spar) {
      n <- dim(x)[1]
      # ## Observation number follows Possion Distribution
      # n.obs <- rpois(1, obs_spar * n)
      # T.obs <- sample(n, size = n.obs)

      ## Observation follows Bimonial Distribution
      T.obs <- replicate(n, sample(c(0,1), 1, prob = c(1 - obs_spar, obs_spar)))
      T.obs <- which(T.obs == 1)
      x.obs <- x[T.obs, ]
      return(x.obs)
    }, obs_spar = obs_spar)
    Z_t.obs <- plyr::ldply(Z_t.obs, data.frame, .id = "Subject_ID")
  }
  # <<< Compositioanal Data output: Z_t.obs >>>

  # >>> Output <<<
  if(m == 0) {
    Z_control <- NULL
    beta_c <- NULL}

  data <- list(y = Y.obs, Comp = Z_t.obs, Zc = Z_control, intercept = intercept)
  beta <- c(beta_C, beta_c, ifelse(intercept, beta0, 0))
  data.raw <- list(Z_t.full = Z_t.full, Z_ITG = Z_ITG,
                   Y.tru = Y.tru, X = X
                   #W = W
                   )
  basis.info <- cbind(sseq, beta.basis)
  parameter <- list(p = p, n = n, m = m,
                    k = df_beta, n_T = n_T, basis_beta = basis_beta,
                    rho_X = rho_X, rho_T = rho_T, Corr_X = Corr_X, Corr_T = Corr_T, sigma = sigma,
                    degree_beta = degree_beta,
                    SNR = SNR,  sigma_e = sigma_e, sd_ratio = sd_ratio, y_sd = y_sd,
                    theta.add = add.on, obs_spar = obs_spar)

  output <- list(data = data, beta = beta, basis.info = basis.info, data.raw = data.raw,
                 parameter = parameter)
  output$call <- this.call
  # <<< Output >>>

  return(output)
}




#' @title Simulation for log-contrast model.
#'
#' @description Simulate data for log-contrast model with a single set of compositional data.
#'
#' @usage comp_Model(n, p, rho = 0.2, sigma = 0.5, gamma = 0.5, add.on = 1:5,
#'            beta = c(c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2), rep(0, times = p - 8)),
#'            beta0  = 1, intercept = TRUE)
#'
#' @param n sample size
#'
#' @param p number of components in the compositional data
#'
#' @param rho parameter used to generate the \eqn{p \times p}{p*p} autocorrelation
#'            matrix for correlations among the components.
#'            Default is 0.2.
#'
#' @param sigma standard deviation for the noise terms, which are iid normal with mean 0.
#'              Default is 0.5.
#'
#' @param add.on an index vector with value(s) in \code{[1,p]}, specifying which component(s)
#'               of compositions is of high level mean. Default is \code{1:5}.
#'
#' @param gamma a scaler. For the high level mean component(s), \code{log(p * gamma)} is added to
#'              the "non-normalized" data \code{w_i} before the data are
#'              converted to compositional.
#'
#' @param beta coefficients for the compositional variables.
#'
#' @param beta0 coefficient for the intercept. Default is 1.
#'
#' @param intercept whether to include an intercept. Default is \code{FALSE}.
#'
#' @details
#' The setup of this simulation follows
#' Lin, W., Shi, P., Peng, R. and Li, H. (2014) \emph{Variable selection in regression with compositional covariates},
#' \href{https://academic.oup.com/biomet/article/101/4/785/1775476}{https://academic.oup.com/biomet/article/101/4/785/1775476}.
#' \cr
#' Specifically, we first generate the correlation matrix among the components
#' \code{X.Sigma} by \code{rho} with an autoregressive correlation structure. we then
#' generate the "non-normalized" data \eqn{w_i} for each subject from
#' multivariate normal distribution with covariance \code{X.Sigma} and mean
#' determined by \code{add.on} and \code{gamma}. Each \eqn{w_i} is a vector of
#' length \code{p}.
#' Finally, the compositional covariates are obtained as
#' \deqn{
#' x_{ij}=exp(w_{ij})/\sum_{k=1}^{p}exp(w_{ik}),
#' }
#' for each subject \eqn{i=1,...,n} and component \eqn{j=1,...,p}.
#' \cr
#'
#' @return A list containing:
#' \item{y}{a n-vector of the simulated response}
#' \item{X.comp}{a matrix of the simulated compositional predictors of dimension
#'               \eqn{n \times p}{n*p}}
#' \item{Z}{a matrix of the log-transformed compositional predictors}
#' \item{Zc}{a matrix of the simulated covariates}
#' \item{intercept}{whether an intercept is included}
#' \item{beta}{the true coefficient vector}
#'
#' @examples
#' p = 30
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c( beta, rep(0, times = p - length(beta)) )
#' Data = comp_Model(n = 50, p = p, intercept = FALSE,
#'                   rho = 0.2, sigma = 0.5, gamma  = 0.5, add.on = 1:5,
#'                   beta = beta)
#'
#' @inherit coef.compCL author references
#'
#' @export



comp_Model <- function(n, p, rho = 0.2, sigma = 0.5, gamma = 0.5, add.on = 1:5,
                       beta = c( c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2), rep(0, times = p - 8) ),
                       beta0  = 1, intercept = TRUE ) {
  this.call <- match.call()
  SIGMA <- CorrAR(p, rho)
  theta <- rep(0, times = p)
  theta[add.on] <- log(gamma * p)
  W <- mvrnorm(n, theta, SIGMA)
  X <- exp(W)
  X.comp <- X / rowSums(X)
  Z <- log(X.comp)
  Y.tru <- Z %*% beta + as.integer(intercept) * beta0
  y <- Y.tru + rnorm(n, 0, sigma)
  data <- list(y = y, Z = Z, Zc = NULL, intercept = intercept, beta = beta, X.comp = X.comp)
  data$call <- this.call
  return(data)
}



