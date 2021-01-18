#' Dynamic spatio-temporal model with EOFs
#' @useDynLib ideq
#' @importFrom Rcpp sourceCpp
#' @importFrom stats cov
#' @description 
#' Fits a dynamic spatio-temporal model using empirical orthogonal functions
#' (EOFs).
#' The model does not require the spatial locations because the process model
#' is based on the principal components of the data matrix.
#' Three broad model types are supported:
#' 
#' 1. RW: A random walk model for which the process matrix is the identity.
#' 
#' 2. AR: An auto-regressive model for which the process matrix is diagonal
#'    and its elements are estimated.
#'    
#' 3. Dense: A model in which the process matrix is a dense, estimated
#'    matrix.
#'    
#' For each broad model type,
#' users can specify a variety of options including
#' the size of the state space,
#' the form of the process error,
#' and whether to sample the observation error.
#' Users can specify prior distributions for all sampled quantities using the
#' `params` argument.
#' 
#' Each model type mentioned above is a dynamic linear model (DLM),
#' and the state vectors can be estimated using the forward filtering
#' backward sampling algorithm.
#' The other parameters are estimated with conditionally conjugate updates.
#' 
#' @param Y 
#' (numeric matrix) S by T data matrix containing the response variable
#' at S spatial locations and T time points.
#' The t-th column (NOT row) corresponds to the t-th observation vector.
#' @param proc_model 
#' (character string) Process model: one of "RW" (identity process matrix),
#' "AR" (diagonal process matrix), or "Dense" (dense process matrix).
#' @param P 
#' (integer) Number of EOFs or, in other words, the state space size.
#' @param proc_error 
#' (character string) Process error: 
#' "IW" (inverse-Wishart) or "Discount" (discount factor).
#' @param n_samples (numeric scalar) Number of samples to draw
#' @param sample_sigma2 
#' (logical) Whether to sample the variance of the iid observation error.
#' @param verbose 
#' (logical) Whether to print additional information;
#' e.g., iteration in sampling algorithm.
#' @param params 
#' (list) List of hyperparameter values; see details.
#' 
#' @details 
#' This section explains how to specify custom hyperparameters using the `params` argument.
#' For each distribution referenced below,
#' we use the scale parameterization found on the distribution's Wikipedia page.
#' You may specify the following as named elements of the `params` list:
#' 
#' m_0: (numeric vector) The prior mean of the state vector at time zero
#'  (\eqn{\theta_0}).
#'
#' C_0: (numeric matrix) The prior variance-covariance matrix of the state
#' vector at time zero (\eqn{\theta_0}).
#'
#' alpha_sigma2, beta_sigma2: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on the observation error 
#' (\eqn{\sigma^2}).
#' 
#' sigma2: (numeric scalar) The value to use for the observation error 
#' (\eqn{\sigma^2}) if `sample_sigma2` = FALSE.
#' 
#' mu_G: (numeric matrix) The prior mean for the process matrix G.
#' If `proc_model` = "AR", then `mu_G` must be a diagonal matrix.
#' If `proc_model` = "Dense", then `mu_G` has no constraints.
#' 
#' Sigma_G: (numeric matrix) The prior variance-covariance matrix for the
#' process matrix. If proc_model = "AR", then Sigma_G should be P by 
#' P and is the variance-covariance matrix for diag(G).
#' If proc_model = "Dense", then Sigma_G should be P^2 by P^2 and is the 
#' variance-covariance matrix for vec(G).
#' 
#' alpha_lambda, beta_lambda: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on 
#' \eqn{\lambda = (1 - \delta) / \delta},
#' where \eqn{\delta} is the discount factor.
#' 
#' scale_W: (numeric matrix) The scale matrix for the inverse-Wishart prior
#' distribution on the variance-covariance matrix of the process error (`W`).
#' 
#' df_W: (numeric scalar) The degees of freedom for the inverse-Wishart prior
#' distribution on the variance-covariance matrix of the process error (`W`).
#' 
#' @seealso [dstm_ide]
#' 
#' @references 
#' Cressie, N., and Wikle, C. K. (2011), Statistics for spatio-temporal data,
#' John Wiley and Sons, New York, ISBN:978-0471692744.
#' 
#' Fruhwirth-Schnatter, S. (1994), 
#' \dQuote{Data Augmentation and Dynamic Linear Models,}
#' Journal of Time Series Analysis, 15, 183–202,
#' <doi:10.1111/j.1467-9892.1994.tb00184.x>.
#' 
#' Petris, G., Petrone, S., and Campagnoli, P. (2009), 
#' Dynamic Linear Models with R, useR!, Springer-Verlag, New York,
#' ISBN:978-0387772370, <doi:10.1007/b135794>.
#' 
#' @examples
#' # Load example data
#' data("ide_standard")
#' 
#' # Illustrate methods
#' rw_model <- dstm_eof(ide_standard, proc_model="RW", verbose=TRUE)
#' summary(rw_model)
#' predict(rw_model) 
#' 
#' # Other model types
#' dstm_eof(ide_standard, proc_model="AR") # Diagonal process matrix
#' dstm_eof(ide_standard, proc_model="Dense") # Dense process matrix
#' 
#' # Specify hyperparameters
#' P <- 4
#' dstm_eof(ide_standard, P=P,
#'          params=list(m_0=rep(1, P), C_0=diag(0.01, P),
#'                      scale_W=diag(P), df_W=100))
#' @export
dstm_eof <- function(Y, proc_model = "Dense", P = 4L, proc_error = "IW",
                     n_samples = 1L, sample_sigma2 = TRUE, verbose = FALSE,
                     params = NULL) {
  
  # Observation Model; creates F, m_0, C_0
  e <- eigen(cov(t(Y)))
  F_ <-e$vectors[, 1:P]
  
  # Process Model; creates mu_G and Sigma_G
  mu_G <- Sigma_G <- Sigma_G_inv <- matrix()
  
  # mu_G
  mu_G <- if ( is.null(params[["mu_G"]]) ) diag(P) else params[["mu_G"]]
  if (proc_model=="AR" && !matrixcalc::is.diagonal.matrix(mu_G))
    stop("mu_G must be diagonal for proc_model=\"AR\"")
  check.numeric.matrix(mu_G, P)
  
  # Sigma_G
  if (proc_model != "RW") {
    if (proc_model == "AR") {
      Sigma_G <- params[["Sigma_G"]] %else% diag(1e-2, P)
      check.dim(Sigma_G, P, "Sigma_G")
      check.cov.matrix(Sigma_G, P)
    } else if (proc_model == "Dense") {
      Sigma_G <- params[["Sigma_G"]] %else% diag(1e-2, P^2)
      check.dim(Sigma_G, P^2, "Sigma_G", "P^2")
      check.cov.matrix(Sigma_G, P^2)
    }
    # Sigma_G_inv
    Sigma_G_inv <- chol_inv(Sigma_G)
  }
    
  # Process remaining params
  new_params <- process_common_params(params, proc_error, P, sample_sigma2)
  scalar_params <- c(alpha_lambda=new_params[["alpha_lambda"]], 
                     beta_lambda=new_params[["beta_lambda"]],
                     alpha_sigma2=new_params[["alpha_sigma2"]], 
                     beta_sigma2=new_params[["beta_sigma2"]],
                     sigma2=new_params[["sigma2"]], 
                     df_W=new_params[["df_W"]])

  # Fit the model
  results <- eof(Y, F_, mu_G, Sigma_G_inv, new_params[["m_0"]], 
                 new_params[["C_0"]], new_params[["scale_W"]],
                 scalar_params, proc_model, n_samples, verbose)
  
  # Process results
  if ("lambda" %in% names(results))
    results[["lambda"]] <- as.numeric(results[["lambda"]])
  if ("sigma2" %in% names(results))
    results[["sigma2"]] <- as.numeric(results[["sigma2"]])
  
  # Parameters
  results[["scalar_params"]] <- as.list(
    c(scalar_params, 
      n_samples=n_samples, 
      verbose=verbose)
  )
  results[["other_params"]] <- list(m_0=new_params[["m_0"]], 
                                    C_0=new_params[["C_0"]],
                                    scale_W=new_params[["scale_W"]])

  # Set class and attributes
  class(results) <- c("dstm_eof", "dstm", "list")
  attr(results, "proc_model") <- proc_model
  attr(results, "proc_error") <- proc_error
  attr(results, "sample_sigma2") <- sample_sigma2

  return(results)
}

#' Integrodifference equation (IDE) model
#' @useDynLib ideq
#' @importFrom Rcpp sourceCpp
#' @description
#' dstm_ide fits a type of dynamic spatio-temporal model called
#' an integrodifference equation (IDE) model.
#' It estimates a redistribution kernel---a
#' probability distribution controlling diffusion across time and space.
#' Currently, only Gaussian redistribution kernels are supported.
#' 
#' The process model is decomposed with an orthonormal basis function expansion
#' (a Fourier series).
#' It can then be estimated as a special case of a dynamic linear model (DLM),
#' using the forward filtering backward sampling algorithm to estimate the
#' state vector.
#' The kernel parameters are estimated with a random walk Metropolis-Hastings
#' update.
#' The other parameters are estimated with conditionally conjugate updates.
#' 
#' @param Y 
#' (numeric matrix) S by T data matrix containing response variable at S spatial
#' locations and T time points.
#' The t-th column (NOT row) corresponds to the t-th observation vector.
#' @param locs
#' (numeric matrix)
#' S by 2 matrix containing the spatial locations of the observed data.
#' The rows of `locs` correspond with the rows of `Y`.
#' @param knot_locs 
#' (integer or numeric matrix) Knot locations for the spatially varying IDE model.
#' The kernel parameters are estimated at these locations and then mapped to the
#' spatial locations of the observed data via process convolution.
#' If an integer is provided, then the knots are located on an equally spaced
#' grid with dimension (`knot_locs`, `knot_locs`).
#' If a matrix is provided, 
#' then each row of the matrix corresponds to a knot location.
#' If NULL, then the standard (spatially constant) IDE is fit.
#' @param proc_error 
#' (character string) Process error:
#' "IW" (inverse-Wishart) or "Discount" (discount factor).
#' "IW" is recommended because it is more computationally stable.
#' @param J 
#' (integer) Extent of the Fourier approximation.
#' The size of the state space is (2*`J` + 1)^2.
#' @param n_samples 
#' (integer) Number of posterior samples to draw.
#' @param sample_sigma2 
#' (logical) Whether to sample the variance of the iid observation error.
#' @param verbose 
#' (logical) Whether to print additional information;
#' e.g., iteration in sampling algorithm.
#' @param params 
#' (list) List of hyperparameter values; see details.
#' 
#' @details
#' This section explains how to specify custom hyperparameters using the `params` argument.
#' For each distribution referenced below,
#' we use the scale parameterization found on the distribution's Wikipedia page.
#' You may specify the following as named elements of the `params` list:
#' 
#' m_0: (numeric vector) The prior mean of the state vector at time zero
#'  (\eqn{\theta_0}).
#'
#' C_0: (numeric matrix) The prior variance-covariance matrix of the state
#' vector at time zero (\eqn{\theta_0}).
#'
#' alpha_sigma2, beta_sigma2: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on the observation error 
#' (\eqn{\sigma^2}).
#' 
#' sigma2: (numeric scalar) The value to use for the observation error 
#' (\eqn{\sigma^2}) if `sample_sigma2` = FALSE.
#' 
#' alpha_lambda, beta_lambda: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on 
#' \eqn{\lambda = (1 - \delta) / \delta},
#' where \eqn{\delta} is the discount factor.
#' 
#' scale_W: (numeric matrix) The scale matrix for the inverse-Wishart prior
#' distribution on the variance-covariance matrix of the process error (`W`).
#' 
#' df_W: (numeric scalar) The degees of freedom for the inverse-Wishart prior
#' distribution on the variance-covariance matrix of the process error (`W`).
#' 
#' L: (numeric scalar) The period of the Fourier series approximation.
#' The spatial locations and knot locations are rescaled
#' to range from -`L`/4 to `L`/4 because the Fourier decomposition assumes that
#' the spatial surface is periodic.
#' Regardless of the value of `L`, 
#' kernel parameter estimates are back-transformed to the original scale.
#' 
#' smoothing: (numeric scalar) Controls the degree of smoothing in the 
#' process convolution for models with spatially varying kernel parameters.
#' The values in the process convolution matrix are proportional to
#' exp(d/`smoothing`) where d is the distance between spatial locations
#' before rescaling with `L`
#' 
#' mean_mu_kernel: (numeric vector) The mean of the normal prior distribution
#' on `mu_kernel`, the mean of the redistribution kernel.
#' In the spatially varying case, the prior distribution for `mu_kernel`
#' is assumed to be the same at every knot location.
#' 
#' var_mu_kernel: (numeric matrix) The variance of the normal prior distribution
#' on `mu_kernel`, the mean of the redistribution kernel.
#' 
#' scale_Sigma_kernel: (numeric matrix) The scale matrix for the 
#' inverse-Wishart prior distribution on `Sigma_kernel`,
#' the variance-covariance matrix of the redistribution kernel.
#' 
#' df_Sigma_kernel: (numeric scalar) The degrees of freedom for the 
#' inverse-Wishart prior distribution on `Sigma_kernel`,
#' the variance-covariance matrix of the redistribution kernel.
#' 
#' proposal_factor_mu: (numeric scalar) Controls the variance of the proposal distribution for
#' `mu_kernel`. The proposals have a variance of `proposal_factor_mu`^2 * `var_mu_kernel`.
#' `proposal_factor_mu` must generally be set lower for spatially varying models.
#' 
#' proposal_factor_Sigma: (numeric scalar) Controls the variance of the proposal distribution
#' for `Sigma_kernel`. As is the case with `proposal_factor_mu`, a higher value
#' corresponds to a higher variance.
#' The degrees of freedom for the proposal distribution for `Sigma_kernel` is 
#' ncol(`locs`) + `df_Sigma_kernel` / `proposal_factor_Sigma`.
#' `proposal_factor_Sigma` must generally be set lower for spatially varying 
#' models.
#' 
#' kernel_samples_per_iter: (numeric scalar) Number of times to update the kernel
#' parameters per iteration of the sampling loop.
#' 
#' @seealso [dstm_eof]
#' 
#' @references 
#' Cressie, N., and Wikle, C. K. (2011), Statistics for spatio-temporal data,
#' John Wiley and Sons, New York, ISBN:978-0471692744.
#' 
#' Fruhwirth-Schnatter, S. (1994), 
#' \dQuote{Data Augmentation and Dynamic Linear Models,}
#' Journal of Time Series Analysis, 15, 183–202,
#' <doi:10.1111/j.1467-9892.1994.tb00184.x>.
#' 
#' Petris, G., Petrone, S., and Campagnoli, P. (2009), 
#' Dynamic Linear Models with R, useR!, Springer-Verlag, New York,
#' ISBN:978-0387772370, <doi:10.1007/b135794>.
#' 
#' Wikle, C. K., and Cressie, N. (1999), 
#' \dQuote{A dimension-reduced approach to space-time Kalman filtering,} 
#' Biometrika, 86, 815–829, <https://www.jstor.org/stable/2673587>.
#' 
#' Wikle, C. K. (2002), 
#' \dQuote{A kernel-based spectral model for non-Gaussian spatio-temporal processes,}
#' Statistical Modelling, 2, 299–314,
#' <doi:10.1191/1471082x02st036oa>.
#' 
#' @examples
#' # Load example data
#' data("ide_standard", "ide_spatially_varying", "ide_locations")
#' 
#' # Basic IDE model with one kernel
#' mod <- dstm_ide(ide_standard, ide_locations)
#' predict(mod)
#' summary(mod)
#' 
#' # IDE model with spatially varying kernel
#' dstm_ide(ide_spatially_varying, ide_locations, knot_locs=4)
#' 
#' # Fix sigma2
#' dstm_ide(ide_standard, ide_locations, 
#'          sample_sigma2=FALSE, params=list(sigma2=1))
#' 
#' # Set proposal scaling factors, number of kernel updates per iteration,
#' # and prior distribution on kernel parameters
#' dstm_ide(ide_standard, ide_locations, 
#'          params=list(proposal_factor_mu=2, proposal_factor_Sigma=3,
#'                      kernel_updates_per_iter=2,
#'                      scale_Sigma_kernel=diag(2), df_Sigma_kernel=100,
#'                      mean_mu_kernel=c(0.2, 0.4), var_mu_kernel=diag(2)))
#'                      
#' # Set priors on state vector, process error, and observation error
#' J <- 1
#' P <- (2*J + 1)^2
#' dstm_ide(ide_standard, ide_locations,
#'          params=list(m_0=rep(1, P), C_0=diag(0.01, P),
#'                      alpha_sigma2=20, beta_sigma2=20,
#'                      scale_W=diag(P), df_W=100))
#' @export
dstm_ide <- function(Y, locs=NULL, knot_locs=NULL, proc_error = "IW", J=1L,
                     n_samples = 1L, sample_sigma2 = TRUE,
                     verbose = FALSE, params = NULL) {
  
  # Error checking for J, L, locs, knot_locs
  if (is.null(locs)) stop("locs must be specified for IDE models")
  if ("data.frame" %in% class(locs)) locs <- as.matrix(locs)
  if (! "matrix" %in% class(locs)) stop("locs must be data.frame or matrix")
  if (ncol(locs) != 2) stop("locs must have 2 columns")

  J <- as.integer(J)
  if (is.null(J) || is.na(J) || J<1) stop("J must be an integer > 0")
  P <- (2*J + 1)^2

  L <- params[["L"]] %else% 2
  check.numeric.scalar(L)
  transformed_range <- L/2
  
  smoothing <- params[["smoothing"]] %else% 1
  check.numeric.scalar(smoothing)
  
  SV <- is.numeric(knot_locs)
  if (!(SV || is.null(knot_locs))) stop("knot_locs must be numeric or null")
  
  # Center spatial locations to have range of transformed_range
  locs_ranges <- apply(locs, 2, function(x) diff(range(x)))
  locs_offsets <- apply(locs, 2, function(x) max(x) - diff(range(x))/2)
  locs_scaled <- scale_all(locs, transformed_range)
  
  # Process Model; creates kernel parameters
  locs_dim <- ncol(locs)
  
  # kernel_samples_per_iter
  kernel_samples_per_iter <- params[["kernel_samples_per_iter"]] %else% 1
  check.numeric.scalar(kernel_samples_per_iter, x_min=1, strict_inequality=FALSE)
  kernel_samples_per_iter <- as.integer(kernel_samples_per_iter)
  
  # mean_mu_kernel
  mean_mu_kernel <- params[["mean_mu_kernel"]] %else% rep(0, locs_dim)
  check.numeric.vector(mean_mu_kernel, locs_dim, dim_name="ncol(locs)")
  mean_mu_kernel <- as.matrix(mean_mu_kernel)
  mean_mu_kernel <- mean_mu_kernel * (transformed_range) / locs_ranges
  
  # var_mu_kernel
  var_mu_kernel <- params[["var_mu_kernel"]] %else% diag(L/4, locs_dim)
  check.cov.matrix(var_mu_kernel, locs_dim, dim_name="ncol(locs)")
  D <- diag(1/locs_ranges)
  var_mu_kernel <- (transformed_range)^2 * D %*% var_mu_kernel %*% D
  
  # proposal_factor_mu
  proposal_factor_mu <- params[["proposal_factor_mu"]] %else% 1
  check.numeric.scalar(proposal_factor_mu)
  
  # df_Sigma_kernel
  k <- 10
  df_Sigma_kernel <- params[["df_Sigma_kernel"]] %else% (k * locs_dim)
  check.numeric.scalar(df_Sigma_kernel, x_min=locs_dim-1)
  
  # scale_Sigma_kernel
  Sigma_kernel_mean <- transformed_range/10
  scale_Sigma_kernel <- params[["scale_Sigma_kernel"]] %else% 
                          diag((df_Sigma_kernel-locs_dim-1)*Sigma_kernel_mean, locs_dim)
  check.cov.matrix(scale_Sigma_kernel, locs_dim, dim_name="ncol(locs)")
  scale_Sigma_kernel <- (transformed_range)^2 * D %*% scale_Sigma_kernel %*% D
  
  # proposal_factor_Sigma
  proposal_factor_Sigma <- params[["proposal_factor_Sigma"]] %else%
                           ifelse(SV, 1/12, 1)
  check.numeric.scalar(proposal_factor_Sigma)
  
  # Error checking for knot_locs
  # And adjustment to above quantities in spatially varying case
  n_knots <- 1
  K  <- array(0, dim=c(1, 1, 1))
  
  if (SV) {
    # prepare knot_locs
    if (length(knot_locs) > 1) {
      knot_locs_scaled <- scale_all(knot_locs, transformed_range, locs_ranges)
    } else {
      knot_locs_scaled <- gen_grid(as.integer(knot_locs), transformed_range)
      knot_locs <- knot_locs_scaled
      for (i in seq(2)) {
        knot_locs[,i] <- knot_locs[,i] * locs_ranges[i] + locs_offsets[i]
      }
    }
    
    # Modify kernel parameters
    n_knots <- nrow(knot_locs)
    mean_mu_kernel <- rep(1, n_knots) %x% matrix(mean_mu_kernel, nrow=1)
    var_mu_kernel <- var_mu_kernel %x% diag(n_knots)
    
    # Create K matrix
    K <- pdist::pdist(knot_locs, locs)
    K <- as.matrix(K)
    K <- exp(-K / smoothing)
    K <- apply(K, 2, function(x) x / sum(x))
    K <- t(K) # Makes K of dimension n_locs by n_knots
    K <- array(K, dim=c(dim(K), 1))
  }
  
  scale_Sigma_kernel <- array(1, dim=c(1, 1, n_knots)) %x%
                        scale_Sigma_kernel
  
  
  # Process remaining params
  new_params <- process_common_params(params, proc_error, P, sample_sigma2)
  scalar_params <- c(alpha_lambda=new_params[["alpha_lambda"]], 
                     beta_lambda=new_params[["beta_lambda"]],
                     alpha_sigma2=new_params[["alpha_sigma2"]], 
                     beta_sigma2=new_params[["beta_sigma2"]],
                     sigma2=new_params[["sigma2"]], 
                     df_W=new_params[["df_W"]], 
                     kernel_samples_per_iter=kernel_samples_per_iter,
                     proposal_factor_mu=proposal_factor_mu,
                     proposal_factor_Sigma=proposal_factor_Sigma,
                     df_Sigma_kernel=df_Sigma_kernel,
                     SV=SV, J=J, L=L)

  # Fit the model
  results <- ide(Y, locs_scaled, new_params[["m_0"]], new_params[["C_0"]],
                 mean_mu_kernel, var_mu_kernel, K, scale_Sigma_kernel, 
                 new_params[["scale_W"]], scalar_params, n_samples, verbose)
  
  # Process results
  if ("lambda" %in% names(results)) {
    results[["lambda"]] <- as.numeric(results[["lambda"]])
  }
  if ("sigma2" %in% names(results)) {
    results[["sigma2"]] <- as.numeric(results[["sigma2"]])
  }
  if (length(results[["Sigma_kernel"]]) > 1) {
    results[["Sigma_kernel"]] <- results[["Sigma_kernel"]][-1]
  }
  
  # Back transform based on locs_ranges
  for (i in seq(2)) {
    if (SV) {
      for (item_name in c("mu_kernel", "mu_kernel_knots")) {
        results[[item_name]][,i,] <- results[[item_name]][,i,] *
                                     (locs_ranges[i] / (transformed_range))
      }
    } else {
      results[["mu_kernel"]][i,,] <- results[["mu_kernel"]][i,,] *
                                     (locs_ranges[i] / (transformed_range))
    }
  }
  
  results[["Sigma_kernel"]] <- simplify2array(results[["Sigma_kernel"]])
  if (SV) {
    results[["Sigma_kernel_knots"]] <- simplify2array(results[["Sigma_kernel_knots"]])
  }
  for (i in seq(2)) { for (j in seq(2)) {
      results[["Sigma_kernel"]][i,j,,] <- 
        results[["Sigma_kernel"]][i,j,,] *
        (locs_ranges[i] * locs_ranges[j] / (transformed_range)^2)
      if (SV) {
        results[["Sigma_kernel_knots"]][i,j,,] <- 
          results[["Sigma_kernel_knots"]][i,j,,] * 
          (locs_ranges[i] * locs_ranges[j] / (transformed_range)^2)
      }
  }}
  
  # Eliminate extra index when possible
  if (SV) {
    if (dim(results[["process_convolution_map"]])[3] == 1) {
      results[["process_convolution_map"]] <- results[["process_convolution_map"]][,,1]
    }
  } else {
    results[["mu_kernel"]] <- t(results[["mu_kernel"]][,1,])
    results[["Sigma_kernel"]] <- results[["Sigma_kernel"]][,,1,]
  }
  
  # Information from sampling algorithm
  kernel_updates <- n_samples * kernel_samples_per_iter
  results[["mu_acceptance_rate"]] <- results[["mu_acceptances"]] / kernel_updates
  results[["Sigma_acceptance_rate"]] <- results[["Sigma_acceptances"]] / kernel_updates
  results[["mu_acceptances"]] <- NULL
  results[["Sigma_acceptances"]] <- NULL
  
  # Knots locations
  results[["locs"]] <- locs
  results[["locs_scaled"]] <- locs_scaled
  if (SV) {
    results[["knot_locs"]] <- knot_locs
    results[["knot_locs_scaled"]] <- knot_locs_scaled
  }
  
  # Parameters
  results[["scalar_params"]] <- as.list(
    c(scalar_params, 
      n_samples=n_samples, 
      verbose=verbose)
  )
  results[["other_params"]] <- list(m_0=new_params[["m_0"]], 
                                    C_0=new_params[["C_0"]],
                                    mean_mu_kernel=mean_mu_kernel, 
                                    var_mu_kernel=var_mu_kernel, 
                                    scale_Sigma_kernel=scale_Sigma_kernel, 
                                    scale_W=new_params[["scale_W"]])
  
  # Set class and attributes
  class(results) <- c("dstm_ide" , "dstm", "list")
  attr(results, "proc_model") <- "ide"
  attr(results, "proc_error") <- proc_error
  attr(results, "sample_sigma2") <- sample_sigma2

  return(results)
}

.onUnload <- function (libpath) {
  library.dynam.unload("ideq", libpath)
}