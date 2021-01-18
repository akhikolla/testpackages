source('R/init.R')
source('R/ref.R')
source('R/c.R')

# x list of data blocks
# k rank of approximation
# p list of view dimensions
# n number of views
# m number of data blocks

#' Multiview principal component analysis
#'
#' Analyzes several related matrices of data.
#'
#' @param x List of matrices to analyze
#' @param inds Matrix containing view indices. The matrix should have two
#'   columns and the same number of rows as the length of \code{x}. The first
#'   (second) column contains the view index of the rows (columns) of the
#'   corresponding matrix.
#' @param k Integer giving the maximum rank of the analysis, i.e. the maximum
#'   number of principal components for each view.
#' @param lambda Vector or matrix of lambda values. The length (or width if it
#'   is a matrix) depends on the number of active penalties (2, 3 or 4). If it
#'   is a matrix, try different lambda values (one try for each row). Default: a
#'   matrix where each column is the sequence \code{exp(seq(-6, 0)))}.
#' @param trace Integer selecting the amount of log messages. 0 (default): no
#'   output, 3: all output.
#' @param init_theta NULL, functions or numeric. NULL (default) use initial
#'   values based on ordinary SVD. If init_theta is a list of three functions
#'   (\code{CMF}, \code{matrix_to_triplets} and \code{getCMFopts} from package
#'   \code{CMF}) use the supplied functions to find initial values with
#'   collaborative matrix factorization (CMF). If init_theta is a numeric vector
#'   it is used as initial value.
#' @param cachepath Character vector with path to directory to store
#'   intermediate results. If NULL (default) intermediate results are not
#'   stored. For caching to work it is required that the random number
#'   generation seed is constant between calls to mmpca, so \code{set.seed}
#'   needs to be called before mmpca.
#' @param enable_rank_selection Boolean deciding if the second penalty that
#'   imposes a low rank model should be enabled.
#' @param enable_sparsity Boolean deciding if the third penalty that imposes
#'   sparsity in V should be enabled.
#' @param enable_variable_selection Boolean deciding if the fourth penalty that
#'   increases the tendence for sparsity structure of different V columns to be
#'   similar. Defaults to FALSE meaning this penalty is not used.
#' @param parallel Boolean deciding if computations should be run on multiple
#'   cores simultaneously.
#'
#' @return A list with components
#'   \item{initial}{initial values used in optimization}
#'   \item{cmf}{solution found with CMF (if init_theta == c(CMF,
#'     matrix_to_triplets, getCMFopts))}
#'   \item{training}{solutions for different values of lambda}
#'   \item{solution}{solution for optimal lambda value}
#'
#' @examples
#' x <- list(matrix(rnorm(110), 10, 11), matrix(rnorm(120), 10, 12))
#' inds <- matrix(c(1, 1, 2, 3), 2, 2)
#' result <- mmpca(x, inds, 3, parallel=FALSE)
#'
#' @author Jonatan Kallus, \email{kallus@@chalmers.se}
#' @keywords pca models multivariate
#'
#' @export
mmpca <- function(x, inds, k, lambda=NULL, trace=0, init_theta=NULL,
    cachepath=NULL, enable_rank_selection=TRUE, enable_sparsity=TRUE,
    enable_variable_selection=FALSE, parallel=TRUE) {

  if (enable_variable_selection && !enable_sparsity) {
    stop('Variable selection can not be enabled when sparisty is not enabled.')
  }

  nparam <- 2 + enable_sparsity + enable_variable_selection

  if (is.null(lambda)) {
    lambda <- matrix(rep(exp(seq(-6, 0)), nparam), ncol=nparam)
    if (!enable_rank_selection) {
      lambda[, 2] <- 0
    }
  }

  if (is.null(dim(lambda))) {
    lambda <- matrix(lambda, 1)
  }

  validate_inds(x, inds)

  theta <- init_theta
  p <- init_view_dimensions(x, inds)
  n <- length(p)

  # create cachepath directory unless is exists
  if (!is.null(cachepath)) {
    if (file.exists(cachepath) && !dir.exists(cachepath)) {
      stop('cachepath is not a directory')
    }
    if (!dir.exists(cachepath)) {
      dir.create(cachepath)
    }
  }

  # missing values, training set and test set
  missing_masks <- lapply(x, function(x) !is.na(x))
  for (i in 1:length(x)) x[[i]][is.na(x[[i]])] <- 0
  train_masks <- lapply(x,
    function(x) 0.9 > matrix(stats::runif(prod(dim(x))), nrow(x), ncol(x)))
  test_masks <- lapply(train_masks, function(x) !x)
  train_masks <- lapply(1:length(x),
    function(i) 1 * (train_masks[[i]] & missing_masks[[i]]))
  test_masks <- lapply(1:length(x),
    function(i) 1 * (test_masks[[i]] & missing_masks[[i]]))
  missing_masks <- lapply(1:length(x), function(i) 1 * missing_masks[[i]])

  # rescale data so that biggest sing val is pi^2
  # factor is used to scale found D before returning
  max_sing_val <- max(sapply(x, function(x) svd(x)$d[1]))
  data_scale_factor <- pi^2 / max_sing_val
  x <- lapply(x, function(x) x * data_scale_factor)

  # find scaling factor for lambda
  lambda_factor <- mean(sapply(1:length(x),
    function(i) sqrt(sum((x[[i]] * train_masks[[i]])^2))))^(3/2)

  # find initial values
  cmf_result <- NULL
  if (length(theta) == 3 && is.function(theta[[1]]) && is.function(theta[[2]])
      && is.function(theta[[3]])) { # find initial values xi and D with CMF
    cmf_fun <- theta
    if (trace) msg('Finding initial values... ')
    if (trace > 1) msg('\n')
    # set NA in data for CMF for missing data and test data
    x_cmf <- lapply(1:length(x), function(i) {
      xx <- x[[i]]
      xx[test_masks[[i]] | !missing_masks[[i]]] <- NA
      return(xx)
    })
    if (trace > 1) msg('CMF:\n')
    cmf_result <- cmf_cached(x_cmf, inds, k, trace, cachepath, cmf_fun)
    # calculate training and test error of CMF
    cmf_result$loss <- sum(sapply(1:length(x), function(i) {
      U <- cmf_result$U[inds[i, ]]
      sum(train_masks[[i]] * (x[[i]] - U[[1]] %*% t(U[[2]]))^2)
    }))
    cmf_result$test_loss <- sum(sapply(1:length(x), function(i) {
      U <- cmf_result$U[inds[i, ]]
      sum(test_masks[[i]] * (x[[i]] - U[[1]] %*% t(U[[2]]))^2)
    }))
    normalize <- function(M) {
      for (i in 1:ncol(M)) {
        M[, i] <- M[, i] / sqrt(sum(M[, i]^2))
      }
      return(M)
    }
    V <- lapply(cmf_result$U, normalize)
    # set 0 for test data when initializing singular values
    x_init <- lapply(1:length(x), function(i) {
      xx <- x[[i]]
      xx[test_masks[[i]]] <- 0
      return(xx)
    })
    singular_values <- init_singular_values(V, x_init, inds)
    if (trace > 1) msg('Inverse Euler transformation...')
    if (parallel) {
      xi <- parallel::mclapply(V, function(v) init_inv_v_cached(v, cachepath),
        mc.preschedule=FALSE)
    } else {
      xi <- lapply(V, function(v) init_inv_v_cached(v, cachepath))
    }
    if (trace > 1) msg('done\n')
    D <- init_d(singular_values)
    theta <- ref_vectorize(xi, D)
    if (trace) {
      msg('done\nCMF training loss: ')
      msg(round(cmf_result$loss, 4))
      msg('\nCMF test loss: ')
      msg(round(cmf_result$test_loss, 4))
      msg('\nOptimizing hyper parameters...\n')
    }
  } else if (is.null(theta)) {
    V <- init_v(x, inds, k)
    singular_values <- init_singular_values(V, x, inds)
    D <- init_d(singular_values)
    if (trace > 1) msg('Inverse Euler transformation...')
    if (parallel) {
      xi <- parallel::mclapply(V, function(v) init_inv_v_cached(v, cachepath),
        mc.preschedule=FALSE)
    } else {
      xi <- lapply(V, function(v) init_inv_v_cached(v, cachepath))
    }
    theta <- ref_vectorize(xi, D)
    #theta <- stats::runif(sum(p * k) + k * n, -pi, pi)
  }

  result <- list()
  xiD <- ref_unvectorize(theta, k, max(inds), p)
  result$initial <- list(xi=xiD$xi, D=xiD$D, theta=theta)
  result$cmf <- cmf_result

  if (dim(lambda)[1] > 1) {
  results <- list()

  L <- function(lambda) {
    if (trace > 1) {
      msg('lambda: ', lambda, '\n')
    }
    c_init_parallel()
    res <- optim_mmpca_cached(theta, x, train_masks, inds, k, p, lambda,
      lambda_factor, trace > 2, cachepath, 1)
    # change small values to exact zeros in D
    theta <- res[[1]]
    ix <- (1+sum(p*k)):length(theta)
    theta[ix][abs(theta[ix]) < 1e-5] <- 0
    res[[1]] <- theta
    loss <- c_objective(theta, x, test_masks, inds, k, p,
      rep(0, length(lambda)))
    if (trace > 0) {
      msg('lambda: ', lambda, ' -> test loss: ', loss, '\n')
    }
    return(list(value=loss,
      extra=list(lambda=lambda, theta=theta, loss=loss, res=res)))
  }

  if (parallel) {
    c_init_parallel()
    tmp_res <- parallel::mclapply(1:nrow(lambda), function(i) L(lambda[i, ]),
      mc.preschedule=FALSE)
  } else {
    tmp_res <- lapply(1:nrow(lambda), function(i) L(lambda[i, ]))
  }

  for (i in 1:length(tmp_res)) {
    results[[i]] <- tmp_res[[i]]$extra
  }

  if (trace) msg('Postprocessing solutions... ')
  solutions <- list()
  losses <- rep(NA, length(results))
  for (i in 1:length(results)) {
    r <- results[[i]]
    losses[i] <- r$loss
    res <- r$res

    solutions[[i]] <- list(lambda=r$lambda, theta=r$theta,
      iterations=res[[2]], status=res[[3]],
      cv_error=r$loss, unpenalized_objective_value=res[[4]],
      stepsize=res[[5]], message=res[[6]])
  }

  # order solutions from worst to best
  result$training <- solutions[order(-losses)]
  theta <- result$training[[length(result$training)]]$theta
  lambda <- result$training[[length(result$training)]]$lambda

  if (trace) msg('done\n')
  }

  # calculate final solution using all data
  if (trace) msg('Calculating final solution using all data... ')
  if (trace > 1) msg('\n')
  res <- optim_mmpca_cached(theta, x, missing_masks, inds, k, p, lambda,
    lambda_factor, trace > 1, cachepath,
    ifelse(parallel, parallel::detectCores(), 1))
  # change small values to exact zeros in D
  theta <- res[[1]]
  ix <- (1+sum(p*k)):length(theta)
  theta[ix][abs(theta[ix]) < 1e-5] <- 0
  res[[1]] <- theta

  xiD <- ref_unvectorize(theta, k, n, p)

  # rescale D
  xiD$D <- xiD$D / sqrt(data_scale_factor)

  component_importances <- calculate_component_importances(xiD$D, inds,
    lapply(x, function(x) x / data_scale_factor))

  # order components by total importance and calculate xhat
  ix <- order(-component_importances$total)
  V <- lapply(xiD$xi, function(xi) ref_Vxi(xi)[, ix, drop=FALSE])
  # set exact zeros in V
  for (i in 1:n) V[[i]][abs(V[[i]]) < 1e-5] <- 0
  xiD$D <- xiD$D[ix, ]
  component_importances$block <- component_importances$block[ix, ]
  component_importances$total <- component_importances$total[ix]
  xhat <- list()
  for (j in 1:length(x)) {
    row <- inds[j, 1]
    col <- inds[j, 2]
    xhat[[j]] <- V[[row]] %*% diag(xiD$D[, row] * xiD$D[, col]) %*%
      t(V[[col]])
  }

  result$solution <- list(V=V, D=xiD$D, xhat=xhat, lambda=lambda,
    R2_blockwise=component_importances$block,
    R2_total=component_importances$total, test_masks=test_masks,
    iterations=res[[2]], status=res[[3]], stepsize=res[[5]], message=res[[6]])
  if (trace) msg('done\n')

  return(result)
}

optim_mmpca_cached <- function(theta, x, masks, inds, k, p, lambda,
    lambda_factor, trace, path, nparallel) {
  hash <- digest::digest(list(x, masks, inds, k, p, lambda_factor))
  lambda_str <- paste(lambda, collapse='_')
  if (!is.null(path)) {
    filename <- file.path(path, paste(hash, lambda_str, sep=':'))
    if (file.exists(filename)) {
      return(readRDS(filename))
    }
  }
  # make sparsity factor have same magnitude
  factor <- rep(1, length(lambda))
  if (length(lambda) > 2) {
    factor[3] <- 1/length(p)
    if (length(lambda) > 3) {
      factor[4] <- 1/length(p)
    }
  }
  res <- c_optim_mmpca(theta, x, masks, inds, k, p,
    lambda * lambda_factor * factor, trace, nparallel)
  if (!is.null(path)) {
    saveRDS(res, filename)
  }
  return(res)
}

cmf_cached <- function(data, views, K, trace, path, cmf_fun) {
  hash <- digest::digest(list(data, views, K))
  if (!is.null(path)) {
    filename <- file.path(path, paste(hash, 'cmf', sep=':'))
    if (file.exists(filename)) {
      return(readRDS(filename))
    }
  }
  # make sparsity factor have same magnitude
  res <- cmf(data, views, K, trace, cmf_fun)
  if (!is.null(path)) {
    saveRDS(res, filename)
  }
  return(res)
}

init_inv_v_cached <- function(v, path) {
  hash <- digest::digest(v)
  if (!is.null(path)) {
    filename <- file.path(path, paste(hash, 'init_inv_v', sep=':'))
    if (file.exists(filename)) {
      return(readRDS(filename))
    }
  }
  # make sparsity factor have same magnitude
  res <- init_inv_v(v)
  if (!is.null(path)) {
    saveRDS(res, filename)
  }
  return(res)
}

cmf <- function(data, views, K, trace, cmf_fun) {
  D <- rep(NA, max(views))
  for (i in 1:nrow(views)) {
    D[views[i, 1]] <- nrow(data[[i]])
    D[views[i, 2]] <- ncol(data[[i]])
  }
  data <- lapply(data, cmf_fun[[2]])
  opts <- cmf_fun[[3]]()
  opts$verbose <- trace - 1
  cmf_fun[[1]](data, views, K, rep('gaussian', length(data)), D, opts=opts)
}

mmpca_lambda1 <- function(x, inds, k, lambda, nparallel, init=FALSE,
    trace=FALSE) {
  result <- list()
  p <- init_view_dimensions(x, inds)

  # handle missing values in x
  masks <- lapply(x, function(x) 1 * !is.na(x))
  for (i in 1:length(x)) {
    x[[i]][is.na(x[[i]])] <- 0
  }

  if (init) { # find initial values xi and D
    V <- init_v(x, inds, k)
    xi <- lapply(V, init_inv_v)
    singular_values <- init_singular_values(V, x, inds)
    D <- init_d(singular_values)
    theta <- ref_vectorize(xi, D)
  } else { # not using the initialization!
    theta <- seq(0, 1, length=sum(p * k) + k * length(p))
  }
  xiD <- ref_unvectorize(theta, k, max(inds), p)
  solutions <- list()
  result$initial <- list(xi=xiD$xi, D=xiD$D, theta=theta)
  for (i in 1:length(lambda)) {
    if (trace) {
      msg('lambda: ', lambda[i], '\n')
    }

    res <- c_optim_mmpca(theta, x, masks, inds, k, p, lambda[i], trace,
      nparallel)
    #res <- ref_optim_mmpca(theta, x, inds, k, p, lambda[i], trace)
    theta <- res[[1]]
    xiD <- ref_unvectorize(theta, k, length(p), p)
    component_importances <- calculate_component_importances(xiD$D, inds, x)

    if (i == 1) { # reorder components according to importances
      ix <- order(component_importances$total, decreasing=TRUE)
      for (j in 1:length(xiD$xi)) {
        V <- ref_Vxi(xiD$xi[[j]])
        xiD$xi[[j]] <- init_inv_v(V[, ix])
      }
      xiD$D <- xiD$D[ix, ]
      theta <- ref_vectorize(xiD$xi, xiD$D)
      component_importances$block <- component_importances$block[ix, ]
      component_importances$total <- component_importances$total[ix]
    }

    if (trace) {
      msg('iter: ', res[[2]], ' status: ', res[[3]])
      msg(' upval: ', res[[4]], ' step: ', res[[5]])
      msg(' msg: ', res[[6]], '\n')
    }

    V <- lapply(xiD$xi, function(xi) ref_Vxi(xi))
    xhat <- list()
    for (j in 1:length(x)) {
      row <- inds[j, 1]
      col <- inds[j, 2]
      xhat[[j]] <- V[[row]] %*% diag(xiD$D[, row] * xiD$D[, col]) %*%
        t(V[[col]])
    }

    solutions[[i]] <- list(xi=xiD$xi, D=xiD$D, theta=theta, V=V,
      xhat=xhat, iterations=res[[2]], status=res[[3]],
      unpenalized_value=res[[4]], stepsize=res[[5]], message=res[[6]],
      component_importances_blockwise=component_importances$block,
      component_importances_total=component_importances$total)
  }
  result$solutions <- solutions
  result$lambda <- lambda
  return(result)
}

# calculate R^2, the proportion of variance explained
calculate_component_importances <- function(D, inds, x) {
  res <- matrix(NA, nrow(D), nrow(inds))
  normalizer <- 0
  for (i in 1:nrow(inds)) {
    res[, i] <- (D[, inds[i, 1]] * D[, inds[i, 2]])^2
    normalizer <- normalizer + sum(x[[i]]^2)
  }
  list(block=res/normalizer, total=rowSums(res)/normalizer)
}

validate_inds <- function(x, inds) {
  for (i in unique(c(inds))) {
    p <- NA
    for (j in which(inds[, 1] == i)) {
      if (!is.na(p)) {
        if (p != nrow(x[[j]])) {
          stop('Some view has inconsistent size')
        }
      }
      p <- nrow(x[[j]])
    }
    for (j in which(inds[, 2] == i)) {
      if (!is.na(p)) {
        if (p != ncol(x[[j]])) {
          stop('Some view has inconsistent size')
        }
      }
      p <- ncol(x[[j]])
    }
  }
}

msg <- function(...) {
  message(..., appendLF=F)
}
