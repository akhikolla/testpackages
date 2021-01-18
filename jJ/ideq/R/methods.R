# Methods for dstm objects

#' Predict Method for DSTM Fits
#' @importFrom stats rnorm
#' @description 
#' Generates samples from the posterior predictive distribution 
#' at future time points for 
#' (1) the observation vector and (2) the state vector.
#' @param object
#' A `dstm` object
#' 
#' @param K 
#' (integer scalar) The number of future time periods for which to generate
#' predictions
#' 
#' @param only_K 
#' (logical scalar) Whether to return predictions for time period T+K only
#' (as opposed to T+1, T+2, ..., T+K)
#' 
#' @param return_ys
#' (logical scalar) Whether to return samples from the posterior predictive
#' distribution of the observation vector (ys)
#' 
#' @param return_thetas
#' (logical scalar) Whether to return samples from the posterior predictive
#' distribution of the state vector (thetas)
#' 
#' @param burnin
#' (integer scalar) The number of samples to discard as burn-in.
#' If object$burnin exists, this argument will override it.
#' 
#' @param ... Arguments passed to other methods (necessary for 
#' S3 generic compatibility)
#' 
#' @details
#' The posterior predictive samples are returned in a matrix or 3-D array,
#' depending on whether samples from multiple time points are requested.
#' The dimensions are always in the following order:
#' 
#' 1. The index of the value within the state or observation vector.
#' 
#' 2. The time period
#' 
#' 3. The sample number
#' @examples
#  Load example data
#' data("ide_standard", "ide_locations")
#' 
#' # IDE example
#' mod_ide <- dstm_ide(ide_standard, ide_locations)
#' predict(mod_ide)
#' predict(mod_ide, K=4, return_thetas=TRUE)
#' 
#' # EOF example
#' mod_eof <- dstm_eof(ide_standard, n_samples=2)
#' predict(mod_eof, K=2, only_K=TRUE, burnin=1)
#' @export
predict.dstm <- function(object, K = 1, only_K = FALSE, return_ys = TRUE,
                         return_thetas = FALSE, burnin = NULL, ...) {
  # Argument burnin is used if provided
  # If not, we use object$burnin if available
  # Otherwise, we use no burnin
  if (is.null(burnin)) {
    burnin <- if ( is.null(object$burnin) ) 0 else object$burnin
  }
  if (burnin < 0) stop("burnin must be non-negative")
  n_samples_total <- dim(object[["theta"]])[3]
  if (burnin >= n_samples_total) stop("burnin must be less than n_samples")
  if (!return_ys && !return_thetas)
    stop("return_ys and/or return_thetas must be true")
  
  # Model meta information
  Discount <- attr(object, "proc_error") == "Discount"
  RW <- attr(object, "proc_model") == "RW"
  
  # Dimensions
  S  <- nrow(object[["F"]])
  P  <- dim(object[["theta"]])[1]
  Tp1 <- dim(object[["theta"]])[2] # T plus 1
  n_samples <- dim(object[["theta"]])[3] - burnin
  idx_post_burnin <- seq(n_samples) + burnin
  
  # Create copies of objects needed for sampling
  thetas_prev <- object[["theta"]][,Tp1,idx_post_burnin]
  thetas_prev <- as.matrix(thetas_prev) # In case it's a vector
  sigma2 <- object[["sigma2"]][idx_post_burnin]
  if (RW) {
    G <- diag(P)
  } else {
    G <- object[["G"]][,,idx_post_burnin]
  }
  # In case G is a matrix instead of an array
  if (length(dim(G)) == 2) G <- array(1, dim=c(1, 1, 1)) %x% G
  
  if (Discount) {
    lambda <- object[["lambda"]][idx_post_burnin]
    C_T <- object[["C_T"]][,,idx_post_burnin]
    if (length(dim(C_T)) == 2) C_T <- array(1, dim=c(1, 1, 1)) %x% C_T
  } else {
    W <- object[["W"]][,,idx_post_burnin]
    if (length(dim(W)) == 2) W <- array(1, dim=c(1, 1, 1)) %x% W
  }
  
  # Step 1: Sample thetas from posterior predictive distribution
  thetas <- array(NA, dim = c(P, K, n_samples))
  
  # Create W for discount models
  if (Discount) {
    C_Tpk <- update_C(C_T, G)
    W <- calc_W(lambda, C_Tpk)
  }
  
  thetas[, 1,] <- next_thetas(thetas_prev, G, W)
  if (K > 1) {
    for (k in seq(2, K)) {
      if (Discount) {
        C_Tpk <- update_C(C_Tpk, G)
        W <- calc_W(lambda, C_Tpk)
      }
      thetas[,k,] <- next_thetas(as.matrix(thetas[,k-1,]), G, W)
    }
  }

  # Calculate ys
  if (return_ys) {
    # Calculate standard deviation for observation model
    sample_sigma2 <- is.logical(attr(object, "sample_sigma2")) &&
                     attr(object, "sample_sigma2")
    if (sample_sigma2) my_sd <- rep(sqrt(sigma2), each=P)
    else my_sd <- sqrt(object[["sigma2"]])

    # Function to get predicted y values for a given time period
    get_preds <- function(k) object[["F"]] %*% thetas[,k,] + rnorm(n_samples*S, sd=my_sd)

    # Get predicted y values for requested time period
    if (only_K || K < 2) {
      ys <- get_preds(K)
    }
    else { # Get predicted y values for all time period <= T+K
      ys <- array(NA, dim = c(S, n_samples, K))
      for (k in seq(K)) ys[,,k] <- get_preds(k)
    }

    # Create output list depending on whether user wants thetas
    if (return_thetas) {
      if (only_K || K < 2) thetas <- thetas[,K,]
      results <- list(ys = ys, thetas = thetas)
    } else results <- ys

  # Create output for case when user does not want ys
  } else {
    if (only_K || K < 2) results <- thetas[,K,]
    else results <- thetas
  }

  return(results)
}

#' Summary Method for DSTM Fits
#' @importFrom stats var quantile
#' @description 
#' Prints summary information for `dstm` objects.
#' @param object A `dstm` object
#' @param object_name The name to be printed in the summary (if desired)
#' @param ... Arguments passed to other methods (necessary for 
#' S3 generic compatibility)
#' @examples
#' # Load example data
#' data("ide_standard", "ide_locations")
#' mod_ide <- dstm_ide(ide_standard, ide_locations)
#' summary(mod_ide)
#' @seealso [print.summary()]
#' @export
summary.dstm <- function(object, object_name = deparse(substitute(object)), ...) {
  cat("Summary for dstm object \`", object_name, "\`\n", sep = "")
  cat("Process model: `", attr(object, "proc_model"), "`\n", sep="")
  cat("Process error: `", attr(object, "proc_error"), "`\n", sep="")
  if (attr(object, "sample_sigma2")) {
    cat("sigma2 was sampled\n\n")
  } else {
    cat("sigma2 was fixed\n\n")
  }

  cat("List elements (in order) are as follows:\n")
  cat(names(object), "\n")
  
  # Numeric Scalars
  numeric_bool <- sapply(object, function(y) is.vector(y) && is.numeric(y))
  scalar_bool <- sapply(object, function(y) length(y) == 1)
  scalar_idx <- which(numeric_bool &  scalar_bool)
  vector_idx <- which(numeric_bool & !scalar_bool)
  
  if ( length(scalar_idx) > 0 ) {
    cat("\nScalar Objects:\n")
    scalars <- numeric()
    for (i in scalar_idx) scalars <- c(scalars, object[[i]])
    names(scalars) <- names(object)[scalar_idx]
    print(scalars)
  }

  # Numeric Vectors
  if ( length(vector_idx) > 0 ) {
    cat("\nVector Objects:\n")
    vector_summary <- matrix(NA, nrow = length(vector_idx), ncol = 8)
    my_probs = c(0.0, 0.25, 0.5, 0.75, 1.0)
    counter <- 1
    for (i in vector_idx) {
      vector_summary[counter, 1]   <- length(object[[i]])
      vector_summary[counter, 2]   <- mean(object[[i]])
      vector_summary[counter, 3]   <- var(object[[i]])
      vector_summary[counter, 4:8] <- quantile(object[[i]], probs = my_probs)
      counter <- counter + 1
    }
    colnames(vector_summary) <- c("Length", "Mean", "Var",
                                 paste0(as.character(my_probs * 100), "%"))
    rownames(vector_summary) <- names(object)[vector_idx]
    print(vector_summary)
  }

  # Matrices/Arrays
  mat_arr_idx <- which(sapply(object, function(y) is.matrix(y) || is.array(y)))
  if ( length(mat_arr_idx) > 0 ) {
    cat("\n\nMatrices/Arrays:\n")
    mat_arr_summary <- matrix(NA, nrow = length(mat_arr_idx), ncol = 5)
    counter <- 1
    for (i in mat_arr_idx) {
      mat_arr_summary[counter, 1] <- class(object[[i]])[1]
      dims_i <- dim(object[[i]])
      mat_arr_summary[counter, 2:(1 + length(dims_i))] <- dims_i
      counter <- counter + 1
    }
    colnames(mat_arr_summary) <- c("class", "dim 1", "dim 2", "dim 3", "dim 4")
    rownames(mat_arr_summary) <- names(object)[mat_arr_idx]
    mat_arr_summary <- as.data.frame(mat_arr_summary)
    
    # Remove dim 4 if not needed
    if (all(is.na(mat_arr_summary[,5]))) {
      mat_arr_summary <- mat_arr_summary[,-5]
    }
    print(mat_arr_summary)
  }
  
  # Lists
  list_idx <- which(sapply(object, function(y) is.list(y)))
  if ( length(list_idx) > 0 ) {
    cat("\n\nLists:\n")
    list_summary <- matrix(sapply(object[list_idx], length))
    colnames(list_summary) <- c("length")
    rownames(list_summary) <- names(object)[list_idx]
    print(list_summary)
  }
}

#' Print Method for DSTM Fits
#' @description 
#' Prints a summary for a `dstm` object by calling summary.dstm().
#' @param x A `dstm` object
#' @param x_name (optional) Object name to display
#' @param ... Arguments passed to other methods (necessary for 
#' S3 generic compatibility)
#' @seealso [summary.dstm()]
#' @export
print.dstm <- function(x, x_name = deparse(substitute(x)), ...) {
  summary.dstm(x, object_name = x_name)
}
