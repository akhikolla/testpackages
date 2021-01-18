myscale <- function(x, center = TRUE, scale = TRUE){
  # My own scale function
  # Treat specifically the case where columns have zero sds
  n <- nrow(x)
  x <- scale(x, center = center, scale = scale)
  # for those columns with 0 sds
  # add small Gaussian noise to it
  mean_vec <- attr(x = x, which = "scaled:center")
  sd_vec <- attr(x = x, which = "scaled:scale")
  ind <- which(sd_vec == 0)
  if(length(ind) > 0){
    warning("Found columns with zero sd! Replaced this column with a column of all zeros!")
    submat <- matrix(0, nrow = n, ncol = length(ind))
    x[, ind] <- submat
    attr(x = x, which = "scaled:center")[ind] <- 0
    attr(x = x, which = "scaled:scale")[ind] <- 1
  }
  return(x)
}

get_lambda <- function(x, y,
                       nlam = 100,
                       lam_min_ratio = ifelse(nrow(x) < ncol(x),
                                              0.01, 1e-04)){
  # get a path of tuning parameters
  n <- nrow(x)
  lam_max <- max(abs(crossprod(x, y))) / n
  if(lam_max == 0)
    stop('the calculated maximum of tuning parameter is 0! Please check your response and data matrix.')
  return(lam_max * exp(seq(0, log(lam_min_ratio), length = nlam)))
}

extract_inter_indices <- function(iT, p){
  # iT is a n-by-2 matrix, i.e., t(combn(p, 2))
  # Extract position of 2nd-order terms
  # from the index representation iT
  # the column numbers of (x_1, ..., x_p, x_1^2, ..., x_p^2, x_1*x_2, x_1*x_3, ..., x_{p-1}*x_p)
  # by the ordering 1-2, ..., 1-p, 2-3, ..., 2-p, ..., (p-1)-p
  # P(i,j) = 2 * p + (i-1)*p - i*(i-1)/2 + j - i
  # where i,j are row and column indices
  # in case of of ordering 1-1, ..., 1-p, 2-2, 2-3, ..., p-p
  # P(i,j) = (i-1)*(p+1) - i*(i-1)/2 + j - i + 1
  # where i,j are row and column indices
  if(length(iT) == 0)
    return(integer(0))
  else{
    apply(iT, 1, function(x, p){
      if(x[1] == 0)
        x[2]
      else if(x[1] == x[2])
        p + x[1]
      else
        2 * p + (x[1] - 1) * p - x[1] * (x[1] - 1) / 2 + x[2] - x[1]
    }, p)
  }
}

run_path <- function(x_tr, y_tr, x_te, y_te, lambda, intercept = TRUE){
  # if intercept, we need to model intercept
  # otherwise, no need to model intercept
  if(intercept)
    mean_y <- mean(y_tr)
  else
    mean_y <- 0
  # fit a lasso path of lambda on (x_tr, y_tr)
  fit <- glmnet::glmnet(x = x_tr, y = y_tr - mean_y,
                        lambda = lambda,
                        intercept = FALSE,
                        standardize = FALSE)
  # get prediction on x_te
  fitted <- mean_y + x_te %*% fit$beta
  # compute the test mean-squared error on (x_te, y_te)
  mse <- sqrt(colMeans((matrix(y_te, nrow = nrow(x_te),
                          ncol = length(lambda)) - fitted)^2))
  return(as.numeric(mse))
}
