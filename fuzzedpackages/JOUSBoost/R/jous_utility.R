
#' Return indices to be used for jittered data in oversampling
#'
#' @param ix_pos Indices for positive examples in data.
#' @param ix_neg Indices for negative examples in data.
#' @param q Quantiles for which to construct tilted datasets.
#' @return returns a list, each of element of which gives indices to be used on
#'         a particular cut (note: will be of length delta - 1)
index_over = function(ix_pos, ix_neg, q){

  if(length(ix_neg) == 0 | length(ix_pos) == 0 )
    stop("Must have at least one positive and one negative example!")

  if(!any(abs(q-0.5) < 1e-6))
    stop("Vector of quantiles must contain median!")

  ncut = length(q)

  ix_neg_cut = vector(mode="list", length=ncut)
  ix_pos_cut = vector(mode="list", length=ncut)

  for(i in 2:ncut){
    ix_neg_cut[[i]] = c(ix_neg_cut[[i-1]], sample(ix_neg, replace=T))
    ix_pos_cut[[ncut - i + 1]] = c(ix_pos_cut[[ncut - i + 2]],
                                   sample(ix_pos, replace=T))
  }

  # Keep the original data for the median classifier
  # gymnastics to assign NULL ...
  median_loc = which(q == 0.5)
  ix_neg_cut[median_loc] = list(NULL)
  ix_pos_cut[median_loc] = list(NULL)

  out = list(ix_neg_cut = ix_neg_cut, ix_pos_cut = ix_pos_cut)
  out
}

#' Return indices to be used in original data for undersampling
#'
#' (note: sampling is done without replacement)
#'
#' @param ix_pos Indices for positive examples in data.
#' @param ix_neg Indices for negative examples in data.
#' @param q Quantiles for which to construct tilted datasets.
#' @param delta Number of quantiles.
#' @return returns a list, each of element of which gives indices to be used on
#'         a particular cut (note: will be of length delta - 1)
index_under = function(ix_pos, ix_neg, q, delta){

  if(length(ix_neg) == 0 | length(ix_pos) == 0 )
    stop("Must have at least one positive and one negative example!")

  if(!any(abs(q - 0.5) < 1e-6))
    stop("Vector of quantiles must contain median!")

  ncut = length(q)
  neach_pos = floor(1/delta * length(ix_pos))
  neach_neg = floor(1/delta * length(ix_neg))

  if(neach_pos == 0 | neach_neg == 0)
    stop("Must have at least floor(1/delta * n_sign) observations!")

  ix_neg_cut = vector(mode="list", length=ncut)
  ix_pos_cut = vector(mode="list", length=ncut)
  ix_neg_cut[[1]] = sample(ix_neg, size=neach_neg)
  ix_pos_cut[[ncut]] = sample(ix_pos, size=neach_pos)

  for(i in 2:ncut){

    ix_neg_cut[[i]] = c(ix_neg_cut[[i-1]],
                        sample(setdiff(ix_neg, ix_neg_cut[[i-1]]),
                               size=neach_neg))
    ix_pos_cut[[ncut - i + 1]] = c(ix_pos_cut[[ncut - i + 2]],
                                   sample(setdiff(ix_pos,
                                                  ix_pos_cut[[ncut - i + 2]]),
                                          size=neach_pos))
  }

  # Keep the original data for the median classifier
  median_loc = which(q == 0.5)
  ix_neg_cut[[median_loc]] = ix_neg
  ix_pos_cut[[median_loc]] = ix_pos

  out = list(ix_neg_cut = ix_neg_cut, ix_pos_cut = ix_pos_cut)
  out
}


#' @useDynLib JOUSBoost, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

