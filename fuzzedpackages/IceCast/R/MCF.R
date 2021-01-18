#' Computes conditional probability of the observed event under some probability
#' model
#' @title Compute conditional probability of observed event
#' @param obs matrix with binary observations indicating if sea ice
#'            concentration of at least 15$\%$ was observed. Dimension is lon x lat.
#' @param mod matrix with estimated sea ice probability from a model. Dimension
#'            is lon x lat.
cond_prob <- function(obs, mod) {
  stopifnot(dim(obs) == dim(mod))
  cond <- array(dim = dim(obs), data = NA)
  stopifnot(sum(is.na(mod)) == sum(is.na(obs)))
  cond[which(obs == 1)] <- mod[which(obs == 1)]
  cond[which(obs == 0)] <- 1 - mod[which(obs == 0)]
  return(cond)
}


#' Compute weighting between two models based on accuracy in predicting a
#' set of observations. Computation is via the Expectation-Maximization algorithm.
#' @title Compute weighting between two models
#' @param mod1 array with estimated sea ice probability from model 1. Dimensions
#'             are nuumber of training years x lon x lat.
#' @param mod2 array with estimated sea ice probability from model 2. Dimensions
#'             are nuumber of training years x lon x lat.
#' @param obs array with observations of sea ice presence (1) and absence (0).
#'            Dimensions are nuumber of training years x lon x lat.
#' @param prop_area matrix that gives the proportion of area in each grid box.
#'                  Should sum to 1. Dimensions are lon x lat.
#' @param w_ini initial value of all w, defaults to 0.5.
#' @param z_ini initial value of all z, defaults to 0.5.
#' @param eps tolerance for EM algorithm to reach convergence, defaults to 0.01.
#' @return value between 0 and 1 giving the weight on the first model
#' @export
#' @examples
#' \dontrun{
#' weight <- fit_weights(mod1 = clim_9_2005_2007, mod2 = ppe_9_2005_2007,
#' obs = obs_9_2005_2007, prop_area = prop_area)
#' }
fit_weights <- function(mod1, mod2, obs, prop_area, w_ini = 0.5, z_ini = 0.5,
                        eps = 0.01) {
  #set up and initialization
  stopifnot(dim(mod1) == dim(mod2))
  stopifnot(dim(obs) == dim(mod1))
  z <- array(dim = dim(mod1), data = z_ini)
  w <- w_last <- w_ini
  diff = eps + .01 #sets diff so that will enter while loop

  #conditional probabilities of having observed what was observed with both models
  cond_prob_mod1 <- cond_prob(obs, mod1)
  cond_prob_mod2 <- cond_prob(obs, mod2)
  stopifnot(dim(cond_prob_mod1) == dim(cond_prob_mod2))

  #make array of area weights
  n_years <- dim(cond_prob_mod1)[1]
  prop_area_all <- array(dim = dim(cond_prob_mod1))
  for (i in 1:n_years) {
    prop_area_all[i,,] <- prop_area
  }

  #only evaluate points where the two models differ, if both cond_prob_mod1 ==
  #cond_prob_mod2 == 0, z cannot be evaluated
  eval <- which((mod1 != mod2))
  while(diff > eps) {
    #E-step
    z_temp <- w*prop_area_all*cond_prob_mod1/(w*prop_area_all*cond_prob_mod1 +
                                                (1 - w)*prop_area_all*cond_prob_mod2)
    z[eval] <- z_temp[eval]

    #M-step
    w = sum(prop_area_all[eval]*z[eval])/sum(prop_area_all[eval])
    diff <- abs(w - w_last)
    w_last <- w
  }
  return(w)
}

#' Function to weight two models
#' @param w weight on model 1
#' @param mod1 array with estimated sea ice probability from model 1. Dimensions
#'             are nuumber of training years x lon x lat.
#' @param mod2 array with estimated sea ice probability from model 1. Dimensions
#'             are nuumber of training years x lon x lat.
#' @export
#' @examples
#' \dontrun{
#' weight <- fit_weights(mod1 = clim_9_2005_2007, mod2 = ppe_9_2005_2007,
#' obs = obs_9_2005_2007, prop_area = prop_area)
#' wght_mod(w = weight, mod1 = clim_9_2008, mod2 = ppe_9_2008)
#' }
wght_mod <- function(w, mod1, mod2) {
  w*mod1 + (1 - w)*mod2
}


