######################################################################################################*
######################################################################################################*
#'@title Calculates Relative Efficiency for Locally Optimal Designs
#' @description
#' Given a vector of initial estimates for the parameters, this function calculates the D-and PA- efficiency of a design \eqn{\xi_1} with respect to a design \eqn{\xi_2}.
#' Usually, \eqn{\xi_2} is an  optimal design.
#'
#'
#' @details For a known \eqn{\theta_0}, relative D-efficiency is
#' \deqn{exp(\frac{log|M(\xi_1, \theta_0)| - log|M(\xi_2, \theta_0)|}{npar})}{
#' exp\{(log|M(\xi_1, \theta_0)| - log|M(\xi_2, \theta_0)|)/npar\}.}
#' The relative P-efficiency is
#' \deqn{\exp(\log(\sum_{i=1}^k w_{1i}p(x_{1i}, \theta_0) - \log(\sum_{i=1}^k w_{2i}p(x{2_i}, \theta_0))}{
#' exp(log (\sum w1_i p(x1_i, \theta_0) - log(\sum w2_i p(x2_i, \theta_0)),
#' }
#' where \eqn{x_2}{x1} and \eqn{w_2}{w1} are usually the support points and the corresponding weights of the optimal design, respectively.
#'
#' The argument  \code{x1} is the vector of design points.
#'  For design points with more than one dimension (the models with more than one predictors),
#'    it is a concatenation of the design points, but \strong{dimension-wise}.
#'    For example, let the model has three predictors   \eqn{(I, S, Z)}.
#'     Then,  a two-point optimal design has the following points:
#'    \eqn{\{\mbox{point1} = (I_1, S_1, Z_1), \mbox{point2} = (I_2, S_2, Z_2)\}}{{point1 = (I1, S1, Z1), point2 = (I2, S2, Z2)}}.
#'     Then, the argument \code{x1} is equal to
#'     \code{x = c(I1, I2, S1, S2, Z1, Z2)}.
#'
#' @export
#' @inheritParams locallycomp
#' @param x1 Vector of design (support) points of \eqn{\xi_1}. See 'Details' of \code{\link{leff}}.
#' @param w1 Vector of corresponding design weights for \code{x1}.
#' @param x2 Vector of design (support) points of the optimal design (\eqn{\xi_2}). Similar to \code{x1}.
#' @param w2 Vector of corresponding design weights for \code{x2}.
#' @param type A character. \code{"D"} denotes the D-efficiency and \code{"PA"} denotes the average P-efficiency.
#' @return A value between 0 and 1.
#' @references McGree, J. M., Eccleston, J. A., and Duffull, S. B. (2008). Compound optimal design criteria for nonlinear models. Journal of Biopharmaceutical Statistics, 18(4), 646-661.
#' @example inst/examples/leff_examples.R
leff <- function(formula,
                 predvars,
                 parvars,
                 family = gaussian(),
                 inipars,
                 type = c("D", "PA"),
                 fimfunc = NULL,
                 x2, w2, x1, w1,
                 npar = length(inipars),
                 prob = NULL){
  ## bayesian D-efficiency
  ### relative efficieny of x with respect to x2

  if (!(type[1] %in% c("D", "PA")))
    stop("'type' must be 'D' or 'PA'")
  npred <- length(x1)/length(w1)

  fimfunc_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                       predvars = predvars, parvars = parvars,
                                       family = family, lx =rep(0, npred), ux = rep(1, npred), iter = 1, k = 1,
                                       paramvectorized = FALSE,
                                       prior = NULL, x = NULL,
                                       user_crtfunc = NULL,
                                       user_sensfunc = NULL)


  if (!missing(formula)){
    if (length(inipars) != length(parvars))
      stop("length of 'inipars' is not equal to the length of 'parvars'")
  }

  if(missing(formula)){
    # to handle ...
    fimfunc2 <- function(x, w, param)
      fimfunc(x = x, w = w, param = param)
  } else{
    fimfunc2 <- fimfunc_formula$fimfunc_formula ## can be vectorized with respect to parameters!
  }
  if (type[1] == "D")
    releff <- (det2(fimfunc2(x = x1, w = w1, param = inipars))/
                 det2(fimfunc2(x = x2, w = w2, param = inipars)))^(1/npar)

  if (type[1] == "PA"){

    if (is.null(prob))
      stop("'prob' must be given for P-optimality")

    if (is.formula(prob)){
      prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
    }else{
      if (!is.function(prob))
        stop("'prob' must be either a function or a formula")
      if (!formalArgs(prob) %in% c("x", "param"))
        stop("arguments of 'prob' must be 'x' and 'param'")
    }
    releff <- sum(w1 * prob(x = x1, param = inipars))/sum(w2 * prob(x = x2, param = inipars))
  }
  return(releff)
}

######################################################################################################*
######################################################################################################*
#' @title  Calculates Relative Efficiency for Bayesian Optimal Designs
#' @description
#' Given a prior distribution for the parameters, this function calculates the Bayesian D-and PA- efficiency of a design \eqn{\xi_1} with respect to a design \eqn{\xi_2}.
#' Usually, \eqn{\xi_2} is an optimal design.
#' This function is especially useful for investigating the robustness of Bayesian optimal designs under different prior distributions (See 'Examples').
#'
#' @inheritParams leff
#' @inheritParams bayes
#' @details
#' See Masoudi et al. (2018) for formula details (the paper is under review and will be updated as soon as accepted).
#'
#' The argument  \code{x1} is the vector of design points.
#'  For design points with more than one dimension (the models with more than one predictors),
#'    it is a concatenation of the design points, but \strong{dimension-wise}.
#'    For example, let the model has three predictors   \eqn{(I, S, Z)}.
#'     Then,  a two-point optimal design has the following points:
#'    \eqn{\{\mbox{point1} = (I_1, S_1, Z_1), \mbox{point2} = (I_2, S_2, Z_2)\}}{{point1 = (I1, S1, Z1), point2 = (I2, S2, Z2)}}.
#'     Then, the argument \code{x} is equal to
#'     \code{x = c(I1, I2, S1, S2, Z1, Z2)}.
#'
#'
#' @export
#' @example inst/examples/beff_examples.R
beff <- function(formula,
                 predvars,
                 parvars,
                 family = gaussian(),
                 prior,
                 fimfunc = NULL,
                 x2, w2, x1, w1,
                 crt.bayes.control = list(),
                 npar = NULL,
                 type = c("D", "PA"),
                 prob = NULL){
  ## bayesian D-efficiency
  ### relative efficieny of x with respect to x2
  if (!(type[1] %in% c("D", "PA")))
    stop("'type' must be 'D' or 'PA'")
  if (is.null(npar)){
    if (!missing(formula))
      npar <- length(parvars)
    else
      npar <- prior$npar
  }

  ## only to pass the check_common_eargs
  npred <- length(x1)/length(w1)
  fimfunc_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                       predvars = predvars, parvars = parvars,
                                       family = family, lx =rep(0, npred), ux = rep(1,npred),
                                       iter = 1, k = length(w1),
                                       paramvectorized = FALSE, prior = prior,
                                       x = NULL,   user_crtfunc = NULL,
                                       user_sensfunc = NULL)
  if(missing(formula)){
    # to handle ...
    fimfunc2 <- function(x, w, param){
      fimfunc(x = x, w = w, param = param)
      #fimfunc(x = x, w = w, param = param,...)
    }
  } else{

    #if (length(prior$lower) != length(parvars))
    if (length(prior$lower) != fimfunc_formula$num_unknown_param)
      stop("length of 'prior$lower' is not equal to the number of unknown (not fixed) parameters")
    # fim_localdes <- fimfunc_formula$fimfunc_formula
    fimfunc2 <- fimfunc_formula$fimfunc_formula ## can be vectorized with respect to parameters!
  }

  # do not change its position
  if (!is.null(crt.bayes.control$method))
    if (crt.bayes.control$method == "quadrature")
      warning("The 'quadrature' method is not available for 'beff'. 'method' will be coerced to 'cubature'")

  crt.bayes.control <- do.call("crt.bayes.control", crt.bayes.control)
  if (type[1] == "D")
    cr_integrand <- function(param, x, w){
      # bcrfunc1 <- apply(param, 2,
      #                   FUN = function(col_par)-det2(fimfunc2(x = x, w = w, param = col_par), logarithm = TRUE)) * prior$fn(t(param))
      bcrfunc1 <- apply(param, 2,
                        FUN = function(col_par)-det2(fimfunc2(x = x, w = w, param = col_par), logarithm = TRUE)) * prior$fn(t(param))

      dim(bcrfunc1) <- c(1, length(bcrfunc1))
      return(bcrfunc1)
    }
  if (type[1] == "PA"){
    if (is.formula(prob)){
      prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
    }else{
      if (!is.function(prob))
        stop("'prob' must be either a function or a formula when type is 'PA'")
      if (!formalArgs(prob) %in% c("x", "param"))
        stop("arguments of 'prob' must be 'x' and 'param'")
    }
    cr_integrand <- function(param, x, w){
      bcrfunc1 <- (apply(param, 2, function(col_par) -log(sum(w * prob(x = x, param = col_par))))) * prior$fn(t(param))
      dim(bcrfunc1) <- c(1, length(bcrfunc1))
      return(bcrfunc1)
    }

  }
  crfunc_bayesian  <- function(x, w, maxEval, tol) {
    out <- hcubature(f = cr_integrand, lowerLimit = prior$lower,
                     upperLimit = prior$upper,
                     vectorInterface = TRUE,
                     x = x, w = w, tol = tol, maxEval = maxEval)

    val <- out$integral
    return(list(val = val, fneval = out$functionEvaluations))
  }

  if (type[1] == "D")
    releff<- exp((crfunc_bayesian(x = x2, w = w2, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val-
                    crfunc_bayesian(x = x1, w = w1, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val)/npar)
  # if (type[1] == "PA")
  #   releff<-  (crfunc_bayesian(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val/
  #                   crfunc_bayesian(x = x2, w = w2, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val)
  if (type[1] == "PA")
    releff <- exp(crfunc_bayesian(x = x2, w = w2, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val - crfunc_bayesian(x = x1, w = w1, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val)

  # releff<- crfunc_bayesian_D(x = x2, w = w2, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val/
  #                crfunc_bayesian_D(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val

  return(releff)
}
######################################################################################################*
######################################################################################################*
######################################################################################################*
######################################################################################################*
#'@title Calculates Relative Efficiency for Minimax Optimal Designs
#' @description
#' Given a parameter space for the unknown parameters, this function calculates the D-efficiency of a design \eqn{\xi_1} with respect to a design \eqn{\xi_2}.
#' Usually, \eqn{\xi_2} is an  optimal design.
#'
#'
#' @details
#' See Masoudi et al. (2017) for formula details.
#'
#' The argument  \code{x1} is the vector of design points.
#'  For design points with more than one dimension (the models with more than one predictors),
#'    it is a concatenation of the design points, but \strong{dimension-wise}.
#'    For example, let the model has three predictors   \eqn{(I, S, Z)}.
#'     Then,  a two-point optimal design has the following points:
#'    \eqn{\{\mbox{point1} = (I_1, S_1, Z_1), \mbox{point2} = (I_2, S_2, Z_2)\}}{{point1 = (I1, S1, Z1), point2 = (I2, S2, Z2)}}.
#'     Then, the argument \code{x} is equal to
#'     \code{x = c(I1, I2, S1, S2, Z1, Z2)}.
#'
#' @export
#' @inheritParams minimax
#' @param x1 Vector of design (support) points of \eqn{\xi_1}. See 'Details' of \code{\link{leff}}.
#' @param w1 Vector of corresponding design weights for \code{x}.
#' @param x2 Vector of design (support) points of the optimal design (\eqn{\xi_2}). Similar to \code{x}.
#' @param w2 Vector of corresponding design weights for \code{x2}.
#' @return A value between 0 and 1.
#' @example inst/examples/meff_examples.R
meff <- function(formula,
                 predvars,
                 parvars,
                 family = gaussian(),
                 lp, up,
                 fimfunc = NULL,
                 x2, w2, x1, w1,
                 standardized = FALSE,
                 localdes = NULL,
                 crt.minimax.control = list(),
                 npar = length(lp)){
  ## minimax D-efficiency
  ### relative efficieny of x with respect to x2
  if (!is.logical(standardized))
    stop("'standardized' must be logical")


  npred <- length(x1)/length(w1)
  if (standardized){
    if (!missing(formula))
      localdes_par <- create_localdes(parvars = parvars, localdes = localdes) else
        localdes_par <- localdes
  }else
    localdes_par <- NULL

  fimfunc_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                       predvars = predvars, parvars = parvars,
                                       family = family, lx =rep(0, npred), ux = rep(1, npred),
                                       iter = 1, k = 1,
                                       paramvectorized = FALSE,
                                       prior = NULL, x = NULL,
                                       user_crtfunc = NULL,
                                       user_sensfunc = NULL)

  if (!missing(formula)){
    if (length(lp) != length(parvars))
      stop("length of 'lp' is not equal to the length of 'parvars'")
  }

  if(missing(formula)){
    # to handle ...
    fimfunc2 <- function(x, w, param)
      fimfunc(x = x, w = w, param = param)
    #fimfunc(x = x, w = w, param = param,...)
  } else{
    fimfunc2 <- fimfunc_formula$fimfunc_formula ## can be vectorized with respect to parameters!
  }
  crt.minimax.control <- do.call("crt.minimax.control", crt.minimax.control)
  if (!standardized){
    crfunc <- function(x, w, param){
      minimax_crfunc <-  -det2(fimfunc2(x = x, w = w, param = param), logarithm = TRUE) + 5000 * (sum(w) - 1) ^ 2
      return(minimax_crfunc)}
  }else{
    crfunc <- function(x, w, param) {
      localdes_res <- localdes_par(param = param)
      denominator <- det2(fimfunc2(x = localdes_res$x, w = localdes_res$w, param = param), logarithm = FALSE)
      numerator <- det2(fimfunc2(x = x, w = w, param = param), logarithm = FALSE)
      eff <- (numerator/denominator)
      if (is.nan(eff))
        stop("The criterion (D-efficiency) value is 'NaN' for ", paste("iniparam = c(", paste(round(param, 5), collapse = ","), ").",sep = ""),
             "\nCheck the Fisher information matrix, number of design points, lx, ux and .... for possible unmatched design parameters.")
      if (round(eff, 7) > 1)
        stop("Efficiency value ", eff,
             " is greater than one. This results in a wrong conclusion that the non-optimal design is better than the true optimal design when ",
             paste("iniparam = c(", paste(round(param, 5), collapse = ","), ").",sep = ""),
             "\n  Probably 'localdes' does not return the true locally optimal designfor iniparam.")

      if (npar %% 2 != 0) {
        eff <- (eff)^(1/npar)
      }else{
        eff <-  ifelse(eff < 0, 0,(eff)^(1/npar))
      }
      return(-eff)
    }
  }
  any_fixed <- sapply(1:length(lp), function(i) lp [i] == up[i]) # is a vector
  if (any(any_fixed)){
    is_fixed <- TRUE
    fixedpar_id <- which(any_fixed)
    fixedpar <- lp[fixedpar_id]
    lp_nofixed <- lp[-fixedpar_id]
    up_nofixed <- up[-fixedpar_id]
    ## warning: lp and up are channged here if 'any_fix == TRUE'
  }else{
    fixedpar <- NA
    fixedpar_id <- NA
    is_fixed <- FALSE
    lp_nofixed <- lp
    up_nofixed <- up
  }
  if (is_fixed){
    crfunc2 <- function(param, x, w, fixedpar = NA, fixedpar_id = NA){
      param1 <- insert_fixed(param = param, fixedpar = fixedpar, fixedpar_id = fixedpar_id)
      out <- crfunc(param = param1,  x = x, w = w)
      return(out)
    }
  }else{
    crfunc2 <- function(param, x, w, fixedpar = NA, fixedpar_id = NA){
      out <- crfunc(param = param, x = x, w = w)
      return(out)
    }
  }

  vertices_inner <- make_vertices(lower = lp_nofixed, upper = up_nofixed)
  optim_nloptr <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id, standardized){
    if (is.null(crt.minimax.control$x0))
      x0 <- (lower + upper)/2 else
        x0 <- crt.minimax.control$x0
      out_nloptr <- nloptr::nloptr(x0= x0, eval_f = fn, lb = lower, ub = upper,
                                   opts = crt.minimax.control$optslist,
                                   x = x, w = w, fixedpar = fixedpar,
                                   fixedpar_id = fixedpar_id)
      out <- find_on_points(fn = fn,
                            points = vertices_inner,
                            x = x, w = w,
                            fixedpar = fixedpar,
                            fixedpar_id = fixedpar_id)
      minima <- out$minima
      minima_nloptr <- c(out_nloptr$solution, out_nloptr$objective)
      #minima_nloptr <- c(out_nloptr$par, out_nloptr$value)
      minima <- rbind(minima, minima_nloptr)
      if (standardized)
        crt_val <- -min(minima[, dim(minima)[2]]) else
          crt_val <- max(minima[, dim(minima)[2]])
      return(crt_val)
  }

  crt_x1 <- optim_nloptr(fn = crfunc2, lower = lp_nofixed,
                        upper = up_nofixed, x = x1, w = w1,
                        fixedpar_id = fixedpar_id, fixedpar = fixedpar,
                        standardized = standardized)
  # optimal
  crt_x2 <- optim_nloptr(fn = crfunc2, lower = lp_nofixed,
                           upper = up_nofixed, x = x2, w = w2,
                           fixedpar_id = fixedpar_id, fixedpar = fixedpar,
                           standardized = standardized)
  # (crt_x1/crt_x2)^(1/npar) when log = FALSE
  #0.011237
  if (!standardized)
    releff <- exp((crt_x2 - crt_x1)/npar) else
      releff <- crt_x1/crt_x2

  return(releff)
}




######################################################################################################*
######################################################################################################*

######################################################################################################*
######################################################################################################*
