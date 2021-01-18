#' A function to predict N and P excretion, CNP egestion, CNP ingestion rate, using MCMC and stan
#'
#' This function combines MTE and stoichiometric theory in order to predict nescessary ingestion and excretion processes. A probability distribution is obtained by including uncertainty of parameters and using MCMC sampling with stan.
#'
#' @param TL Total length(s) in cm
#' @param param List of all parameter means (add "_m") and standard deviations (add "_sd") Default parameters are set with very low sd's.
#' parameters:
#' \itemize{
#' \item{Qc_m, Qc_sd:} percentage C of dry mass fish
#' \item{Qn_m, Qn_sd:} percentage N of dry mass fish
#' \item{Qp_m, Qp_sd:} percentage P of dry mass fish
#' \item{Dc_m, Dc_sd:} percentage C of dry mass food
#' \item{Dn_m, Dn_sd:} percentage N of dry mass food
#' \item{Dp_m, Dp_sd:} percentage P of dry mass food
#' \item{ac_m, ac_sd:} C-specific assimilation efficiency
#' \item{an_m, an_sd:} N-specific assimilation efficiency
#' \item{ap_m, ap_sd:} P-specific assimilation efficiency
#' \item{linf_m, linf_sd:} Von Bertalanffy Growth parameter, theoretical maximum size in TL (cm)
#' \item{k_m, k_sd:} Von Bertalanffy Growth parameter, growth rate (yr^-1)
#' \item{t0_m, tO_sd:} Von Bertalanffy Growth parameter (yr)
#' \item{lwa_m, lwa_sd:} Parameter length-weight relationship (g cm^-1)
#' \item{lwb_m, lwb_sd:} Parameter length-weight relationship
#' \item{mdw_m, wprop_sd:} Ratio between dry weight and wet weight of fish
#' \item{F0nz_m, F0nz_sd:} N-specific turnover rate
#' \item{F0pz_m, F0pz_sd:} P-specific turnover rate
#' \item{f0_m, f0_sd:} Metabolic normalisation constant independent of body mass (g C g^-alpha d^-1)
#' \item{alpha_m, alpha_sd:} Metabolic rate mass-scaling exponent
#' \item{theta_m, theta_sd:} Activity scope
#' \item{r_m, r_sd:} Aspect ratio of caudal fin
#' \item{h_m, h_sd:} Trophic level
#' \item{v_m, v_sd:} Environmental temperature (degrees celcius)
#' }
#' @param cor     A list of correlations between certain parameters: ro_Qc_Qn, ro_Qc_Qp, ro_Qn_Qp,
#' ro_Dc_Dn, ro_Dc_Dp, ro_Dn_Dp, ro_lwa_lwb, ro_alpha_f0
#' @param iter    A positive integer specifying the number of iterations. The default is 2000.
#' @param ... Additional arguments rstan::sampling, see ?rstan:sampling
#'
#' @details Returns a list with two objects: A stanfit object and a data.frame with a summary of all model components.
#'  See \code{\link{extract}} to extract a summary of predicted variables and
#'  \code{\link{limitation}} to get information on the limiting element.
#' 
#' @keywords fish stoichiometry excretion mcmc
#' @importFrom stats median quantile sd
#' @importFrom parallel mclapply
#' @importFrom plyr ldply
#'
#' @examples
#' library(fishflux)
#' model <- cnp_model_mcmc(TL = 10, param = list(
#' Qc_m = 40, Qn_m = 10, Qp_m = 4, theta_m = 3))
#' 
#' @export
cnp_model_mcmc <- function(TL, param, iter = 1000,
                           cor = list(ro_Qc_Qn = 0.5, ro_Qc_Qp = -0.3, ro_Qn_Qp = -0.2,
                                      ro_Dc_Dn = 0.2, ro_Dc_Dp = -0.1, ro_Dn_Dp = -0.1,
                                      ro_lwa_lwb = 0.9, ro_alpha_f0 = 0.9), ...) {

  ##standard parameters, all sd's are quite low here!
  params_st <- list(lt_m = 10,
                    ac_m = 0.8,
                    an_m = 0.8,
                    ap_m = 0.7,
                    Dc_m = 2.5,
                    Dn_m = 0.3,
                    Dp_m = 0.1,
                    linf_m = 20,
                    k_m = 0.4,
                    t0_m = 0,
                    theta_m = 2,
                    r_m = 2,
                    h_m = 2,
                    lwa_m = 0.0137,
                    lwb_m = 3.083,
                    mdw_m = 0.309,
                    v_m = 27,
                    F0nz_m = 0.01,
                    F0pz_m = 0.0007,
                    Qc_m = 40,
                    Qn_m = 10,
                    Qp_m = 4,
                    alpha_m = 0.8,
                    f0_m = 0.002,

                    lt_sd = 0.0000000001,
                    ac_sd = 0.0000000001,
                    an_sd = 0.0000000001,
                    ap_sd = 0.0000000001,
                    Dc_sd = 0.0000000001,
                    Dn_sd = 0.0000000001,
                    Dp_sd = 0.0000000001,
                    linf_sd = 0.0000000001,
                    k_sd = 0.0000000001,
                    t0_sd = 0.0000000001,
                    theta_sd = 0.0000000001,
                    r_sd = 0.0000000001,
                    h_sd = 0.0000000001,
                    lwa_sd = 0.0000000001,
                    lwb_sd = 0.0000000001,
                    mdw_sd = 0.0000000001,
                    v_sd = 0.0000000001,
                    F0nz_sd = 0.0000000001,
                    F0pz_sd = 0.0000000001,
                    Qc_sd = 0.0000000001,
                    Qn_sd = 0.0000000001,
                    Qp_sd = 0.0000000001,
                    alpha_sd = 0.0000000001,
                    f0_sd = 0.0000000001
  )

  #check input variable names
  if (TRUE %in% (!names(param) %in% names(params_st))) {
    wrong <- names(param)[!(names(param) %in% names(params_st))]
    error <- paste("The following input parameters do not exist: ", paste(wrong, collapse = ", "),
                   "  Check ?cnp_model_mcmc for a description of valid input parameters")
    stop(error)
  }

  if (missing(TL)) {
    stop("please provide TL: total length")
  }

  if ("Qc_m" %in% names(param)) {
    if (param$Qc_m <= 0 | param$Qc_m >= 100) {
      stop("Qc_m should be between 0 and 100")
   }
  }

  if ("Qn_m" %in% names(param)) {
    if (param$Qn_m <= 0 | param$Qn_m >= 100) {
      stop("Qn_m should be between 0 and 100")
    }
  }

  if ("Qp_m" %in% names(param)) {
    if (param$Qp_m <= 0 | param$Qp_m >= 100) {
      stop("Qp_m should be between 0 and 100")
    }
  }

  if ("Dc_m" %in% names(param)) {
    if (param$Dc_m <= 0 | param$Dc_m >= 100) {
      stop("Dc_m should be between 0 and 100")
    }
  }

  if ("Dn_m" %in% names(param)) {
    if (param$Dn_m <= 0 | param$Dn_m >= 100) {
      stop("Dn_m should be between 0 and 100")
    }
  }

  if ("Dp_m" %in% names(param)) {
    if (param$Dp_m <= 0 | param$Dp_m >= 100) {
      stop("Dp_m should be between 0 and 100")
    }
  }

  if ("ac_m" %in% names(param)) {
    if (param$ac_m <= 0 | param$ac_m >= 1) {
      stop("ac_m should be between 0 and 1")
    }
  }

  if ("Dp_m" %in% names(param)) {
    if (param$an_m <= 0 | param$an_m >= 1) {
      stop("an_m should be between 0 and 1")
    }
  }

  if ("ap_m" %in% names(param)) {
    if (param$ap_m <= 0 | param$ap_m >= 1) {
      stop("ap_m should be between 0 and 1")
    }
  }

  if (length(TL) == 1) { ## option for only one length ##
    result <- cnp_mcmc(TL, param, iter, params_st, cor, ...)
    list(stanfit = result[[1]], summary = result[[2]])
  } else { ## option for vector of lengths ##
    result <- mclapply(TL, FUN = cnp_mcmc, param = param,
                       iter = iter, params_st = params_st,
                       cor = cor, ...)
    stanfit <- lapply(result, FUN = function(x) {x[[1]]})
    summary <- lapply(result, FUN = function(x) {x[[2]]})
    summary <- ldply(summary)
    list(stanfit = stanfit, summary = summary)
  }
}

#' cnp_mcmc
#' 
#' @inheritParams cnp_model_mcmc
#' 
#' @param params_st Standard parameters.
#' @param ... Additional arguments rstan::sampling, see ?rstan:sampling
#' 
#' @importFrom rstan sampling extract
#' @importFrom plyr ldply
cnp_mcmc <- function(TL, param, iter, params_st, cor, ...) {
  ## add TL to parameter list
  param[["lt_m"]] <- TL

  p_given <- names(param)
  p_all <- names(params_st)
  unknown <- p_all[!p_all %in% p_given]

  if (length(unknown > 0)) {
    warning("not inputting certain parameters may give wrong results")
    for (v in unknown) {
      warning("adding standard values for ", v)
    }
  }

  params_missing <- params_st[which(p_all %in% unknown)]
  param <- append(param, params_missing)
  param <- append(param, cor)

  ##add others plus replace Qcnp
  stanfit <-  sampling(stanmodels$cnpmodelmcmc, data = param,
                       iter = iter, algorithm = "Fixed_param", chains = 1,
                       ...)

  ee <- rstan::extract(stanfit)
  par <- names(ee)
  
  result <-
    ldply(lapply(par, function(x, ee) {
        data.frame(          
          TL = mean(ee[["lt"]]),
          variable = x,
          mean = mean(ee[[x]]),
          median = median(ee[[x]]),
          se = sd(ee[[x]]) / sqrt(length(ee[[x]])),
          sd = sd(ee[[x]]),
          Q_2.5 = quantile(ee[[x]], 0.025),
          Q_97.5 = quantile(ee[[x]], 0.975),
          Q_25 = quantile(ee[[x]], 0.25),
          Q_75 = quantile(ee[[x]], 0.75))
    }, ee = ee))
  list(stanfit, summary = result)
}
