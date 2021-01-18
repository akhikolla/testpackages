#' @title Bayesian Discount Prior: Survival Analysis
#' @description \code{bdpsurvival} is used to estimate the survival probability
#'   (single arm trial; OPC) or hazard ratio (two-arm trial; RCT) for
#'   right-censored data using the survival analysis implementation of the
#'   Bayesian discount prior. In the current implementation, a two-arm analysis
#'   requires all of current treatment, current control, historical treatment,
#'   and historical control data. This code is modeled after
#'   the methodologies developed in Haddad et al. (2017).
#' @param formula an object of class "formula." Must have a survival object on
#'   the left side and at most one input on the right side, treatment. See
#'   "Details" for more information.
#' @param data a data frame containing the current data variables in the model.
#'   Columns denoting 'time' and 'status' must be present. See "Details" for required
#'   structure.
#' @param data0 optional. A data frame containing the historical data variables in the model.
#'   If present, the column labels of data and data0 must match.
#' @param breaks vector. Breaks (interval starts) used to compose the breaks of the
#'   piecewise exponential model. Do not include zero. Default breaks are the
#'   quantiles of the input times.
#' @param a0 scalar. Prior value for the gamma shape of the piecewise
#'   exponential hazards. Default is 0.1.
#' @param b0 scalar. Prior value for the gamma rate of the piecewise
#'   exponential hazards. Default is 0.1.
#' @param surv_time scalar. Survival time of interest for computing the
#'   probability of survival for a single arm (OPC) trial. Default is
#'   overall, i.e., current+historical, median survival time.
#' @param discount_function character. Specify the discount function to use.
#'   Currently supports \code{weibull}, \code{scaledweibull}, and
#'   \code{identity}. The discount function \code{scaledweibull} scales
#'   the output of the Weibull CDF to have a max value of 1. The \code{identity}
#'   discount function uses the posterior probability directly as the discount
#'   weight. Default value is "\code{identity}".
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param fix_alpha logical. Fix alpha at alpha_max? Default value is FALSE.
#' @param number_mcmc scalar. Number of Monte Carlo simulations. Default is 10000.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. For a two-arm trial, users may specify a vector of two
#'   values where the first value is used to estimate the weight of the
#'   historical treatment group and the second value is used to estimate the
#'   weight of the historical control group.
#' @param method character. Analysis method with respect to estimation of the weight
#'   paramter alpha. Default method "\code{mc}" estimates alpha for each
#'   Monte Carlo iteration. Alternate value "\code{fixed}" estimates alpha once
#'   and holds it fixed throughout the analysis.  See the the
#'   \code{bdpsurvival} vignette \cr
#'   \code{vignette("bdpsurvival-vignette", package="bayesDP")} for more details.
#' @param compare logical. Should a comparison object be included in the fit?
#'   For a one-arm analysis, the comparison object is simply the posterior
#'   chain of the treatment group parameter. For a two-arm analysis, the comparison
#'   object is the posterior chain of the treatment effect that compares treatment and
#'   control. If \code{compare=TRUE}, the comparison object is accessible in the
#'   \code{final} slot, else the \code{final} slot is \code{NULL}. Default is
#'   \code{TRUE}.
#' @details \code{bdpsurvival} uses a two-stage approach for determining the
#'   strength of historical data in estimation of a survival probability outcome.
#'   In the first stage, a \emph{discount function} is used that
#'   that defines the maximum strength of the
#'   historical data and discounts based on disagreement with the current data.
#'   Disagreement between current and historical data is determined by stochastically
#'   comparing the respective posterior distributions under noninformative priors.
#'   With a single arm survival data analysis, the comparison is the
#'   probability (\code{p}) that the current survival is less than the historical
#'   survival. For a two-arm survival data, analysis the comparison is the
#'   probability that the hazard ratio comparing treatment and control is
#'   different from zero. The comparison metric \code{p} is then
#'   input into the discount function and the final strength of the
#'   historical data is returned (alpha).
#'
#'   In the second stage, posterior estimation is performed where the discount
#'   function parameter, \code{alpha}, is used incorporated in all posterior
#'   estimation procedures.
#'
#'   To carry out a single arm (OPC) analysis, data for the current and
#'   historical treatments are specified in separate data frames, data and data0,
#'   respectively. The data frames must have matching columns denoting time and status.
#'   The 'time' column is the survival (censor) time of the event and the 'status' column
#'   is the event indicator. The results are then based on the posterior probability of
#'   survival at \code{surv_time} for the current data augmented by the historical data.
#'
#'   Two-arm (RCT) analyses are specified similarly to a single arm trial. Again
#'   the input data frames must have columns denoting time and status, but now
#'   an additional column named 'treatment' is required to denote treatment and control
#'   data. The 'treatment' column must use 0 to indicate the control group. The current data
#'   are augmented by historical data (if present) and the results are then based
#'   on the posterior distribution of the hazard ratio between the treatment
#'   and control groups.
#'
#'   For more details, see the \code{bdpsurvival} vignette: \cr
#'   \code{vignette("bdpsurvival-vignette", package="bayesDP")}
#'
#' @return \code{bdpsurvival} returns an object of class "bdpsurvival".
#' The functions \code{\link[=summary,bdpsurvival-method]{summary}} and \code{\link[=print,bdpsurvival-method]{print}} are used to obtain and
#' print a summary of the results, including user inputs. The \code{\link[=plot,bdpsurvival-method]{plot}}
#' function displays visual outputs as well.
#'
#' An object of class "\code{bdpsurvival}" is a list containing at least
#' the following components:
#' \describe{
#'  \item{\code{posterior_treatment}}{
#'    list. Entries contain values related to the treatment group:}
#'    \itemize{
#'      \item{\code{alpha_discount}}{
#'        numeric. Alpha value, the weighting parameter of the historical data.}
#'      \item{\code{p_hat}}{
#'        numeric. The posterior probability of the stochastic comparison
#'        between the current and historical data.}
#'      \item{\code{posterior_survival}}{
#'        vector. If one-arm trial, a vector of length \code{number_mcmc}
#'        containing the posterior probability draws of survival at
#'        \code{surv_time}.}
#'      \item{\code{posterior_flat_survival}}{
#'        vector. If one-arm trial, a vector of length \code{number_mcmc}
#'        containing the probability draws of survival at \code{surv_time}
#'        for the current treatment not augmented by historical treatment.}
#'      \item{\code{prior_survival}}{
#'        vector. If one-arm trial, a vector of length \code{number_mcmc}
#'        containing the probability draws of survival at \code{surv_time}
#'        for the historical treatment.}
#'      \item{\code{posterior_hazard}}{
#'        matrix. A matrix with \code{number_mcmc} rows and \code{length(breaks)}
#'        columns containing the posterior draws of the piecewise hazards
#'        for each interval break point.}
#'      \item{\code{posterior_flat_hazard}}{
#'        matrix. A matrix with \code{number_mcmc} rows and \code{length(breaks)}
#'        columns containing the draws of piecewise hazards for each interval
#'        break point for current treatment not augmented by historical treatment.}
#'      \item{\code{prior_hazard}}{
#'        matrix. A matrix with \code{number_mcmc} rows and \code{length(breaks)}
#'        columns containing the draws of piecewise hazards for each interval break point
#'        for historical treatment.}
#'   }
#'  \item{\code{posterior_control}}{
#'    list. If two-arm trial, contains values related to the control group
#'    analagous to the \code{posterior_treatment} output.}
#'
#'  \item{\code{final}}{
#'    list. Contains the final comparison object, dependent on the analysis type:}
#'    \itemize{
#'      \item{One-arm analysis:}{
#'        vector. Posterior chain of survival probability at requested time.}
#'      \item{Two-arm analysis:}{
#'        vector. Posterior chain of log-hazard rate comparing treatment and control groups.}
#'   }
#'
#'  \item{\code{args1}}{
#'    list. Entries contain user inputs. In addition, the following elements
#'    are ouput:}
#'    \itemize{
#'      \item{\code{S_t}, \code{S_c}, \code{S0_t}, \code{S0_c}}{
#'        survival objects. Used internally to pass survival data between
#'        functions.}
#'      \item{\code{arm2}}{
#'        logical. Used internally to indicate one-arm or two-arm analysis.}
#'   }
#' }
#'
#' @seealso \code{\link[=summary,bdpsurvival-method]{summary}},
#'   \code{\link[=print,bdpsurvival-method]{print}},
#'   and \code{\link[=plot,bdpsurvival-method]{plot}} for details of each of the
#'   supported methods.
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' # One-arm trial (OPC) example - survival probability at 5 years
#'
#' # Collect data into data frames
#' df_ <- data.frame(status = rexp(50, rate=1/30),
#'                   time   = rexp(50, rate=1/20))
#' df_$status <- ifelse(df_$time < df_$status, 1, 0)
#'
#' df0 <- data.frame(status = rexp(50, rate=1/30),
#'                   time   = rexp(50, rate=1/10))
#' df0$status <- ifelse(df0$time < df0$status, 1, 0)
#'
#'
#' fit1 <- bdpsurvival(Surv(time, status) ~ 1,
#'                     data  = df_,
#'                     data0 = df0,
#'                     surv_time = 5,
#'                     method = "fixed")
#'
#' print(fit1)
#' \dontrun{
#' plot(fit1)
#' }
#'
#' # Two-arm trial example
#' # Collect data into data frames
#' df_ <- data.frame(time = c(rexp(50, rate=1/20),  # Current treatment
#'                            rexp(50, rate=1/10)), # Current control
#'                   status = rexp(100, rate=1/40),
#'                   treatment = c(rep(1,50), rep(0,50)))
#' df_$status <- ifelse(df_$time < df_$status, 1, 0)
#'
#' df0 <- data.frame(time = c(rexp(50, rate=1/30),  # Historical treatment
#'                            rexp(50, rate=1/5)),  # Historical control
#'                   status =  rexp(100, rate=1/40),
#'                   treatment = c(rep(1,50), rep(0,50)))
#' df0$status <- ifelse(df0$time < df0$status, 1, 0)
#'
#' fit2 <- bdpsurvival(Surv(time, status) ~ treatment,
#'                     data = df_,
#'                     data0 = df0,
#'                     method = "fixed")
#' summary(fit2)
#'
#' ### Fix alpha at 1
#' fit2_1 <- bdpsurvival(Surv(time, status) ~ treatment,
#'                       data = df_,
#'                       data0 = df0,
#'                       fix_alpha = TRUE,
#'                       method = "fixed")
#' summary(fit2_1)
#'
#'
#' @rdname bdpsurvival
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset
#'   model.response pweibull quantile rbeta rgamma rnorm var vcov update
#' @importFrom survival Surv survSplit
#' @aliases bdpsurvival,ANY-method
#' @export bdpsurvival
bdpsurvival <- setClass("bdpsurvival", slots = c(posterior_treatment = "list",
                                                 posterior_control = "list",
                                                 final = "list",
                                                 args1 = "list"))
setGeneric("bdpsurvival",
  function(formula           = formula,
           data              = data,
           data0             = NULL,
           breaks            = NULL,
           a0                = 0.1,
           b0                = 0.1,
           surv_time         = NULL,
           discount_function = "identity",
           alpha_max         = 1,
           fix_alpha         = FALSE,
           number_mcmc       = 10000,
           weibull_scale     = 0.135,
           weibull_shape     = 3,
           method            = "mc",
           compare           = TRUE){
             standardGeneric("bdpsurvival")
           })

setMethod("bdpsurvival",
  signature(),
  function(formula           = formula,
           data              = data,
           data0             = NULL,
           breaks            = NULL,
           a0                = 0.1,
           b0                = 0.1,
           surv_time         = NULL,
           discount_function = "identity",
           alpha_max         = 1,
           fix_alpha         = FALSE,
           number_mcmc       = 10000,
           weibull_scale     = 0.135,
           weibull_shape     = 3,
           method            = "mc",
           compare           = TRUE){


  ### Check validity of data input
  call <- match.call()
  if (missing(data)) {
    stop("Current data not input correctly.")
  }


  ##############################################################################
  ### Parse current data
  ##############################################################################
  ### Check data frame and ensure it has the correct column names
  mf <- mf0 <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$na.action <- NULL
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) {
      names(Y) <- nm
    }
  }

  if(class(Y) != "Surv") stop("Data input incorrectly.")



  ##############################################################################
  ### Parse historical data
  ##############################################################################
  ### Parse historical data
  if(!is.null(data0)){
    m00 <- match("data0", names(mf0), 0L)
    m01 <- match("data", names(mf0), 0L)
    mf0[m01] <- mf0[m00]
    m0  <- match(c("formula", "data"), names(mf0), 0L)
    mf0 <- mf0[c(1L, m0)]
    mf0$na.action <- NULL
    mf0[[1L]] <- quote(stats::model.frame)
    mf0 <- eval(mf0, parent.frame())
    mt0 <- attr(mf0, "terms")
    Y0 <- model.response(mf0, "any")
    if (length(dim(Y0)) == 1L) {
      nm0 <- rownames(Y0)
      dim(Y0) <- NULL
      if (!is.null(nm0)) {
        names(Y0) <- nm0
      }
    }
  } else{
    Y0 <- NULL
  }


  # Check that discount_function is input correctly
  all_functions <- c("weibull", "scaledweibull", "identity")
  function_match <- match(discount_function, all_functions)
  if(is.na(function_match)) {
    stop("discount_function input incorrectly.")
  }


  historical <- NULL
  treatment <- NULL



  ##############################################################################
  # Quick check, if alpha_max, weibull_scale, or weibull_shape have length 1,
  # repeat input twice
  ##############################################################################
  if(length(alpha_max)==1){
    alpha_max <- rep(alpha_max, 2)
  }

  if(length(weibull_scale)==1){
    weibull_scale <- rep(weibull_scale, 2)
  }

  if(length(weibull_shape)==1){
    weibull_shape <- rep(weibull_shape, 2)
  }


  ##############################################################################
  # Format input data
  ##############################################################################
  ### If no breaks input, create intervals along quantiles
  if(is.null(breaks)){
     breaks <- quantile(c(Y[,1], Y0[,1]),probs=c(0.2,0.4,0.6,0.8))
  }


  ### If zero is present in breaks, remove and give warning
  if(any(breaks==0)){
    breaks <- breaks[!(breaks==0)]
    warning("Breaks vector included 0. The zero value was removed.")
  }


  ### Combine current and historical data, and format for analysis
  dataCurrent <- data
  dataCurrent$historical <- 0

  if(!is.null(data0)){
    dataHistorical            <- data0
    dataHistorical$historical <- 1
  } else{
    dataHistorical <- NULL
  }

  dataALL <- rbind(dataCurrent, dataHistorical)


  ### Update formula
  formula <- update(formula, ~ . + historical)

  ### Split the data on the breaks
  dataSplit <- survSplit(formula,
                         cut     = breaks,
                         start   = "start",
                         episode = "interval",
                         data    = dataALL)

  ### Change time and status column names
  vars <- as.character(attr(mt, "variables"))[2]
  vars <- strsplit(vars, "Surv\\(|, |\\)")[[1]]
  var_time   <- vars[2]
  var_status <- vars[3]

  names(dataSplit)[match(var_time, names(dataSplit))] <- "time"
  names(dataSplit)[match(var_status, names(dataSplit))] <- "status"


  # Grab fu-time column
  id_time <- names(dataSplit[ncol(dataSplit)-2])

  # Look for treatment column, if missing, add it and set to 1
  if(!any(names(dataSplit) == "treatment")){
    dataSplit$treatment <- 1
  }

  ### Compute exposure time within each interval
  dataSplit$exposure <- dataSplit$time - dataSplit$start

  ### Create new labels for the intervals
  maxTime  <- max(dataSplit$time)
  labels_t <- unique(c(breaks, maxTime))
  dataSplit$interval <- factor(dataSplit$interval,
                               labels = labels_t)

  ### Parse out the historical and current data
  S_t  <- subset(dataSplit, historical==0 & treatment == 1)
  S_c  <- subset(dataSplit, historical==0 & treatment == 0)
  S0_t <- subset(dataSplit, historical==1 & treatment == 1)
  S0_c <- subset(dataSplit, historical==1 & treatment == 0)

  if(nrow(S0_t) == 0) S0_t <- NULL
  if(nrow(S_c) == 0)  S_c  <- NULL
  if(nrow(S0_c) == 0) S0_c <- NULL


  ### Compute arm2, internal indicator of a two-arm trial
  if(is.null(S_c) & is.null(S0_c)){
    arm2 <- FALSE
  } else{
    arm2 <- TRUE
  }

  #if(arm2) stop("Two arm trials are not currently supported.")

  ### If surv_time is null, replace with median time
  if(is.null(surv_time) & !arm2){
    surv_time <- median(c(Y[,1], Y[,0]))
  }


  ### Check inputs
  if(!arm2){
    if(nrow(S_t) == 0) stop("Current treatment data missing or input incorrectly.")
    if(is.null(S0_t)) warning("Historical treatment data missing or input incorrectly.")
  } else if(arm2){
    if(nrow(S_t) == 0) stop("Current treatment data missing or input incorrectly.")
    if(is.null(S_c)) warning("Current control data missing or input incorrectly.")
    if(is.null(S0_t) & is.null(S0_c)) warning("Historical data input incorrectly.")
  }

  posterior_treatment <- posterior_survival(
    S                 = S_t,
    S0                = S0_t,
    surv_time         = surv_time,
    discount_function = discount_function,
    alpha_max         = alpha_max[1],
    fix_alpha         = fix_alpha,
    a0                = a0,
    b0                = b0,
    number_mcmc       = number_mcmc,
    weibull_shape     = weibull_shape[1],
    weibull_scale     = weibull_scale[1],
    breaks            = breaks,
    arm2              = arm2,
    method            = method)

  if(arm2){
    posterior_control <- posterior_survival(
      S                 = S_c,
      S0                = S0_c,
      discount_function = discount_function,
      alpha_max         = alpha_max[2],
      fix_alpha         = fix_alpha,
      a0                = a0,
      b0                = b0,
      surv_time         = surv_time,
      number_mcmc       = number_mcmc,
      weibull_shape     = weibull_shape[2],
      weibull_scale     = weibull_scale[2],
      breaks            = breaks,
      arm2              = arm2,
      method            = method)
  } else{
    posterior_control <- NULL
  }


  args1 <- list(S_t               = S_t,
                S_c               = S_c,
                S0_t              = S0_t,
                S0_c              = S0_c,
                discount_function = discount_function,
                alpha_max         = alpha_max,
                fix_alpha         = fix_alpha,
                a0                = a0,
                b0                = b0,
                surv_time         = surv_time,
                number_mcmc       = number_mcmc,
                weibull_scale     = weibull_scale,
                weibull_shape     = weibull_shape,
                method            = method,
                arm2              = arm2,
                breaks            = breaks,
                data              = dataSplit,
                data_current      = data)


  ##############################################################################
  ### Create final (comparison) object
  ##############################################################################
  if(!compare){
    final <- NULL
  } else{
    if(arm2){
      R0      <- log(posterior_treatment$posterior_hazard)-log(posterior_control$posterior_hazard)
      V0      <- 1/apply(R0,2,var)
      logHR0  <- R0%*%V0/sum(V0)
      final   <- list()
      final$posterior_loghazard <- logHR0
    } else{
      final                    <- list()
      final$posterior_survival <- posterior_treatment$posterior_survival
    }
  }


  me <- list(posterior_treatment = posterior_treatment,
             posterior_control   = posterior_control,
             final               = final,
             args1               = args1)

  class(me) <- "bdpsurvival"

  return(me)
})




################################################################################
# Survival posterior estimation
# 1) Estimate the discount function (if current+historical data both present)
# 2) Estimate the posterior of the augmented data
################################################################################
### Combine  loss function and posterior estimation into one function
posterior_survival <- function(S, S0, surv_time, discount_function,
                               alpha_max, fix_alpha, a0, b0,
                               number_mcmc, weibull_shape, weibull_scale,
                               breaks, arm2, method){


  ### Extract intervals and count number of intervals
  ### - It should be that S_int equals S0_int
  if(!is.null(S)){
    S_int  <- levels(S$interval)
    nInt   <- length(S_int)
  }

  if(!is.null(S0)){
    S0_int  <- levels(S0$interval)
    nInt    <- length(S0_int)
  }

  interval <- NULL


  ##############################################################################
  # Discount function
  # - Comparison is made only if both S and S0 are present
  ##############################################################################
  # Compute hazards for historical and current data efficiently
  if(!is.null(S) & !is.null(S0)){
    ### Compute posterior of interval hazards
    a_post  <- b_post  <- numeric(nInt)
    a_post0 <- b_post0 <- numeric(nInt)

    if(!is.null(S)){
      posterior_flat_hazard <- prior_hazard <- matrix(NA, number_mcmc, nInt)
    } else{
      posterior_flat_hazard <- prior_hazard <- matrix(NA, number_mcmc, nInt)
    }


    ### Compute posterior values
    for(i in 1:nInt){
      a_post[i] <- a0 + sum(subset(S, interval==S_int[i])$status)
      b_post[i] <- b0 + sum(subset(S, interval==S_int[i])$exposure)

      a_post0[i] <- a0 + sum(subset(S0, interval==S0_int[i])$status)
      b_post0[i] <- b0 + sum(subset(S0, interval==S0_int[i])$exposure)

      ### Interval hazards - add on a very small value to avoid underflow
      posterior_flat_hazard[,i]  <- rgamma(number_mcmc, a_post[i],  b_post[i])+1e-12
      prior_hazard[,i]           <- rgamma(number_mcmc, a_post0[i], b_post0[i])+1e-12
    }
  } else if(!is.null(S) & is.null(S0)){
    ### Compute posterior of interval hazards
    a_post  <- b_post  <- numeric(nInt)

    posterior_flat_hazard <- matrix(NA, number_mcmc, nInt)
    prior_hazard          <- NULL

    ### Compute posterior values
    for(i in 1:nInt){
      a_post[i] <- a0 + sum(subset(S, interval==S_int[i])$status)
      b_post[i] <- b0 + sum(subset(S, interval==S_int[i])$exposure)

      ### Interval hazards - add on a very small value to avoid underflow
      posterior_flat_hazard[,i]  <- rgamma(number_mcmc, a_post[i],  b_post[i])+1e-12
    }
  } else if(is.null(S) & !is.null(S0)) {
    ### Compute posterior of interval hazards
    a_post0 <- b_post0 <- numeric(nInt)

    prior_hazard          <- matrix(NA, number_mcmc, nInt)
    posterior_flat_hazard <- NULL

    ### Compute posterior values
    for(i in 1:nInt){
      a_post0[i] <- a0 + sum(subset(S0, interval==S0_int[i])$status)
      b_post0[i] <- b0 + sum(subset(S0, interval==S0_int[i])$exposure)

      ### Interval hazards - add on a very small value to avoid underflow
      prior_hazard[,i]           <- rgamma(number_mcmc, a_post0[i], b_post0[i])+1e-12
    }
  }


  ### If only one of S or S0 is present, return related hazard and (if !arm2), return survival
  if(!is.null(S) & is.null(S0)){
    posterior_hazard <- posterior_flat_hazard

    if(!arm2){
      posterior_survival <- posterior_flat_survival <- 1 - ppexp(q=surv_time,
                                                                 x=posterior_hazard,
                                                                 cuts = c(0,breaks))
      prior_survival      <- NULL
    } else{
      posterior_survival <- posterior_flat_survival <- prior_survival <- NULL
    }
  } else if(is.null(S) & !is.null(S0)){
    posterior_hazard <- prior_hazard

    if(!arm2){
      posterior_survival <- prior_survival <- 1 - ppexp(q=surv_time,
                                                        x=posterior_hazard,
                                                        cuts = c(0,breaks))
      posterior_flat_survival <- NULL
    } else{
      posterior_survival  <- prior_survival <- posterior_flat_survival <- NULL
    }
  }

  if(!(!is.null(S) & !is.null(S0))){
    return(list(alpha_discount          = NULL,
                p_hat                   = NULL,
                posterior_survival      = posterior_survival,
                posterior_flat_survival = posterior_flat_survival,
                prior_survival          = prior_survival,
                posterior_hazard        = posterior_hazard,
                posterior_flat_hazard   = posterior_flat_hazard,
                prior_hazard            = prior_hazard))
  }


  ### If both S and S0 are present, carry out the comparison and compute alpha
  if(!arm2){
    ### Posterior survival probability
    posterior_flat_survival  <- 1 - ppexp(q=surv_time, x=posterior_flat_hazard, cuts = c(0,breaks))
    prior_survival           <- 1 - ppexp(q=surv_time, x=prior_hazard, cuts = c(0,breaks))

    ### Compute probability that survival is greater for current vs historical
    if(method == "mc"){
      logS  <- log(posterior_flat_survival)
      logS0 <- log(prior_survival)

      ### Variance of log survival, computed via delta method of hazards
      nIntervals <- sum(surv_time > c(0,breaks))
      IntLengths <- c(c(0,breaks)[1:nIntervals], surv_time)
      surv_times <- diff(IntLengths)

      v          <- as.matrix(posterior_flat_hazard[,1:nIntervals]^2) %*% (surv_times^2/a_post[1:nIntervals])
      v0         <- as.matrix(prior_hazard[,1:nIntervals]^2) %*% (surv_times^2/a_post0[1:nIntervals])
      Z          <- abs(logS - logS0) / (v+v0)
      p_hat      <- 2*(1-pnorm(Z))

    } else if(method == "fixed"){
      p_hat <- mean(posterior_flat_survival > prior_survival)   # higher is better survival
    } else{
      stop("Unrecognized method. Use one of 'fixed' or 'mc'")
    }

    if(fix_alpha){
      alpha_discount <- alpha_max
    } else{
        p_hat          <- 2*ifelse(p_hat > 0.5, 1 - p_hat, p_hat)

        # Compute alpha discount based on distribution
        if(discount_function == "weibull"){
          alpha_discount <- pweibull(p_hat, shape=weibull_shape,
                                     scale=weibull_scale)*alpha_max
        } else if(discount_function == "scaledweibull"){
          max_p <- pweibull(1, shape=weibull_shape, scale=weibull_scale)

          alpha_discount <- pweibull(p_hat, shape=weibull_shape,
                                     scale=weibull_scale)*alpha_max/max_p
        } else if(discount_function == "identity"){
          alpha_discount <- p_hat*alpha_max
        }
    }
  } else{
    ### Weight historical data via (approximate) hazard ratio comparing
    ### current vs historical
    if(method == "mc"){
      R0    <- log(prior_hazard)-log(posterior_flat_hazard)
      v     <- 1/(a_post)
      v0    <- 1/(a_post0)
      R     <- rowSums(as.matrix(R0/(v+v0) / sum(1/(v+v0))))
      V0    <- 1 / sum(1/(v+v0))
      Z     <- abs(R) / V0
      p_hat <- 2*(1-pnorm(Z))
    } else if(method == "fixed"){
      R0     <- log(prior_hazard)-log(posterior_flat_hazard)
      V0     <- 1/apply(R0,2,var)
      logHR0 <- R0%*%V0/sum(V0)    #weighted average  of SE^2
      p_hat  <- mean(logHR0 > 0)   #larger is higher failure
      p_hat  <- 2*ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
    } else{
      stop("Unrecognized method. Use one of 'fixed' or 'mc'")
    }

    if(fix_alpha){
      alpha_discount <- alpha_max
    } else{

        # Compute alpha discount based on distribution
        if(discount_function == "weibull"){
          alpha_discount <- pweibull(p_hat, shape=weibull_shape,
                                     scale=weibull_scale)*alpha_max
        } else if(discount_function == "scaledweibull"){
          max_p <- pweibull(1, shape=weibull_shape, scale=weibull_scale)

          alpha_discount <- pweibull(p_hat, shape=weibull_shape,
                                     scale=weibull_scale)*alpha_max/max_p
        } else if(discount_function == "identity"){
          alpha_discount <- p_hat*alpha_max
        }
    }
  }



  ##############################################################################
  # Posterior augmentation via the interval hazards
  # - If current or historical data are missing, this will not augment(see above)
  ##############################################################################
  posterior_hazard <- matrix(NA, number_mcmc, nInt)

  for(i in 1:nInt){
    a_post0[i] <- a0 + sum(subset(S0, interval==S0_int[i])$status)
    b_post0[i] <- b0 + sum(subset(S0, interval==S0_int[i])$exposure)

    a_post[i] <- a0 + sum(subset(S, interval==S_int[i])$status)
    b_post[i] <- b0 + sum(subset(S, interval==S_int[i])$exposure)

    ### Add on a very small value to avoid underflow
    posterior_hazard[,i]  <- rgamma(number_mcmc,
                                    a_post[i]+a_post0[i]*alpha_discount,
                                    b_post[i]+b_post0[i]*alpha_discount) + 1e-12
  }

  ### Posterior of survival time (if !arm2)
  if(!arm2){
    posterior_survival <- 1-ppexp(q=surv_time, x=posterior_hazard, cuts=c(0,breaks))
  } else{
    posterior_survival <- posterior_flat_survival <- prior_survival <- NULL
  }


  return(list(alpha_discount          = alpha_discount,
              p_hat                   = p_hat,
              posterior_survival      = posterior_survival,
              posterior_flat_survival = posterior_flat_survival,
              prior_survival          = prior_survival,
              posterior_hazard        = posterior_hazard,
              posterior_flat_hazard   = posterior_flat_hazard,
              prior_hazard            = prior_hazard))
}

