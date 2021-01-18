#' @title Bayesian Discount Prior: Two-Arm Logistic Regression
#' @description \code{bdplogit} is used to estimate the treatment effect
#'   in the presence of covariates using the logistic regression analysis
#'   implementation of the Bayesian discount prior. \code{summary} and
#'   \code{print} methods are supported. Currently, the function only supports
#'   a two-arm clinical trial where all of current treatment, current control,
#'   historical treatment, and historical control data are present
#' @param formula an object of class "formula." See "Details" for
#'   more information, including specification of treatment
#'   data indicators.
#' @param data a data frame containing the current data variables in the model.
#'   A column named \code{treatment} must be present; \code{treatment} must
#'   be binary and indicate treatment group vs. control group.
#' @param data0 a data frame containing the historical data variables in the model.
#'   The column labels of data and data0 must match.
#' @param prior_treatment_effect scalar. The historical adjusted treatment effect.
#'   If left \code{NULL}, value is estimated from the historical data.
#' @param prior_control_effect scalar. The historical adjusted control effect.
#'   If left \code{NULL}, value is estimated from the historical data.
#' @param prior_treatment_sd scalar. The standard deviation of the historical
#'   adjusted treatment effect. If left \code{NULL}, value is estimated from
#'   the historical data.
#' @param prior_control_sd scalar. The standard deviation of the historical
#'   adjusted control effect. If left \code{NULL}, value is estimated from
#'   the historical data.
#' @param prior_covariate_effect vector. The prior mean(s) of the covariate
#'   effect(s). Default value is zero. If a single value is input, the
#'   the scalar is repeated to the length of the input covariates. Otherwise,
#'   care must be taken to ensure the length of the input matches the number of
#'   covariates.
#' @param prior_covariate_sd vector. The prior standard deviation(s) of the
#'   covariate effect(s). Default value is 1e4. If a single value is input, the
#'   the scalar is repeated to the length of the input covariates. Otherwise,
#'   care must be taken to ensure the length of the input matches the number of
#'   covariates.
#' @param discount_function character. Specify the discount function to use.
#'   Currently supports \code{weibull}, \code{scaledweibull}, and
#'   \code{identity}. The discount function \code{scaledweibull} scales
#'   the output of the Weibull CDF to have a max value of 1. The \code{identity}
#'   discount function uses the posterior probability directly as the discount
#'   weight. Default value is "\code{identity}".
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1. Users may specify a vector of two values where the first
#'   value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param fix_alpha logical. Fix alpha at alpha_max? Default value is FALSE.
#' @param number_mcmc_alpha scalar. Number of Monte Carlo
#'   simulations for estimating the historical data weight. Default is 5000.
#' @param number_mcmc_beta scalar. Number of Monte Carlo simulations for
#'   estimating beta, the vector of regression coefficients. Default is 10000.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. Users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group. Not used when \code{discount_function} = "identity".
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. Users may specify a vector of two
#'   values where the first value is used to estimate the weight of the
#'   historical treatment group and the second value is used to estimate the
#'   weight of the historical control group.
#'   Not used when \code{discount_function} = "identity".
#' @param method character. Analysis method with respect to estimation of the weight
#'   paramter alpha. Default method "\code{mc}" estimates alpha for each
#'   Monte Carlo iteration. Alternate value "\code{fixed}" estimates alpha once and
#'   holds it fixed throughout the analysis.  See the the \code{bdplm} vignette \cr
#'   \code{vignette("bdplm-vignette", package="bayesDP")} for more details.
#' @details \code{bdplogit} uses a two-stage approach for determining the
#'   strength of historical data in estimation of an adjusted mean or covariate
#'   effect. In the first stage, a \emph{discount function} is used that
#'   that defines the maximum strength of the
#'   historical data and discounts based on disagreement with the current data.
#'   Disagreement between current and historical data is determined by stochastically
#'   comparing the respective posterior distributions under noninformative priors.
#'   Here with a two-arm regression analysis, the comparison is the
#'   proability (\code{p}) that the covariate effect of an historical data indicator is
#'   significantly different from zero. The comparison metric \code{p} is then
#'   input into the discount function and the final strength of the
#'   historical data is returned (\code{alpha}).
#'
#'   In the second stage, posterior estimation is performed where the discount
#'   function parameter, \code{alpha}, is used to weight the historical data
#'   effects.
#'
#'   The formula must include an intercept (i.e., do not use \code{-1} in
#'   the formula) and both data and data0 must be present.
#'   The column names of data and data0 must match. See \code{examples} below for
#'   example usage.
#'
#'   The underlying model uses the \code{MCMClogit} function of the MCMCpack
#'   package to carryout posterior estimation. Add more.
#'
#' @return \code{bdplogit} returns an object of class "bdplogit".
#'
#' An object of class "\code{bdplogit}" is a list containing at least
#' the following components:
#' \describe{
#'   \item{\code{posterior}}{
#'     data frame. The posterior draws of the covariates, the intercept, and
#'     the treatment effect. The grid of sigma values are included.}
#'   \item{\code{alpha_discount}}{
#'     vector. The posterior probability of the stochastic comparison
#'     between the current and historical data for each of the treatment
#'     and control arms. If \code{method="mc"}, the result is a matrix of
#'     estimates, otherwise for \code{method="fixed"}, the result is a vector
#'     of estimates.}
#'   \item{\code{estimates}}{
#'     list. The posterior means and standard errors of the intercept,
#'     treatment effect, covariate effect(s) and error variance.}
#' }
#'
#' @examples
#' # Set sample sizes
#' n_t  <- 30     # Current treatment sample size
#' n_c  <- 30     # Current control sample size
#' n_t0 <- 80     # Historical treatment sample size
#' n_c0 <- 80     # Historical control sample size
#'
#' # Treatment group vectors for current and historical data
#' treatment   <- c(rep(1,n_t), rep(0,n_c))
#' treatment0  <- c(rep(1,n_t0), rep(0,n_c0))
#'
#' # Simulate a covariate effect for current and historical data
#' x  <- rnorm(n_t+n_c, 1, 5)
#' x0 <- rnorm(n_t0+n_c0, 1, 5)
#'
#' # Simulate outcome:
#' # - Intercept of 10 for current and historical data
#' # - Treatment effect of 31 for current data
#' # - Treatment effect of 30 for historical data
#' # - Covariate effect of 3 for current and historical data
#' Y  <- 10 + 31*treatment  + x*3 + rnorm(n_t+n_c,0,5)
#' Y0 <- 10 + 30*treatment0 + x0*3 + rnorm(n_t0+n_c0,0,5)
#'
#' # Place data into separate treatment and control data frames and
#' # assign historical = 0 (current) or historical = 1 (historical)
#' df_ <- data.frame(Y=Y, treatment=treatment, x=x)
#' df0 <- data.frame(Y=Y0, treatment=treatment0, x=x0)
#'
#' # Fit model using default settings
#' fit <- bdplm(formula=Y ~ treatment+x, data=df_, data0=df0,
#'              method="fixed")
#'
#' # Look at estimates and discount weight
#' summary(fit)
#' print(fit)
#'
#' @rdname bdplogit
#' @import methods
#' @importFrom stats density is.empty.model median model.offset model.response
#'   pweibull pnorm quantile rbeta rgamma rnorm var vcov contrasts<- dt gaussian
#'   lm.fit model.frame model.matrix.default offset terms terms.formula
#'   coefficients lm qgamma runif glm binomial
#' @importFrom MCMCpack MCMClogit 
#' @aliases bdplogit,ANY-method
#' @useDynLib bayesDP
#' @export bdplogit
bdplogit <- setClass("bdplogit", slots = c(posterior_treatment = "list",
                                        posterior_control = "list",
                                        args1 = "list"))
setGeneric("bdplogit",
           function(formula                   = formula,
                    data                      = data,
                    data0                     = NULL,
                    prior_treatment_effect    = NULL, #0,
                    prior_control_effect      = NULL, #0,
                    prior_treatment_sd        = NULL, #1e4,
                    prior_control_sd          = NULL, #1e4,
                    prior_covariate_effect    = 0,
                    prior_covariate_sd        = 1e4,
                    number_mcmc_alpha         = 5000,
                    number_mcmc_beta          = 10000,
                    discount_function         = "identity",
                    alpha_max                 = 1,
                    fix_alpha                 = FALSE,
                    weibull_scale             = 0.135,
                    weibull_shape             = 3,
                    method                    = "mc"){
             standardGeneric("bdplogit")
           })

setMethod("bdplogit",
          signature(),
          function(formula                   = formula,
                   data                      = data,
                   data0                     = NULL,
                   prior_treatment_effect    = NULL, #0,
                   prior_control_effect      = NULL, #0,
                   prior_treatment_sd        = NULL, #1e4,
                   prior_control_sd          = NULL, #1e4,
                   prior_covariate_effect    = 0,
                   prior_covariate_sd        = 1e4,
                   number_mcmc_alpha         = 5000,
                   number_mcmc_beta          = 10000,
                   discount_function         = "identity",
                   alpha_max                 = 1,
                   fix_alpha                 = FALSE,
                   weibull_scale             = 0.135,
                   weibull_shape             = 3,
                   method                    = "mc"){
            
            ### Check validity of data input
            call <- match.call()
            if (missing(data)) {
              stop("Current data not input correctly.")
            }
            
            if (is.null(data0)) {
              stop("Historical data not input correctly.")
            }
            
            if(method == "mc"){
              stop("Method 'mc' not currently supported.")
            }
            
            ### Parse current data
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
            
            X <- if (!is.empty.model(mt)) {
              model.matrixBayes(object = mt, data = data, contrasts.arg = NULL,
                                keep.order = TRUE, drop.baseline = TRUE)
            } else {
              matrix(, NROW(Y), 0L)
            }
            
            ### Alter intercept column name, if present
            imatch <- match("(Intercept)", colnames(X))
            if(!is.na(imatch)) colnames(X)[imatch] <- "intercept"
            
            # Create indicator of whether each column is present
            cmatch <- match(c("intercept","treatment") , colnames(X))
            cmatch <- !is.na(cmatch)  ### == TRUE == present
            
            if(!all(cmatch)){
              stop("Current data is input incorrectly. Intercept and treatment must be present.")
            }
            
            ### Count levels of treatment data and sure 0 and 1 are present
            trt_levels <- levels(as.factor(X[,"treatment"]))
            
            if(!(all(trt_levels %in% c(0,1)))){
              stop("Treatment input has wrong levels. Values should be 0 and 1.")
            }
            
            ### Parse historical data
            if(missing(data0)){
              stop("Historical data is required.")
            } else{
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
              
              X0 <- if (!is.empty.model(mt0)) {
                model.matrixBayes(object = mt0, data = data0, contrasts.arg = NULL,
                                  keep.order = TRUE, drop.baseline = TRUE)
              } else {
                matrix(, NROW(Y0), 0L)
              }
              
              ### Alter intercept column name, if present
              imatch <- match("(Intercept)", colnames(X0))
              if(!is.na(imatch)) colnames(X0)[imatch] <- "intercept"
              
              # Create indicator of whether each column is present
              cmatch <- match(c("intercept","treatment") , colnames(X0))
              cmatch <- !is.na(cmatch)  ### == TRUE == present
              
              if(!all(cmatch)){
                stop("Historical data is input incorrectly. Intercept and treatment must be present.")
              }
              
              cnames  <- colnames(X)
              cnames0 <- colnames(X0)
              if(!all(cnames==cnames0)){
                stop("Column names in the historical data must match the current data.")
              }
              
              ### Count levels of treatment data and ensure 0 and 1 are present
              trt_levels0 <- levels(as.factor(X0[,"treatment"]))
              
              # Check that the levels are 0 and/or 1
              if(!(all(trt_levels0 %in% c(0,1)))){
                stop("Treatment input has wrong levels. Values should be 0 and 1.")
              }
              
            }
            
            ### If method == "mc", fix_alpha must be FALSE
            if(method=="mc" & fix_alpha){
              stop("mc method not possible with fix_alpha == TRUE. Set fix_alpha == FALSE." )
            }
            
            
            # Check that discount_function is input correctly
            all_functions <- c("weibull", "scaledweibull", "identity")
            function_match <- match(discount_function, all_functions)
            if(is.na(function_match)) {
              stop("discount_function input incorrectly.")
            }
            
            historical <- NULL
            treatment  <- NULL
            intercept  <- NULL
            
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
            ### Create a master dataframe
            df  <- data.frame(Y=Y, data.frame(X), historical=0)
            df0 <- data.frame(Y=Y0, data.frame(X0), historical=1)
            
            df <- rbind(df, df0)
            
            ### Split data into separate treatment & control dataframes
            df_t <- subset(df, treatment == 1)
            df_c <- subset(df, treatment == 0)
            
            ### Also create current data dataframe
            df_current <- subset(df, historical == 0, select=-intercept)
            
            ### Count number of covariates
            names_df <- names(df)
            covs_df  <- names_df[!(names_df %in% c("Y", "intercept", "treatment", "historical"))]
            n_covs   <- length(covs_df)
            
            
            ##############################################################################
            # Estimate discount weights for each of the treatment and control arms
            ##############################################################################
            discount_treatment <- discount_logit(df                 = df_t,
                                              discount_function  = discount_function,
                                              alpha_max          = alpha_max[1],
                                              fix_alpha          = fix_alpha,
                                              number_mcmc_alpha  = number_mcmc_alpha,
                                              weibull_shape      = weibull_shape[1],
                                              weibull_scale      = weibull_scale[1],
                                              method             = method)
            
            discount_control <- discount_logit(df                 = df_c,
                                            discount_function  = discount_function,
                                            alpha_max          = alpha_max[2],
                                            fix_alpha          = fix_alpha,
                                            number_mcmc_alpha  = number_mcmc_alpha,
                                            weibull_shape      = weibull_shape[2],
                                            weibull_scale      = weibull_scale[2],
                                            method             = method)
            
            ##############################################################################
            # Estimate historical adjusted treatment and control effects to use as the
            # priors for the current data
            ##############################################################################
            if(is.null(prior_treatment_effect) | is.null(prior_control_effect) |
               is.null(prior_treatment_sd)     | is.null(prior_control_sd)){
              
              df0$control <- 1-df0$treatment
              
              cnames0 <- names(df0)
              cnames0 <- cnames0[!(cnames0 %in% c("Y", "intercept", "historical","treatment","control"))]
              cnames0 <- c("treatment","control", cnames0)
              f0      <- paste0("Y~ -1+",paste0(cnames0,collapse="+"))
              
              fit_0   <- glm(f0, data=df0, family=binomial)
              summ_0  <- summary(fit_0)
              
              if(is.null(prior_treatment_effect)) prior_treatment_effect <- summ_0$coef[1,1]
              if(is.null(prior_control_effect))   prior_control_effect   <- summ_0$coef[2,1]
              
              if(is.null(prior_treatment_sd))     prior_treatment_sd     <- summ_0$coef[1,2]
              if(is.null(prior_control_sd))       prior_control_sd       <- summ_0$coef[2,2]
            }
            
            
            ### Prior covariate effects
            if(length(prior_covariate_effect) == 1){
              prior_covariate_effect <- rep(prior_covariate_effect,n_covs)
            }
            
            if(length(prior_covariate_sd) == 1){
              prior_covariate_sd <- rep(prior_covariate_sd,n_covs)
            }
            
            ##############################################################################
            # Estimate augmented treatment effect
            ##############################################################################
            ### Compute prior terms
            ### - Covarance/variance needs to be parameterized as precision
            tau2   <- c(1/prior_treatment_sd, 1/prior_control_sd, 1/prior_covariate_sd)^2
            mu0    <- c(prior_treatment_effect, prior_control_effect, prior_covariate_effect)
            
            
            ### Calculate constants from current data
            df_current$control <- 1-df_current$treatment
          
            # Create analysis formula
            df_analysis  <- df_current[,c("treatment", "control",covs_df)]
            cnames       <- names(df_analysis)
            f            <- paste0("Y~ -1 +",paste0(cnames,collapse="+"))
            


            ### Estimation scheme differs conditional on discounting method
            if(method == "fixed"){
              
              ### Extract alpha0, append "zero" weight for the covariate effect(s)
              alpha0 <- c(discount_treatment$alpha_discount + 1e-12,
                          discount_control$alpha_discount + 1e-12,
                          rep(1e-12, n_covs))
              
              ### Create precision matrix
              B0 <- diag(alpha0/tau2)
              
              
              ### Get Bayesian fit
              fit       <- MCMCpack::MCMClogit(f, data=df_analysis, 
                                               mcmc = number_mcmc_beta,
                                               b0 = mu0,
                                               B0 = diag(tau2))

              beta_samples <- data.frame(fit)
              
              
            } else if(method == "mc"){
              stop("Method 'mc' not currently supported.")
            }

            ### Estimate posterior of intercept and treatment effect
            beta_samples$intercept <- beta_samples$control
            beta_samples$treatment <- beta_samples$treatment-beta_samples$control

            ### Format posterior table
            m            <- ncol(beta_samples)
            beta_samples <- beta_samples[,c(m,1,3:(m-1))]
            
            ### Format alpha_discount values
            alpha_discount <- data.frame(treatment = discount_treatment$alpha_discount,
                                         control   = discount_control$alpha_discount)
            
            ### Format estimates
            estimates <- list()
            estimates$coefficients <- data.frame(t(colMeans(beta_samples)))
            estimates$coefficients <- estimates$coefficients[c("intercept", "treatment", covs_df)]
            estimates$se           <- data.frame(t(apply(beta_samples,2,sd)))
            estimates$se           <- estimates$se[c("intercept", "treatment", covs_df)]
            
            
            ### Format input arguments
            args1 <- list(call  = call,
                          data  = df)
            
            ### Format output
            me <- list(posterior      = beta_samples,
                       alpha_discount = alpha_discount,
                       estimates      = estimates,
                       args1          = args1)
            
            class(me) <- "bdplm"
            return(me)
          })



################################################################################
# Logistic Regression discount weight estimation
# 1) Estimate the discount function
#    - Test that the historical effect is different from zero
#    - Estimate alpha based on the above comparison
################################################################################
discount_logit <- function(df, discount_function, alpha_max, fix_alpha,
                           number_mcmc_alpha,
                           weibull_shape, weibull_scale,
                           method){
  
  # Create formula
  cnames <- names(df)
  cnames <- cnames[!(cnames %in% c("Y", "intercept", "historical","treatment"))]
  cnames <- c("historical", cnames)
  f      <- paste0("Y~",paste0(cnames,collapse="+"))
  
  ### Get Bayesian fit
  fit       <- MCMCpack::MCMClogit(f, data=df, mcmc = number_mcmc_alpha)
  posterior <- data.frame(fit)
  
  ### Monte Carlo simulations of historical effect
  beta      <- posterior$historical
  
  ### Compute posterior comparison
  if(method == "fixed"){
    p_hat  <- mean(beta>0)
    p_hat  <- 2*ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
  } else if(method == "mc"){
    stop("Method 'mc' not currently supported.")
    #Z     <- abs(beta)/sigma2_beta
    #p_hat <- 2*(1-pnorm(Z))
  }
  
  ### Compute alpha, the discount weight
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
  
  res <- list(p_hat          = p_hat,
              alpha_discount = alpha_discount)
  res
}



model.matrixBayes <- function(object, data=environment(object), contrasts.arg=NULL,
                              xlev=NULL, keep.order=FALSE, drop.baseline=FALSE,...){
  
  t <- if( missing( data ) ) {
    terms( object )
  } else{
    terms.formula(object, data = data, keep.order=keep.order)
  }
  
  attr(t, "intercept") <- attr(object, "intercept")
  if (is.null(attr(data, "terms"))){
    data <- model.frame(object, data, xlev=xlev)
  } else {
    reorder <- match(sapply(attr(t,"variables"), deparse, width.cutoff=500)[-1], names(data))
    if (anyNA(reorder)) {
      stop( "model frame and formula mismatch in model.matrix()" )
    }
    if(!identical(reorder, seq_len(ncol(data)))) {
      data <- data[,reorder, drop = FALSE]
    }
  }
  
  int <- attr(t, "response")
  if(length(data)) {      # otherwise no rhs terms, so skip all this
    if (drop.baseline){
      contr.funs <- as.character(getOption("contrasts"))
    }else{
      contr.funs <- as.character(list("contr.bayes.unordered", "contr.bayes.ordered"))
    }
    
    namD <- names(data)
    
    ## turn  character columns into factors
    for(i in namD)
      if(is.character( data[[i]] ) ) {
        data[[i]] <- factor(data[[i]])
        warning( gettextf( "variable '%s' converted to a factor", i ), domain = NA)
      }
    isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)
    isF[int] <- FALSE
    isOF <- vapply(data, is.ordered, NA)
    for( nn in namD[isF] )            # drop response
      if( is.null( attr( data[[nn]], "contrasts" ) ) ) {
        contrasts( data[[nn]] ) <- contr.funs[1 + isOF[nn]]
      }
    
    ## it might be safer to have numerical contrasts:
    ##    get(contr.funs[1 + isOF[nn]])(nlevels(data[[nn]]))
    if ( !is.null( contrasts.arg ) && is.list( contrasts.arg ) ) {
      if ( is.null( namC <- names( contrasts.arg ) ) ) {
        stop( "invalid 'contrasts.arg' argument" )
      }
      for (nn in namC) {
        if ( is.na( ni <- match( nn, namD ) ) ) {
          warning( gettextf( "variable '%s' is absent, its contrast will be ignored", nn ), domain = NA )
        }
        else {
          ca <- contrasts.arg[[nn]]
          if( is.matrix( ca ) ) {
            contrasts( data[[ni]], ncol( ca ) ) <- ca
          }
          else {
            contrasts( data[[ni]] ) <- contrasts.arg[[nn]]
          }
        }
      }
    }
  } else {
    isF  <-  FALSE
    data <- data.frame(x=rep(0, nrow(data)))
  }
  ans  <- model.matrix.default(object=t, data=data)
  cons <- if(any(isF)){
    lapply( data[isF], function(x) attr( x,  "contrasts") )
  }else { NULL }
  attr(ans, "contrasts" ) <- cons
  ans
}





# library(MCMCpack)
# df0 <- data.frame(y = c(1,1,1,rep(0,138)),
#                   x = rnorm(141,5,0.5),
#                   historical = 1)
# 
# df <- data.frame(y = rbinom(100,1,0.03),
#                  x = rnorm(100,3,1),
#                  historical = 0)
# 
# df_ <- rbind(df0,df)
# 
# fit <- MCMClogit(y~x+historical, data=df_,
#                  mcmc = 1e4)
# 
# posterior <- data.frame(fit)
# p_hat <- mean(posterior$historical > 0)
# p_hat  <- 2*ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
# 
# 
# 
# # b0=c(0,0,1), 
# # B0=B0)

