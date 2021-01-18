#' Strategy Estimation Function
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param data a \code{stratEst.data} object or \code{data.frame}.
#' @param strategies a list of strategies. Each element if the list must be an object of class \code{stratEst.strategy}.
#' @param shares a numeric vector of strategy shares. The order of the elements corresponds to the order in \code{strategies}. Elements which are \code{NA} are estimated from the data. Use a list of numeric vectors if data has more than one sample and shares are sample specific.
#' @param coefficients a matrix of latent class regression coefficients.
#' @param covariates a character vector with the names of the covariates of the latent class regression model in the data. The covariates cannot have missing values.
#' @param sample.id a character string indicating the name of the variable which identifies the samples in data. Individual observations must be nested in samples.
#' @param response a character string which is either \code{"pure"} or \code{"mixed"}. If \code{"pure"} the estimated choice probabilities are either zero or one. If \code{"mixed"} the estimated choice probabilities are mixed parameters. The default is \code{"mixed"}.
#' @param sample.specific a character vector, Defines the model parameters that are sample specific. Can contain the character strings \code{"shares"} (\code{"probs"},  \code{"trembles"}. If the vector contains \code{"shares"} (\code{"probs"}, \code{"trembles"}), the estimation function estimates a set of shares (choice probabilities, trembles) for each sample in the data.
#' @param r.probs a character string. Options are \code{"no"}, \code{"strategies"}, \code{"states"} or \code{"global"}. Option \code{"no"} yields one vector of choice probabilities per strategy and state. Option \code{"strategies"} yields one vector of choice probabilities per strategy. Option \code{"states"} yields one vector of choice probabilities per state. Option \code{"global"} yields a single vector of choice probabilities. Default is \code{"no"}.
#' @param r.trembles a character string. Options are \code{"no"}, \code{"strategies"}, \code{"states"} or \code{"global"}. Option \code{"no"} yields one tremble probability per strategy and state. Option \code{"strategies"} yields one tremble probability per strategy. Option \code{"states"} yields one tremble probability per state. Option \code{"global"} yields a single tremble probability. Default is \code{"no"}.
#' @param select a character vector. Indicates the classes of model parameters that are selected. Can contain the character strings \code{"strategies"}, (\code{"probs"}, and \code{"trembles"}.  If the vector contains\code{"strategies"} (\code{"probs"}, \code{"trembles"}), the number of strategies (choice probabilities, trembles) is selected based on the selection criterion specified in \code{"crit"}. The selection can be restricted with the arguments \code{r.probs} and \code{r.trembles}. Default is \code{NULL}.
#' @param min.strategies an integer. The minimum number of strategies in case of strategy selection. The strategy selection procedure stops if the minimum is reached.
#' @param crit a character string. Defines the information criterion used for model selection. Options are \code{"bic"} (Bayesian information criterion), \code{"aic"} (Akaike information criterion) or \code{"icl"} (Integrated Classification Likelihood). Default is \code{"bic"}.
#' @param se a string. Defines how standard errors are obtained. Options are \code{"analytic"} or \code{"bootstrap"}. Default is \code{"analytic"}.
#' @param outer.runs an integer. The number of outer runs of the solver. Default is 1.
#' @param outer.tol a number close to zero. The tolerance of the stopping condition of the outer runs. The iterative algorithm stops if the relative decrease of the log likelihood is smaller than this number. Default is 1e-10.
#' @param outer.max an integer. The maximum number of iterations of the outer runs of the solver. The iterative algorithm stops after \code{"outer.max"} iterations if it does not converge. Default is 1000.
#' @param inner.runs  an integer. The number of inner runs of the solver. Default is 10.
#' @param inner.tol a number close to zero. The tolerance of the stopping condition of the inner runs. The iterative algorithm stops if the relative decrease of the log likelihood is smaller than this number. Default is 1e-5.
#' @param inner.max an integer. The maximum number of iterations of the outer runs of the solver. The iterative algorithm stops after \code{"inner.max"} iterations if it does not converge. Default is 10.
#' @param lcr.runs an integer. The number of latent class regression runs of the solver. Default is 100.
#' @param lcr.tol a number close to zero. The tolerance of the stopping condition of the latent class regression runs. The iterative algorithm stops if the relative decrease of the log likelihood is smaller than this number. Default is 1e-10.
#' @param lcr.max an integer. The maximum number of iterations of the latent class regression runs of the solver. The iterative algorithm stops after \code{"lcr.max"} iterations if it does not converge. Default is 1000.
#' @param bs.samples an integer. The number of bootstrap samples.
#' @param quantiles a numeric vector. The quantiles of the sampling distribution of the estimated parameters. Depending on the option of \code{se}, the quantiles are either estimated based on a t-distribution with \code{res.degrees} of freedom and the analytic standard errors or based the bootstrap.
#' @param step.size a number between zero and one. The step size of the Fisher scoring step which updates the coefficients. Values smaller than one slow down the convergence of the algorithm and prevent overshooting. Default is one.
#' @param penalty a logical. If \code{TRUE} the Firth penalty is used to estimate the coefficients of the latent class regression model. Default is \code{FALSE}.
#' @param verbose a logical. If \code{TRUE} information about the estimation process are printed to the console. Default is \code{FALSE}.
#' @note Strategy estimation was introduced by Dal Bo and Frechette (2011) to estimate the maximum likelihood frequencies of a set of candidate strategies in the repeated prisoner's dilemma. Breitmoser (2015) introduces model parameters for the choice probabilities of individual strategies to the strategy estimation model. Dvorak and Fehrler (2018) extend the basic strategy estimation model by individual level covariates to explain the selection of strategies by individuals. The estimation function of the package obtains maximum likelihood estimates for the model parameters based on expectation maximization (Dempster, Laird, and Rubin, 1977) and Newton-Raphson algorithms. To decrease the computation time, the package integrates \proglang{C++} and \proglang{R} with the help of the \proglang{R} packages \pkg{Rcpp} (Eddelbuettel and Francois, 2011) and the open source linear algebra library for the C++ language \pkg{RppArmadillo} (Sanderson and Curtin, 2016).
#' @return An object of class \code{stratEst.model}. A list with the following elements.
#' \item{strategies}{the fitted strategies.}
#' \item{shares}{the strategy shares.}
#' \item{probs}{the choice probabilities of the strategies.}
#' \item{trembles}{the tremble probabilities of the strategies.}
#' \item{gammas}{the gamma parameters of the strategies.}
#' \item{coefficients}{ the coefficients of the covariates.}
#' \item{shares.par}{the estimated strategy share parameters.}
#' \item{probs.par}{the estimated choice probability parameters.}
#' \item{trembles.par}{the estimated tremble parameters.}
#' \item{gammas.par}{the estimated gamma parameters.}
#' \item{coefficients.par}{the estimated coefficient parameters of the covariates.}
#' \item{shares.indices}{the parameter indices of the strategy shares.}
#' \item{probs.indices}{the parameter indices of the choice probabilities.}
#' \item{trembles.indices}{the parameter indices of the tremble probabilities.}
#' \item{coefficients.indices}{the parameter indices of the coefficients.}
#' \item{loglike}{the log likelihood of the model.}
#' \item{num.ids}{the number of individuals.}
#' \item{num.obs}{the number of observations.}
#' \item{num.par}{the total number of model parameters.}
#' \item{free.par}{the number of free model parameters.}
#' \item{res.degrees}{the residual degrees of freedom.}
#' \item{aic}{the Akaike information criterion.}
#' \item{bic}{the Bayesian information criterion.}
#' \item{icl}{The integrated classification likelihood.}
#' \item{crit.val}{the value of the selection criterion defined by the argument \code{crit}.}
#' \item{eval}{the total number of iterations of the solver.}
#' \item{tol.val}{the relative decrease of the log likelihood in the last iteration of the algorithm.}
#' \item{convergence}{the maximum of the absolute scores of the estimated model parameters.}
#' \item{entropy.model}{the entropy of the model.}
#' \item{entropy.assignments}{the entropy of the posterior probability assignments of individuals to strategies.}
#' \item{chi.global}{the chi square statistic for global model fit.}
#' \item{chi.local}{the chi square statistics for local model fit.}
#' \item{state.obs}{the weighted observations for each strategy state.}
#' \item{post.assignments}{the posterior probability assignments of individuals to strategies.}
#' \item{prior.assignments}{the prior probability of each individual to use a strategy as predicted by the individual covariates.}
#' \item{shares.se}{the standard errors of the estimated share parameters.}
#' \item{probs.se}{the standard errors of the estimated choice probability parameters.}
#' \item{trembles.se}{the standard errors of the estimated tremble probability parameters.}
#' \item{gammas.se}{the standard errors of the estimated gamma parameters.}
#' \item{coefficients.se}{the standard errors of the estimated coefficients.}
#' \item{shares.quantiles}{the quantiles of the estimated population shares.}
#' \item{probs.quantiles}{the quantiles of the estimated choice probabilities.}
#' \item{trembles.quantiles}{the quantiles of the estimated trembles.}
#' \item{coefficients.quantiles}{the quantiles of the estimated coefficients.}
#' \item{shares.score}{the scores of the estimated share parameters.}
#' \item{probs.score}{the score of the estimated choice probabilities.}
#' \item{trembles.score}{the score of the estimated tremble probabilities.}
#' \item{coefficients.score}{the score of the estimated coefficient.}
#' \item{shares.fisher}{the Fisher information matrix of the estimated shares.}
#' \item{probs.fisher}{the Fisher information matrix of the estimated choice probabilities.}
#' \item{trembles.fisher}{the Fisher information matrix of the estimated trembles.}
#' \item{coefficients.fisher}{the fisher information matrix of the estimated coefficients.}
#' \item{fit.args}{the input objects of the function call.}
#' @description The estimation function of the package.
#' @details The estimation function of the package obtains maximum likelihood estimates for the model parameters based on expectation maximization and Newton-Raphson algorithms.
#' @references
#' Breitmoser Y (2015). "Cooperation, but no Reciprocity: Individual Strategies in the Repeated Prisoner’s Dilemma." \emph{American Economic Review}, 105(9), 2882-2910.
#'
#' Dal Bo P, Frechette GR (2011). "The Evolution of Cooperation in Infinitely Repeated Games: Experimental Evidence." \emph{American Economic Review}, 101(1), 411-429.
#'
#' Dempster A, Laird N, Rubin DB (1977). "Maximum Likelihood from Incomplete Data via the EM Algorithm." \emph{Journal of the Royal Statistical Society Series B}, 39(1), 1-38.
#'
#' Dvorak F, Fehrler S (2018). "Negotiating Cooperation under Uncertainty: Communication in Noisy, Indefinitely Repeated Interactions." \emph{IZA Working Paper}, No. 11897.
#'
#' Dvorak F, Fischbacher U, Schmelz K (2020). "Incentives for Conformity and Anticonformity." \emph{TWI Working Paper Series}.
#'
#' Eddelbuettel D, Francois R (2011). "Rcpp: Seamless R and C++ Integration." \emph{Journal of Statistical Software}, 40(8), 1-18.
#'
#' Fudenberg D, Rand DG, Dreber A (2012). "Slow to Anger and Fast to Forgive: Cooperation in an Uncertain World." \emph{American Economic Review}, 102(2), 720-749.
#'
#' Sanderson C, Curtin R (2016). "Armadillo: A Template-Based C++ Library for Linear Algebra." \emph{Journal of Open Source Software}, 1, 26.
#'
#' Wang Z, Xu B, Zhou HJ (2014). "Social Cycling and Conditional Responses in the Rock-Paper-Scissors Game." \emph{Scientific Reports}, 4(1), 2045-2322.
#'
#'#' @examples
#' ## Strategy model for rock-paper-scissors data of Wang, Xu, and Zhou (2014).
#' ## Fit a mixture of the Nash strategy and a strategy that imitates the last choice.
#' strategies.mixture = list("nash" = strategies.RPS$nash, "imitate" = strategies.RPS$imitate)
#' model.mixture <- stratEst.model(data.WXZ2014,strategies.mixture)
#'
#' ## Replication of Dal Bo and Frechette (2011), Table 7, page 424
#' model.DF2011 <- stratEst.model(data.DF2011, strategies.DF2011, sample.id = "treatment")
#'
#' ## Replication of Dvorak, Fischbacher, and Schmelz (2020)
#' covs <- c("intercept", "conformity.score")
#' model.DFS2020 <- stratEst.model(data.DFS2020, strategies.DFS2020, covariates = covs)
#'
#' @export
stratEst.model <- function( data, strategies, shares = NULL, coefficients = NULL, covariates = NULL, sample.id = NULL,  response = "mixed", sample.specific = c("shares","probs","trembles"), r.probs = "no", r.trembles = "global", select = NULL, min.strategies = 1, crit = "bic", se = "analytic", outer.runs = 1, outer.tol = 1e-10, outer.max = 1000, inner.runs = 10, inner.tol = 1e-5, inner.max = 10, lcr.runs = 100, lcr.tol = 1e-10, lcr.max = 1000 , bs.samples = 1000, quantiles = c(0.05,0.5,0.95), step.size = 1, penalty = F, verbose = FALSE ){

  #fit args
  fit.args = list()

  # check data
  if( missing(data) ) {
    stop("stratEst error: Mandatory input object 'data' missing.")
  }
  else{
    fit.args[[1]] = data
    stratEst.return.data <- stratEst.check.data( data )
    data <- stratEst.return.data$data
    id <- stratEst.return.data$id
    game <- stratEst.return.data$game
    period <- stratEst.return.data$period
    input <- stratEst.return.data$input
    output <- stratEst.return.data$output
    input_factor <- stratEst.return.data$input.factor
    output_factor <- stratEst.return.data$output.factor
    levels_input <- stratEst.return.data$levels.input
    levels_output <- stratEst.return.data$levels.output
    input.is.null <- stratEst.return.data$input.is.null
  }

  # check sample.id
  if( is.null(sample.id) ){
    sample <- rep(1,length(id))
    num_samples <- 1
    sample_is_factor = FALSE
    sample_levels = FALSE
  }
  else{
    stratEst.check.sample.id.return <- stratEst.check.sample.id( data , sample.id )
    sample <- stratEst.check.sample.id.return$sample
    sample_factor <- stratEst.check.sample.id.return$sample.factor
    num_samples <- stratEst.check.sample.id.return$num.sample
    sample_levels <- stratEst.check.sample.id.return$sample.levels
    data$sample <- sample
    sample_is_factor = TRUE
  }

  # check cluster.id
  #if( T ){
    cluster <- matrix(0,1,1)
  #}
  # else{
  #   stratEst.check.cluster.id.return <- stratEst.check.cluster.id( data , cluster.id )
  #   cluster <- stratEst.check.cluster.id.return$cluster
  #   cluster_factor <- stratEst.check.cluster.id.return$cluster.factor
  # }

  # check covariates
  if( is.null(covariates) ){
    covariate_mat <- matrix(0,1,1)
    LCR = FALSE
  }
  else{
    covariate_mat <- stratEst.check.covariates( data , covariates )
    LCR = TRUE
  }

  # check other inputs
  stratEst.check.other.return <- stratEst.check.other( response , sample.specific , r.probs , r.trembles , select , min.strategies , crit , se , outer.runs , outer.tol , outer.max , inner.runs , inner.tol , inner.max , lcr.runs , lcr.tol , lcr.max , bs.samples , step.size , penalty , verbose , quantiles )
  select_strategies = stratEst.check.other.return$select.strategies
  select_responses = stratEst.check.other.return$select.responses
  select_trembles = stratEst.check.other.return$select.trembles
  specific_shares = stratEst.check.other.return$specific.shares
  specific_responses = stratEst.check.other.return$specific.responses
  specific_trembles = stratEst.check.other.return$specific.trembles
  specific_coefficients = stratEst.check.other.return$specific.coefficients
  quantile_vec = stratEst.check.other.return$quantile.vec
  print.messages = stratEst.check.other.return$print.messages
  print.summary = stratEst.check.other.return$print.summary

  # check strategies
  if( missing(strategies) ) {
    stop("stratEst error: Mandatory input object 'strategies' is missing. Specify an integer or create a data.frame object.")
  }else{
    fit.args[[2]] = strategies
    stratEst.check.strategies.return <- stratEst.check.strategies( strategies , input_factor , output_factor , input , output , select_strategies )
    strategies <- stratEst.check.strategies.return$strategies
    strategies_matrix <- stratEst.check.strategies.return$strategies.matrix
    trembles <- stratEst.check.strategies.return$trembles
    num_strats <- stratEst.check.strategies.return$num_strats
    unique_inputs <- stratEst.check.strategies.return$unique.inputs
    unique_outputs <- stratEst.check.strategies.return$unique.outputs
    num_unique_inputs <- stratEst.check.strategies.return$num.unique.inputs
    num_unique_outputs <- stratEst.check.strategies.return$num.unique.outputs
    sid <- stratEst.check.strategies.return$sid
    integer_strategies <- stratEst.check.strategies.return$integer.strategies
    response_mat_col_index <- stratEst.check.strategies.return$response.mat.col.index
    names_strategies <- stratEst.check.strategies.return$names.strategies
  }

  #check shares
  if( is.null(shares) ) {
    shares_matrix = matrix( NA , num_strats , num_samples )
  }else{
    shares_matrix <- stratEst.check.shares( shares , LCR , specific_shares , num_samples , num_strats , sample.id , sample_levels , select_strategies )
  }

  #check coefficients
  if( is.null(coefficients) ) {
    if( LCR ){
      coefficient_mat = matrix(NA,length(covariates),num_strats,1)
      coefficient_mat[,1] = 0
    }else{
      coefficient_mat = matrix(0,1,1)
    }
  }else{
    coefficient_mat <- stratEst.check.coefficients( coefficients , covariates , num_strats , names_strategies )
    if( LCR == F ){
      warning("stratEst warning: No covariates specified. The input object 'coefficients' is ignored.");
    }
  }

  #########################################################################################################
  # PREPARE AND ORDER DATA
  #########################################################################################################

  # prepare data
  data_matrix <- cbind(id,game,period,input,output,sample)

  # sort data
  data_matrix <- data_matrix[order(data_matrix[,3]), ]
  data_matrix <- data_matrix[order(data_matrix[,2]), ]
  data_matrix <- data_matrix[order(data_matrix[,1]), ]

  # unique ids
  num_obs <- nrow(data_matrix)
  unique_ids <- sort(unique(id))
  num_unique_ids <- length(unique_ids)

  # names of fit arguments
  fit.args[[3]] = shares
  fit.args[[4]] = coefficients
  fit.args[[5]] = covariates
  fit.args[[6]] = sample.id
  fit.args[[7]] = response
  fit.args[[8]] = sample.specific
  fit.args[[9]] = r.probs
  fit.args[[10]] = r.trembles
  fit.args[[11]] = select
  fit.args[[12]] = min.strategies
  fit.args[[13]] = crit
  fit.args[[14]] = se
  fit.args[[15]] = outer.runs
  fit.args[[16]] = outer.tol
  fit.args[[17]] = outer.max
  fit.args[[18]] = inner.runs
  fit.args[[19]] = inner.tol
  fit.args[[20]] = inner.max
  fit.args[[21]] = lcr.runs
  fit.args[[22]] = lcr.tol
  fit.args[[23]] = lcr.max
  fit.args[[24]] = bs.samples
  fit.args[[25]] = quantiles
  fit.args[[26]] = step.size
  fit.args[[27]] = penalty
  fit.args[[28]] = verbose
  names(fit.args) = c("data","strategies","shares","coefficients","covariates","sample.id","response","sample.specific","r.probs","r.trembles","select","min.strategies","crit","se","outer.runs","outer.tol","outer.max","inner.runs","inner.tol","inner.max","lcr.runs","lcr.tol","lcr.max","bs.samples","quantiles","step.size","penalty","verbose")

  #########################################################################################################
  # ESTIMATION
  #########################################################################################################

  # create cpp output
  cpp.output <- stratEst_cpp( data_matrix, strategies_matrix, sid , shares_matrix , trembles , coefficient_mat, covariate_mat, LCR, cluster, quantile_vec , response, specific_shares , specific_responses , specific_trembles , specific_coefficients , r.probs, r.trembles, select_strategies , select_responses , select_trembles , min.strategies, crit, se, outer.runs, outer.tol, outer.max, inner.runs, inner.tol, inner.max, lcr.runs, lcr.tol, lcr.max, bs.samples , print.messages, integer_strategies, step.size , penalty )

  # stratEst.return object
  stratEst.return <- list("strategies" = cpp.output$strategies, "shares" = cpp.output$shares, "coefficients" = cpp.output$coefficients, "probs" = cpp.output$responses, "trembles" = cpp.output$trembles, "gammas" = 0, "shares.par" = cpp.output$shares.list$shares.par, "probs.par" = cpp.output$responses.list$responses.par, "trembles.par" = cpp.output$trembles.list$trembles.par,"gammas.par" = 0,  "coefficients.par" =  cpp.output$coefficients.list$coefficients.par,  "shares.indices" = cpp.output$shares.list$shares.indices, "probs.indices" = cpp.output$responses.list$responses.indices, "trembles.indices" = cpp.output$trembles.list$trembles.indices, "coefficients.indices" = cpp.output$coefficients.list$coefficients.indices, "loglike" = cpp.output$fit[1,1], "num.ids" = NA, "num.obs" = NA, "num.par" = NA, "free.par" = NA, "res.degrees" = NA, "aic" = NA, "bic" = NA, "icl" = NA, "entropy.model" = NA, "entropy.assignments" = NA, "chi.square" = NULL, "fit" = NA, "crit.val" = cpp.output$fit[1,2], "eval" = cpp.output$solver[1,1], "tol.val" = cpp.output$solver[1,2], "convergence" = cpp.output$convergence, "state.obs" = cpp.output$state.obs, "post.assignment" = cpp.output$assignments, "prior.assignment" = cpp.output$priors, "shares.se" = cpp.output$shares.list$shares.se, "probs.se" = cpp.output$responses.list$responses.se, "trembles.se" = cpp.output$trembles.list$trembles.se, "gammas.se" = 0, "coefficients.se" = cpp.output$coefficients.list$coefficients.se, "shares.quantiles" = NA, "probs.quantiles" = NA, "trembles.quantiles" = NA, "coefficients.quantiles" = NA, "shares.covar" = cpp.output$shares.list$shares.covar, "shares.score" =  cpp.output$shares.list$shares.score, "shares.fisher" = cpp.output$shares.list$shares.fisher, "probs.covar" = cpp.output$responses.list$responses.covar, "probs.score" = cpp.output$responses.list$responses.score, "probs.fisher" = cpp.output$responses.list$responses.fisher, "trembles.covar" = cpp.output$trembles.list$trembles.covar, "trembles.score" = cpp.output$trembles.list$trembles.score, "trembles.fisher" = cpp.output$trembles.list$trembles.fisher, "coefficients.covar" = cpp.output$coefficients.list$coefficients.covar, "coefficients.score" = cpp.output$coefficients.list$coefficients.score, "coefficients.fisher" = cpp.output$coefficients.list$coefficients.fisher, "covariate.mat" = cpp.output$coefficients.list$covariate.mat  );

  # post-processing
  stratEst.return <- stratEst.post( data , cpp.output , stratEst.return , strategies , covariates , response , unique_ids , num_unique_ids , input , output , unique_inputs , unique_outputs , num_unique_inputs , num_unique_outputs , sample , sample.id , sample_factor , num_samples , specific_shares , specific_responses , specific_trembles , sample_is_factor , integer_strategies , cpp.output$lcr , response_mat_col_index , crit , num_obs , se , quantile_vec, output_factor, input.is.null )

  stratEst.return$fit.args = fit.args

  class(stratEst.return) <- c("stratEst.model","list")

  # print summary
  if( print.summary ){
    summary( stratEst.return )
  }

  # return result
  return(stratEst.return)

}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
