#' Strategy Estimation Function
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param data A \code{stratEst.data} object or \code{data.frame}. Must contain the variables \code{choice},  \code{input}, \code{id}, \code{game}, \code{period}. The variable \code{id} identifies observations of the same individual across games and periods. The factor \code{input} indicates the discrete information observed by the individual before making a choice. The factor \code{choice} indicates the choice of the individual.
#' @param strategies A list of strategies. Each strategy is a data.frame of class \code{stratEst.strategy}. Each row of the data.frame represents one state of the strategy. The first row defines the initial state which is entered if the variable input is NA. Column names which start with the string 'output.' indicate the columns which contain the multinomial choice probabilities of the strategy. For example, a column labeled 'output.x' contains the probability to observe the output 'x'. The column 'tremble' contains a tremble probability for pure strategies. Column names which start with the string 'input.' indicate the columns which contain the deterministic state transition of the strategy. For example, a column with name 'input.x' indicates the state transition after observing input 'x'.
#' @param shares A vector of strategy shares. The elements to the order of strategies in the list \code{strategies}. Shares which are \code{NA} are estimated from the data. With more than one sample and sample specific shares, a list of column vectors is required.
#' @param coefficients Column vector which contains the latent class regression coefficients. The elements correspond to the vector of estimates.
#' @param covariates A character vector indicating the names of the variables in data that are the covariates of the latent class regression model. Rows with the same id must have the values of covariates. Missing value are not allowed.
#' @param sample.id A character indicating the name of the variable which identifies the samples. Individual observations must be nested in samples. The same must be true for clusters if specified. If more than one sample exists, shares are estimated for each sample. All other parameters are estimated for the data of all samples. If the object is not supplied, it is assumed that the data contains only one sample.
#' @param response A string which can be set to \code{"pure"} or \code{"mixed"}. If set to \code{"pure"} all estimated choice probabilities are pure, i.e. either zero or one. If set to \code{"mixed"} all estimated choice probabilities are mixed. The default is \code{"mixed"}.
#' @param sample.specific A character vector defining which model parameters are sample specific. If the vector contains the character \code{"shares"} (\code{"probs"}, \code{"trembles"}), the estimation function estimates a set of shares (choice probabilities, trembles) for each sample in the data. If the vector does not contains the character \code{"shares"} (\code{"probs"}, \code{"trembles"}) one set of shares (choice probabilities, trembles) is estimated for the pooled data of all samples. Default is \code{c("shares","probs","trembles")}.
#' @param r.probs A string which can be set to \code{"no"}, \code{"strategies"}, \code{"states"} or \code{"global"}. If set to \code{"strategies"}, the estimation function estimates strategies with one strategy specific vector of choice probabilities in every state of the strategy. If set to \code{"states"}, one state specific vector of choice probabilities is estimated for each state. If set to \code{"global"}, a single vector of probabilities is estimated which applies in every state of each strategy. Default is \code{"no"}.
#' @param r.trembles A string which can be set to \code{"no"}, \code{"strategies"}, \code{"states"} or \code{"global"}. If set to \code{"strategies"}, the estimation unction estimates strategies with one strategy specific tremble probability. If set to  \code{"states"}, one state specific tremble probability is estimated for each state. If set to \code{"global"}, a single tremble probability is estimated which globally. Default is \code{"global"}.
#' @param select A character vector indicating which model parameters are selected. If the vector contains the character \code{"strategies"} (\code{"probs"}, \code{"trembles"}), the number of strategies (choice probabilities, trembles) is selected based on the selection criterion specified in \code{"crit"}. The selection of choice probabilities and trembles occurs obeying the restriction specified in \code{r.probs} and \code{r.trembles}. (E.g. if \code{r.probs} is set to \code{"strategies"}, \code{select = "probs"} will select the sets of choice probabilities within each strategy). Default is \code{NULL}.
#' @param min.strategies An integer which specifies the minimum number of strategies in case of strategy selection. The strategy selection procedure stops if the minimum is reached.
#' @param crit A string which can be set to \code{"bic"}, \code{"aic"} or \code{"icl"}. If set to \code{"bic"}, model selection based on the Bayesian Information criterion is performed. If set to \code{"aic"}, the Akaike Information criterion is used. If set to \code{"icl"} the Integrated Classification Likelihood criterion is used. Default is \code{"bic"}.
#' @param se A string which can be set to \code{"analytic"} or \code{"bootstrap"}. If set to \code{"bootstrap"}, bootstrapped standard errors are reported. Default is \code{"analytic"}.
#' @param outer.runs A positive integer which stets the number of outer runs of the solver. Default is 1.
#' @param outer.tol A positive number which stets the tolerance of the continuation condition of the outer runs. The iterative algorithm stops if the relative decrease of the log-likelihood is smaller than \code{outer.tol}. Default is 0.
#' @param outer.max A positive integer which stets the maximum number of iterations of the outer runs of the solver. The iterative algorithm stops if it did not converge after \code{"outer.max"} iterations. Default is 1000.
#' @param inner.runs  A positive integer which stets the number of inner runs of the solver. Default is 10.
#' @param inner.tol A positive number which stets the tolerance of the continuation condition of the inner EM runs. The iterative algorithm stops if the relative decrease of the log-likelihood is smaller than \code{inner.tol}. Default is 0.
#' @param inner.max A positive integer which stets the maximum number of iterations of the inner EM runs. The iterative algorithm stops if it did not converge after \code{inner.max} iterations. Default is 10.
#' @param lcr.runs A positive integer which stets the number of estimation runs for latent class regression. Default is 100.
#' @param lcr.tol A positive number which stets the tolerance of the continuation condition of the Latent Class Regression runs. The iterative algorithm stops if the relative decrease of the log-likelihood is smaller than \code{lcr.tol}. Default is 0.
#' @param lcr.max A positive integer which stets the maximum number of iterations of the Latent Class Regression EM runs. The iterative algorithm stops if it did not converge after \code{lcr.max} iterations. Default is 1000.
#' @param bs.samples A positive integer which sets the number of bootstrap samples drawn with replacement.
#' @param quantiles A numeric vector indicating the quantiles of the sampling distribution of the estimated parameters. The quantiles are identified based on the standard error or based on bootstrapping the sampling distribution of the parameter.
#' @param stepsize A positive number which sets the stepsize of the Fisher scoring algorithm used to estimate the coefficients of the latent class regression model. Default is one. Values smaller than one slow down the convergence of the algorithm.
#' @param penalty A logical indicating if the Firth penalty is used to estimate the coefficients of the latent class regression model. Default is \code{FALSE}. Irrespective of the value specified here, the penalty is used in the case of a bootstrap of the standard errors of latent class regression coefficients.
#' @param verbose A logical, if \code{TRUE} messages of the estimation process and a summary of the estimated model is printed to the console. Default is \code{TRUE}.
#' @note The strategy estimation method was introduced by (Dal Bo & Frechette 2011) to estimate the relative frequency of a fixed set of pure strategies in the indefinitely repeated prisoner's dilemma. Breitmoser (2015) extended the method to the estimation of behavior strategies. The \pkg{stratEst} package uses the EM algorithm (Dempster, Laird & Rubin 1977) and the Newton-Raphson method to obtain maximum-likelihood estimates for the population shares and choice probabilities of a set of candidate strategies. The package builds on other software contributions of the R community. To increase speed the estimation procedures, the package uses integration of C++ and R achieved by the Rcpp package (Eddelbuettel & Francois 2011) and the open source linear algebra library for the C++ language RppArmadillo (Sanderson & Curtin 2016).
#' @return An object of class \code{stratEst}. A list with the following elements.
#' \item{strategies}{A list of fitted strategies.}
#' \item{shares}{Matrix of strategy shares. The order of rows corresponds to the order of strategies defined in the input object \code{strategies}.}
#' \item{probs}{Matrix of choice probabilities. The value \code{NA} indicates that the probability could not be estimated since data does not contain observations the model assigns to the corresponding state.}
#' \item{trembles}{Matrix of tremble probabilities of the strategies. The value \code{NA} indicates that the corresponding probability could not be estimated since data does not contain observations the model assigns to the corresponding state.}
#' \item{coefficients}{Matrix of latent class regression coefficients for strategies.}
#' \item{shares.par}{Estimated strategy shares.}
#' \item{probs.par}{Estimated choice probabilities.}
#' \item{trembles.par}{Estimated tremble probabilities.}
#' \item{coefficients.par}{Estimated latent class regression coefficients.}
#' \item{shares.indices}{Indices of strategy shares.}
#' \item{probs.indices}{Indices of choice probabilities.}
#' \item{trembles.indices}{Indices of tremble probabilities.}
#' \item{coefficients.indices}{Indices of latent class regression coefficients.}
#' \item{loglike}{The log-likelihood of the model. Larger values indicate a better fit of the model to the data.}
#' \item{crit.val}{The value of the selection criterion defined under \code{crit}. Larger values indicate a better fit of the model.}
#' \item{eval}{Number of iterations of the solver. The reported number is the sum of iterations performed in the inner and the outer run which produced the reported estimates.}
#' \item{tol.val}{The relative decrease of the log-likelihood in the last iteration of the algorithm. }
#' \item{convergence}{Maximum absolute score of the model parameters. Small values indicate convergence of the algorithm to a (local) maximum of the negative log likelihood.}
#' \item{entropy}{Entropy of the posterior probability assignments of individuals to strategies.}
#' \item{state.obs}{A column vector with the number of weighted observations for each strategy state corresponding to the rows of \code{strategies}.}
#' \item{posterior.assignments}{Posterior probability of each individual to use a strategy.}
#' \item{prior.assignments}{Prior probability of each individual to use a strategy as predicted by the individual covariates.}
#' \item{shares.se}{Standard errors of the estimated shares.}
#' \item{probs.se}{Standard errors of the estimated choice probabilities.}
#' \item{trembles.se}{Standard errors of the estimated trembles.}
#' \item{coefficients.se}{Standard errors of the estimated coefficients.}
#' \item{shares.score}{Score of the estimated shares.}
#' \item{probs.score}{Score of the reported choice probabilities.}
#' \item{trembles.score}{Score of the reported trembles.}
#' \item{coefficients.score}{Score of the reported coefficients.}
#' \item{shares.fisher}{Fisher information of the estimated shares.}
#' \item{probs.fisher}{Fisher information of the reported choice probabilities.}
#' \item{trembles.fisher}{Fisher information of the reported trembles.}
#' \item{coefficients.fisher}{Fisher information of the reported coefficients.}
#' \item{num.obs}{Number of observations.}
#' \item{num.ids}{Number of individuals.}
#' \item{num.par}{Total number of model parameters.}
#' \item{free.par}{Total number of free model parameters.}
#' \item{res.degrees}{Residual degrees of freedom (num.ids - free.par).}
#' \item{shares.quantiles}{Quantiles of the estimated shares.}
#' \item{probs.quantiles}{Quantiles of the estimated choice probabilities.}
#' \item{trembles.quantiles}{Quantiles of the estimated tremble probabilities.}
#' \item{coefficients.quantiles}{Quantiles of the estimated latent class regression coefficients.}
#' \item{gammas}{Gamma parameter of the model.}
#' \item{gammas.par}{Estimated gamma parameters.}
#' \item{gammas.se}{Standard errors of the gamma parameters.}#
#' \item{aic}{Akaike information criterion.}
#' \item{bic}{Bayesian information criterion.}
#' \item{icl}{Integrated classification likelihood information criteria.}
#' @description Performs variants of the strategy estimation method.
#' @details The estimation function \code{stratEst()} returns maximum-likelihood estimates for the population shares and choice probabilities of a set of candidate strategies given some data from an economic experiment. Candidate strategies can be supplied by the user in the form of deterministic finite-state automata. The number and the complexity of strategies can be restricted by the user or selected based on information criteria. stratEst also features latent class regression to assess the influence of covariates on strategy choice.
#' @references
#' Breitmoser, Y. (2015): Cooperation, but no reciprocity: Individual strategies in the repeated prisoner's dilemma, \emph{American Economic Review}, 105, 2882-2910.
#'
#' Dal Bo, P. and G. R. Frechette (2011): The evolution of cooperation in infinitely repeated games: Experimental evidence, \emph{American Economic Review}, 101, 411-429.
#'
#' Dempster, A., N. Laird, and D. B. Rubin (1977): Maximum likelihood from incomplete data via the EM algorithm," \emph{Journal of the Royal Statistical Society Series B}, 39, 1-38.
#'
#' Eddelbuettel, D. and R. Francois (2011): Rcpp: Seamless R and C++ Integration, \emph{Journal of Statistical Software}, 40, 1-18.
#'
#' Sanderson, C. and R. Curtin (2016): Armadillo: a template-based C++ library for linear algebra. \emph{Journal of Open Source Software}, 1-26.
#' @export
stratEst <- function( data , strategies , shares , coefficients , covariates , sample.id ,  response = "mixed" , sample.specific = c("shares","probs","trembles") , r.probs = "no" , r.trembles = "global" , select = NULL , min.strategies = 1 , crit = "bic" , se = "analytic" , outer.runs = 1 , outer.tol = 1e-10 , outer.max = 1000 , inner.runs = 10 , inner.tol = 1e-5 , inner.max = 10 , lcr.runs = 100 , lcr.tol = 1e-10 , lcr.max = 1000 , bs.samples = 1000 , quantiles = c(0.01,0.05,0.5,0.95,0.99) , stepsize = 1 , penalty = F , verbose = TRUE ){

  .Deprecated("stratEst.model")

  # check data
  if( missing(data) ) {
    stop("stratEst error: Mandatory input object 'data' missing.")
  }
  else{
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
  }

  # check sample.id
  if( missing(sample.id) ){
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
  if( missing(covariates) ){
    covariate_mat <- matrix(0,1,1)
    LCR = FALSE
  }
  else{
    covariate_mat <- stratEst.check.covariates( data , covariates )
    LCR = TRUE
  }

  # check other inputs
  stratEst.check.other.return <- stratEst.check.other( response , sample.specific , r.probs , r.trembles , select , min.strategies , crit , se , outer.runs , outer.tol , outer.max , inner.runs , inner.tol , inner.max , lcr.runs , lcr.tol , lcr.max , bs.samples , stepsize , penalty , verbose , quantiles )
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
  if( missing(shares) ) {
    shares = matrix( NA , num_strats , num_samples )
  }else{
    shares <- stratEst.check.shares( shares , LCR , specific_shares , num_samples , num_strats , sample.id , sample_levels , select_strategies )
  }

  #check coefficients
  if( missing(coefficients) ) {
    if( LCR ){
      coefficient_mat = matrix(NA,length(covariates),num_strats,1)
      coefficient_mat[,1] = 0
    }else{
      coefficient_mat = matrix(0,1,1)
    }
  }else{
    if( is.null( coefficients ) ){
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

  #########################################################################################################
  # ESTIMATION
  #########################################################################################################

  # create cpp output
  cpp.output <- stratEst_cpp( data_matrix, strategies_matrix, sid , shares , trembles , coefficient_mat, covariate_mat, LCR, cluster, quantile_vec , response, specific_shares , specific_responses , specific_trembles , specific_coefficients , r.probs, r.trembles, select_strategies , select_responses , select_trembles , min.strategies, crit, se, outer.runs, outer.tol, outer.max, inner.runs, inner.tol, inner.max, lcr.runs, lcr.tol, lcr.max, bs.samples , print.messages, integer_strategies, stepsize , penalty )

  # stratEst.return object
  stratEst.return <- list("strategies" = cpp.output$strategies, "shares" = cpp.output$shares, "coefficients" = cpp.output$coefficients, "probs" = cpp.output$responses, "trembles" = cpp.output$trembles, "shares.par" = cpp.output$shares.list$shares.par, "probs.par" = cpp.output$responses.list$responses.par, "trembles.par" = cpp.output$trembles.list$trembles.par, "coefficients.par" =  cpp.output$coefficients.list$coefficients.par,  "shares.indices" = cpp.output$shares.list$shares.indices, "probs.indices" = cpp.output$responses.list$responses.indices, "trembles.indices" = cpp.output$trembles.list$trembles.indices, "coefficients.indices" = cpp.output$coefficients.list$coefficients.indices, "loglike" = cpp.output$fit[1,1], "crit.val" = cpp.output$fit[1,2], "eval" = cpp.output$solver[1,1], "tol.val" = cpp.output$solver[1,2], "convergence" = cpp.output$convergence, "entropy" = cpp.output$fit[1,3], "state.obs" = cpp.output$state.obs, "posterior.assignment" = cpp.output$assignments, "prior.assignment" = cpp.output$priors, "shares.se" = cpp.output$shares.list$shares.se, "probs.se" = cpp.output$responses.list$responses.se, "trembles.se" = cpp.output$trembles.list$trembles.se, "coefficients.se" = cpp.output$coefficients.list$coefficients.se, "shares.covar" = cpp.output$shares.list$shares.covar, "shares.score" =  cpp.output$shares.list$shares.score, "shares.fisher" = cpp.output$shares.list$shares.fisher, "probs.covar" = cpp.output$responses.list$responses.covar, "probs.score" = cpp.output$responses.list$responses.score, "probs.fisher" = cpp.output$responses.list$responses.fisher, "trembles.covar" = cpp.output$trembles.list$trembles.covar, "trembles.score" = cpp.output$trembles.list$trembles.score, "trembles.fisher" = cpp.output$trembles.list$trembles.fisher, "coefficients.covar" = cpp.output$coefficients.list$coefficients.covar, "coefficients.score" = cpp.output$coefficients.list$coefficients.score, "coefficients.fisher" = cpp.output$coefficients.list$coefficients.fisher );

  # post-processing
  stratEst.return <- stratEst.post( cpp.output , stratEst.return , strategies , covariates , response , unique_ids , num_unique_ids , input , output , unique_inputs , unique_outputs , num_unique_inputs , num_unique_outputs , sample , sample.id , sample_factor , num_samples , specific_shares , specific_responses , specific_trembles , sample_is_factor , integer_strategies , cpp.output$lcr , response_mat_col_index , crit , num_obs , se , quantile_vec )

  class(stratEst.return) <- c("stratEst.model","list")

  # print summary
  if( print.summary ){
    summary( stratEst.return )
  }

  # return result
  return(stratEst.return)

}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
