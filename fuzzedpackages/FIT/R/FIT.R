#' FIT: a statistical modeling tool for transcriptome dynamics under fluctuating field conditions
#'
#'
#' Provides functionality for constructing statistical models of transcriptomic dynamics in field 
#' conditions. It further offers the function to predict expression of a gene given the attributes 
#' of samples and meteorological data. Nagano, A. J., Sato, Y., Mihara, M., Antonio, B. A., 
#' Motoyama, R., Itoh, H., Naganuma, Y., and Izawa, T. (2012). <doi:10.1016/j.cell.2012.10.048>. 
#' Iwayama, K., Aisaka, Y., Kutsuna, N., and Nagano, A. J. (2017). 
#' <doi:10.1093/bioinformatics/btx049>.
#' 
#' @references
#' [Nagano et al.] A.J.~Nagano, et al.
#' ``Deciphering and prediction of transcriptome dynamics under fluctuating field conditions,''
#' Cell~151, 6, 1358--69 (2012)
#'
#' [Iwayama] K.~Iwayama, et al. 
#' ``FIT: statistical modeling tool for transcriptome dynamics under fluctuating field conditions,''
#' Bioinformatics, btx049 (2017) 
#'
#' @section Overview:
#' The \pkg{FIT} package is an \code{R} implementation
#' of a class of transcriptomic models that
#' relates gene expressions of plants and weather conditions to which
#' the plants are exposed.
#' (The reader is referred to [Nagano et al.] for the detail of
#' the class of models concerned.)
#'
#' By providing
#' (a) gene expression profiles of plants brought up in a field condition,
#' and (b) the relevant weather history (temperature etc.) of the said field,
#' the user of the package is able to
#' (1) construct optimized models (one for each gene) for their expressions,
#' and
#' (2) use them to predict the expressions for another weather history
#' (possibly in a different field).
#'
#' Below, we briefly explain
#' the construction of the optimized models (``training phase'')
#' and the way to use them to make predictions (``prediction phase'').
#'
#' \subsection{Model training phase}{
#'
#' The model of [Nagano et al.] belongs to the class of statistical models
#' called ``linear models''
#' and are specified by a set of ``parameters'' and
#' ``(linear regression) coefficients''.
#' The former are used to convert weather conditions to
#' the ``input variables'' for a regression, and the latter are then
#' multiplied to the input variables to form the expectation values
#' for the gene expressions.
#' The reader is referred to the original article [Nagano et al.]
#' for the formulas for the input variables.
#' (See also [Iwayama] for a review.)
#'
#' The training phase consists of three stages:
#' \enumerate{
#' \item \code{Init}: fixes the initial model parameters
#' \item \code{Optim}: optimizes the model parameters
#' \item \code{Fit}: fixes the linear regression coefficients
#' }
#' The user can configure the training phase
#' through a custom data structure (``recipe''),
#' which can be constructed by using the utility function
#' \code{FIT::make.recipe()}.
#'
#' The role of the first stage \code{Init} is to fix the initial values
#' for the model parameters from which the parameter optimization is performed.
#' At the moment two methods, \code{'manual'} and \code{'gridsearch'},
#' are implemented.
#' With the \code{'manual'} method the user can simply specify the set of
#' initial values that he thinks is promising.
#' For the \code{'gridsearch'} method the user discretizes
#' the parameter space to a grid by providing
#' a finite number of candidate values for each parameter.
#' \pkg{FIT} then performs a search over the grid
#' for the ``best'' combinations of the initial parameters.
#' % In both cases relevant data are passed through \code{init.data}.
#'
#' The second stage \code{Optim} is the main step of the model training,
#' and \pkg{FIT} tries to gradually improve the model parameters
#' using the Nelder-Mead method.
#'
#' This stage could be run one or more times where each can be run
#' using the method \code{'none'}, \code{'lm'} or \code{'lasso'}.
#' The \code{'none'} method passes the given parameter as-is
#' to the next method in the \code{Optim} pipeline or to the next stage \code{Fit}.
#' (Basically, the method is there so that the user can skip the entire
#' \code{Optim} stage, but the method could be used for slightly warming-up the CPU as well.)
#'
#' The \code{'lm'} method uses the a simple (weighted) linear regression to
#' guide the parameter optimization. That is, \pkg{FIT}
#' first computes the ``input variables'' from the current parameters and
#' associated weather data, and then finds the set of linear coefficients
#' that best explains the ``output variables'' (gene expressions).
#' Finally, the quadratic residual is used as the measure for the
#' error and is fed back to the Nelder-Mead method.
#'
#' The \code{'lasso'} method is similar to the \code{'lm'} method
#' but uses the (weighted) Lasso regression
#' (``linear'' regression with an L1-regularization for the regression coefficients)
#' instead of the simple linear regression.
#' \pkg{FIT} uses the \pkg{glmnet} package to perform
#' the Lasso regression and the strength of the L1-regularization
#' is fixed via a cross validation. (See \code{cv.glmnet()} from the \pkg{glmnet}
#' package.
#' The Lasso regression is said to suppress irrelevant input variables automatically
#' and tends to create models with better prediction ability.
#' On the other hand, \code{'lasso'} runs considerably slower than \code{'lm'}.
#'
#' For example, passing a vector \code{c('lm', 'lasso')} to the
#' argument \code{optim} (of \code{make.recipe()}) creates a recipe
#' that instructs the \code{Optim} stage to
#' (1) first optimize using the \code{'lm'} method,
#' (2) and then fine tunes the parameters using the \code{'lasso'} method.
#'
#' After fixing the model parameters in the \code{Optim} stage,
#' the \code{Fit} stage can be used to fix the linear coefficients
#' of the models.
#' Here, either \code{'fit.lm'} or \code{'fit.lasso'} can be used
#' to find the ``best'' coefficients, the main difference being that
#' the coefficients are penalized by an L1-norm for the latter.
#' Note that it is perfectly okay to use \code{'fit.lasso'} for
#' the parameters optimized using \code{'lm'}.
#'
#' In order to prepare for the possibly huge variations
#' of expression data as measured by RNA-seq,
#' \pkg{FIT} provides a way to weight regression penalties from each sample
#' with different weights as in
#' \code{sum_{s in samples} (weight_s) (error_s)^2}.
#'
#' } % subsection model training
#'
#' \subsection{Prediction phase}{
#' For each gene, the trained model of the previous subsection
#' can be thought of as a black box that maps
#' the field conditions (weather data),
#' to which a plant containing the gene is exposed,
#' to its expected expression.
#' \pkg{FIT} provides a simple function
#' \code{FIT::predict()} that does just this.
#' 
#' \code{FIT::predict()} takes as its argument
#' a list of pretrained models
#' as well as actual/hypothetical plant sample attributes and weather data,
#' and returns the predicted values of gene expressions.
#'
#' When there is a set of actually measured expressions,
#' an associated function \code{FIT::prediction.errors()})
#' can be used to check the validity of the predictions made by
#' the models.
#' } % subsection prediction phase
#'
#' @section Namespece contamination:
#' The \pkg{FIT} package exports fairly ubiquitous names
#' auch as \code{optim}, \code{predict} etc.\ as its API.
#' Users, therefore, are advised to load \pkg{FIT}
#' via \code{requireNamespace('FIT')} and use its API function with
#' a namaspace qualifier (e.g.~\code{FIT::optim()})
#' rather than loading \emph{and} attaching it via \code{library('FIT')}.
#' 
#' @section Basic usage:
#' See vignettes for examples of actual scripts that use \pkg{FIT}.
#'
#' @examples
#' \dontrun{
#' # The following snippet shows the structure of a typical
#' # driver script of the FIT package.
#' # See vignettes for examples of actual scripts that use FIT.
#'
#' ##############
#' ## training ##
#' ##############
#' ## discretized parameter space (for 'gridsearch')
#' grid.coords <- list(
#'   clock.phase = seq(0, 23*60, 1*60),
#'   # :
#'   gate.radiation.amplitude = c(-5, 5)
#' )
#' 
#' ## create a training recipe
#' recipe <- FIT::make.recipe(c('temperature', 'radiation'),
#'                            init  = 'gridsearch',
#'                            init.data = grid.coords,
#'                            optim = c('lm'),
#'                            fit   = 'fit.lasso',
#'                            time.step = 10, 
#'                            opts =
#'                              list(lm    = list(maxit = 900),
#'                              lasso = list(maxit = 1000))
#'                            )
#' 
#' ## names of genes to construct models
#' genes <- c('Os12g0189300', 'Os02g0724000')
#' 
#' }
#' 
#' 
#' \dontrun{
#' ## load training data
#' training.attribute  <- FIT::load.attribute('attribute.2008.txt')
#' training.weather    <- FIT::load.weather('weather.2008.dat', 'weather')
#' training.expression <- FIT::load.expression('expression.2008.dat', 'ex', genes)
#' 
#' ## models will be a list of trained models (length: ngenes)
#' models <- FIT::train(training.expression,
#'                      training.attribute,
#'                      training.weather,
#'                      recipe)
#' 
#' }
#' 
#' ################
#' ## prediction ##
#' ################
#'
#' \dontrun{
#' ## load validation data
#' prediction.attribute  <- FIT::load.attribute('attribute.2009.txt');
#' prediction.weather    <- FIT::load.weather('weather.2009.dat', 'weather')
#' prediction.expression <- FIT::load.expression('expression.2009.dat', 'ex', genes)
#' 
#' ## predict
#' prediction.result <- FIT::predict(models[[1]],
#'                                  prediction.attribute,
#'                                  prediction.weather)
#'
#'
#'}
#'
#' @docType package
#' @name FIT
NULL
#> NULL

################################################################
###
#' @useDynLib FIT
#' @importFrom Rcpp sourceCpp
NULL
#> NULL

# cleanup
.onUnload <- function (path) {
  library.dynam.unload("FIT", path)
}

######################################################################
### Namespaces for submodules

# internal use only
Model <- new.env()
Train <- new.env()
IO    <- new.env()
Norm  <- new.env()
Jma   <- new.env()

######################################################################
### Enduser API functions for model construction (training) and prediction

#' Supported weather factors.
#' @examples
#' length(FIT::weather.entries)
#' @export
weather.entries <- c('wind', 'temperature', 'humidity',
                     'atmosphere', 'precipitation', 'radiation')

#' Creates a recipe for training models.
#'
#' @param envs An array of weather factors to be taken into account
#'     during the construction of models.
#'     At the moment, the array \code{envs} can only contain a single weather factor
#'     from \code{weather.entries}, though there is a plan to remove the restriction
#'     in a future version.
#' @param init A string to specify the method to choose the initial parameters.
#'     (One of \code{'gridsearch'} or \code{'manual'}.)
#' @param optim A string to specify the method to be used for optimizing
#'     the model parameters.
#'     (One of \code{'none'}, \code{'lm'} or \code{'lasso'})
#' @param fit A string to specify the method to be used for fixing
#'     the linear regression coefficients.
#'     (One of \code{'fit.lm'} or \code{'fit.lasso'}.)
#' @param init.data Auxiliary data needed to perform the Init stage
#'     using the method specified by the \code{init} argument.
#'     When \code{init} is \code{'gridsearch'}, it should be a list representing
#'     a discretized parameter space.
#'     When \code{init} is \code{'manual'}, it should be a list of parameter
#'     values that is used as the initial values for the parameters in
#'     the Optim stage.
#' @param time.step An integer to specify the basic unit of time (in minute)
#'     for the transcriptomic models.
#'     Must be a multiple of the time step of weather data.
#' @param gate.open.min The minimum opning length in minutes of the gate function for 
#'     environmental inputs. 
#' @param opts An optional named list that specifies the arguments to be passed
#'     to methods that constitute each stage of the model training.
#'     Each key of the list corresponds to a name of a method.
#'
#'     See examples for the supported options.
#' @return An object representing the procedure to construct models.
#' @examples
#' \dontrun{
#' init.params <- .. # choose them wisely
#' # Defined in Train.R:
#' # default.opts <- list(
#' #   none  = list(),
#' #   lm    = list(maxit=1500, nfolds=-1), # nfolds for lm is simply ignored
#' #   lasso = list(maxit=1000, nfolds=10)
#' # )
#' recipe <- FIT::make.recipe(c('wind', 'temperature'),
#'                            init = 'manual',
#'                            init.data = init.params,
#'                            optim = c('lm', 'none', 'lasso'),
#'                            fit = 'fit.lasso',
#'                            time.step = 10,
#'                            opts =
#'                              list(lm    = list(maxit = 900),
#'                                   lasso = list(maxit = 1000)))
#' }
#'
#' @export
make.recipe <- function(envs, init, optim, fit, init.data, time.step,
                        gate.open.min = 0, opts = NULL) {
  all <- weather.entries
  if (!all(vapply(envs, function(o) o %in% all, TRUE))) stop('some envs are invalid: ', envs)
  if (length(time.step) != 1) stop('multiple time.step forbidden: ', time.step)
  if (init != 'gridsearch' && init != 'manual') stop('invalid init method: ', init)
  if (!all(vapply(optim, function(o) o %in% c('none', 'lm', 'lasso'), TRUE)))
    stop('some optim methods are invalid: ', optim)
  if (fit != 'fit.lm' && fit != 'fit.lasso') stop('invalid fitting method: ', fit)
  if (gate.open.min < 0 || gate.open.min > 1440) stop('invalid gate.open.min')
  if (is.null(opts)) opts <- list()
  
  Model$Recipe(envs=envs, init=init, optim=optim, fit=fit, init.data=init.data,
               time.step=as.integer(time.step), gate.open.min=gate.open.min, opts=opts)
}

### Notation: `env` runs over `envs`; `e` runs over entries of an `env`

#' Constructs models following a recipe.
#'
#' @param expression An object that represents gene expression data.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * ngenes}
#'     using \code{FIT::load.expression()}.
#'     (At the moment it is an instance of a hidden class IO$Expression,
#'     but this may be subject to change.)
#' @param attribute An object that represents the attributes of
#'     microarray/RNA-seq data.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * nattributes}
#'     using \code{FIT::load.attribute()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param weather An object that represents actual or hypothetical weather data
#'     with which the training of models are done.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{ntimepoints * nfactors}
#'     using \code{FIT::load.weather()}.
#'     (At the moment it is an instance of a hidden class IO$Weather,
#'     but this may be subject to change.)
#' @param recipe An object that represents the training protocol of models.
#'     A recipe can be created using \code{FIT::make.recipe()}.
#' @param weight An optional numerical matrix of size \code{nsamples * ngenes}
#'     that during regression penalizes errors from each sample 
#'     using the formula
#'     \code{sum_{s in samples} (weight_s) (error_s)^2}.
#' 
#'     This argument is optional for a historical reason,
#'     and when it is omitted, all samples are equally penalized.
#' @param min.expressed.rate A number used to 
#'   A gene with \code{var(expr) < thres.expr} is regarded as unexpressed,
#'   and \pkg{FIT} sets its model as: \code{expr = log(offset) + 0*inputs}.
#' @return A collection of trained models.
#'
#' @examples
#' \dontrun{
#' # create recipe
#' recipe <- FIT::make.recipe(..)
#'
#' #load training data
#  genes <- c('Os01g0182600', 'Os02g0618200')
#' training.attribute  <- FIT::load.attribute('attribute.2008.txt');
#' training.weather    <- FIT::load.weather('weather.2008.dat', 'weather')
#' training.expression <- FIT::load.expression('expression.2008.dat', 'ex', genes)
#' training.weight     <- FIT::load.weight('weight.2008.dat', 'weight', genes)
#'
#' # train models
#' models <- FIT::train(training.expression,
#'                      training.attribute,
#'                      training.weather,
#'                      recipe,
#'                      training.weight)
#' }
#' @export
train <- function(expression, attribute, weather, recipe, weight = NULL, min.expressed.rate = 0.01) {
  genes <- expression$entries
  exprs <- expression$rawdata
  if (!all(genes == colnames(exprs))) stop('inconsistent expression data')

  samples.n <- nrow(exprs)
  # genes.n   <- ncol(exprs)

  if (nrow(attribute$data) != samples.n) stop('inconsistent attribute data')
  if (is.null(weight)) weight <- IO$trivialWeights(samples.n, genes)
  if (!all(dim(expression$rawdata) == dim(weight$rawdata))) stop('inconsistnet weight data')
  weights <- weight$rawdata

  cat('# * Training..\n')
  if (recipe$time.step %% weather$data.step != 0)
    stop('recipe$time.step (= ', recipe$time.step,
         ') must be an integral multiple of weather$data.step, (= ',
         weather$data.step, ')')
  cat('# ** Prep+Init:\n')
  models <- Train$init(exprs, weights, attribute$data, weather$data,
                       recipe$envs, recipe$init, recipe$init.data,
                       weather$data.step, recipe$time.step)

  cat('# ** Optim ('); cat(recipe$optim, sep=', '); cat('):\n')
  os <- recipe$optim
  for (o in os)
    models <- Train$optim(exprs, weights, attribute$data, weather$data, models, o,
                          weather$data.step, recipe$time.step,
                          recipe$opts[[o]]$maxit, recipe$opts[[o]]$nfolds,
                          min.expressed.rate, 
                          recipe$gate.open.min)

  cat('# ** Creating optimized models\n')
  models <- Train$fit(exprs, weights, attribute$data, weather$data, models, recipe$fit,
                      weather$data.step, recipe$time.step)
  cat('# Done (training)\n')

  names(models) <- genes # should be unnecessary, but for document purpose..
  models
}

################################
## From this layer down,
## (1) weights are non-optional, and
## (2) placed next to exprs
##
## Users are not recommended to use them directly.

#' A raw API for initializing model parameters.
#'
#' Note: use \code{train()} unless the user is willing to
#' accept breaking API changes in the future.
#'
#' @param expression An object that represents gene expression data.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * ngenes}
#'     using \code{FIT::load.expression()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param weight A matrix of size \code{nsamples * ngenes}
#'     that during regression penalizes errors from each sample 
#'     using the formula
#'     \code{sum_{s in samples} (weight_s) (error_s)^2}.
#'
#'     Note that, unlike for \code{FIT::train()}, this argument
#'     is NOT optional.
#' @param attribute An object that represents the attributes of a
#'     microarray/RNA-seq data.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * nattributes}
#'     using \code{FIT::load.attribute()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param weather An object that represents actual or hypothetical weather data
#'     with which the training of models are done.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{ntimepoints * nfactors}
#'     using \code{FIT::load.weather()}.
#'     (At the moment it is an instance of a hidden class IO$Weather,
#'     but this may be subject to change.)
#' @param recipe An object that represents the training protocol of models.
#'     A recipe can be created using \code{FIT::make.recipe()}.
#' @return A collection of models whose parameters are
#'     set by using the \code{'init'} method in the argument \code{recipe}.
#' @export
init <- function(expression, weight, attribute, weather, recipe) {
  Train$init(expression$rawdata, weight$rawdata, attribute$data, weather$data,
             recipe$envs, recipe$init, recipe$init.data,
             weather$data.step, recipe$time.step)
}

#' A raw API for optimizing model parameters.
#'
#' Note: use \code{train()} unless the user is willing to
#' accept breaking API changes in the future.
#'
#' @param expression An object that represents gene expression data.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * ngenes}
#'     using \code{FIT::load.expression()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param weight A matrix of size \code{nsamples * ngenes}
#'     that during regression penalizes errors from each sample 
#'     using the formula
#'     \code{sum_{s in samples} (weight_s) (error_s)^2}.
#'
#'     Note that, unlike for \code{FIT::train()}, this argument
#'     is NOT optional.
#' @param attribute An object that represents the attributes of
#'     microarray/RNA-seq data.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * nattributes}
#'     using \code{FIT::load.attribute()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param weather An object that represents actual or hypothetical weather data
#'     with which the training of models are done.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{ntimepoints * nfactors}
#'     using \code{FIT::load.weather()}.
#'     (At the moment it is an instance of a hidden class IO$Weather,
#'     but this may be subject to change.)
#' @param recipe An object that represents the training protocol of models.
#'     A recipe can be created using \code{FIT::make.recipe()}.
#' @param models A collection of models being trained as is returnd by
#'     \code{FIT::init()}.
#'
#'     At this moment, it must be a list (genes) of a list (envs) of models
#'     and must contain at least one model.
#'     (THIS MIGHT CHANGE IN A FUTURE.)
#' @param maxit An optional number that specifies
#'     the maximal number of times that the parameter optimization is performed.
#'
#'     The user can control this parameter by using the \code{opts} argument
#'     for \code{FIT::train()}.
#' @param nfolds An optional number that specifies the order of
#'     cross validation when \code{optim} method is \code{'lasso'}.
#'     This is simply ignored when \code{optim} method is \code{'lm'}.
#' @return A collection of models whose parameters are
#'     optimized by using the \code{'optim'} pipeline
#'     in the argument \code{recipe}.
#' @export
optim <- function(expression, weight, attribute, weather, recipe,  models,
                  maxit = NULL, nfolds = NULL) {
  # DOC: 'models' must be a list (genes) of a list (envs) of models
  # and must contain at least one model.
  os  <- recipe$optim
  log <- models[[1]][[1]][['log']][-1]
  log <- log[log != 'none']
  if (length(log) < length(os)) {
    Train$optim(expression$rawdata, weight$rawdata, attribute$data, weather$data,
                models, os[[length(log)+1]],
                weather$data.step, recipe$time.step,
                maxit, nfolds, 
                gate.open.min = recipe$gate.open.min)
  } else {
    models
  }
}

#' A raw API for fixing linear regression coefficients.
#'
#' Note: use \code{train()} unless the user is willing to
#' accept breaking API changes in the future.
#'
#' @param expression An object that represents gene expression data.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * ngenes}
#'     using \code{FIT::load.expression()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param weight A matrix of size \code{nsamples * ngenes}
#'     that during regression penalizes errors from each sample 
#'     using the formula
#'     \code{sum_{s in samples} (weight_s) (error_s)^2}.
#'
#'     Note that, unlike for \code{FIT::train()}, this argument
#'     is NOT optional.
#' @param attribute An object that represents the attributes of
#'     microarray/RNA-seq data.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * nattributes}
#'     using \code{FIT::load.attribute()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param weather An object that represents actual or hypothetical weather data
#'     with which the training of models are done.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{ntimepoints * nfactors}
#'     using \code{FIT::load.weather()}.
#'     (At the moment it is an instance of a hidden class IO$Weather,
#'     but this may be subject to change.)
#' @param recipe An object that represents the training protocol of models.
#'     A recipe can be created using \code{FIT::make.recipe()}.
#' @param models A collection of models being trained as is returnd by
#'     \code{FIT::optim()}.
#' @return A collection of models whose parameters and regression coeffients
#'     are optimized.
#' @export
fit.models <- function(expression, weight, attribute, weather, recipe, models) {
  Train$fit(expression$rawdata, weight$rawdata, attribute$data, weather$data,
            models, recipe$fit,
            weather$data.step, recipe$time.step)
}

################################################################
#' Predicts gene expressions using pretrained models.
#'
#' @param models A list of trained models for the genes of interest.
#' 
#'     At the moment the collection of trained models returned
#'     by \code{FIT::train()} cannot be directly passed to \code{FIT::predict()}:
#'     the user has to explicitly convert it to an appropriate format by using
#'     \code{FIT::train.to.predict.adaptor()}.
#'     (This restriction might be removed in a future.) 
#' @param attribute An object that represents the attributes of
#'     microarray/RNA-seq data.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * nattributes}
#'     using \code{FIT::load.attribute()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param weather An object that represents actual or hypothetical weather data
#'     with which predictions of gene expressions are made.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{ntimepoints * nfactors}
#'     using \code{FIT::load.weather()}.
#'     (At the moment it is an instance of a hidden class IO$Weather,
#'     but this may be subject to change.)
#' @return A list of prediction results as returned by the models.
#'
#' @examples
#' \dontrun{
#' # prepare models
#' # NOTE: FIT::train() returns a nested list of models
#' #   so we have to flatten it using FIT::train.to.predict.adaptor()
#' #   before passing it to FIT::predict().
#' models <- FIT::train(..)
#' models.flattened <- FIT::train.to.predict.adaptor(models)
#'
#' # load data used for prediction
#' prediction.attribute  <- FIT::load.attribute('attribute.2009.txt')
#' prediction.weather    <- FIT::load.weather('weather.2009.dat', 'weather')
#' prediction.expression <- FIT::load.expression('expression.2009.dat', 'ex', genes)
#'
#' prediction.results <- FIT::predict(models.flattened,
#'                                    prediction.attribute,
#'                                    prediction.weather)
#' }
#' @export
predict <- function(models, attribute, weather) {
  lapply(models, function(m) m$predict(attribute$data, weather$data, weather$data.step))
}

#' Computes the prediction errors using the trained models.
#'
#' @param models A list of trained models for the genes of interest.
#' 
#'     At the moment the collection of trained models returned
#'     by \code{FIT::train()} cannot be directly passed to \code{FIT::predict()}:
#'     the user has to explicitly convert it to an appropriate format by using
#'     \code{FIT::train.to.predict.adaptor()}.
#'     (This restriction might be removed in a future.) 
#' @param expression An object that represents the actual measured data of
#'     gene expressions.
#'     The object can be created from a dumped/saved dataframe
#'     of size \code{nsamples * ngenes}
#'     using \code{FIT::load.expression()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param attribute An object that represents the attributes of
#'     microarray/RNA-seq data.
#'     The object can be created from a dumped/saved dataframe
#'     using \code{FIT::load.attribute()}.
#'     (At the moment it is an instance of a hidden class IO$Attribute,
#'     but this may be subject to change.)
#' @param weather An object that represents actual or hypothetical weather data
#'     with which predictions of gene expressions are made.
#'     The object can be created from a dumped/saved dataframe
#'     using \code{FIT::load.weather()}.
#'     (At the moment it is an instance of a hidden class IO$Weather,
#'     but this may be subject to change.)
#' @return A list of deviance (a measure of validity of predictions,
#'     as is defined by each model) between the prediction results
#'     and the measured results (as is provided by the user through
#'     \code{expression} argument).
#' @examples
#' \dontrun{
#' # see the usage of FIT::predict()
#' }
#' @export
prediction.errors <- function(models, expression, attribute, weather) {
  lapply(models, function(m) m$deviance(expression$rawdata, attribute$data, weather$data, weather$data.step))
}

################################################################
### reexport some IO stuff

#' Converts attribute data from a dataframe into an object. 
#' 
#' @param data A dataframe of the attributes of microarray/RNA-seq data.
#' @param sample An optional numeric array that designates
#'     the samples, that is rows, of the dataframe to be loaded.
#' @return An object that represents the attributes of
#'     microarray/RNA-seq data.
#'     Internally, the object holds a dataframe whose number of entries
#'     (rows) equals that of the samples.
#' @export
convert.attribute <- function(data, sample = NULL)
  IO$Attribute$new(data, sample)

#' Loads attribute data.
#'
#' @param path A path of a file that contains attribute data to be loaded.
#'     When the file is a loadable \code{.Rdata},
#'     \code{name} of the dataframe object in the \code{.Rdata}
#'     (that actually contains the relevant data)
#'     has to be specified as well.
#' @param variable An optional string that designates the name of a
#'     dataframe object that has been saved in an \code{.Rdata}.
#'     (See the description of \code{path}.)
#' @param sample An optional numeric array that designates
#'     the samples, that is rows, of the dataframe to be loaded.
#' @return An object that represents the attributes of
#'     microarray/RNA-seq data.
#'     Internally, the object holds a dataframe whose number of entries
#'     (rows) equals that of the samples.
#'
#' @export
load.attribute <- function(path, variable = NULL, sample = NULL){
  file.data <- IO$slurp(path, variable)
  IO$Attribute$new(file.data, sample)
}

#' converts expression data from a dataframe into an object.
#'
#' @param data A dataframe of expression data to be loaded.
#' @param entries An optional string array that designates
#'     the entries of the dataframe to be loaded.
#' @return An object that represents the expression data of microarray/RNA-seq.
#'     Internally, the object holds a matrix of size
#'     \code{nsamples * ngenes}.
#' @export
convert.expression <- function(data, entries = NULL)
  IO$Expression$new(data, entries)
#' Loads expression data.
#'
#' @param path A path of a file that contains attribute data to be loaded.
#'     When the file is a loadable \code{.Rdata},
#'     \code{name} of the dataframe object in the \code{.Rdata}
#'     (that actually contains the relevant data)
#'     has to be specified as well.
#' @param variable An optional string that designates the name of a
#'     dataframe object that has been saved in an \code{.Rdata}.
#'     (See the description of \code{path}.)
#' @param entries An optional string array that designates
#'     the entries of the dataframe to be loaded.
#' @return An object that represents the expression data of microarray/RNA-seq.
#'     Internally, the object holds a matrix of size
#'     \code{nsamples * ngenes}.
#' @export
load.expression <- function(path, variable = NULL, entries = NULL){
  file.data <- IO$slurp(path, variable)
  IO$Expression$new(file.data, entries)
}

#' Converts weather data from a dataframe into an object.
#'
#' @param data A dataframe of weather data to be converted.
#' @param entries An optional string array that designates
#'     the entries of the dataframe to be loaded.
#' @return An object that reprents the timeseries data of weather factors.
#'     Internally, the object holds a dataframe of size
#'     \code{ntimepoints * nfactors}.
#' @export
convert.weather <- function(data, entries = IO$weather.entries)
  IO$Weather$new(data, entries)
  
#' Loads weather data.
#'
#' @param path A path of a file that contains weather data to be loaded.
#'     When the file is a loadable \code{.Rdata},
#'     \code{name} of the dataframe object in the \code{.Rdata}
#'     (that actually contains the relevant data)
#'     has to be specified as well.
#' @param variable An optional string that designates the name of a
#'     dataframe object that has been saved in an \code{.Rdata}.
#'     (See the description of \code{path}.)
#' @param entries An optional string array that designates
#'     the entries of the dataframe to be loaded.
#' @return An object that reprents the timeseries data of weather factors.
#'     Internally, the object holds a dataframe of size
#'     \code{ntimepoints * nfactors}.
#' @export
load.weather <- function(path, variable = NULL, entries = IO$weather.entries){
  file.data <- IO$slurp(path, variable)
  IO$Weather$new(file.data, entries)
}

#' Converts regression weight data from a dataframe into an object.
#'
#' @param data A dataframe that contains weight data to be loaded.
#' @param entries An optional string array that designates
#'     the entries of the dataframe to be loaded.
#' @return An object that represents the weights 
#'     Internally, the object holds a matrix of size
#'     \code{nsamples * ngenes}.
#' @export
convert.weight <- function(data, entries = NULL)
  IO$Weights$new(data, entries)
#' Loads regression weight data.
#'
#' @param path A path of a file that contains weight data to be loaded.
#'     When the file is a loadable \code{.Rdata},
#'     \code{name} of the dataframe object in the \code{.Rdata}
#'     (that actually contains the relevant data)
#'     has to be specified as well.
#' @param variable An optional string that designates the name of a
#'     dataframe object that has been saved in an \code{.Rdata}.
#'     (See the description of \code{path}.)
#' @param entries An optional string array that designates
#'     the entries of the dataframe to be loaded.
#' @return An object that represents the weights 
#'     Internally, the object holds a matrix of size
#'     \code{nsamples * ngenes}.
#' @export
load.weight <- function(path, variable = NULL, entries = NULL){
  file.data <- IO$slurp(path, variable)
  IO$Weights$new(file.data, entries)
}
#' Makes trivial weight data
#' 
#' @param samples.n A number of samples. 
#' @param genes A list of genes. 
#' @return An object that represens the trivial weights. 
#'        Internally, the object holds an identity matrix of size 
#'        \code{nsamples * ngenes}.
#' @export
make.trivial.weights <- function(samples.n, genes)
  IO$trivialWeights(samples.n, genes)
