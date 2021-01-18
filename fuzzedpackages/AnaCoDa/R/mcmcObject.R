#' Initialize MCMC 
#' 
#' @param samples Number of samples to be produced when running the 
#' MCMC algorithm. No default value.
#' 
#' @param thinning The thinning interval between consecutive observations. If set to 
#' 1, every step will be saved as a sample. Default value is 1.
#' 
#' @param adaptive.width Number that determines how often the acceptance/rejection
#' window should be altered. Default value is 100 samples.
#' Proportion of MCMC steps where the proposal distribution is adaptive can be set using \code{mcmc$setStepsToAdapt}. The default parameter passed in as -1 uses the full iterations.
#' 
#' @param est.expression Boolean that tells whether or not synthesis rate values
#' should be estimated in the MCMC algorithm run. Default value is TRUE.
#' 
#' @param est.csp Boolean that tells whether or not codon specific values
#' should be estimated in the MCMC algorithm run. Default value is TRUE.
#' 
#' @param est.hyper Boolean that tells whether or not hyper parameters
#' should be estimated in the MCMC algorithm run. Default value is TRUE.
#' Setting for expression noise parameter sepsilon can be overridden by setting \code{fix.observation.noise} in \code{initializeModelObject()}
#' 
#' @param est.mix Boolean that tells whether or not the genes' mixture element
#' should be estimated in the MCMC algorithm run. Default value is TRUE.
#' 
#' @return mcmc Returns an intialized MCMC object. 
#' 
#' @description \code{initializeMCMCObject} initializes a MCMC object to 
#' perform a model fitting for a parameter and model object.
#' 
#' @details \code{initializeMCMCObject} sets up the MCMC object 
#' (monte carlo markov chain) and returns the object so a model fitting can be done.
#' It is important to note that est.expression and est.hyper will affect one another
#' negatively if their values differ.
#' 
#' @examples 
#' 
#' ## initializing an object of type mcmc
#' 
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' 
#' ## estimate all parameter types
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width=adaptiveWidth, 
#'                              est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#'                              
#' ## do not estimate expression values, initial conditions will remain constant
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width=adaptiveWidth, 
#'                              est.expression=FALSE, est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#'                              
#' ## do not estimate hyper parameters, initial conditions will remain constant
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width=adaptiveWidth, 
#'                              est.expression=TRUE, est.csp=TRUE, est.hyper=FALSE, est.mix = TRUE) 
#' 
initializeMCMCObject <- function(samples, thinning=1, adaptive.width=100, 
                                 est.expression=TRUE, est.csp=TRUE, 
                                 est.hyper=TRUE, est.mix=TRUE){
  
  # error check given values.
  if (!is.numeric(samples) || samples < 1 || !all(samples == as.integer(samples))) {
    stop("samples must be a positive integer\n")
  }
  if (!is.numeric(thinning) || thinning < 1 || !all(thinning == as.integer(thinning))) {
    stop("thinning must be a positive integer\n")
  }
  if (!is.numeric(adaptive.width) || adaptive.width < 1 || 
      !all(adaptive.width == as.integer(adaptive.width))) {
    stop("adaptive.width must be a positive integer\n")
  }
  if (!identical(est.expression, TRUE) && !identical(est.expression, FALSE)) {
    stop("est.expression must be a boolean value\n")
  }
  if (!identical(est.csp, TRUE) && !identical(est.csp, FALSE)) {
    stop("est.csp must be a boolean value\n")
  }
  if (!identical(est.hyper, TRUE) && !identical(est.hyper, FALSE)) {
    stop("est.hyper must be a boolean value\n")
  }
  if (!identical(est.mix, TRUE) && !identical(est.mix, FALSE)) {
    stop("est.mix must be a boolean value\n")
  }

  mcmc <- new(MCMCAlgorithm, samples, thinning, adaptive.width, est.expression, 
              est.csp, est.hyper)
  mcmc$setEstimateMixtureAssignment(est.mix)
  return(mcmc)
}


#' Run MCMC 
#' 
#' @param mcmc MCMC object that will run the model fitting algorithm.
#' 
#' @param genome Genome that the model fitting will run on. Should be 
#' the same genome associated with the parameter and model objects.
#' 
#' @param model Model to run the fitting on. Should be associated with
#' the given genome.
#' 
#' @param ncores Number of cores to perform the model fitting with. Default
#' value is 1.
#' 
#' @param divergence.iteration Number of steps that the initial conditions
#' can diverge from the original conditions given. Default value is 0.
#' 
#' @return This function has no return value.
#' 
#' @description \code{runMCMC} will run a monte carlo markov chain algorithm
#' for the given mcmc, genome, and model objects to perform a model fitting.
#' 
#' @details \code{runMCMC} will run for the number of samples times the number
#' thinning given when the mcmc object is initialized. Updates are provided every 100
#' steps, and the state of the chain is saved every thinning steps.
#' 
#' @examples 
#' 
#' #fitting a model to a genome using the runMCMC function
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' model <- initializeModelObject(parameter = parameter, model = "ROC")
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
#'                              adaptive.width=adaptiveWidth, est.expression=TRUE, 
#'                              est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#' divergence.iteration <- 10
#' \dontrun{
#' runMCMC(mcmc = mcmc, genome = genome, model = model, 
#'         ncores = 4, divergence.iteration = divergence.iteration)
#' }
#' 
runMCMC <- function(mcmc, genome, model, ncores = 1, divergence.iteration = 0){
  if(class(mcmc) != "Rcpp_MCMCAlgorithm") stop("mcmc is not of class Rcpp_MCMCAlgorithm")
  
  if (ncores < 1 || !all(ncores == as.integer(ncores))) {
    stop("ncores must be a positive integer\n")
  }
  if(class(model) == "Rcpp_PANSEModel")
  {
    mcmc$run_PANSE(genome, model, ncores, divergence.iteration)
  } else
  {
   mcmc$run(genome, model, ncores, divergence.iteration)
  }
}


#' Set Restart Settings 
#' 
#' @param mcmc MCMC object that will run the model fitting algorithm.
#' 
#' @param filename Filename for the restart files to be written.
#' 
#' @param samples Number of samples that should occur before a file is written.
#' 
#' @param write.multiple Boolean that determines if multiple restart files
#' are written. Default value is TRUE.
#' 
#' @return This function has no return value.
#' 
#' @description \code{setRestartSettings} sets the needed information (what the file 
#' is called, how often the file should be written) to write
#' information to restart the MCMC algorithm from a given point.
#' 
#' @details \code{setRestartSettings} writes a restart file every set amount of samples
#' that occur. Also, if write.multiple is true, instead of overwriting the previous restart
#' file, the sample number is prepended onto the file name and multiple rerstart files
#' are generated for a run.
#' 
#' @examples 
#' 
#' ## set restart settings for checkpointing
#' 
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' 
#' ## estimate all parameter types
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
#'                              adaptive.width=adaptiveWidth, est.expression=TRUE, 
#'                              est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#'                              
#' # prompts the mcmc to write a restart file every 100 samples during the run.
#' setRestartSettings(mcmc = mcmc, filename = "test_restart", samples = 100)
#' 
#' # prompts the mcmc to write a restart file every 100 samples during the run, 
#' # but will overwrite it each time.
#' setRestartSettings(mcmc = mcmc, filename = "test_restart", samples = 100, 
#'                    write.multiple = FALSE)
#'            
setRestartSettings <- function(mcmc, filename, samples, write.multiple=TRUE){
  if(class(mcmc) != "Rcpp_MCMCAlgorithm") stop("mcmc is not of class Rcpp_MCMCAlgorithm")
  mcmc$setRestartFileSettings(filename, samples, write.multiple)
}


#' Convergence Test
#' 
#' @param object an object of either class Trace or MCMC
#' 
#' @param samples number of samples at the end of the trace used to determine convergence (< length of trace). 
#' Will use as starting point of convergence test. If the MCMC trace
#' is of length x, then starting point for convergence test will be x - samples.
#' 
#' @param frac1 fraction to use from beginning of samples
#' 
#' @param frac2 fraction to use from end of samples
#' 
#' @param thin the thinning interval between consecutive observations, which is used in creating a coda::mcmc object (according to the Coda documentation, users should specify if a MCMC chain has already been thinned using a the thin parameter). This does not further thin the data.
#' 
#' @param plot (logical) plot result instead of returning an object
#' 
#' @param what (for Trace Object only) which parameter to calculate convergence.test -- current options are Selection, Mutation, MixtureProbability, Sphi, Mphi, ExpectedPhi, and AcceptanceCSP
#' 
#' @param mixture (for Trace Object only) mixture for which to calculate convergence.test
#' 
#' @details  Be aware that convergence.test for Trace objects works primarily for Trace objects from the ROC parameter class. Future updates will adapt this function to work for parameters from other models and expression traces
#' 
#' @return Geweke score object evaluating whether means of two fractions (frac1 and frac2) differ.  Convergence occurs when they don't differ significantly, i.e. pnorm(abs(convergence.test(mcmcObj)$a, ,lower.tail=FALSE)*2 > 0.05
#' 
#' @examples 
#' 
#' ## check for convergence after a run:
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
#'                              adaptive.width=adaptiveWidth, est.expression=TRUE, 
#'                              est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#' divergence.iteration <- 10
#' \dontrun{
#' runMCMC(mcmc = mcmc, genome = genome, model = model, 
#'         ncores = 4, divergence.iteration = divergence.iteration)
#' # check if posterior trace has converged
#' convergence.test(object = mcmc, samples = 500, plot = TRUE)
#' 
#' trace <- getTrace(parameter)
#' # check if Mutation trace has converged
#' convergence.test(object = trace, samples = 500, plot = TRUE, what = "Mutation")
#' # check if Sphi trace has converged
#' convergence.test(object = trace, samples = 500, plot = TRUE, what = "Sphi")
#' # check if ExpectedPhi trace has converged
#' convergence.test(object = trace, samples = 500, plot = TRUE, what = "ExpectedPhi")
#' }
convergence.test <- function(object, samples = 10, frac1 = 0.1, frac2 = 0.5, 
                    thin = 1, plot = FALSE, what = "Mutation", mixture = 1){
  UseMethod("convergence.test", object)
}

convergence.test.Rcpp_MCMCAlgorithm <- function(object, samples = 10, frac1 = 0.1, 
                                       frac2 = 0.5, thin = 1, plot = FALSE, what = "Mutation", mixture = 1){
  # TODO: extend to work with multiple chains once we have that capability.
  trace <- object$getLogPosteriorTrace()
  loglik.trace <- object$getLogPosteriorTrace()
  trace.length <- length(loglik.trace)
  start <- max(1, trace.length - samples)
  
  # the start and end parameter do NOT work, using subsetting to achieve goal
  mcmcobj <- coda::mcmc(data=loglik.trace[start:trace.length], thin = thin)
  diag <- coda::geweke.diag(mcmcobj, frac1=frac1, frac2=frac2)
  if(plot){ 
    coda::geweke.plot(mcmcobj, frac1=frac1, frac2=frac2)
  }else{
    return(diag)
  }
}


#' Write MCMC Object 
#' 
#' @param mcmc MCMC object that has run the model fitting algorithm.
#' 
#' @param file A filename where the data will be stored.
#' 
#' @return This function has no return value.
#' 
#' @description \code{writeMCMCObject} stores the MCMC information from the 
#' model fitting run in a file.
#' 
#' @examples
#'
#' ## saving the MCMC object after model fitting
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
#'                              adaptive.width=adaptiveWidth, est.expression=TRUE, 
#'                              est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#' divergence.iteration <- 10
#' \dontrun{
#' runMCMC(mcmc = mcmc, genome = genome, model = model, 
#'         ncores = 4, divergence.iteration = divergence.iteration)
#' writeMCMCObject(mcmc = mcmc, file = file.path(tempdir(), "file.Rda"))
#' 
#' }
writeMCMCObject <- function(mcmc, file){
  logPostTrace <- mcmc$getLogPosteriorTrace()
  logLikeTrace <- mcmc$getLogLikelihoodTrace()
  samples <- mcmc$getSamples()
  thinning <- mcmc$getThinning()
  adaptiveWidth <- mcmc$getAdaptiveWidth()
  save(list = c("logPostTrace","logLikeTrace", "samples", "thinning", "adaptiveWidth"), file=file)
}


#' Load MCMC Object
#' 
#' @param files The filenames where the data will be stored.
#' 
#' @return This function has no return value.
#' 
#' @description \code{loadMCMCObject} creates a new MCMC object and fills it with
#' the information in the file given.
#' 
#' @details This MCMC object is not intended to be used to do another model fitting, only
#' to graph the stored results.
#' 
#' @examples
#' 
#' ## loading mcmc objects from the filesystem
#' \dontrun{
#' # load one mcmc object
#' mcmc <- loadMCMCObject(files = "mcmc.Rda")
#' 
#' # load and combine multiple mcmc objects. Useful when using checkpointing
#' mcmc <- loadMCMCObject(files = c("mcmc1.Rda", "mcmc2.Rda"))
#' }
loadMCMCObject <- function(files){
  mcmc <- new(MCMCAlgorithm)
  samples <- 0
  logPostTrace <- numeric(0)
  logLikeTrace <- numeric(0)
  for (i in 1:length(files)){
    tempEnv <- new.env();
    load(file = files[i], envir = tempEnv)
    samples <- samples + tempEnv$samples
    max <- tempEnv$samples + 1
    curLogPostTrace <- tempEnv$logPostTrace
    curLoglikelihoodTrace <- tempEnv$logLikeTrace
    logPostTrace <- c(logPostTrace, curLogPostTrace[2:max])
    logLikeTrace <- c(logLikeTrace, curLoglikelihoodTrace[2:max])
   }
    mcmc$setSamples(samples)
    mcmc$setThinning(tempEnv$thinning) #not needed?
    mcmc$setAdaptiveWidth(tempEnv$adaptiveWidth) #not needed?
    mcmc$setLogPosteriorTrace(logPostTrace)
    mcmc$setLogLikelihoodTrace(logLikeTrace)

  return(mcmc)
}


