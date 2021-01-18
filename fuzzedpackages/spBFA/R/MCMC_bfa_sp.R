#' Spatial factor analysis using a Bayesian hierarchical model.
#'
#' \code{bfa_sp} is a Markov chain Monte Carlo (MCMC) sampler for a Bayesian spatial factor analysis model. The spatial component is 
#' introduced using a Probit stick-breaking process prior on the factor loadings. The model is implemented using a Bayesian hierarchical framework.
#'
#' @param formula A \code{formula} object, corresponding to the spatial factor analysis model. The response must be on the left of a \code{~} operator, and the terms on the right 
#'                 must indicate the covariates to be included in the fixed effects. If no covariates are desired a zero should be used, \code{~ 0}. 
#'                 
#' @param data A required \code{data.frame} containing the variables in the model. The data frame must contain \code{M x O x Nu} rows.
#'  Here, \code{M} represents the number of spatial locations, \code{O} the number of different observation types
#'  and \code{Nu} the number of temporal visits. The observations must be first be
#'  ordered spatially, second by observation type and then temporally. This means that the first \code{M x O} observations come from the first time point and
#'  the first \code{M} observations come the first spatial observation type.
#'
#' @param family Character string indicating the distribution of the observed data. Options
#'  include: \code{"normal"}, \code{"probit"}, \code{"tobit"}, and \code{"binomial"}. \code{family} must have either \code{O} or
#'  \code{1} dimension(s) (the one populates the rest). Any combination of likelihoods can be used.
#'  
#' @param dist A \code{M x M} dimensional distance matrix. For a \code{discrete} spatial process the matrix contains binary adjacencies that dictate the
#'  spatial neighborhood structure and for \code{continuous} spatial processes the matrix should be a continuous distance matrix (e.g., Euclidean).
#'  
#' @param time A \code{Nu} dimensional vector containing the observed time points
#'  in increasing order.
#'
#' @param K A scalar that indicates the dimension (i.e., quantity) of latent factors.
#'  
#' @param L The number of latent clusters. If finite, a scalar indicating the number of clusters for each column of the factor loadings matrix. By default \code{L} is set at \code{Inf}
#'  so that the Probit stick-breaking process becomes an infinite mixture model.
#'  
#' @param trials A variable in \code{data} that contains the number of trials for each of the binomial observations. If there is no count data, \code{trials} should be left missing.
#'  
#' @param starting Either \code{NULL} or a \code{list} containing starting values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the starting values may be specified.
#'
#'  When \code{NULL} is chosen then default starting values are automatically generated.
#'  Otherwise a \code{list} must be provided with names \code{Beta}, \code{Delta}, \code{Sigma2}, \code{Kappa}, \code{Rho}, \code{Upsilon} or
#'  \code{Psi} containing appropriate objects. \code{Beta} (or \code{Delta}) must either be a \code{P} (or \code{K}) dimensional
#'  vector or a scalar (the scalar populates the entire vector). \code{Sigma2} must be either a \code{M x (O - C)} matrix or a scalar.
#'  \code{Kappa} must be a \code{O x O} dimensional matrix, \code{Rho} a scalar, \code{Upsilon} a \code{K x K} matrix, and \code{Psi} a scalar.
#'
#' @param hypers Either \code{NULL} or a \code{list} containing hyperparameter values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the hyperparameter values may be specified.
#'
#'  When \code{NULL} is chosen then default hyperparameter values are automatically
#'  generated. These default hyperparameters are described in detail in (Berchuck et al.).
#'  Otherwise a \code{list} must be provided with names \code{Beta}, \code{Delta}, \code{Sigma2}, \code{Kappa}, \code{Rho}, \code{Upsilon} or
#'  \code{Psi} containing further hyperparameter information. These objects are themselves
#'  \code{lists} and may be constructed as follows.
#'
#'  \code{Beta} is a \code{list} with two objects, \code{MuBeta} and \code{SigmaBeta}. These values represent the prior mean and variance 
#'  parameters for the multivariate normal prior.
#'  
#'  \code{Delta} is a \code{list} with two objects, \code{A1} and \code{A2}. These values represent the prior shape 
#'  parameters for the multiplicative Gamma shrinkage prior.
#'  
#'  \code{Sigma2} is a \code{list} with two objects, \code{A} and \code{B}. These values represent the shape and scale for the variance parameters.
#'  
#'  \code{Kappa} is a \code{list} with two objects,
#'  \code{SmallUpsilon} and \code{BigTheta}. \code{SmallUpsilon} represents the degrees of freedom parameter for the
#'  inverse-Wishart hyperprior and must be a real number scalar, while \code{BigTheta} represents
#'  the scale matrix and must be a \code{O x O} dimensional positive definite matrix.
#'  
#'  \code{Rho} is a \code{list} with two objects, \code{ARho} and \code{BRho}. \code{ARho}
#'  represents the lower bound for the uniform hyperprior, while \code{BRho} represents
#'  the upper bound. The bounds must be specified carefully. This is only specified for continuous spatial processes.
#'  
#'  \code{Upsilon} is a \code{list} with two objects,
#'  \code{Zeta} and \code{Omega}. \code{Zeta} represents the degrees of freedom parameter for the
#'  inverse-Wishart hyperprior and must be a real number scalar, while \code{Omega} represents
#'  the scale matrix and must be a \code{K x K} dimensional positive definite matrix.
#'
#'  \code{Psi} is a \code{list} with two objects, dependent on if the temporal kernel is \code{exponential} or \code{ar1}.
#'  For \code{exponential}, the two objects are \code{APsi} and \code{BPsi}. \code{APsi}
#'  represents the lower bound for the uniform hyperprior, while \code{BPsi} represents
#'  the upper bound. The bounds must be specified carefully. For \code{ar1}, the two objects are \code{Beta} and \code{Gamma}, which are the 
#'  two shape parameters of a Beta distribution shifted to have domain in (-1, 1). 
#'  
#' @param tuning Either \code{NULL} or a \code{list} containing tuning values
#'  to be specified for the MCMC Metropolis steps. If \code{NULL} is not chosen then all
#'  of the tuning values must be specified.
#'
#'  When \code{NULL} is chosen then default tuning values are automatically generated to
#'  \code{1}. Otherwise a \code{list} must be provided with names \code{Psi}, 
#'  or \code{Rho}. Each of these entries must be scalars containing tuning variances for their corresponding Metropolis updates.
#'  \code{Rho} is only specified for continuous spatial processes.
#'
#' @param mcmc Either \code{NULL} or a \code{list} containing input values to be used
#'  for implementing the MCMC sampler. If \code{NULL} is not chosen then all
#'  of the MCMC input values must be specified.
#'
#'  \code{NBurn}: The number of sampler scans included in the burn-in phase. (default =
#'  \code{10,000})
#'
#'  \code{NSims}: The number of post-burn-in scans for which to perform the
#'   sampler. (default = \code{10,000})
#'
#'  \code{NThin}: Value such that during the post-burn-in phase, only every
#'  \code{NThin}-th scan is recorded for use in posterior inference (For return values
#'  we define, NKeep = NSims / NThin (default = \code{1}).
#'
#'  \code{NPilot}: The number of times during the burn-in phase that pilot adaptation
#'  is performed (default = \code{20})
#'
#' @param temporal.structure Character string indicating the temporal kernel. Options include:
#'  \code{"exponential"} and \code{"ar1"}.
#'
#' @param spatial.structure Character string indicating the type of spatial process. Options include:
#'  \code{"continuous"} (i.e., Gaussian process with exponential kernel) and \code{"discrete"} (i.e., proper CAR).
#'
#' @param seed An integer value used to set the seed for the random number generator
#'  (default = 54).
#'  
#' @param gamma.shrinkage A logical indicating whether a gamma shrinkage process prior is used for the variances of the factor loadings columns. If FALSE,
#'  the hyperparameters (A1 and A2) indicate the shape and rate for a gamma prior on the precisions. Default is TRUE.
#'
#' @param include.space A logical indicating whether a spatial process should be included. Default is TRUE, however if FALSE the spatial correlation matrix 
#'  is fixed as an identity matrix. This specification overrides the \code{spatial.structure} input.
#'  
#' @param clustering A logical indicating whether the Bayesian non-parametric process should be used, default is TRUE. If FALSE is specified
#'  each column is instead modeled with an independent spatial process.
#'  
#' @details Details of the underlying statistical model proposed by
#'  Berchuck et al. 2019. are forthcoming.
#'
#' @return \code{bfa_sp} returns a list containing the following objects
#'
#'   \describe{
#'
#'   \item{\code{lambda}}{\code{NKeep x (M x O x K)} \code{matrix} of posterior samples for factor loadings matrix \code{lambda}.
#'   The labels for each column are Lambda_O_M_K.}
#'
#'   \item{\code{eta}}{\code{NKeep x (Nu x K)} \code{matrix} of posterior samples for the latent factors \code{eta}.
#'   The labels for each column are Eta_Nu_K.}
#'
#'   \item{\code{beta}}{\code{NKeep x P} \code{matrix} of posterior samples for \code{beta}.}
#'
#'   \item{\code{sigma2}}{\code{NKeep x (M * (O - C))} \code{matrix} of posterior samples for the variances \code{sigma2}.
#'   The labels for each column are Sigma2_O_M.}
#'
#'   \item{\code{kappa}}{\code{NKeep x ((O * (O + 1)) / 2)} \code{matrix} of posterior samples for \code{kappa}. The
#'   columns have names that describe the samples within them. The row is listed first, e.g.,
#'   \code{Kappa3_2} refers to the entry in row \code{3}, column \code{2}.}
#'
#'   \item{\code{delta}}{\code{NKeep x K} \code{matrix} of posterior samples for \code{delta}.}
#'
#'   \item{\code{tau}}{\code{NKeep x K} \code{matrix} of posterior samples for \code{tau}.}
#'
#'   \item{\code{upsilon}}{\code{NKeep x ((K * (K + 1)) / 2)} \code{matrix} of posterior samples for \code{Upsilon}. The
#'   columns have names that describe the samples within them. The row is listed first, e.g.,
#'   \code{Upsilon3_2} refers to the entry in row \code{3}, column \code{2}.}
#'
#'   \item{\code{psi}}{\code{NKeep x 1} \code{matrix} of posterior samples for \code{psi}.}
#'
#'   \item{\code{xi}}{\code{NKeep x (M x O x K)} \code{matrix} of posterior samples for factor loadings cluster labels \code{xi}.
#'   The labels for each column are Xi_O_M_K.}
#'   
#'   \item{\code{rho}}{\code{NKeep x 1} \code{matrix} of posterior samples for \code{rho}.}
#'
#'   \item{\code{metropolis}}{\code{2 (or 1) x 3} \code{matrix} of metropolis
#'   acceptance rates, updated tuners, and original tuners that result from the pilot
#'   adaptation.}
#'
#'   \item{\code{runtime}}{A \code{character} string giving the runtime of the MCMC sampler.}
#'
#'   \item{\code{datobj}}{A \code{list} of data objects that are used in future \code{bfa_sp} functions
#'   and should be ignored by the user.}
#'
#'   \item{\code{dataug}}{A \code{list} of data augmentation objects that are used in future
#'   \code{bfa_sp} functions and should be ignored by the user.}
#'
#'   }
#'
#' @examples
#' \donttest{
#' ###Load womblR for example visual field data
#' library(womblR)
#' 
#' ###Format data for MCMC sampler
#' blind_spot <- c(26, 35) # define blind spot
#' VFSeries <- VFSeries[order(VFSeries$Location), ] # sort by location
#' VFSeries <- VFSeries[order(VFSeries$Visit), ] # sort by visit
#' VFSeries <- VFSeries[!VFSeries$Location %in% blind_spot, ] # remove blind spot locations
#' dat <- data.frame(Y = VFSeries$DLS / 10) # create data frame with scaled data
#' Time <- unique(VFSeries$Time) / 365 # years since baseline visit
#' W <- HFAII_Queen[-blind_spot, -blind_spot] # visual field adjacency matrix (data object from womblR)
#' M <- dim(W)[1] # number of locations
#' 
#' ###Prior bounds for psi
#' TimeDist <- as.matrix(dist(Time))
#' BPsi <- log(0.025) / -min(TimeDist[TimeDist > 0])
#' APsi <- log(0.975) / -max(TimeDist)
#' 
#' ###MCMC options
#' K <- 10 # number of latent factors
#' O <- 1 # number of spatial observation types
#' Hypers <- list(Sigma2 = list(A = 0.001, B = 0.001),
#'                Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
#'                Delta = list(A1 = 1, A2 = 20),
#'                Psi = list(APsi = APsi, BPsi = BPsi),
#'                Upsilon = list(Zeta = K + 1, Omega = diag(K)))
#' Starting <- list(Sigma2 = 1,
#'                  Kappa = diag(O),
#'                  Delta = 2 * (1:K),
#'                  Psi = (APsi + BPsi) / 2,
#'                  Upsilon = diag(K))
#' Tuning <- list(Psi = 1)
#' MCMC <- list(NBurn = 1000, NSims = 1000, NThin = 2, NPilot = 5)
#' 
#' ###Fit MCMC Sampler
#' reg.bfa_sp <- bfa_sp(Y ~ 0, data = dat, dist = W, time = Time,  K = 10, 
#'                      starting = Starting, hypers = Hypers, tuning = Tuning, mcmc = MCMC,
#'                      L = Inf,
#'                      family = "tobit",
#'                      trials = NULL,
#'                      temporal.structure = "exponential",
#'                      spatial.structure = "discrete",
#'                      seed = 54, 
#'                      gamma.shrinkage = TRUE,
#'                      include.space = TRUE,
#'                      clustering = TRUE)
#' 
#' ###Note that this code produces the pre-computed data object "reg.bfa_sp"
#' ###More details can be found in the vignette.
#' 
#' }
#'
# @author Samuel I. Berchuck
#' @references Reference for Berchuck et al. 2019 is forthcoming.
#' @export
bfa_sp <- function(formula, data, dist, time, K, L = Inf, trials = NULL,
                   family = "normal", temporal.structure = "exponential", spatial.structure = "discrete",
                   starting = NULL, hypers = NULL, tuning = NULL, mcmc = NULL, seed = 54,
                   gamma.shrinkage = TRUE, include.space = TRUE, clustering = TRUE) {
  
  ###Function Inputs
  # formula = sens ~ age
  # data = dat
  # dist = W
  # time = Time
  # trials = NULL
  # starting = Starting
  # hypers = Hypers
  # tuning = Tuning
  # mcmc = MCMC
  # family = "tobit"
  # temporal.structure = "exponential"
  # spatial.structure = "discrete"
  # seed = 54
  # K = K
  # L = Inf
  # gamma.shrinkage = TRUE
  # include.space = TRUE
  # clustering = TRUE

  ###Check for missing objects
  if (missing(formula)) stop("formula: missing")
  if (missing(data)) stop("data: missing")
  if (missing(dist)) stop("dist: missing")
  if (missing(time)) stop("time: missing")
  if (missing(K)) stop("K: missing")

  ###Check model inputs
  CheckInputs(formula, data, dist, time, K, L, trials, starting, hypers, tuning, mcmc, family, temporal.structure, spatial.structure, gamma.shrinkage, include.space, clustering)

  ####Set seed for reproducibility
  set.seed(seed)

  ###Check to see if the job is interactive
  Interactive <- interactive()

  ###Create objects for use in sampler
  DatObj <- CreateDatObj(formula, data, dist, time, trials, K, L, family, temporal.structure, spatial.structure, gamma.shrinkage, include.space, clustering)
  HyPara <- CreateHyPara(hypers, DatObj) 
  MetrObj <- CreateMetrObj(tuning, DatObj)
  Para <- CreatePara(starting, DatObj, HyPara)
  DatAug <- CreateDatAug(DatObj)
  McmcObj <- CreateMcmc(mcmc, DatObj)
  RawSamples <- CreateStorage(DatObj, McmcObj)

  ###Time MCMC sampler
  BeginTime <- Sys.time()

  ###Run MCMC sampler in Rcpp
  RegObj <- bfa_sp_Rcpp(DatObj, HyPara, MetrObj, Para, DatAug, McmcObj, RawSamples, Interactive)

  ###Set regression objects
  RawSamples <- RegObj$rawsamples
  MetropRcpp <- RegObj$metropolis

  ###End time
  FinishTime <- Sys.time()
  RunTime <- FinishTime - BeginTime

  ###Collect output to be returned
  DatObjOut <- OutputDatObj(DatObj)
  DatAugOut <- OutputDatAug(DatAug)
  Metropolis <- SummarizeMetropolis(DatObj, MetrObj, MetropRcpp, McmcObj)
  Samples <- FormatSamples(DatObj, RawSamples)

  ###Return spBFA object
  spBFA <- list(lambda = Samples$Lambda,
                eta = Samples$Eta,
                beta = Samples$Beta,
                sigma2 = Samples$Sigma2,
                kappa = Samples$Kappa,
                delta = Samples$Delta,
                tau = Samples$Tau,
                upsilon = Samples$Upsilon,
                psi = Samples$Psi,
                xi = Samples$Xi,
                rho = Samples$Rho,
                metropolis = Metropolis,
                datobj = DatObjOut,
                dataug = DatAugOut,
                runtime = paste0("Model runtime: ",round(RunTime, 2), " ", attr(RunTime, "units")))
  spBFA <- structure(spBFA, class = "spBFA")
  return(spBFA)

###End sampler
}
