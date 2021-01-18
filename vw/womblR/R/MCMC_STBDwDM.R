#' MCMC sampler for spatiotemporal boundary detection with dissimilarity metric.
#'
#' \code{STBDwDM} is a Markov chain Monte Carlo (MCMC) sampler for a spatiotemporal
#'  boundary detection model using the Bayesian hierarchical framework.
#'
#' @param Y An \code{N} dimensional vector containing the observed outcome data.
#'  Here, \code{N = M * Nu}, where \code{M} represents the number of spatial locations
#'  and \code{Nu} the number of temporal visits. The observations in \code{Y} must be first
#'  ordered spatially and then temporally, meaning the first \code{M} observations
#'  in \code{Y} should come from the initial time point.
#'
#' @param DM An \code{M} dimensional vector containing a dissimilarity metric
#'  for each spatial location. The order of the spatial locations must match the order from
#'  \code{Y}.
#'
#' @param W An \code{M x M} dimensional binary adjacency matrix for dictating the
#'  spatial neigborhood structure.
#'
#' @param Time A \code{Nu} dimensional vector containing the observed time points for each
#'  vector of outcomes in increasing order.
#'
#' @param Starting Either \code{NULL} or a \code{list} containing starting values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the starting values may be specified.
#'
#'  When \code{NULL} is chosen then default starting values are automatically generated.
#'  Otherwise a \code{list} must be provided with names \code{Delta}, \code{T} or
#'  \code{Phi} containing appropriate objects. \code{Delta} must be a \code{3} dimensional
#'  vector, \code{T} a \code{3 x 3} dimensional matrix and \code{Phi} a scalar.
#'
#' @param Hypers Either \code{NULL} or a \code{list} containing hyperparameter values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the hyperparameter values may be specified.
#'
#'  When \code{NULL} is chosen then default hyperparameter values are automatically
#'  generated. These default hyperparameters are described in detail in (Berchuck et al.).
#'  Otherwise a \code{list} must be provided with names \code{Delta}, \code{T} or
#'  \code{Phi} containing further hyperparameter information. These objects are themselves
#'  \code{lists} and may be constructed as follows.
#'
#'  \code{Delta} is a \code{list} with two objects, \code{MuDelta} and \code{OmegaDelta}.
#'  \code{MuDelta} represents the mean component of the multivariate normal hyperprior and
#'  must be a \code{3} dimensional vector, while \code{OmegaDelta} represents the covariance
#'  and must be a \code{3 x 3} dimensional matrix.
#'
#'  \code{T} is a \code{list} with two objects, \code{Xi} and \code{Psi}. \code{Xi}
#'  represents the degrees of freedom parameter for the inverse-Wishart hyperprior and
#'  must be a real number scalar, while \code{Psi} represents the scale matrix
#'  and must be a \code{3 x 3} dimensional positive definite matrix.
#'
#'  \code{Phi} is a \code{list} with two objects, \code{APhi} and \code{BPhi}. \code{APhi}
#'  represents the lower bound for the uniform hyperprior, while \code{BPhi} represents
#'  the upper bound. The bounds must be specified carefully. For example, if the exponential
#'  temporal correlation structure is chosen both bounds must be restricted to be non-negative.
#'
#' @param Tuning Either \code{NULL} or a \code{list} containing tuning values
#'  to be specified for the MCMC Metropolis steps. If \code{NULL} is not chosen then all
#'  of the tuning values must be specified.
#'
#'  When \code{NULL} is chosen then default tuning values are automatically generated to
#'  \code{1}. Otherwise a \code{list} must be provided with names \code{Theta2},
#'  \code{Theta3} and \code{Phi}. \code{Theta2} and \code{Theta3} must be
#'  \code{Nu} dimensional vectors and \code{Phi} a scalar. Each containing tuning variances
#'  for their corresponding Metropolis updates.
#'
#' @param MCMC Either \code{NULL} or a \code{list} containing input values to be used
#'  for implementing the MCMC sampler. If \code{NULL} is not chosen then all
#'  of the MCMC input values must be specified.
#'
#'  \code{NBurn}: The number of sampler scans included in the burn-in phase. (default =
#'  \code{10,000})
#'
#'  \code{NSims}: The number of post-burn-in scans for which to perform the
#'   sampler. (default = \code{100,000})
#'
#'  \code{NThin}: Value such that during the post-burn-in phase, only every
#'  \code{NThin}-th scan is recorded for use in posterior inference (For return values
#'  we define, NKeep = NSims / NThin (default = \code{10}).
#'
#'  \code{NPilot}: The number of times during the burn-in phase that pilot adaptation
#'  is performed (default = \code{20})
#'
#' @param Family Character string indicating the distribution of the observed data. Options
#'  include: \code{"normal"}, \code{"probit"}, \code{"tobit"}.
#'
#' @param TemporalStructure Character string indicating the temporal structure of the
#'  time observations. Options include: \code{"exponential"} and \code{"ar1"}.
#'
#' @param Distance Character string indicating the distance metric for computing the
#'  dissimilarity metric. Options include: \code{"euclidean"} and \code{"circumference"}.
#'
#' @param Weights Character string indicating the type of weight used. Options include:
#'  \code{"continuous"} and \code{"binary"}.
#'
#' @param Rho A scalar in \code{(0,1)} that dictates the magnitude of local spatial sharing.
#'  By default it is fixed at \code{0.99} as suggested by Lee and Mitchell (2012).
#'
#' @param ScaleY A positive scalar used for scaling the observed data, \code{Y}. This is
#'  used to aid numerically for MCMC convergence, as scaling large observations often
#'  stabilizes chains. By default it is fixed at \code{10}.
#'
#' @param ScaleDM A positive scalar used for scaling the dissimilarity metric distances,
#'  \code{DM}. This is used to aid numerically for MCMC convergence. as scaling spatial
#'  distances is often used for improved MCMC convergence. By default it is fixed at \code{100}.
#'
#' @param Seed An integer value used to set the seed for the random number generator
#'  (default = 54).
#'
#' @details Details of the underlying statistical model can be found in the article by
#' Berchuck et al. (2018), "Diagnosing Glaucoma Progression with Visual Field Data Using
#' a Spatiotemporal Boundary Detection Method", <arXiv:1805.11636>.
#'
#' @return \code{STBDwDM} returns a list containing the following objects
#'
#'   \describe{
#'
#'   \item{\code{mu}}{\code{NKeep x Nu} \code{matrix} of posterior samples for \code{mu}. The
#'   t-th column contains posterior samples from the the t-th time point.}
#'
#'   \item{\code{tau2}}{\code{NKeep x Nu} \code{matrix} of posterior samples for \code{tau2}.
#'   The t-th column contains posterior samples from the the t-th time point.}
#'
#'   \item{\code{alpha}}{\code{NKeep x Nu} \code{matrix} of posterior samples for \code{alpha}.
#'   The t-th column contains posterior samples from the the t-th time point.}
#'
#'   \item{\code{delta}}{\code{NKeep x 3} \code{matrix} of posterior samples for \code{delta}.
#'   The columns have names that describe the samples within them.}
#'
#'   \item{\code{T}}{\code{NKeep x 6} \code{matrix} of posterior samples for \code{T}. The
#'   columns have names that describe the samples within them. The row is listed first, e.g.,
#'   \code{t32} refers to the entry in row \code{3}, column \code{2}.}
#'
#'   \item{\code{phi}}{\code{NKeep x 1} \code{matrix} of posterior samples for \code{phi}.}
#'
#'   \item{\code{metropolis}}{\code{(2 * Nu + 1) x 2} \code{matrix} of metropolis
#'   acceptance rates and tuners that result from the pilot adaptation. The first \code{Nu}
#'   correspond to the \code{Theta2} (i.e. \code{tau2}) parameters, the next \code{Nu} correspond to
#'   the \code{Theta3} (i.e. \code{alpha}) parameters and the last row give the \code{phi} values.}
#'
#'   \item{\code{runtime}}{A \code{character} string giving the runtime of the MCMC sampler.}
#'
#'   \item{\code{datobj}}{A \code{list} of data objects that are used in future \code{STBDwDM} functions
#'   and should be ignored by the user.}
#'
#'   \item{\code{dataug}}{A \code{list} of data augmentation objects that are used in future
#'   \code{STBDwDM} functions and should be ignored by the user.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @references Berchuck et al. (2018), "Diagnosing Glaucoma Progression with Visual Field Data Using
#' a Spatiotemporal Boundary Detection Method", <arXiv:1805.11636>.
#' @export
STBDwDM <- function(Y, DM, W, Time, Starting = NULL, Hypers = NULL, Tuning = NULL,
			 		          MCMC = NULL, Family = "tobit", TemporalStructure = "exponential",
			 		          Distance = "circumference", Weights = "continuous", Rho = 0.99,
			 		          ScaleY = 10, ScaleDM = 100, Seed = 54) {


  ###Function inputs
  # Y = Y
  # DM = DM
  # W = W
  # Time = Time
  # Starting = Starting
  # Hypers = Hypers
  # Tuning = Tuning
  # MCMC = MCMC
  # Family = "tobit"
  # TemporalStructure = "exponential"
  # Distance = "circumference"
  # Weights = "continuous"
  # Rho = 0.99
  # ScaleY = 10
  # ScaleDM = 100
  # Seed = 54

  ###Check for missing objects
  if (missing(Y)) stop("Y: missing")
  if (missing(DM)) stop("DM: missing")
  if (missing(W)) stop("W: missing")
  if (missing(Time)) stop("Time: missing")

  ###Check model inputs
  CheckInputs(Y, DM, W, Time, Starting, Hypers, Tuning, MCMC, Family, TemporalStructure, Distance, Weights, Rho, ScaleY, ScaleDM)

  ####Set seed for reproducibility
  set.seed(Seed)

  ###Check to see if the job is interactive
  Interactive <- interactive()

  ###Create objects for use in sampler
  DatObj <- CreateDatObj(Y, DM, W, Time, Rho, ScaleY, ScaleDM, TemporalStructure, Family, Distance, Weights)
  HyPara <- CreateHyPara(Hypers, DatObj)
  MetrObj <- CreateMetrObj(Tuning, DatObj)
  Para <- CreatePara(Starting, DatObj, HyPara)
  DatAug <- CreateDatAug(DatObj)
  McmcObj <- CreateMcmc(MCMC, DatObj)
  RawSamples <- CreateStorage(DatObj, McmcObj)

  ###Time MCMC sampler
  BeginTime <- Sys.time()

  ###Run MCMC sampler in Rcpp
  RegObj <- STBDwDM_Rcpp(DatObj, HyPara, MetrObj, Para, DatAug, McmcObj, RawSamples, Interactive)

  ###Set regression objects
  RawSamples <- RegObj$rawsamples
  MetrObj <- RegObj$metropolis

  ###End time
  FinishTime <- Sys.time()
  RunTime <- FinishTime - BeginTime

  ###Collect output to be returned
  DatObjOut <- OutputDatObj(DatObj)
  DatAugOut <- OutputDatAug(DatAug)
  Metropolis <- SummarizeMetropolis(DatObj, MetrObj, McmcObj)
  Samples <- FormatSamples(DatObj, RawSamples)

  ###Return spBDwDM object
  STBDwDM <- list(mu = Samples$Mu,
                  tau2 = Samples$Tau2,
                  alpha = Samples$Alpha,
                  delta = Samples$Delta,
                  T = Samples$T,
                  phi = Samples$Phi,
                  metropolis = Metropolis,
                  datobj = DatObjOut,
                  dataug = DatAugOut,
                  runtime = paste0("Model runtime: ",round(RunTime, 2), " ",attr(RunTime, "units")))
  STBDwDM <- structure(STBDwDM, class = "STBDwDM")
  return(STBDwDM)

###End STBDwDM function
}
