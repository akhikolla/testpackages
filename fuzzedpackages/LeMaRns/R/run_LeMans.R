#' An S4 class representing the outputs of the LeMans model
#'
#' @slot N An array with dimensions \code{nsc}, \code{nfish} and \code{tot_time} representing the number of individuals in each length class at each time step of the model.
#' @slot Catch An array with dimensions \code{nsc}, \code{nfish} and \code{tot_time} representing the biomass of each species that is caught in each length class at each time step of the model.
#' @slot M2 An array with dimensions \code{nsc}, \code{nfish} and \code{tot_time} representing the predation mortality of each species in each length class at each time step of the model.
#' @slot R A matrix with dimensions \code{tot_time} and \code{nfish} representing the number of recruits of each species at each time step of the model.
#' @export
setClass("LeMans_outputs", slots=
           c(N="array",
             Catch="array",
             M2="array",
             R="matrix"))

#' Project the LeMans model
#'
#' @description Project the LeMans model forward in time.
#' @param params A \linkS4class{LeMans_param} object containing the parameter values of the current LeMans model.
#' @param years A numeric value representing the number of years that the model is run for. The default is 10.
#' @param N0 A matrix with dimensions \code{nsc} and \code{nfish} representing the number of individuals in each length class when the model is initialised.
#' @param effort A matrix with dimensions \code{years} and the number of fishing gears, representing fishing effort in each year for each gear. This parameter is required only if \code{Fs} is missing.
#' @param Fs An array with dimensions \code{nsc}, \code{nfish} and \code{tot_time} representing the fishing mortality of each species in each length class at each time step.
#' @param intercept A numeric value representing the number of individuals in the first length class. This parameter is only required if \code{N0} is missing. The default is \code{1e10}.
#' @param slope A numeric value representing the slope of the community size spectrum. This parameter is only required if \code{N0} is missing. The default is -5.
#' @param tot_time A numeric value representing the number of time steps to run the model for.
#' @param nsc A numeric value representing the number of length classes in the model.
#' @param nfish A numeric value representing the number of fish species in the model.
#' @param phi_min A numeric value representing the time step of the model.
#' @param mature A matrix with dimensions \code{nsc} and \code{nfish} and elements in the range 0-1 representing the proportion of individuals that are mature for each species and length class.
#' @param sc_Linf A numeric vector of length \code{nsc} representing the length class at which each species reaches its asymptotic length.
#' @param wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class.
#' @param phi A matrix with dimensions \code{nsc} and \code{nfish} representing the proportion of individuals that leave each length class.
#' @param ration A matrix with dimensions \code{nsc} and \code{nfish} representing the amount of food required for fish of a given species and length class to grow according to the von Bertalanffy growth curve in a time step.
#' @param other A numeric value representing the amount of other food (g) available from prey that is not explicitly represented in the model.
#' @param M1 A matrix of dimensions \code{nsc} and \code{nfish} representing the natural mortality of each species for each length class.
#' @param suit_M2 A list object of length \code{nfish}. Each element in the list is an array of dimensions \code{nsc}, \code{nsc} and \code{nfish} containing a value between zero and 1 representing prey preference and prey suitability for each species and length class.
#' @param stored_rec_funs A list object of length \code{nfish} where each element includes the stock recruitment function for each species. If an invalid recruitment function is selected, \code{NULL} is returned and a warning message is shown.
#' @param recruit_params A list object of length \code{nfish} specifying the parameters for the recruitment function.
#' @param eps A numeric value specifying a numerical offset. The default value is \code{1e-5}.
#' @param ... Additional arguments.
#' @return An object of class \code{\linkS4class{LeMans_outputs}}.
#' @seealso \code{\linkS4class{LeMans_outputs}}, \code{\linkS4class{LeMans_param}}, \code{\link{LeMansParam}}
#' @examples
#' # Run the model with all inputs specified explicitly:
#' # Set up the inputs to the function - species-independent parameters
#' nfish <- nrow(NS_par)
#' nsc <- 32
#' maxsize <- max(NS_par$Linf)*1.01 # the biggest size is 1% bigger than the largest Linf
#' l_bound <- seq(0, maxsize, maxsize/nsc); l_bound <- l_bound[-length(l_bound)]
#' u_bound <- seq(maxsize/nsc, maxsize, maxsize/nsc)
#' mid <- l_bound+(u_bound-l_bound)/2
#'
#' # Set up the inputs to the function - species-specific parameters
#' Linf <- NS_par$Linf # the von-Bertalanffy asymptotic length of each species (cm).
#' W_a <- NS_par$W_a # length-weight conversion parameter.
#' W_b <- NS_par$W_b # length-weight conversion parameter.
#' k <- NS_par$k # the von-Bertalnaffy growth parameter.
#' Lmat <- NS_par$Lmat # the length at which 50\% of individuals are mature (cm).
#'
#' # Get phi_min
#' tmp <- calc_phi(k, Linf, nsc, nfish, u_bound, l_bound, calc_phi_min=FALSE,
#'                   phi_min=0.1) # fixed phi_min
#' phi <- tmp$phi
#' phi_min <- tmp$phi_min
#'
#' # Calculate growth increments
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' ration <- tmp$ration
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#' g_eff <- tmp$g_eff
#'
#' # Calculate maturity
#' mature <- calc_mature(Lmat, nfish, mid, kappa=rep(10, nfish), sc_Linf)
#'
#' # Create recruitment functions
#' stored_rec_funs <- get_rec_fun(rep("hockey-stick", nfish))
#' recruit_params <- do.call("Map", c(c, list(a=NS_par$a, b=NS_par$b)))
#'
#' # Calculate background mortality
#' M1 <- calc_M1(nsc, sc_Linf, phi_min)
#'
#' # Calculate predator-prey size preferences
#' prefs <- calc_prefs(pred_mu=-2.25, pred_sigma=0.5, wgt, sc_Linf)
#'
#' # Calculate prey preference and prey suitability
#' suit_M2 <- calc_suit_vect(nsc, nfish, sc_Linf, prefs, NS_tau)
#'
#' # Calculate catchability
#' Qs <- calc_Q(curve=rep("logistic", nfish), species=NS_par$species_names,
#'              max_catchability=rep(1, nfish), gear_name=NS_par$species_names,
#'              nsc=nsc, nfish=nfish, mid=mid, l_bound=l_bound, u_bound=u_bound,
#'              species_names=NS_par$species_names, eta=rep(0.25, nfish), L50=Lmat)
#'
#' # Get an initial population
#' N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
#' years <- 10 # run the model for 10 years
#' tot_time <- years*phi_min # total number of time steps
#'
#' # Define fishing effort to be 0.5 for all species
#' effort <- matrix(0.5, tot_time, dim(Qs)[3])
#'
#' # Calculate F
#' Fs <- array(0, dim=c(nsc, nfish, tot_time))
#' for (j in 1:ncol(effort)) {
#'   for (ts in 1:tot_time) {
#'     Fs[,,ts] <- Fs[,,ts]+effort[ts, j]*Qs[,,j]
#'   }
#' }
#'
#' # Run the model
#' model_run <- run_LeMans(N0=N0, tot_time=tot_time, Fs=Fs, nsc=nsc, nfish=nfish,
#'                         phi_min=phi_min, mature=mature, sc_Linf=sc_Linf, wgt=wgt,
#'                         phi=phi, ration=ration, other=NS_other, M1=M1, suit_M2=suit_M2,
#'                         stored_rec_funs=stored_rec_funs, recruit_params=recruit_params,
#'                         eps=1e-05)
#'
#' ##############################################
#' # Alternatively:
#' NS_params <- LeMansParam(NS_par,tau=NS_tau,eta=rep(0.25,21),L50=NS_par$Lmat,other=NS_other)
#'
#' # Define fishing effort
#' effort <- matrix(0.5, 10, dim(NS_params@Qs)[3])
#'
#' # Run the model
#' model_run <- run_LeMans(NS_params, years=10, effort=effort)
#' @export
setGeneric("run_LeMans", function(params, ...)
  standardGeneric('run_LeMans'))

# Script to be used when signature(params = "missing")
#' @rdname run_LeMans
setMethod("run_LeMans", signature(params = "missing"),
function(N0, Fs, tot_time, nsc, nfish, phi_min, mature, sc_Linf, wgt, phi, ration, other, M1, suit_M2, stored_rec_funs, recruit_params, eps=1e-5){
  M2 <- Catch <- N <- array(0, dim=c(nsc, nfish, tot_time+1))
  R <- matrix(0, tot_time+1, nfish)
  SSB <- rep(0, nfish)
  N[,,1] <- N0
  for (ts in 2:(tot_time+1)) {
    N[,,ts] = N[,,ts-1]
    # Recruitment
    if (ceiling((ts-1)*phi_min)!=ceiling((ts-2)*phi_min)) { # check whether it is the start of the year
      SSB <- calc_SSB(mature, N[,,ts], wgt)
      R[ts,] <- calc_recruits(SSB, stored_rec_funs, recruit_params)
      N[1,,ts] <- N[1,,ts] + R[ts, ]
    }

    # Mortality
    M2[,,ts] <- calc_M2(N[,,ts], ration, wgt, nfish, nsc, other, sc_Linf, suit_M2)
    Z <- Fs[,,ts-1]*phi_min+M1+M2[,,ts]+eps
    if (max(Fs[,,ts-1])>0) {
      Catch[,,ts] <- (Fs[,,ts-1]*phi_min/Z)*N[,,ts]*(1-exp(-Z))*wgt
    }
    N[,,ts] <- N[,,ts]*exp(-Z)

    # Growth
    N[,,ts] <- calc_growth(N[,,ts], phi, nfish, nsc)
  }
  return(new("LeMans_outputs",
             N=N,
             Catch=Catch,
             M2=M2,
             R=R))
})

# Script to be used when signature(params = "LeMans_param")
#' @rdname run_LeMans
setMethod("run_LeMans", signature(params = "LeMans_param"),
function(params, years=10, N0=NULL, effort=matrix(0, years, dim(params@Qs)[3]), Fs, intercept=1e10, slope=-5, tot_time) {
  if (is.null(N0)) {
    N0 <- get_N0(params@nsc, params@nfish, params@mid, params@wgt, params@sc_Linf, intercept, slope)
  }
  if (years<=0) {
    stop("years must be positive")
  }
  if (missing(tot_time)){
    tot_time <- floor(years/params@phi_min)
  }
  if ((1/params@phi_min)%%1!=0) {
    warning("The time step is not a multiple of a year")
  }
  if (missing(Fs)) {
    Fs <- array(0, dim=c(params@nsc, params@nfish, tot_time))
    if (ncol(effort)!=dim(params@Qs)[3]) {
      stop("Incorrect number of gears in effort")
    }
    year <- ceiling(1:tot_time*params@phi_min)
    if (is.null(colnames(effort))) {
      for (j in 1:ncol(effort)) {
        for (ts in 1:tot_time) {
          Fs[,,ts] <- Fs[,,ts]+effort[year[ts], j]*params@Qs[,,j]
        }
      }
    } else {
      if (identical(colnames(effort), dimnames(params@Qs)[[3]])) {
        for (j in colnames(effort)){
          for (ts in 1:tot_time){
            Fs[,,ts] <- Fs[,,ts]+effort[year[ts], j]*params@Qs[,,j]
          }
        }
      } else {
        stop("Incorrect gear names in effort")
      }
    }
  } else {
    if (identical(as.numeric(dim(Fs)), c(params@nsc,params@nfish, tot_time))==FALSE) {
      stop("Fs must be of dimensions nsc, nfish and tot_time")
    }
  }
  return(run_LeMans(N0=N0, tot_time=tot_time, Fs=Fs, nsc=params@nsc, nfish=params@nfish, phi_min=params@phi_min, mature=params@mature, sc_Linf=params@sc_Linf, wgt=params@wgt, phi=params@phi, ration=params@ration, other=params@other, M1=params@M1, suit_M2=params@suit_M2, stored_rec_funs=params@stored_rec_funs, recruit_params=params@recruit_params, eps=params@eps))
})

#' Generate a starting value for N
#'
#' @description Generate a starting value for \code{N}, which represents the number of individuals in each length class for each species.
#' @param nsc A numeric value representing the number of length classes in the model.
#' @param nfish A numeric value representing the number of fish species in the model.
#' @param mid A numeric vector of length \code{nfish} representing the mid-point of the length classes.
#' @param wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class.
#' @param sc_Linf A numeric vector of length \code{nsc} representing the length class at which each species reaches its asymptotic length.
#' @param intercept A numeric value representing the number of individuals in the first length class. The default is \code{1e10}.
#' @param slope A numeric value representing the slope of the community size spectrum. The default is -5.
#' @details The total number of individuals in the community in each length class is equal to \code{intercept*mid^slope}. Within each length class, the number of individuals of each species is determined using the proportion of each species' biomass that is found in that length class.
#' @return A matrix with dimensions \code{nsc} and \code{nfish} representing the number of individuals in each length class.
#' @seealso \code{\link{run_LeMans}}
#' @references Andersen, K.H., Blanchard, J.L., Fulton, E.A., Gislason, H., Jacobsen, N.S., van Kooten, T. (2016). Assumptions behind size-based ecosystem models are realistic. \emph{ICES Journal of Marine Science}, 73(6):1651-1655.
#' @examples
#' # Set up the inputs to the function
#' nfish <- nrow(NS_par)
#' nsc <- 32
#' maxsize <- max(NS_par$Linf)*1.01 # the biggest size is 1% bigger than the largest Linf
#' l_bound <- seq(0, maxsize, maxsize/nsc); l_bound <- l_bound[-length(l_bound)]
#' u_bound <- seq(maxsize/nsc, maxsize, maxsize/nsc)
#' mid <- l_bound+(u_bound-l_bound)/2
#'
#' # Set up the inputs to the function - species-specific parameters
#' Linf <- NS_par$Linf # the von-Bertalanffy asymptotic length of each species (cm).
#' W_a <- NS_par$W_a # length-weight conversion parameter.
#' W_b <- NS_par$W_b # length-weight conversion parameter.
#' k <- NS_par$k # the von-Bertalnaffy growth parameter.
#' Lmat <- NS_par$Lmat # the length at which 50\% of individuals are mature (cm).
#'
#' # Get phi_min
#' tmp <- calc_phi(k, Linf, nsc, nfish, u_bound, l_bound, calc_phi_min=FALSE,
#'                   phi_min=0.1) # fixed phi_min
#' phi_min <- tmp$phi_min
#'
#' # Calculate growth increments
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#'
#' # Get an initial population
#' get_N0(nsc, nfish, mid, wgt, sc_Linf)
#' @export
get_N0 <- function(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5) {
  tot_n <- intercept*mid^slope
  tmp <- matrix(tot_n, nsc, nfish)
  prop_l <- sapply(sc_Linf, function(x, tot_n, nsc){tmp <- (tot_n)[1:x];c(tmp/sum(tmp), rep(0, nsc-x))}, tot_n=tot_n, nsc=nsc)
  return(tmp*t(apply(prop_l*matrix(1/colSums(prop_l*wgt), nsc, nfish, byrow = T), 1, function(x){x/sum(x)})))
}

#' Combine two \code{LeMans_outputs} objects
#'
#' @description Combines two \linkS4class{LeMans_outputs} objects.
#' @param LeMans_run_x A \linkS4class{LeMans_outputs} object.
#' @param LeMans_run_y A \linkS4class{LeMans_outputs} object.
#' @param cont A logical statement indicating whether or not \code{LeMans_run_y} is a continuation of \code{LeMans_run_x}. The default is \code{TRUE}.
#' @details If \code{cont==T}, the first years output from \code{LeMans_run_y} is removed as this will be the same as the last year of \code{LeMans_run_x}.
#' @return A \code{LeMans_outputs} object.
#' @seealso \linkS4class{LeMans_outputs}, \code{\link{run_LeMans}}
#' @examples
#' # Set up the inputs to the model
#' NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=NS_other)
#'
#' # Define fishing effort
#' effort <- matrix(0.5, 10, dim(NS_params@Qs)[3])
#'
#' # Run the model for the first time
#' model_run1 <- run_LeMans(NS_params, years=10, effort=effort)
#'
#' # Run the model for another 12 years
#' effort1 <- matrix(0.5, 12, dim(NS_params@Qs)[3])
#' model_run2 <- run_LeMans(N=model_run1@N[,,101], NS_params, years=12, effort=effort1)
#'
#' # Combine the two model runs into a single run
#' out <- comb_LeMans_run(model_run1, model_run2, cont=TRUE)
#' @export
comb_LeMans_run <- function(LeMans_run_x, LeMans_run_y, cont=TRUE) {
  if (class(LeMans_run_x)!=class(LeMans_run_y)|class(LeMans_run_x)!="LeMans_outputs") {
    stop("LeMans_run_x and LeMans_run_y must be of class LeMans_outputs")
  }
  if (cont==T) {
    N <- abind(LeMans_run_x@N, LeMans_run_y@N[,,-1])
    Catch <- abind(LeMans_run_x@Catch, LeMans_run_y@Catch[,,-1])
    M2 <- abind(LeMans_run_x@M2, LeMans_run_y@M2[,,-1])
    R <- rbind(LeMans_run_x@R, LeMans_run_y@R[-1, ])
  } else {
    N <- abind(LeMans_run_x@N, LeMans_run_y@N)
    Catch <- abind(LeMans_run_x@Catch, LeMans_run_y@Catch)
    M2 <- abind(LeMans_run_x@M2, LeMans_run_y@M2)
    R <- rbind(LeMans_run_x@R, LeMans_run_y@R)
  }
  return(new("LeMans_outputs", N=N, Catch=Catch, M2=M2, R=R))
}
