#' Calculate the number of recruits
#'
#' @description Calculates the number of recruits of each species in a year.
#' @param SSB A numeric vector of length \code{nfish} representing the Spawning Stock Biomass (SSB) of each species.
#' @param stored_rec_funs A list object of length \code{nfish} that is used to store the generated stock recruitment functions.
#' @param rec_args A list object of length \code{nfish} that is used to store the parameters of the stock recruitment functions.
#' @return A numeric vector of length \code{nfish} representing the number of new recruits of each species.
#' @seealso \code{\link{get_rec_fun}}, \code{\link{make_rec_fun}}, \code{\link{rec_BH}}, \code{\link{rec_Ricker}}, \code{\link{rec_hockey}}, \code{\link{rec_const}}, \code{\link{rec_linear}} and \code{\link{calc_SSB}}
#' @references Barrowman, N.J., Myers, R.A. (2000).  Still more spawner-recruit curves: the hockey stick and its generalisations. \emph{Canadian Journal of Fisheries and Aquatic Science}, 57:665–676.
#' @references Beverton, R.J.H., Holt, S.J. (1957). On the Dynamics of Exploited Fish Populations, volume 19 of Fisheries Investigations (Series 2). United Kingdom Ministry of Agriculture and Fisheries.
#' @references Hall, S. J., Collie, J. S., Duplisea, D. E., Jennings, S., Bravington, M., & Link, J. (2006). A length-based multispecies model for evaluating community responses to fishing. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 63(6):1344-1359.
#' @references Ogle, D.H. (2016). Introductory Fisheries Analyses with R. CRC Press.
#' @references Ricker, W.E. (1954). Stock and recruitment. \emph{Journal of the Fisheries Research Board of Canada}, 11:559-623.
#' @references Thorpe, R.B., Le Quesne, W.J.F., Luxford, F., Collie, J.S., Jennings, S. (2015). Evaluation and management implications of uncertainty in a multispecies size-structured model of population and community responses to fishing. \emph{Methods in Ecology and Evolution}, 6:49-58.
#' @examples
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
#' phi_min <- tmp$phi_min
#'
#' # Run calc_ration_growthfac()
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#'
#' # Calculate maturity
#' mature <- calc_mature(Lmat, nfish, mid, kappa=rep(10, nfish), sc_Linf)
#'
#' # Create recruitment functions
#' stored_rec_funs <- get_rec_fun(rep("hockey-stick", nfish))
#' recruit_params <- do.call("Map", c(c, list(a=NS_par$a, b=NS_par$b)))
#'
#' # Get an initial population
#' N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
#'
#' # Calculate the SSB
#' SSB <- calc_SSB(mature, N0, wgt)
#'
#' # Calculate the number of recruits
#' R <- calc_recruits(SSB, stored_rec_funs, recruit_params)
#' @export
calc_recruits <- function(SSB, stored_rec_funs, rec_args){
  rec <- mapply(function(stored_rec_funs, SSB, rec_args){stored_rec_funs(SSB, rec_args)}, SSB=SSB, stored_rec_funs=stored_rec_funs, rec_args=rec_args)
  return(rec)
}

#' Collate the stock recruitment functions
#'
#' @description Collates the stock recruitment functions for all of the species in the model.
#' @param rec_fun A character vector representing the stock recruitment function to be applied to each species. The default is \code{"hockey-stick"}, but \code{rec_fun} can also take \code{"Ricker"}, \code{"Beverton-Holt"}, \code{"constant"}, or \code{"linear"} for each species.
#' @details For \code{"Beverton-Holt"}, the stock recruitment function is defined as \code{a*SSB/(1+b*SSB)}; for \code{"Ricker"} it is defined as \code{a*SSB*exp(-b*SSB)}; for \code{"hockey-stick"} it is defined as \code{min(a*SSB, b)}; for \code{"constant"} it is defined as \code{a}, and for \code{"linear"} it is defined as \code{a*SSB}. In all cases, \code{SSB} is the Spawning Stock Biomass in 1000s of tonnes and \code{a} and \code{b} are parameters of the specific stock recruitment functions.
#' @return A list object of length \code{rec_fun} where each element includes the stock recruitment function for a given species. If an invalid recruitment function is selected, \code{NULL} is returned and a warning message is shown.
#' @seealso \code{\link{calc_recruits}}, \code{\link{make_rec_fun}}, \code{\link{rec_BH}}, \code{\link{rec_Ricker}}, \code{\link{rec_hockey}}, \code{\link{rec_const}}, \code{\link{rec_linear}} and \code{\link{calc_SSB}}
#' @references Barrowman, N.J., Myers, R.A. (2000).  Still more spawner-recruit curves: the hockey stick and its generalisations. \emph{Canadian Journal of Fisheries and Aquatic Science}, 57:665–676.
#' @references Beverton, R.J.H., Holt, S.J. (1957). On the Dynamics of Exploited Fish Populations, volume 19 of Fisheries Investigations (Series 2). United Kingdom Ministry of Agriculture and Fisheries.
#' @references Hall, S. J., Collie, J. S., Duplisea, D. E., Jennings, S., Bravington, M., & Link, J. (2006). A length-based multispecies model for evaluating community responses to fishing. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 63(6):1344-1359.
#' @references Ogle, D.H. (2016). Introductory Fisheries Analyses with R. CRC Press.
#' @references Ricker, W.E. (1954). Stock and recruitment. \emph{Journal of the Fisheries Research Board of Canada}, 11:559-623.
#' @references Thorpe, R.B., Le Quesne, W.J.F., Luxford, F., Collie, J.S., Jennings, S. (2015). Evaluation and management implications of uncertainty in a multispecies size-structured model of population and community responses to fishing. \emph{Methods in Ecology and Evolution}, 6:49-58.
#' @examples
#' nfish <- nrow(NS_par)
#' stored_rec_funs <- get_rec_fun(rep("hockey-stick", nfish))
#' @export
get_rec_fun <- function(rec_fun="hockey-stick"){
  return(sapply(rec_fun, make_rec_fun, simplify=FALSE))
}

#' Generate the stock recruitment functions
#'
#' @description Generates the stock recruitment function for a given species.
#' @param rec_fun A character vector representing the stock recruitment function to be applied to each species. The default is \code{"hockey-stick"}, but \code{rec_fun} can also take \code{"Ricker"}, \code{"Beverton-Holt"}, \code{"constant"}, or \code{"linear"} for each species.
#' @details For \code{"Beverton-Holt"}, the stock recruitment function is defined as \code{a*SSB/(1+b*SSB)}; for \code{"Ricker"} it is defined as \code{a*SSB*exp(-b*SSB)}; for \code{"hockey-stick"} it is defined as \code{min(a*SSB, b)}; for \code{"constant"} it is defined as \code{a}, and for \code{"linear"} it is defined as \code{a*SSB}. In all cases, \code{SSB} is the Spawning Stock Biomass in 1000s of tonnes and \code{a} and \code{b} are parameters of the specific stock recruitment functions.
#' @return The stock recruitment function for a given species. If an invalid recruitment function is selected, \code{NULL} is returned and a warning message is shown.
#' @seealso \code{\link{calc_recruits}}, \code{\link{get_rec_fun}}, \code{\link{rec_BH}}, \code{\link{rec_Ricker}}, \code{\link{rec_hockey}}, \code{\link{rec_const}}, \code{\link{rec_linear}} and \code{\link{calc_SSB}}
#' @references Barrowman, N.J., Myers, R.A. (2000).  Still more spawner-recruit curves: the hockey stick and its generalisations. \emph{Canadian Journal of Fisheries and Aquatic Science}, 57:665–676.
#' @references Beverton, R.J.H., Holt, S.J. (1957). On the Dynamics of Exploited Fish Populations, volume 19 of Fisheries Investigations (Series 2). United Kingdom Ministry of Agriculture and Fisheries.
#' @references Hall, S. J., Collie, J. S., Duplisea, D. E., Jennings, S., Bravington, M., & Link, J. (2006). A length-based multispecies model for evaluating community responses to fishing. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 63(6):1344-1359.
#' @references Ogle, D.H. (2016). Introductory Fisheries Analyses with R. CRC Press.
#' @references Ricker, W.E. (1954). Stock and recruitment. \emph{Journal of the Fisheries Research Board of Canada}, 11:559-623.
#' @references Thorpe, R.B., Le Quesne, W.J.F., Luxford, F., Collie, J.S., Jennings, S. (2015). Evaluation and management implications of uncertainty in a multispecies size-structured model of population and community responses to fishing. \emph{Methods in Ecology and Evolution}, 6:49-58.
#' @examples
#' # Run the function
#' make_rec_fun("Ricker")
#' @export
make_rec_fun <- function(rec_fun="hockey-stick"){
  if (rec_fun == "hockey-stick"){
    return(rec_hockey)
  }
  if (rec_fun == "Ricker"){
    return(rec_Ricker)
  }
  if (rec_fun == "Beverton-Holt"){
    return(rec_BH)
  }
  if (rec_fun == "constant"){
    return(rec_const)
  }
  if (rec_fun == "linear"){
    return(rec_linear)
  }
  warning("Invalid recruitment function. Please see the help file (?get_rec_fun) for valid options.")
  return(NULL)
}

#' The Beverton-Holt stock recruitment function
#'
#' @description Calculates the number of recruits as given by the Beverton-Holt stock recruitment function.
#' @param SSB A numeric value representing the Spawning Stock Biomass (SSB) of a given species (g).
#' @param rec_args A list object of length \code{nfish}, with each element in the list including a value of \code{a} and \code{b} for each species. \code{a} is a positive numeric value, often referred to as the \emph{density-independent} parameter. The default is 1. \code{b} is a positive numeric value, often referred to as the \emph{density-dependent} parameter. The default is 0.001.
#' @details The Beverton-Holt stock recruitment function is defined as \code{a*(SSB/1e9))/(1+b*(SSB/1e9))}.
#' @return A numeric value representing the number of recruits of a given species.
#' @seealso \code{\link{calc_recruits}}, \code{\link{make_rec_fun}}, \code{\link{get_rec_fun}}, \code{\link{rec_Ricker}}, \code{\link{rec_hockey}}, \code{\link{rec_const}}, \code{\link{rec_linear}} and \code{\link{calc_SSB}}
#' @references Beverton, R.J.H., Holt, S.J. (1957). On the Dynamics of Exploited Fish Populations, volume 19 of Fisheries Investigations (Series 2). United Kingdom Ministry of Agriculture and Fisheries.
#' @examples
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
#' phi_min <- tmp$phi_min
#'
#' # Run calc_ration_growthfac()
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#'
#' # Calculate maturity
#' mature <- calc_mature(Lmat, nfish, mid, kappa=rep(10, nfish), sc_Linf)
#'
#' # Create recruitment functions
#' stored_rec_funs <- get_rec_fun(rep("hockey-stick", nfish))
#' recruit_params <- do.call("Map", c(c, list(a=NS_par$a, b=NS_par$b)))
#'
#' # Get an initial population
#' N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
#'
#' # Calculate the SSB
#' SSB <- calc_SSB(mature, N0, wgt)
#'
#' # Run the function
#' rec_BH(SSB[1], recruit_params[[1]])
#' @export
rec_BH <- function(SSB, rec_args){
  SSB <- SSB/1e9 # SSB in tonnes x 10^3
  return(rec_args["a"]*SSB/(1+rec_args["b"]*SSB))
}

#' The Ricker stock recruitment function
#'
#' @description Calculates the number of recruits as given by the Ricker stock recruitment function.
#' @param SSB A numeric value representing the Spawning Stock Biomass (SSB) of a given species (g).
#' @param rec_args  A list object of length \code{nfish}, with each element in the list including a value of \code{a} and \code{b} for each species. \code{a} is a positive numeric value, often referred to as the \emph{density-independent} parameter. The default is 1. \code{b} is a positive numeric value, often referred to as the \emph{density-dependent} parameter. The default is 0.001.
#' @details The Ricker stock recruitment function is defined as \code{a*(SSB/1e9))*exp(-b*(SSB/1e9))}.
#' @return A numeric value representing the number of recruits of a given species.
#' @seealso \code{\link{calc_recruits}}, \code{\link{make_rec_fun}}, \code{\link{get_rec_fun}}, \code{\link{rec_BH}}, \code{\link{rec_hockey}}, \code{\link{rec_const}}, \code{\link{rec_linear}} and \code{\link{calc_SSB}}
#' @references Hall, S. J., Collie, J. S., Duplisea, D. E., Jennings, S., Bravington, M., & Link, J. (2006). A length-based multispecies model for evaluating community responses to fishing. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 63(6):1344-1359.
#' @references Ricker, W.E. (1954). Stock and recruitment. \emph{Journal of the Fisheries Research Board of Canada}, 11:559-623.
#' @examples
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
#' phi_min <- tmp$phi_min
#'
#' # Run calc_ration_growthfac()
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#'
#' # Calculate maturity
#' mature <- calc_mature(Lmat, nfish, mid, kappa=rep(10, nfish), sc_Linf)
#'
#' # Create recruitment functions
#' stored_rec_funs <- get_rec_fun(rep("hockey-stick", nfish))
#' recruit_params <- do.call("Map", c(c, list(a=NS_par$a, b=NS_par$b)))
#'
#' # Get an initial population
#' N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
#'
#' # Calculate the SSB
#' SSB <- calc_SSB(mature, N0, wgt)
#'
#' rec_Ricker(SSB[1], recruit_params[[1]])
#' @export
rec_Ricker <- function(SSB, rec_args){
  SSB <- SSB/1e9 # SSB in tonnes x 10^3
  return(rec_args["a"]*SSB*exp(-rec_args["b"]*SSB)*1e6)
}

#' The hockey-stick stock recruitment function
#'
#' @description Calculates the number of recruits as given by the hockey-stick stock recruitment function.
#' @param SSB A numeric value representing the Spawning Stock Biomass (SSB) of a given species (g).
#' @param rec_args A list object of length \code{nfish}, with each element in the list including a value of \code{a} and \code{b} for each species. \code{a} is a positive numeric value, often referred to as the \emph{density-independent} parameter. The default is 1. \code{b} is a positive numeric value, often referred to as the \emph{density-dependent} parameter. The default is 0.001.
#' @details The stock recruitment function is defined as \code{min(a*(SSB/1e9), b)}.
#' @return A numeric value representing the number of recruits of a given species.
#' @seealso \code{\link{calc_recruits}}, \code{\link{make_rec_fun}}, \code{\link{get_rec_fun}}, \code{\link{rec_BH}}, \code{\link{rec_Ricker}}, \code{\link{rec_const}}, \code{\link{rec_linear}} and \code{\link{calc_SSB}}
#' @references Barrowman, N.J., Myers, R.A. (2000).  Still more spawner-recruit curves: the hockey stick and its generalisations. \emph{Canadian Journal of Fisheries and Aquatic Science}, 57:665–676.
#' @references Thorpe, R.B., Le Quesne, W.J.F., Luxford, F., Collie, J.S., Jennings, S. (2015). Evaluation and management implications of uncertainty in a multispecies size-structured model of population and community responses to fishing. \emph{Methods in Ecology and Evolution}, 6:49-58.
#' @examples
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
#' phi_min <- tmp$phi_min
#'
#' # Run calc_ration_growthfac()
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#'
#' # Calculate maturity
#' mature <- calc_mature(Lmat, nfish, mid, kappa=rep(10, nfish), sc_Linf)
#'
#' # Create recruitment functions
#' stored_rec_funs <- get_rec_fun(rep("hockey-stick", nfish))
#' recruit_params <- do.call("Map", c(c, list(a=NS_par$a, b=NS_par$b)))
#'
#' # Get an initial population
#' N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
#'
#' # Calculate the SSB
#' SSB <- calc_SSB(mature, N0, wgt)
#'
#' rec_hockey(SSB[1], recruit_params[[1]])
#' @export
rec_hockey <- function(SSB, rec_args){
  SSB <- SSB/1e9 # SSB in tonnes x 10^3
  return(min(rec_args["a"]*SSB, rec_args["b"])*1e6)
}

#' The constant stock recruitment function
#'
#' @description Calculates the number of recruits as given by the constant stock recruitment function.
#' @param SSB A numeric value representing the Spawning Stock Biomass (SSB) of a given species (g).
#' @param rec_args A list object of length \code{nfish}, with each element in the list including a value of \code{a} for each species. \code{a} is a positive numeric value, often referred to as the \emph{density-independent} parameter. The default is 1.
#' @details The number of recruits is defined as the value of \code{a}.
#' @return A numeric value representing the number of recruits of a given species.
#' @seealso \code{\link{calc_recruits}}, \code{\link{make_rec_fun}}, \code{\link{get_rec_fun}}, \code{\link{rec_BH}}, \code{\link{rec_Ricker}}, \code{\link{rec_hockey}}, \code{\link{rec_linear}} and \code{\link{calc_SSB}}
#' @examples
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
#' phi_min <- tmp$phi_min
#'
#' # Run calc_ration_growthfac()
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#'
#' # Calculate maturity
#' mature <- calc_mature(Lmat, nfish, mid, kappa=rep(10, nfish), sc_Linf)
#'
#' # Create recruitment functions
#' stored_rec_funs <- get_rec_fun(rep("hockey-stick", nfish))
#' recruit_params <- do.call("Map", c(c, list(a=NS_par$a, b=NS_par$b)))
#'
#' # Get an initial population
#' N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
#'
#' # Calculate the SSB
#' SSB <- calc_SSB(mature, N0, wgt)
#'
#' rec_const(SSB[1], recruit_params[[1]])
#' @export
rec_const <- function(SSB, rec_args){
  return(rec_args["a"]*1e6)
}

#' The density-independent stock recruitment function
#'
#' @description Calculates the number of recruits as given by the density-independent stock recruitment function.
#' @param SSB A numeric value representing the Spawning Stock Biomass (SSB) of a given species (g).
#' @param rec_args A list object of length \code{nfish}, with each element in the list including a value of \code{a} for each species. \code{a} is a positive numeric value, often referred to as the \emph{density-independent} parameter. The default is 1.
#' @details The number of recruits is defined as \code{a*(SSB/1e9)}.
#' @return A numeric value representing the number of recruits of a given species.
#' @seealso \code{\link{calc_recruits}}, \code{\link{make_rec_fun}}, \code{\link{get_rec_fun}}, \code{\link{rec_BH}}, \code{\link{rec_Ricker}}, \code{\link{rec_hockey}}, \code{\link{rec_const}} and \code{\link{calc_SSB}}
#' @references Ogle, D.H. (2016). Introductory Fisheries Analyses with R. CRC Press.
#' @examples
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
#' phi_min <- tmp$phi_min
#'
#' # Run calc_ration_growthfac()
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#'
#' # Calculate maturity
#' mature <- calc_mature(Lmat, nfish, mid, kappa=rep(10, nfish), sc_Linf)
#'
#' # Create recruitment functions
#' stored_rec_funs <- get_rec_fun(rep("hockey-stick", nfish))
#' recruit_params <- do.call("Map", c(c, list(a=NS_par$a, b=NS_par$b)))
#'
#' # Get an initial population
#' N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
#'
#' # Calculate the SSB
#' SSB <- calc_SSB(mature, N0, wgt)
#'
#' rec_linear(SSB[1], recruit_params[[1]])
#' @export
rec_linear <- function(SSB, rec_args){
  SSB <- SSB/1e9 # SSB in tonnes x 10^3
  return(rec_args["a"]*SSB*1e6)
}

