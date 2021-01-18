#' Calculate gear catchability
#'
#' @description Calculates the catchability for each fishing gear.
#' @param curve A character vector of almost any length describing the type of curve to be used to determine the catchability of each species by the fishing gear. By default, \code{curve} is a character vector of length \code{nfish} that takes the value \code{"logistic"} in each element, but it can also take \code{"log-gaussian"} or \code{"knife-edge"}. If a custom curve is required for a particular species and/or fishing gear, the curve must be specified in \code{custom}. See 'Details' for more information.
#' @param species A numeric or character vector of the same length as \code{curve} describing the species that each element of \code{curve} relates to. By default, \code{species} is a numeric vector of length \code{curve} that takes the values \code{(0:(length(curve)-1))\%\%(nfish)+1}.
#' @param max_catchability A numeric vector of length \code{curve} describing the maximum catchability for each catchability curve.
#' @param gear_name A character vector of the same length as \code{curve} and \code{species} describing the fishing gear that each element of \code{curve} and \code{species} relates to. By default, \code{gear_name} is a character vector of length \code{curve} that takes the value \code{paste("gear_", 1:length(curve), sep = "")}.
#' @param custom An array of dimensions \code{nsc}, \code{nfish} and the number of custom catchability curves that are required. \code{custom} represents the catchability of each species by the gears specified using custom catchability curves. By default, \code{custom} is set to \code{NULL}.
#' @param nsc A numeric value representing the number of length classes in the model.
#' @param nfish A numeric value representing the number of fish species in the model.
#' @param l_bound A numeric vector of length \code{nsc} representing the lower bounds of the length classes.
#' @param u_bound A numeric vector of length \code{nsc} representing the upper bounds of the length classes.
#' @param mid A numeric vector of length \code{nfish} representing the mid-point of the length classes.
#' @param eps A numeric value specifying a numerical offset. The default is \code{1e-5}.
#' @param species_names A character vector of length \code{nfish} that is the name of the species in the same order as the species parameters.
#' @param ... Vectors of the same length as \code{curve} representing the parameters of the catchability curves. See 'Details' for more information.
#' @param L50 A numeric value representing the length at 50\% of the maximum catchability of the catchability curve. This is used with the \code{logistic_catch} function. The default value is 50.
#' @param eta A numeric value representing the steepness of the slope of the catchability curve. This is used with the \code{logistic_catch} function. The default value is 0.25
#' @param Lmu A numeric value representing the length at the maximum catchability of the catchability curve. This is used with the \code{log_gaussian_catch} function. The default value is 50.
#' @param Lsigma A numeric value representing the standard deviation of the catchability curve. See 'Details' for more information. This is used with the \code{log_gaussian_catch} function. The default value is 1.
#' @param Lmin A numeric value representing the minimum length that is caught by the catchability curve. This is used with the \code{knife_edge_catch} function. The default value is 50.
#' @details This function allows three different models for catchability, all of which are rescaled so that their maximum catchability is equal to one:
#' @details (1) \code{"logistic"} is proportional to
#' @details \code{1/(1+exp(-eta*(mid-L50)))};
#' @details (2) \code{"log-gaussian"} is proportional to
#' @details \code{dlnorm(mid, log(Lmu), Lsigma)};
#' @details and (3) \code{"knife-edge"} is equal to 1 for the length classes indexed by \code{max(which(l_bound<Lmin)):nsc} and 0 for all other length classes. This means that all of the individuals in the length class containing \code{Lmin} will have a catchability of 1.
#' @details In \code{calc_Q} the catchability is rescaled so that the maximum catchability is one.
#' @return \code{calc_Q} returns an array of dimensions \code{nsc}, \code{nfish} and \code{gear} representing the catchability of each species by each of the fishing gears, scaled from 0 to \code{max_catchability}.
#' @return \code{get_Q} returns a catchability curve for a given species and gear scaled from 0 to 1. If an invalid catchability curve is selected, \code{NULL} is returned and a warning message is shown.
#' @return \code{logistic_catch} returns a matrix with dimensions \code{nsc} and \code{nfish} with all elements set to zero excluding one column that represents the logistic catchability curve for the selected species scaled from 0 to 1.
#' @return \code{log_gaussian_catch} returns a matrix with dimensions \code{nsc} and \code{nfish} with all elements set to zero excluding one column that represents the log-gaussian catchability curve for the selected species scaled from 0 to 1.
#' @return \code{knife_edge_catch} returns a matrix with dimensions \code{nsc} and \code{nfish} with all elements set to zero excluding one column that represents the knife-edge catchability curve for the selected species scaled from 0 to 1.
#' @seealso \code{\link[stats]{dlnorm}}
#' @references Hall, S. J., Collie, J. S., Duplisea, D. E., Jennings, S., Bravington, M., & Link, J. (2006). A length-based multispecies model for evaluating community responses to fishing. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 63(6):1344-1359.
#' @references Thorpe, R.B., Jennings, S., Dolder, P.J. (2017). Risks and benefits of catching pretty good yield in multispecies mixed fisheries. \emph{ICES Journal of Marine Science}, 74(8):2097-2106.
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
#' # Calculate gear catchability
#' Qs <- calc_Q(curve=rep("logistic", nfish), species=NS_par$species_names,
#'              max_catchability=rep(1, nfish), gear_name=NS_par$species_names,
#'              nsc=nsc, nfish=nfish, mid=mid, l_bound=l_bound, u_bound=u_bound,
#'              species_names=NS_par$species_names, eta=rep(0.25, nfish), L50=Lmat)
#'
#' # Calculate logistic catchability for the first species
#' logistic_catch(species=1, nsc, nfish, mid, eps=1e-5,L50=Lmat[1], eta=0.25)
#'
#' # Calculate log-gaussian catchability for the first species
#' log_gaussian_catch(species=1, nsc, nfish, mid, eps=1e-5, Lmu=50, Lsigma=1)
#'
#' # Calculate knife-edge catchability for the first species
#' knife_edge_catch(species=1, nsc, nfish, l_bound, u_bound, Lmin=50)
#' @export
calc_Q <- function(curve=rep("logistic", nfish), species=((0:(length(curve)-1))%%nfish)+1, max_catchability=rep(1, length(curve)), gear_name=paste("gear_", 1:length(curve), sep=""), custom=NULL, nsc, nfish, l_bound, u_bound, mid, eps=1e-5, species_names,...) {
  if (is.character(species)){
    species <- sapply(species,function(species){which(species==species_names)})
  }
  tmp_Q <- mapply(get_Q, curve=curve, species=species, ...=..., MoreArgs = list(nsc=nsc, nfish=nfish, mid=mid, l_bound=l_bound, u_bound=u_bound, eps=eps))
  Qs <- array(0, dim=c(nsc, nfish, length(unique(gear_name))));
  dimnames(Qs)[[3]] <- unique(gear_name)

  for(i in 1:length(curve)){
    Qs[,,as.character(gear_name[i])] <- tmp_Q[, i]*max_catchability[i]+Qs[,,as.character(gear_name[i])]
  }

  Qs <- abind(Qs, custom)
  return(Qs)
}

#' @rdname calc_Q
#' @export
get_Q <- function(curve, species, gear_name, nsc, nfish, l_bound, u_bound, mid, eps, ...) {
  if (curve == "logistic") {
    return(logistic_catch(species=species, nsc=nsc, nfish=nfish, mid=mid, eps=eps, ...))
  }
  if (curve == "log-gaussian") {
    return(log_gaussian_catch(species=species, nsc=nsc, nfish=nfish, mid=mid, eps=eps, ...))
  }
  if (curve == "knife-edge") {
    return(knife_edge_catch(species=species, nsc=nsc, nfish=nfish, l_bound=l_bound, u_bound=u_bound, ...))
  }
  stop("Invalid curve specified. Please see the help file (help(calc_Q)) for valid options.")
}

#' @rdname calc_Q
#' @export
logistic_catch <- function(species, nsc, nfish, mid, eps, L50=50, eta=0.25, ...){
  ret <- matrix(0, nsc, nfish)
  tmp <- 1/(1+exp(-eta*(mid-L50)))
  tmp <- tmp/max(tmp)
  tmp[tmp<eps] <- 0
  ret[, species] <- tmp
  return(ret)
}

#' @rdname calc_Q
#' @export
log_gaussian_catch <- function(species, nsc, nfish, mid, eps, Lmu=50, Lsigma=1, ...){
  ret <- matrix(0, nsc, nfish)
  tmp <- dlnorm(mid, log(Lmu), Lsigma)
  tmp <- tmp/max(tmp)
  tmp[tmp<eps] <- 0
  ret[, species] <- tmp
  return(ret)
}

#' @rdname calc_Q
#' @export
knife_edge_catch <- function(species, nsc, nfish, l_bound, u_bound, Lmin=50, ...){
  ret <- matrix(0, nsc, nfish)
  if (Lmin>max(u_bound)){
    warning("Lmin cannot be greater than the largest size class")
    return(ret)
  }
  ret[max(which(l_bound<Lmin)):nsc, species] <- 1
  return(ret)
}
