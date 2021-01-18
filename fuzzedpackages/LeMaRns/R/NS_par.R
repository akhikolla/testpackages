#' North Sea data
#'
#' Data for the 21 species in the North Sea version of the LeMans model.
#'
#' @format A data frame with 21 rows and 8 variables, including:
#' \describe{
#' \item{\code{species_names}}{A numeric or character vector of length \code{nfish} that denotes the names of the species in the model.}
#'  \item{\code{Linf}}{A numeric vector of length \code{nfish} representing the asymptotic length of each species (cm).}
#'  \item{\code{W_a}}{A numeric vector of length \code{nfish} representing the parameter \code{a} in the length-weight conversion.}
#'  \item{\code{W_b}}{A numeric vector of length \code{nfish} representing the parameter \code{b} in the length-weight conversion.}
#'  \item{\code{k}}{A numeric vector of length \code{nfish} representing the von Bertalanffy growth parameter \code{(1/yr)} for each species.}
#'  \item{\code{Lmat}}{A numeric vector of length \code{nsc} representing the length at which 50\% of the individuals are mature (cm).}
#'  \item{\code{a}}{A numeric value representing the density independent part of the hockey-stick recruitment curve.}
#'  \item{\code{b}}{A numeric value representing the density dependent part of the hockey-stick recruitment curve.}
#' }
#' @references Thorpe, R.B., Le Quesne, W.J.F., Luxford, F., Collie, J.S., Jennings, S. (2015). Evaluation and management implications of uncertainty in a multispecies size-structured model of population and community responses to fishing. \emph{Methods in Ecology and Evolution}, 6:49-58.
"NS_par"

#' Other food for the North Sea
#'
#' Other food for the North Sea dataset, \code{NS_par}.
#'
#' @format A numerical value representing other food. To be used with \code{NS_par}.
"NS_other"

#' North Sea interaction matrix
#'
#' A predator-prey interaction matrix for the 21 species in \code{NS_par}.
#'
#' @format A matrix with 21 rows and 21 columns:
#' \describe{
#'  Row indices represent predators and column indices represent prey. A value of 1 at location \code{i}, \code{j} indicates prey \code{j} is eaten by predator \code{i}.
#' }
#' @references Thorpe, R.B., Le Quesne, W.J.F., Luxford, F., Collie, J.S., Jennings, S. (2015). Evaluation and management implications of uncertainty in a multispecies size-structured model of population and community responses to fishing. \emph{Methods in Ecology and Evolution}, 6:49-58.
"NS_tau"

#' Gear selectivity data frame
#'
#' A gear selectivity data frame for the 21 species in \code{NS_par}.
#'
#' @format A data frame with 21 rows and 4 variables, including:
#' \describe{
#' \item{\code{catch_species}}{A character string describing the species to apply the catchability curve to.}
#'  \item{\code{curve}}{A character vector describing the type of curve to be used to determine the catchability of each species by the fishing gear.}
#'  \item{\code{gear_name}}{A character string describing the name of the gear.}
#'  \item{\code{max_catchability}}{A numeric vector describing the maximum catchability for each catchability curve.}
#' }
#' @references Thorpe, R.B., Dolder, P.J. , Reeves, S., Robinson, P., Jennings, S. (2015). Assessing fishery and ecological consequences of alternative management options for multispecies fisheries \emph{ICES Journal of Marine Science}, 73(6):1503-1512.
"NS_mixed_fish"

#' The steepness of the slope of the catchability curve
#'
#' The steepness of the slope of the catchability curve associated with \code{NS_mixed_fish}.
#'
#' @format A numeric value representing the steepness of the slope of the catchability curve.
#' @references Thorpe, R.B., Dolder, P.J. , Reeves, S., Robinson, P., Jennings, S. (2015). Assessing fishery and ecological consequences of alternative management options for multispecies fisheries \emph{ICES Journal of Marine Science}, 73(6):1503-1512.
"NS_eta"

#' The length at 50\% of the maximum catchability of the catchability curve
#'
#' The length at 50\% of the maximum catchability of the catchability curve associated with \code{NS_mixed_fish}.
#'
#' @format A numeric value representing the length at 50\% of the maximum catchability of the catchability curve.
#' @references Thorpe, R.B., Dolder, P.J. , Reeves, S., Robinson, P., Jennings, S. (2015). Assessing fishery and ecological consequences of alternative management options for multispecies fisheries \emph{ICES Journal of Marine Science}, 73(6):1503-1512.
"NS_L50"
