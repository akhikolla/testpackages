#' Create a vector of scenarios from different dimensions
#' 
#' `rw_scen_gen_names()` creates a vector of full scenario names by combining 
#' multiple dimensions together.
#' 
#' Many RiverWare runs are specified by multiple dimensions (or assumptions), 
#' and RiverSMART creates folder names by combining the dimension names together
#' for a full scenario name. `rw_scen_gen_names()` makes it quick to create all
#' of the full scenario names by passing in the names of the individual 
#' dimensions and creating all possible combinations of all dimensions. 
#' 
#' For example, the RiverWare run might consist of a supply dimension and a 
#' demand dimension, each consisting of two scenarios. This would result in four 
#' total scenarios. 
#' 
#' The function will work with two or more dimensions, as there is no need for 
#' this function if there is only one dimension.
#' 
#' @param dim1 A character vector with all of the first dimension's names.
#' 
#' @param dim2 A character vector with all of the second dimension's names. 
#' 
#' @param ... As many individual character vectors as necessary for the 
#'   remaining dimension's names.
#'   
#' @param sep The character used to separate the different dimension names. 
#'   Defaults to`","`.
#'   
#' @return A character vector of all possible combinations of the dimensions.
#'
#' @examples
#' rw_scen_gen_names("DNF", "CT", c("IG", "NA"), c("MTOM", "24-MS"))
#' rw_scen_gen_names("DNF", "CT", c("IG", "NA"), sep = "_")
#' 
#' @export

rw_scen_gen_names <- function(dim1, dim2, ..., sep = ",") 
{
  scens <- expand.grid(dim1, dim2, ...)
  scens <- apply(scens, 1, paste, collapse = sep)
  scens
}

#' @describeIn rw_scen_gen_names Deprecated version of `rw_scen_gen_names()`
#' @export
makeAllScenNames <- function(dim1, dim2, ..., sep = ",")
{
  .Deprecated("rw_scen_gen_names()")
  rw_scen_gen_names(dim1, dim2, ..., sep = ",")
}
