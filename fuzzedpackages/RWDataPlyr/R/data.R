#' Example rdf file with monthly data.
#'
#' An example of an rdf file that has already been read into R via 
#' [read.rdf()]. This example contains 39 slots, at the monthly 
#' timestep for 11 years and 25 runs. Slots include pool elevation, 
#' flow, and flags. Use this with [rdf_slot_names()] or 
#' [rdf_get_slot()] to use the data.
#'
#' @format A multi level list. `keyRdf$meta` provides a description
#' of the RiverWare run used to generate this data.
#' 
#' @source Bureau of Reclamation, 2016
"keyRdf"

#' Example rdf file with annual data.
#'
#' An example of an rdf file that has already been read into R via 
#' [read.rdf()]. This example contains 23 slots, at the annual 
#' timestep for 11 years and 25 runs. Slots only include  
#' flags. Use this with [rdf_slot_names()] or 
#' [rdf_get_slot()] to use the data.
#'
#' @format A multi level list. `sysRdf$meta` provides a description
#' of the RiverWare run used to generate this data.
#' 
#' @source Bureau of Reclamation, 2016
"sysRdf"

#' Example aggregated scenario data
#' 
#' An example of the tbl_df returned by [rw_scen_aggregate()] containing two 
#' scenarios of data. 
"scen_data"
