#' Create a `rwd_agg` template
#' 
#' `rwd_agg_template()` creates a template csv file to use to create a RiverWare 
#' data aggregator ([rwd_agg]). 
#' 
#' @param file The file name to use for the template
#' @param path The path to create the template at
#' @param examples Boolean; When `FALSE` (default), the template includes only 
#'   headers. When `TRUE`, the template will include several examples of 
#'   specifying how each slot should be summarized. 
#'   
#' @examples 
#' rwd_agg_template(file = "rwa_slots.csv", path = tempdir())
#' rwd_agg_template(file = "rwa_slots.csv", path = tempdir(), examples = TRUE)
#' 
#' @seealso [read_rwd_agg()]
#' 
#' @export

rwd_agg_template <- function(file, path = ".", examples = FALSE)
{
  if (length(file) != 1 || length(path) != 1 || length(examples) != 1)
    stop(
      "`rdf_agg_template()` expects all parameters to have length == 1", 
      call. = FALSE
    )
  
  if (!is.character(file) || !is.character(path))
    stop(
      "`rdf_agg_template()` expects `file` and `path` to be characers", 
      call. = FALSE
    )
  
  if (!is.logical(examples) || (is.logical(examples) && is.na(examples)))
    stop(
      "`rdf_agg_template()` expects `examples` to be TRUE or FALSE",
      call. = FALSE
    )
  
  if (tools::file_ext(file) != "csv")
    stop("`rdf_agg_template()` expects `file` to be a .csv file", call. = FALSE)
  
  if (!dir.exists(path)) 
    stop(path, " does not exist")
  
  if (examples) {
    x <- stats::setNames(data.frame(matrix( 
        c("KeySlots.rdf", "Mead.Pool Elevation", "asis", NA, NA, NA, "meadPe", 
        "KeySlots.rdf", "Powell.Outflow", "wy", "sum", "<", 8230000, "pwylt823",
        "KeySlots.rdf", "Powell.Outflow", "July", NA, NA, 0.001, "pjulrel",
        "KeySlots.rdf", "Powell.Outflow", "wy", "min", NA, NA, "pminrel",
        "KeySlots.rdf", "Powell.Outflow", "djf", "max", ">=", 600000, "pwinmax"
        ),
        nrow = 5, 
        byrow = TRUE
      )), 
      c("file", "slot", "period", "summary" , "eval", "t_s", "variable")
    )
    
  } else {
    x <- stats::setNames(
      data.frame(
        matrix(ncol = 7, nrow = 0),
        stringsAsFactors = FALSE
      ), 
      c("file", "slot", "period", "summary" , "eval", "t_s", "variable")
    )
  }
    
  utils::write.csv(x, file.path(path, file), row.names = FALSE)
}

#' Read in a rwd_agg file
#' 
#' `read_rwd_agg()` reads in a csv file and creates a [rwd_agg] object. 
#' Therefore, if the csv file is not properly formatted to contain the correct
#' information for a [rwd_agg] object, it will fail. [rwd_agg_template()] will
#' create a blank template file for the user to fill in, which has the correct 
#' headers.
#' 
#' @param file The csv file to be read in and converted
#' 
#' @examples
#' read_rwd_agg( 
#'   system.file(
#'     "extdata/rwd_agg_files/passing_aggs.csv", 
#'     package = "RWDataPlyr"
#'   )
#' )
#' 
#' @seealso [rwd_agg_template()]
#' 
#' @export

read_rwd_agg <- function(file)
{
  if (length(file) != 1 || 
      !is.character(file) || 
      tools::file_ext(file) != "csv") {
    stop(
      "`read_rwd_agg()` expects `file` to be a length 1, character, with a .csv extension", 
      call. = FALSE
    )
  }
  
  if (!file.exists(file))
    stop(file, " does not exist.", call. = FALSE)
  
  rwd_agg(utils::read.csv(file, stringsAsFactors = FALSE))
}
