#' Visual field series for one patient.
#'
#' A dataset containing 9 visual field series from a patient of the Vein Pulsation Study
#' Trial in Glaucoma and the Lions Eye Institute trial registry, Perth, Western Australia.
#'
#' @usage data(VFSeries)
#'
#' @format A data frame with 486 rows and 4 variables:
#' \describe{
#'   \item{Visit}{The visual field test visit number, (1, 2, ... , 9).}
#'   \item{DLS}{The observed outcome variable, differential light sensitivity (DLS).}
#'   \item{Time}{The time of the visual field test (in days from baseline).}
#'   \item{Location}{The location on the visual field of a Humphrey Field Analyzer-II
#'  (Carl Zeiss Meditec Inc., Dublin, CA) (1, 2, ... , 54).}
#' }
#' @source \url{https://anzctr.org.au/Trial/Registration/TrialReview.aspx?ACTRN=12608000274370}
"VFSeries"
