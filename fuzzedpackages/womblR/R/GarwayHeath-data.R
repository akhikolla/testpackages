#' Garway-Heath angles for the HFA-II
#'
#' These Garway-Heath angles are used as the dissimilarity metric when implementing the boundary
#' detection model for a longitudinal series of visual fields.
#'
#' @usage data(GarwayHeath)
#'
#' @format A vector with length 54, where each entry represents the angle (in degrees) that the
#'  underlying retinal nerve fiber enters the optic nerve head. The measure ranges from 0-360,
#'  where 0  is designated at the 9-o’clock position (right eye) and angles are counted counter
#'  clockwise. These angles are estimates for the Humphrey Field Analyzer-II (Carl Zeiss Meditec
#'  Inc., Dublin, CA). The 26th and 35th entries are missing as they correspond to a natural
#'  blind spot.
#'
#' @references Garway-Heath, et al. (2000). Ophthalmology 107:10:1809–1815.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/11013178}{PubMed})
"GarwayHeath"
