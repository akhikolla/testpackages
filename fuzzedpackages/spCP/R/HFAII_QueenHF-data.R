#' HFAII Queen Hemisphere Adjacency Matrix
#'
#' Binary adjacency matrix for the Humphrey Field Analyzer-II (Carl Zeiss Meditec Inc., Dublin, CA)
#'
#' @usage data(HFAII_QueenHF)
#'
#' @format This adjacency matrix is formatted using queen neighbor criteria, meaning two locations on
#'  the visual field are only considered neighbors if they share an edge or corner. An additional criterion
#'  is included that locations are not considered neighbors if they fall within different hemispheres on the
#'  visual field. The adjacency matrix is a 54 x 54 dimensional binary object with zeros on the diagonal and
#'  the column and row sums are equal to the number of neighbors.
"HFAII_QueenHF"
