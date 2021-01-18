#' HFAII Queen Adjacency Matrix
#'
#' Binary adjacency matrix for the Humphrey Field Analyzer-II (Carl Zeiss Meditec Inc., Dublin, CA)
#'
#' @usage data(HFAII_Queen)
#'
#' @format This adjacency matrix is formated using queen neighbor criteria, meaning two locations on
#'  the visual field are only considered neighbors if they share an edge or corner. The adjacency matrix is a
#'  54 x 54 dimensional binary object with zeros on the diagonal and the column and row sums are equal
#'  to the number of neighbors.
#'
"HFAII_Queen"
