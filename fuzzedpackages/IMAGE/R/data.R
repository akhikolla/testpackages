
#' Example data set to test the IMAGE method
#'
#' A list containing the data needed to use the \code{image} function.
#'
#'
#' @format A \code{list} three elements:
#' \describe{
#'   \item{data}{A \code{data.frame} with columns r, y, r1, r2, y1, y2}
#'   \item{geno}{A \code{data.frame} with columns hap1 and hap2}
#'   \item{K}{A 100X100 matrix}
#' }
"ExampleData"


#' Example results using the example data set
#'
#' A data set containing the output of \code{image} function call on the
#' included data set \code{ExampleData}.
#'
#' @format A \code{data.frame} with 1000 rows and 9 variables (columns):
#' \describe{
#'   \item{Loc}{ordinal number of SNP-CpG pair being analyzed}
#'   \item{numIDV}{number of observations of SNP-CpG pair being analyzed}
#'   \item{beta}{the fixed effect parameter estimate for the predictor of interest}
#'   \item{se_beta}{the standard deviation of fixed effect}
#'   \item{pvalue}{P-value for the fixed effect, based on the wald test}
#'   \item{h2}{heritability of the transformed rate}
#'   \item{sigma}{total variance component}
#'   \item{converged}{a logical indicator for convergence}
#'   \item{time}{time to converge}
#' }
"example_results"
