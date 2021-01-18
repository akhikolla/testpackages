#' Compute Selected Evaluation Metrics
#' 
#' @description
#' Compute selected evaluation metrics for binary (i.e. two-class) confusion 
#' matrices.
#' 
#' @param mat Binary confusion \code{matrix} (2-by-2; see Examples).
#' @param type Target evaluation metric as \code{character}, defaults to 
#' \code{"accuracy"}. Other available options are \code{"precision"} and 
#' \code{"recall"}.
#' 
#' @return A single \code{numeric}.
#' 
#' @author Florian Detsch
#' 
#' @references 
#' University of Michigan (2017) Applied Machine Learning in Python. Available 
#' online: \url{https://www.coursera.org/learn/python-machine-learning/home/welcome}.
#'
#' @examples 
#' in1 = matrix(c(96, 4, 8, 19), nc = 2L, byrow = TRUE)
#' rownames(in1) = c("Condition Positive", "Condition Negative")
#' colnames(in1) = c("Predicted Positive", "Predicted Negative")
#' 
#' evalMetrics(in1) # default: "accuracy"
#' evalMetrics(in1, "precision")
#' evalMetrics(in1, "recall")
#' 
#' in2 = matrix(c(26, 17, 7, 400), nc = 2, byrow = TRUE)
#' evalMetrics(in2, "precision")
#' evalMetrics(in2, "recall")
#' 
#' @export evalMetrics
#' @name evalMetrics
evalMetrics = function(mat, type = c("accuracy", "precision", "recall")) {
  
  if (!all(dim(mat) == 2))
    stop("More than 2 dimensions not implemented, yet.\n")
  
  tp = mat[1, 1]; tn = mat[2, 2]
  fp = mat[2, 1]; fn = mat[1, 2]
  
  out = if (type[1] == "accuracy") {
    sum(tp, tn) / (sum(tp, tn, fp, fn))
  } else if (type[1] == "precision") {
    tp / sum(tp, fp)
  } else if (type[1] == "recall") {
    tp / sum(tp + fn)
  }
  
  return(out)
}
