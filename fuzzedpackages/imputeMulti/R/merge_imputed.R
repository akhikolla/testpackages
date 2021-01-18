
#' @title Merge imputed data and original dataset
#' @description Merge the imputed dataset from an  \code{imputeMulti} object with the original dataset. 
#' Merging is done by rownames, since imputeMulti maintains row-order during imputation.
#' @param impute_obj An object of class "imputeMulti".
#' @param y The dataset from which the missing data was imputed.
#' @param ... Arguments to be passed to other methods
#' @export
merge_imputed <- function(impute_obj, y, ...) {
  if (!is.imputeMulti(impute_obj)) stop("impute_obj must have class imputeMulti.")
  if (nrow(y) != nrow(impute_obj@data$imputed_data)) 
    stop("target and y do not have equal row counts.")
  
  x <- impute_obj@data$imputed_data
  x$rownames <- rownames(x)
  y$rownames <- rownames(y)
  return(merge(x, y, by= "rownames", ...))
}