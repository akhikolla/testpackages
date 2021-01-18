#' @title DyspiaSig
#' @description Returns the significant summary of DysPIA results.
#'
#' @param DyspiaRes Table with results of running DysPIA().
#' @param fdr Significant threshold of `padj` (a BH-adjusted p-value).
#' @return A list of significant DysPIA results, including correlation gain and correlation loss.
#' @examples
#' data(pathway_list,package="DysPIAData")
#' data(DyspiaRes_p53)
#' summary_p53 <- DyspiaSig(DyspiaRes_p53, 0.05)       # filter with padj<0.05
#' 
DyspiaSig <- function(DyspiaRes, fdr){
  DyspiaRes_gain <- DyspiaRes[DyspiaRes$DysPS >= 0 & DyspiaRes$padj < fdr, ]
  DyspiaRes_loss <- DyspiaRes[DyspiaRes$DysPS < 0 & DyspiaRes$padj < fdr, ]
  
  results <- list(DyspiaRes_gain, DyspiaRes_loss)
  names(results) <- c("DyspiaRes_gain", "DyspiaRes_loss")
  
  return(results)
}
