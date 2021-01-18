################################################################################
### File: data.R
### Description: data object
###
################################################################################

#' Artifical data of 54 subjects
#'
#' An artificial dataset containing data of 54 subjects where where a substance was administered in three different concentrations (1,2 and 3).
#' This data set can be used to show the paradoxical results obtained from rank tests, i.e., the Hettmansperger-Norton test.
#'
#' The columns are as follows:
#' \itemize{
#'   \item conc. Grouping variable specifying which concentration was used. This factor is ordered, i.e., 1 < 2 < 3.
#'   \item score. The response variable.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ParadoxicalRanks
#' @usage data(ParadoxicalRanks)
#' @format A data frame with 54 rows and 2 variables.
#' @references Happ M, Zimmermann G, Brunner E, Bathke AC (2020). Pseudo-Ranks: How to Calculate Them Efficiently in R. Journal of Statistical Software, Code Snippets, *95*(1), 1-22. doi: 10.18637/jss.v095.c01 (URL:https://doi.org/10.18637/jss.v095.c01).
#' @example R/example_paradoxical_results.txt
"ParadoxicalRanks"