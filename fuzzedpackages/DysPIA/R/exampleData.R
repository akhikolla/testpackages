#' Example matrix of gene expression value.
#'
#' A dataset of transcriptional profiles from p53+ and p53 mutant cancer cell lines. 
#' It includes the normalized gene expression for 6385 genes in 50 samples.
#' Rownames are genes, columnnames are samples.
#' @docType data
#' @name gene_expression_p53
#' @usage data(gene_expression_p53)
NULL

#' Example vector of category labels. 
#'
#' The labels for the 50 cell lines in p53 data. Control group's label is 'WT', case group's label is 'MUT'.
#' @docType data
#' @name class.labels_p53
#' @usage data(class.labels_p53)
NULL

#' Example vector of DysGPS in p53 data. 
#'
#' The score vector of 164923 gene pairs from p53 dataset. 
#' It can be loaded from the example datasets of R-package 'DysPIA', 
#' and also can be obtained by running DysGPS(), details see DysGPS.R
#' @docType data
#' @name DysGPS_p53
#' @usage data(DysGPS_p53)
NULL

#' Example list of DysPIA result in p53 data. 
#'
#' The list includes 81 pathway results from 'DisPIA.R' as an example used in 'DyspiaSig.R'. 
#' @docType data
#' @name DyspiaRes_p53
#' @usage data(DyspiaRes_p53)
NULL



#' Example list of gene pair background. 
#'
#' The list of background was used in ''DysGPS.R' and 'calEdgeCorScore_ESEA.R' which is a part of the 'combined_background' in 'DysPIAData'.
#' @name sample_background
#' @usage data(sample_background)
NULL
