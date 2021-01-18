# Project documentation and imports

#' The SCOPer package
#'
#' \code{scoper} is a member of the Immcantation framework and provides computational approaches 
#' for the identification of B cell clones from adaptive immune receptor repertoire sequencing 
#' (AIRR-Seq) datasets. It includes methods for assigning clonal identifiers using
#' sequence identity, hierarchical clustering, and spectral clustering.
#'
#' @section Clonal clustering:
#'
#' \itemize{
#'   \item  \link{identicalClones}:  Clonal assignment using sequence identity partitioning.
#'   \item  \link{hierarchicalClones}:  Hierarchical clustering approach to clonal assignment.
#'   \item  \link{spectralClones}:  Spectral clustering approach to clonal assignment.
#' }
#' 
#' @section Visualization:
#' \itemize{
#'   \item  \link{plotCloneSummary}:  Visualize inter- and intra-clone distances.
#' }
#' 
#' @name        scoper
#' @docType     package
#' @references
#' \enumerate{
#'   \item  Nouri N and Kleinstein SH (2018). A spectral clustering-based method for identifying clones
#'   from high-throughput B cell repertoire sequencing data. Bioinformatics, 34(13):i341-i349.
#'   \item  Nouri N and Kleinstein SH (2019). Somatic hypermutation analysis for improved identification 
#'   of B cell clonal families from next-generation sequencing data. bioRxiv, 10.1101/788620.
#'   \item Gupta NT, et al. (2017). Hierarchical clustering can identify B cell clones with high 
#'   confidence in Ig repertoire sequencing data. The Journal of Immunology, 198(6):2489-2499.
#' }
#'
#' @import      ggplot2
#' @import      graphics
#' @import      methods
#' @import      utils
#' @importFrom  alakazam    pairwiseDist checkColumns getDNAMatrix getAAMatrix
#'                          padSeqEnds progressBar groupGenes baseTheme translateDNA
#' @importFrom  data.table  as.data.table .I
#' @importFrom  doParallel  registerDoParallel
#' @importFrom  dplyr       n %>%
#'                          filter select arrange bind_rows
#'                          group_by ungroup group_indices
#'                          mutate summarize slice 
#' @importFrom  foreach     foreach %dopar% registerDoSEQ
#' @importFrom  rlang       sym syms
#' @importFrom  Rcpp        evalCpp
#' @importFrom  scales      pretty_breaks
#' @importFrom  shazam      consensusSequence
#' @importFrom  stats       density kmeans sd cor
#'                          as.dist hclust cutree
#' @importFrom  stringi     stri_split_fixed stri_length stri_count
#' @importFrom  tidyr       separate
#' @useDynLib   scoper, .registration=TRUE
NULL

# Package loading actions
.onAttach <- function(libname, pkgname) {
    msg <- paste("As of v1.0.0 the AIRR Rearrangement schema is now the default file format.",
                 "A description of the standard is available at https://docs.airr-community.org.",
                 "The legacy Change-O format is supported through arguments to each function",
                 "that allow the input column names to be explicitly defined.",
                 sep="\n")
    packageStartupMessage(msg)
}
