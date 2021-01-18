#### Data ####

#' Example database
#'
#' A small example database subset from Laserson and Vigneault et al, 2014.
#'
#' @format   A data.frame with the following columns:
#'   \itemize{
#'     \item  \code{sequence_id}:                Sequence identifier
#'     \item  \code{sequence_alignment}:         IMGT-gapped observed sequence.
#'     \item  \code{germline_alignment}:         IMGT-gapped germline sequence.
#'     \item  \code{germline_alignment_d_mask}:  IMGT-gapped germline sequence with N, P and
#'                                               D regions masked.
#'     \item  \code{v_call}:                     V region allele assignments.
#'     \item  \code{v_call_genotyped}:           TIgGER corrected V region allele assignment.
#'     \item  \code{d_call}:                     D region allele assignments.
#'     \item  \code{j_call}:                     J region allele assignments.
#'     \item  \code{junction}:                   Junction region sequence.
#'     \item  \code{junction_length}:            Length of the junction region in nucleotides.
#' }
#'
#'
#' @references
#' \enumerate{
#'   \item  Laserson U and Vigneault F, et al. High-resolution antibody dynamics of
#'            vaccine-induced immune responses.
#'            Proc Natl Acad Sci USA. 2014 111:4928-33.
#' }
"ExampleDb"