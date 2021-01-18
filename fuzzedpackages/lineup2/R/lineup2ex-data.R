#' Example dataset for lineup2 package
#'
#' Example dataset for lineup2 package, with gene expression data
#' for a selected set of 200 genes on two tissues on a set of about 500 mice,
#' with 100 genes chosen to be highly correlated between the two tissues
#' and 100 chosen at random.
#'
#' @docType data
#'
#' @usage data(lineup2ex)
#'
#' @format
#' List of two matrices, with gene expression data for gastrocnemius
#' muscle (`gastroc`) and pancreatic islets (`islet`), at a selected
#' set of 200 genes (100 are highly correlated between the two
#' tissues, and 100 others chosen at random). The matrices have
#' samples as rows and genes as columns. The row names are sample
#' identifiers. There are 498 samples for gastroc and 499 samples for
#' islet, with 497 samples in common.
#'
#' @keywords datasets
#'
#' @source <https://phenome.jax.org/projects/Attie1>
#'
#' @references
#' Broman KW, Keller MP, Broman AT, Kendziorski C, Yandell BS, Sen Ś,
#' Attie AD (2015) Identification and correction of sample mix-ups in
#' expression genetic data: A case study. G3 5:2177--2186
#'
#' Tian J, Keller MP, Oler AT, Rabaglia ME, Schueler KL, Stapleton DS,
#' Broman AT, Zhao W, Kendziorski C, Yandell BS, Hagenbuch B, Broman
#' KW, Attie AD (2015) Identification of the bile acid transporter
#' Slco1a6 as a candidate gene that broadly affects gene expression in
#' mouse pancreatic islets. Genetics 201:1253--1262
#'
#' @examples
#' data(lineup2ex)
#' common_ind <- align_matrix_rows(lineup2ex$gastroc, lineup2ex$islet)
"lineup2ex"
