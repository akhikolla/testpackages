#' Convert a table with host and symbiont associations to a matrix
#' @description Converts a table of associations to an association matrix with rows as symbionts and columns as host.
#'
#' @param assoc_table A dataframe with two columns
#'
#' @return A matrix with rows as symbionts and columns as hosts with 1's representing an association.
#'
#' @details Converts a dataframe with first column listing the host individually and the second column as the symbionts.
#' If hosts have more than one symbiont list these with commas.
#' For example, if the table is a tab-delimited file then a row should read:
#' "Hostus_mostus    Symbiont_1, Symbiont_2".
#'
#' @examples
#' file_path <- system.file("extdata",
#'                          "gopher_lice_mapping.txt",
#'                           package = "treeducken")
#' gopher_lice_map <- read.table(file_path,
#'                               stringsAsFactors = FALSE,
#'                               header = TRUE)
#' gopher_lice_assoc_matrix <- convert_assoc_table_to_matrix(gopher_lice_map)
#'
#' @export
convert_assoc_table_to_matrix <- function(assoc_table){
    if(dim(assoc_table)[2] != 2)
        stop("The 'assoc_table' should have two columns.")

    hosts <- unique(assoc_table[,1])
    symbs <- assoc_table[,2]
    num_hosts <- length(hosts)
    num_symbs <- length(symbs)
    strplit_and_unlist <- function(x) {unlist(strsplit(x, split=",", fixed=TRUE))}
    symbs_on_host <- sapply(symbs, strplit_and_unlist)
    names(symbs_on_host) <- NULL
    symbs <- unlist(symbs_on_host)

    assoc_mat <- matrix(0, nrow = length(symbs), ncol = length(hosts))
    rownames(assoc_mat) <- symbs
    colnames(assoc_mat) <- hosts
    for(h in 1:length(hosts)){
        assoc_mat[symbs_on_host[[h]],h] <- 1
    }
    assoc_mat
}