#' Summarize potential problems in a distance matrix
#'
#' For the inviduals represented in a distance matrix, collect the
#' self-self, best, and 2nd best distances, and summarize the results
#' in a data frame.
#'
#' @param d A distance or similarity matrix
#' @param dimension Whether to determine the best distances within rows or columns
#' @param get_min If TRUE, get the minimum (for a distance matrix);
#'     if FALSE, get the maximum (for a similarity matrix)
#' @param subset Whether to return just the rows with potential problems, or all of the rows.
#' @param threshold If `subset="problems"`, the threshold on the difference between the self and best distances.
#'
#' @return A data frame containing individual ID, distance to self,
#'     best distance and corresponding individual, 2nd best distance
#'     and the corresponding individual.
#'
#' @examples
#' # align rows in the provided dataset, lineup2ex
#' aligned <- align_matrix_rows(lineup2ex$gastroc, lineup2ex$islet)
#' # find correlated columns
#' selected_genes <- (corr_betw_matrices(aligned[[1]], aligned[[2]], "paired") > 0.75)
#' # calculate correlation between rows
#' similarity <- corr_betw_matrices(t(lineup2ex$gastroc[,selected_genes]),
#'                                  t(lineup2ex$islet[,selected_genes]), "all")
#' # pull out the problems, looking by row (where best > self + 0.3)
#' problems_byrow <- get_problems(similarity, get_min=FALSE, threshold=0.3)
#'
#' # pull out the problems, looking by column (where best > self + 0.3)
#' problems_bycol <- get_problems(similarity, get_min=FALSE, threshold=0.3,
#'                                dimension="column")
#'
#' @seealso [get_self()], [get_best()], [get_2ndbest()], [which_best()], [get_nonself()]
#'
#' @export
get_problems <-
    function(d, dimension=c("row", "column"), get_min=TRUE,
             subset=c("problems", "all"), threshold=0)

{
    dimension <- match.arg(dimension)

    self <- get_self(d)
    result <- data.frame(ind=names(self),
                         self=self,
                         best=get_best(d, dimension, get_min),
                         which_best=which_best(d, dimension, get_min),
                         next_best=get_2ndbest(d, dimension, get_min),
                         which_next_best=which_2ndbest(d, dimension, get_min),
                         stringsAsFactors=FALSE)
    if(get_min) result <- result[order(result$best - result$self),,drop=FALSE]
    else result <- result[order(result$self - result$best),,drop=FALSE]

    subset <- match.arg(subset)

    if(subset == "problems") {
        if(get_min) {
            result <- result[(is.na(result$self) & !is.na(result$best)) |
                             (!is.na(result$best) & !is.na(result$self) & result$best < result$self - threshold),]
        } else {
            result <- result[(is.na(result$self) & !is.na(result$best)) |
                             (!is.na(result$best) & !is.na(result$self) & result$best > result$self + threshold),]
        }
    }

    rownames(result) <- seq_len(nrow(result))
    result
}
