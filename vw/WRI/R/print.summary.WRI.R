#' print the summary of WRI object
#' @param x a 'summary.WRI' object
#' @param ... further arguments passed to or from other methods.
#' @export
print.summary.WRI <- function(x, ...){

        cat('Call:\n')
        print(x$call)

        cat('\nPartial F test for individual effects:\n\n')
        printCoefmat(x$partial_F_table)

        cat(paste0('\nWasserstein R-squared: ', x$r.square, '\n'))

        cat('F-statistic (by Satterthwaite method): ')

        cat(paste0(x$global_wasserstein_F_stat, ' on ', x$global_wasserstein_F_df,
                   ' DF,', ' p-value: ', formatC(x$global_F_pvalue, format = "e", digits = 3), '\n'))
        invisible(x)
}
