print.summary.etm <- function(x, ...) {
    
    if (!inherits(x, "summary.etm"))
        stop("'x' must be of class 'summary.etm'")

    ## Find out if we have strata
    if (is.data.frame(x[[1]])) {
        ns <- 1
    } else {
        ns <- length(x)
    }

    if (ns == 1) {
    
        time <- x[[1]]$time
        qtime <- quantile(time, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
        ind <- findInterval(qtime, time)
    
        for (i in seq_along(x)) {
            cat(paste("Transition", names(x)[i], "\n", sep = " "))
            print(x[[i]][ind, ], row.names = FALSE)
            cat("\n")
        }

    } else {

        nn <- names(x)
        for (i in seq_len(ns)) {
            cat(nn[i], "\n\n")
            
            print(x[[i]])
        }
    }
        
    invisible()
}
