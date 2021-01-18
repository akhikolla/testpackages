#' Summary of a count variable
#'
#' Summary of a count variable.
#'
#' The function does a similar job to \code{table()} with more flexibility
#' introduced by the argument \code{breaks}. The user can decide how to break
#' the count values and decide to merge some cells if needed.
#'
#' @param count integer, observed count value for every individual in the
#'     sample.
#' @param breaks integer, values at which the breaks should happen. The function
#'     will compute the observed frequency in \code{[breaks[i], breaks[i + 1])}.
#' @param formatChar logical, should the values be converted to character and
#'     formatted?
#' @return \code{matrix} with 2 rows and \code{length(breaks)} columns.  The
#'     column names are the cells names. The rows are the observed frequencies
#'     and relative frequencies (probabilities).
#' @export
count_table <- function(count, breaks, formatChar = FALSE) {
    if (missing(breaks))
        breaks <- sort(unique(count))

    cl <- character()
    nb <- numeric()

    for (i in 1:(length(breaks) - 1)) {
        cl <- c(cl, .getNames_(breaks[i], breaks[i + 1]))
        nb <- c(nb, .getValues_(breaks[i], breaks[i + 1], count))
    }

    ## final count
    cl <- c(cl, paste(">=", breaks[length(breaks)]))
    nb <- c(nb,  sum(count >= breaks[length(breaks)]))

    mat <- matrix(nrow = 2, ncol = length(cl))
    mat[1, ] <- nb
    mat[2, ] <- nb / sum(nb)
    colnames(mat) <- cl
    rownames(mat) <- c("Frequency", "Relative frequency")

    if (formatChar) {
        mat[1, ]  <- formatC(mat[1, ])
        mat[2, ]  <- trimws(formatC(as.numeric(mat[2, ]), digits = 2))
    }

    mat
}

.getNames_ <- function(start, end) {
    if (end - start == 1)
        return(as.character(start))
    else
        return(as.character(paste(start, end - 1, sep = "-")))
}

.getValues_ <- function(start, end, count) {
    sum(count >= start & count < end)
}

.getValues2_ <- function(start, end, count, values) {
    sum(values[which(count >= start & count < end)])
}

.getValuesMat_ <- function(start, end, count, values) {
    if (missing(end)) {
        ## val <- values[, which(count < start), drop = FALSE]
        ## return(1 - rowSums(val))
        val <- values[, which(count >= start), drop = FALSE]
        return(rowSums(val))
    } else {
        val <- values[, which(count >= start & count < end), drop = FALSE]
        return(rowSums(val))
    }
}

.getValues3_ <- function(start, count, values) {
    sum(values[which(count >= start)])
}

.computeCountRange <- function(count, count_range) {
    for (i in 1:length(count_range)) {
        ri <- count_range[i]

        if (grepl("-", ri)) {
            tmp <- unlist(strsplit(ri, "-"))
            start <- as.numeric(tmp[1])
            end <- as.numeric(tmp[2])
            ans <- count >= start & count <= end
        } else if (grepl(">=", ri)) {
            tmp <- unlist(strsplit(ri, ">="))
            ans <- count >= as.numeric(trimws(tmp[2]))
        } else
            ans <- count == as.numeric(ri)

        if (ans) return(i)
    }
}

.adjust_breaks <- function(breaks, res, pij) {
     counts <- character()
     actual <- numeric()
     model <- numeric()
     p <- matrix(nrow = nrow(pij), ncol = length(breaks))

     for (i in 1:(length(breaks) - 1)) {
         counts <- c(counts, .getNames_(breaks[i], breaks[i + 1]))
         actual <- c(actual, .getValues2_(breaks[i], breaks[i + 1],
                                          res$Counts, res$Actual)
                     )

         model <- c(model,
                    .getValues2_(breaks[i], breaks[i + 1],
                                 res$Counts, res$Predicted)
                    )

         p[, i] <- .getValuesMat_(breaks[i], breaks[i + 1],
                                  res$Counts, pij)
     }

     ## produce the final counts
     i <- length(breaks)
     counts <- c(counts, paste(">=", breaks[i]))
     actual <- c(actual, .getValues3_(breaks[i], res$Counts, res$Actual))

     p[, i] <- .getValuesMat_(breaks[i],, res$Counts, pij)

     model <- c(model, .getValues3_(breaks[i], res$Counts, res$Predicted))

     res <- data.frame(Counts = counts, Actual = actual, Predicted = model)
     colnames(p) <- counts
     attr(p, "details") <- res

     p
}

.get_dij_mat <- function(lhs_dat, count_range_) {
    cl <- sapply(lhs_dat, .computeCountRange,
                 count_range = count_range_
                 )

    .ff <- function(ind) {
        res <- rep(0, length(count_range_))
        res[ind] <- 1
        res
    }

    t(sapply(cl, .ff))
}

.run_chisq_reg <- function(si, dij_pij, cols_) {
    ## build the regression matrix: 1_i ~ s_i + dij_pij, j =1,'' J-1
    mat <- data.frame(y = 1, si, dij_pij[, -ncol(dij_pij)])
    mod <- mod <- lm(y~. -1, data = mat)
    stat <- nrow(mat) * summary(mod)$r.squared

    ## format the result
    rval <- matrix(NA, nrow = 1, ncol = 3)
    colnames(rval) <- c("DF", "Chisq", "Pr(>Chisq)")

    rval[, 1] <- length(cols_) - 1 ## number of cells - 1
    rval[, 2] <- stat
    rval[, 3] <- pchisq(stat, rval[, 1], lower.tail = FALSE)

    title <- "chi-square goodness-of-fit test\n"
    topnote <- paste("Cells considered ", paste(cols_, collapse = " "),
    sep = "", collapse = "\n")

    structure(as.data.frame(rval), heading = c(title, topnote),
              class = c("anova", "data.frame"))
}
