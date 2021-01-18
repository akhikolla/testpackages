## This file was originally part of SimInf, a framework for stochastic
## disease spread simulations.
##
## It has been modified by Trevelyan McKinley from the original "mparse.R" file
## in SimInf version 6.3.0.9000, by removing some of the functions in the file. 
## Any errors arising from these changes are therefore my (T. McKinley's) 
## responsibility. Modifications made on 12 August 2019.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2019 Stefan Widgren
##
## SimInf is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SimInf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Split the propensity in order to separate preprocessor and
## punctuator tokens from identifiers, for example:
##
## > tokens(" bR * R ")
## [1] "bR" "*"  "R"
tokens <- function(propensity) {
    ## List of valid preprocessor operator or punctuator tokens.
    operators <- c("<<=", ">>=", "!=", "%=", "##", "&&", "&=", "*=",
                   "++", "+=", "--", "-=", "->", "/=", "<<", "<=", "==",
                   ">=", ">>", "^=", "|=", "||", "!", "~", "%", "&", "(",
                   ")", "*", "+", ",", "-", "/", ":", ";", "<", "=",
                   ">", "?", "[", "]", "^", "{", "|", "}", "#")

    ## Create a matrix (1 x 2) of the propensity, where the first
    ## column is the token and the second column indicates if the
    ## token is one of the operators (indicated with 'op').
    propensity <- cbind(token = propensity, type = "")

    ## Iterate over each operator and try to split each row in the
    ## propensity in smaller pieces.
    for (op in operators) {
        propensity <- lapply(seq_len(nrow(propensity)), function(i) {
            x <- propensity[i, seq_len(ncol(propensity)), drop = FALSE]

            ## Is it a non-operator token that we could split?
            if (nchar(x[1, 2]) == 0) {
                m <- gregexpr(op, x[1, 1], fixed = TRUE)[[1]]
                if (m[1] != -1) {
                    ## The operator exists in the token. Split the
                    ## token in smaller pieces. The cut-points are
                    ## deterimined by the position and length of op
                    ## e.g. "A op B" -> "A", "op", "B".
                    x <- as.character(x[1, 1])
                    j <- 1
                    xx <- NULL
                    for (i in seq_len(length(m))) {
                        if (m[i] > j)
                            xx <- c(xx, substr(x, j, m[i] - 1))
                        j <- m[i] + attr(m, "match.length")[i]
                        xx <- c(xx, substr(x, m[i], j - 1))
                    }

                    ## Make sure last sub-string is copied.
                    if (j <= nchar(x))
                        xx <- c(xx, substr(x, j, nchar(x)))

                    ## Remove leading and trailing whitespace and drop
                    ## empty strings
                    xx <- gsub("(^\\s+)|(\\s+$)", "", xx)
                    xx <- xx[nchar(xx) > 0]

                    ## Create a 2-column matrix from all sub-strings
                    x <- cbind(token = xx, type = ifelse(xx == op, "op", ""))
                }
            }

            x
        })

        propensity <- do.call("rbind", propensity)
    }

    propensity[, 1]
}

## Rewrite propensity
##
## Rewrite the propensity by replacing all compartments by
## \code{u[compartments[j]]} where \code{j} is the numbering in
## compartments. On return, 'depends' contains all compartments upon
## which the propensity depends.
rewrite_propensity <- function(propensity, compartments, ldata_names,
                               gdata_names, v0_names) {
    propensity <- tokens(propensity)
    G_rowname <- paste0(propensity, collapse = "")
    depends <- integer(length(compartments))

    ## Find compartments in propensity
    i <- match(propensity, compartments)
    propensity <- ifelse(is.na(i), propensity, sprintf("u[%i]", i - 1L))
    i <- i[!is.na(i)]
    if (length(i))
        depends[i] <- 1

    ## Find ldata parameters in the propensity
    i <- match(propensity, ldata_names)
    propensity <- ifelse(is.na(i), propensity, sprintf("ldata[%i]", i - 1L))

    ## Find gdata parameters in the propensity
    i <- match(propensity, gdata_names)
    propensity <- ifelse(is.na(i), propensity, sprintf("gdata[%i]", i - 1L))

    ## Find v0 parameters in the propensity
    i <- match(propensity, v0_names)
    propensity <- ifelse(is.na(i), propensity, sprintf("v[%i]", i - 1L))

    list(propensity = paste0(propensity, collapse = ""),
         depends    = depends,
         G_rowname  = G_rowname)
}

## Generate the 'from' or 'dest' labels in the G rownames.
G_label <- function(x) {
  if (length(x) == 0)
    return("@")
  
  ## Prefix compartments if more than one unit, e.g., '2*S'.
  lbl <- ifelse(abs(x) > 1, paste0(abs(x), "*"), "")
  lbl <- paste0(lbl, names(x))
  
  ## Combine all compartments, e.g., 'S + I'
  paste0(lbl, collapse = " + ")
}

parse_compartments <- function(x, compartments) {
    ## Split into 'compartment1 + compartment2 + ..'
    x <- unlist(strsplit(x, "+", fixed = TRUE))

    ## Remove spaces.
    x <- gsub(" ", "", x)

    ## Replace 'n*compartment' with n replicates of 'compartment'
    x <- unlist(sapply(x, function(xx) {
        m <- regexpr("^[[:digit:]]+[*]", xx)
        if (m != 1)
            return(xx)

        ## Determine number of replicates and remove 'n*'
        n <- regmatches(xx, m)
        xx <- sub(n, "", xx, fixed = TRUE)
        n <- as.integer(substr(n, 1, nchar(n) - 1))

        rep(xx, n)
    }))

    ## Check for valid usage of the empty set.
    if (any(x == "@") && length(x) > 1)
        stop("Invalid usage of the empty set '@'.")
    x <- x[x != "@"]

    ## Assign each compartment into its number according to the
    ## ordering in compartments
    i <- match(x, compartments)
    if (any(is.na(i)))
        stop(sprintf("Unknown compartment: '%s'.", x[is.na(i)]), call. = FALSE)

    tabulate(i, length(compartments))
}

parse_transitions <- function(transitions, compartments, ldata_names,
                              gdata_names, v0_names) {
    lapply(strsplit(transitions, "->", fixed = TRUE), function(x) {
        if (length(x) < 3) {
            stop("Invalid transition: '",
                 paste0(x, collapse = "->"),
                 "'.",
                 call. = FALSE)
        }

        ## Remove spaces
        propensity <- gsub(" ", "", x[c(-1, -length(x))])
        propensity <- paste0(propensity, collapse = "->")

        ## Determine the corresponding column in the state change
        ## vector S.
        from <- parse_compartments(x[1], compartments)
        dest <- parse_compartments(x[length(x)], compartments)
        S <- dest - from

        propensity <- rewrite_propensity(propensity, compartments,
                                         ldata_names, gdata_names,
                                         v0_names)

        ## Determine the G rowname
        names(from) <- compartments
        names(dest) <- compartments
        from <- G_label(from[which(from > 0)])
        dest <- G_label(dest[which(dest > 0)])
        G_rowname <- paste(from, "->", propensity$G_rowname, "->", dest)

        list(propensity = propensity$propensity,
             depends    = propensity$depends,
             S          = S,
             G_rowname  = G_rowname)
    })
}
