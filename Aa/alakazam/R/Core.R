# Common input/output and data structure manipulation functions for Alakazam

#### File I/O functions ####

#' Read a Change-O tab-delimited database file
#' 
#' \code{readChangeoDb} reads a tab-delimited database file created by a Change-O tool 
#' into a data.frame.
#'
#' @param    file       tab-delimited database file output by a Change-O tool.
#' @param    select     columns to select from database file.
#' @param    drop       columns to drop from database file.
#' @param    seq_upper  if \code{TRUE} convert sequence columns to upper case;
#'                      if \code{FALSE} do not alter sequence columns. See Value 
#'                      for a list of which columns are effected.

#' 
#' @return   A data.frame of the database file. Columns will be imported as is, except for 
#'           the following columns which will be explicitly converted into character 
#'           values:
#'           \itemize{
#'             \item  SEQUENCE_ID
#'             \item  CLONE
#'             \item  SAMPLE
#'           }
#'           And the following sequence columns which will converted to upper case if
#'           \code{seq_upper=TRUE} (default).
#'           \itemize{
#'             \item  SEQUENCE_INPUT
#'             \item  SEQUENCE_VDJ
#'             \item  SEQUENCE_IMGT
#'             \item  JUNCTION
#'             \item  GERMLINE_IMGT
#'             \item  GERMLINE_IMGT_D_MASK
#'           }
#'                   
#' @seealso  Wraps \link[readr]{read_delim}. 
#'           See \link{writeChangeoDb} for writing to Change-O files.
#' 
#' @examples
#' \dontrun{
#'     # Read all columns in and convert sequence fields to upper case
#'     db <- readChangeoDb("changeo.tsv")
#' 
#'     # Subset columns and convert sequence fields to upper case
#'     db <- readChangeoDb("changeo.tsv", select=c("SEQUENCE_ID", "SEQUENCE_IMGT"))
#' 
#'     # Drop columns and do not alter sequence field case
#'     db <- readChangeoDb("changeo.tsv", drop=c("D_CALL", "DUPCOUNT"), 
#'                         seq_upper=FALSE)
#' }
#' 
#' @export
readChangeoDb <- function(file, select=NULL, drop=NULL, seq_upper=TRUE) {
    # Define column data types
    seq_columns <- c("SEQUENCE_INPUT", "SEQUENCE_VDJ", "SEQUENCE_IMGT", 
                     "JUNCTION", "JUNCTION_AA", "CDR3_IGBLAST_NT", "CDR3_IGBLAST_AA",
                     "GERMLINE_VDJ", "GERMLINE_VDJ_V_REGION", "GERMLINE_VDJ_D_MASK",
                     "GERMLINE_IMGT", "GERMLINE_IMGT_V_REGION", "GERMLINE_IMGT_D_MASK",
                     "FWR1_IMGT", "FWR2_IMGT", "FWR3_IMGT", "FWR4_IMGT",
                     "CDR1_IMGT", "CDR2_IMGT", "CDR3_IMGT")

    # Define types
    header <- names(suppressMessages(readr::read_tsv(file, n_max=1)))
    types <- do.call(readr::cols, CHANGEO[intersect(names(CHANGEO), header)])
    
    # Read file    
    db <- suppressMessages(readr::read_tsv(file, col_types=types, na=c("", "NA", "None")))

    # Select columns
    select_columns <- colnames(db)
    if(!is.null(select)) { select_columns <- intersect(select_columns, select) }
    if(!is.null(drop)) { select_columns <- setdiff(select_columns, drop) }
    db <- select(db, select_columns)
    
    # Convert sequence fields to upper case
    upper_cols <- intersect(seq_columns, names(db))
    if (seq_upper & length(upper_cols) > 0) {
        db <- mutate_at(db, upper_cols, toupper)
    }
    
    return(db)
}


#' Write a Change-O tab-delimited database file
#' 
#' \code{writeChangeoDb} is a simple wrapper around \link[readr]{write_delim} with defaults 
#' appropriate for writing a Change-O tab-delimited database file from a data.frame.
#'
#' @param    data  data.frame of Change-O data.
#' @param    file  output file name.
#' 
#' @return   NULL
#' 
#' @seealso  Wraps \link[readr]{write_delim}. See \link{readChangeoDb} for reading to Change-O files.
#' 
#' @examples
#' \dontrun{
#'   # Write a database
#'   writeChangeoDb(data, "changeo.tsv")
#' }
#' 
#' @export
writeChangeoDb <- function(data, file) {
    write_tsv(data, file, na="NA")
}


#' Create a temporary folder
#'
#' \code{makeTempDir} creates a randomly named temporary folder in the 
#' system temp location.
#' 
#' @param    prefix  prefix name for the folder.
#' 
#' @return   The path to the temporary folder.
#' 
#' @seealso  This is just a wrapper for \link{tempfile} and 
#'           \link{dir.create}.
#' 
#' @examples
#' makeTempDir("Clone50")
#' 
#' @export
makeTempDir <- function(prefix) {
    temp_path <- tempfile(paste0(prefix, "-temp-"))
    dir.create(temp_path)
    
    return(temp_path)
}


#### Data manipulation functions ####

#' Translate a vector of strings
#' 
#' \code{translateStrings} modifies a character vector by substituting one or more 
#' strings with a replacement string.
#'
#' @param    strings      vector of character strings to modify.
#' @param    translation  named character vector or a list of character vectors specifying 
#'                        the strings to replace (values) and their replacements (names).
#' 
#' @return   A modified \code{strings} vector.
#' 
#' @details
#' Does not perform partial replacements. Each translation value must match a complete 
#' \code{strings} value or it will not be replaced.  Values that do not have a replacement
#' named in the \code{translation} parameter will not be modified.
#' 
#' Replacement is accomplished using \link{gsub}.
#' 
#' @seealso  See \link{gsub} for single value replacement in the base package.
#' 
#' @examples
#' # Using a vector translation
#' strings <- LETTERS[1:5]
#' translation <- c("POSITION1"="A", "POSITION5"="E")
#' translateStrings(strings, translation)
#' 
#' # Using a list translation
#' strings <- LETTERS[1:5]
#' translation <- list("1-3"=c("A","B","C"), "4-5"=c("D","E"))
#' translateStrings(strings, translation)
#' 
#' @export
translateStrings <- function(strings, translation) {
    # TODO:  use match instead for complete matching?  Currently regex characters in values will mess up the matching.
    for (n in names(translation)) {
        rep_regex <- paste(translation[[n]], collapse='|')
        strings <- gsub(paste0("^(", rep_regex, ")$"), n, strings)
    }
    
    return(strings)
}


#' Check data.frame for valid columns and issue message if invalid
#'
#' @param   data     data.frame to check.
#' @param   columns  vector of column names to check.
#' @param   logic    one of \code{"all"} or \code{"any"} controlling whether all,
#'                   or at least one, of the columns must be valid, respectively.
#' @return  \code{TRUE} if columns are valid and a string message if not.
# 
#' @examples
#' df <- data.frame(A=1:3, B=4:6, C=rep(NA, 3))
#' checkColumns(df, c("A", "B"), logic="all")
#' checkColumns(df, c("A", "B"), logic="any")
#' checkColumns(df, c("A", "C"), logic="all")
#' checkColumns(df, c("A", "C"), logic="any")
#' checkColumns(df, c("A", "D"), logic="all")
#' checkColumns(df, c("A", "D"), logic="any")
#' 
#' @export
checkColumns <- function(data, columns, logic=c("all", "any")) {
    ## DEBUG
    # data=df; columns=c("A", "D"); logic="any"
    
    # Check arguments
    logic <- match.arg(logic)
    
    data_names <- names(data)
    if (logic == "all") {
        # Check that all columns exist
        for (f in columns) {
            if (!(f %in% data_names)) { 
                msg <- paste("The column", f, "was not found") 
                return(msg)
            }
        }        
        # Check that all values are not NA
        for (f in columns) {
            if (all(is.na(data[[f]]))) { 
                msg <- paste("The column", f, "contains no data") 
                return(msg)
            }
        }
    } else if (logic == "any") {
        # Check that columns exist
        if (!any(columns %in% data_names)) {
            msg <- paste("Input must contain at least one of the columns:", 
                         paste(columns, collapse=", "))
            return(msg)
        }
        # Check that all values are not NA
        
        columns_found <- columns[columns %in% data_names]
        invalid <- sapply(columns_found, function(f) all(is.na(data[[f]])))
        if (all(invalid)) { 
            msg <- paste("None of the columns", paste(columns_found, collapse=", "), 
                         "contain data") 
            return(msg)
        }
    }
    
    # Return TRUE if all checks pass
    return(TRUE)
}

#### Plotting and progress functions ####

#' Standard progress bar
#' 
#' \code{progressBar} defines a common progress bar format.
#' 
#' @param   n  maximum number of ticks
#' 
#' @return  A \link[progress]{progress_bar} object.
#'
#' @export
progressBar <- function(n) {
    pb <- progress::progress_bar$new(format="  PROGRESS> [:bar] :percent :elapsed",
                                     width=40, clear=FALSE, stream=stdout(), force=TRUE,
                                     total=n)
    return(pb)
}


#' Standard ggplot settings
#' 
#' \code{baseTheme} defines common ggplot theme settings for plotting.
#'
#' @param    sizing  defines the style and sizing of the theme. One of 
#'                   \code{c("figure", "window")} where \code{sizing="figure"} is appropriately
#'                   sized for pdf export at 7 to 7.5 inch width, and \code{sizing="window"}
#'                   is sized for an interactive session.
#'
#' @return   A ggplot2 object.
#' 
#' @seealso \link[ggplot2]{theme}.
#' 
#' @export
baseTheme <- function(sizing=c("figure", "window")) {
    # Check arguments
    sizing <- match.arg(sizing)
    
    base_theme <- theme_bw() +
        theme(strip.background=element_blank(),
              plot.background=element_blank(),
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank())
              
    # Define universal plot settings appropriate for PDF figures
    if (sizing == "figure") {
        base_theme <- base_theme + 
            theme(text=element_text(size=8),
                  plot.title=element_text(size=8),
                  strip.text=element_text(size=7, face="bold"),
                  axis.title=element_text(size=8, vjust=0.25),
                  axis.text.x=element_text(size=8, vjust=0.5, hjust=0.5),
                  axis.text.y=element_text(size=8),
                  legend.text=element_text(size=7),
                  legend.title=element_text(size=7),
                  legend.key.height=grid::unit(10, "points"), 
                  legend.key.width=grid::unit(10, "points"))
    } else if (sizing == "window") {
        # Define universal plot settings appropriate for an interactive session
        base_theme <- base_theme + 
            theme(text=element_text(size=14),
                  plot.title=element_text(size=16),
                  strip.text=element_text(size=14, face="bold"),
                  axis.title=element_text(size=16, vjust=0.25),
                  axis.text.x=element_text(size=16, vjust=0.5, hjust=0.5),
                  axis.text.y=element_text(size=16),
                  legend.text=element_text(size=14),
                  legend.title=element_text(size=14),
                  legend.key.height=grid::unit(18, "points"), 
                  legend.key.width=grid::unit(18, "points"))
    }
    
    return(base_theme)
}


#' Plot multiple ggplot objects
#' 
#' Plots multiple ggplot objects in an equally sized grid.
#' 
#' @param   ...    ggplot objects to plot.
#' @param   ncol   number of columns in the plot.
#' @return  NULL
#' 
#' @seealso \link{ggplot}.
#' 
#' @references
#' Modified from:
#' http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
#' 
#' @export
gridPlot <- function(..., ncol=1) {
    p <- list(...)
    n <- length(p)
    layout <- matrix(seq(1, ncol*ceiling(n/ncol)), ncol=ncol, nrow=ceiling(n/ncol))
    
    # Plot
    if (n == 1) {
        plot(p[[1]])
    } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:n) {
            idx <- as.data.frame(which(layout == i, arr.ind=T))
            plot(p[[i]], vp=grid::viewport(layout.pos.row = idx$row, layout.pos.col=idx$col))
        }
    }
}

#### Multiprocessing functions ####

#' Available CPU cores
#'
#' \code{cpuCount} determines the number of CPU cores available.
#' 
#' @return  Count of available cores. Returns 1 if undeterminable.
#'
#' @examples
#' cpuCount()
#' 
#' @export
cpuCount <-function(){
    if (.Platform$OS.type == "windows") {
        nproc <- as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))
    } else if (.Platform$OS.type == "unix") { 
        nproc <- parallel::detectCores()
    } else {
        nproc <- 1    
    }
	
	# in case an NA is returned
	if(is.na(nproc)){nproc <- 1}
    
    return(nproc)
}

#### Generic statistical functions ####

#' Weighted meta-analysis of p-values via Stouffer's method
#'
#' \code{stoufferMeta} combines multiple weighted p-values into a meta-analysis p-value
#' using Stouffer's Z-score method.
#' 
#' @param    p   numeric vector of p-values.
#' @param    w   numeric vector of weights.
#' 
#' @return   A named numeric vector with the combined Z-score and p-value in the form
#'           \code{c(Z, pvalue)}.
#' 
#' @examples
#' # Define p-value and weight vectors
#' p <- c(0.1, 0.05, 0.3)
#' w <- c(5, 10, 1)
#'
#' # Unweighted
#' stoufferMeta(p)
#' 
#' # Weighted
#' stoufferMeta(p, w)
#' 
#' @export
stoufferMeta <- function(p, w=NULL) {
    if (is.null(w)) {
        w <- rep(1, length(p))/length(p)
    } else {
        if (length(w) != length(p)) { stop("Length of p and w must equal.") }
        w <- w/sum(w)
    }
    x <- qnorm(1 - p)
    Z  <- sum(w*x) / sqrt(sum(w^2))
    pvalue <- 1 - pnorm(Z)
    
    return(c(Z=Z, pvalue=pvalue))
}

# Stirling's approximation of the binomial coefficient
# 
# Calculates Stirling's approximation of the binomial coefficient for large numbers.
#
# @param    n  a vector of n.
# @param    k  a vector of k.
#
# @return   The approximation of log(n choose k). For n < 100 \link{lchoose} is used.
lchooseStirling <- function(n, k) {
    if (any(n < k)) {
        stop("n must be >= k")
    }
    
    n_len <- length(n)
    k_len <- length(k)
    nCk <- rep(NA, max(n_len, k_len))
    nCk[n == k] <- 0
    
    # a = index n_small
    # i = index k_small
    # x = index nCk_small
    #
    # b = index n_large
    # j = index k_large
    # y = index nCk_large
    #
    # Check for vector inputs and assign indexing
    if (n_len >= 1 & k_len >= 1 & n_len == k_len) {
        a <- i <- x <- (n < 100 & n != k)
        b <- j <- y <- (n >= 100 & n != k)
    } else if (n_len > 1 & k_len == 1) {
        a <- x <- (n < 100 & n != k)
        b <- y <- (n >= 100 & n != k)
        i <- j <- TRUE
    } else if (n_len == 1 & k_len > 1) {
        a <- (n < 100)
        b <- !a
        i <- j <- (n != k)
        x <- if (n < 100) { i } else { NULL }
        y <- if (n >= 100) { i } else { NULL }
    } else {
        stop("Inputs are wrong. n and k must have the same length or be length one.")
    } 
    
    
    # Small n
    nCk[x] <-  lchoose(n[a], k[i])
    
    # Large n indices
    nCk[y] <- n[b]*log(n[b]) - k[j]*log(k[j]) - (n[b] - k[j])*log(n[b] - k[j]) + 
        0.5*(log(n[b]) - log(k[j]) - log(n[b] - k[j]) - log(2*pi))
    
    #     .nCk <- function(n, k) {
    #         n*log(n) - k*log(k) - (n - k)*log(n - k) + 
    #         0.5*(log(n) - log(k) - log(n - k) - log(2*pi))
    #     }
    
    return(nCk)
}

