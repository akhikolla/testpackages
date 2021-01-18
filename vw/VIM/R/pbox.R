# ----------------------------------------------------------
# Authors: Andreas Alfons, Bernd Prantner, Matthias Templ
#          and Daniel Schopfhauser
#          Vienna University of Technology
# ----------------------------------------------------------



#' Parallel boxplots with information about missing/imputed values
#' 
#' Boxplot of one variable of interest plus information about missing/imputed
#' values in other variables.
#' 
#' This plot consists of several boxplots. First, a standard boxplot of the
#' variable of interest is produced. Second, boxplots grouped by observed and
#' missing/imputed values according to `selection` are produced for the
#' variable of interest.
#' 
#' Additionally, the frequencies of the missing/imputed values can be
#' represented by numbers.  If so, the first line corresponds to the observed
#' values of the variable of interest and their distribution in the different
#' groups, the second line to the missing/imputed values.
#' 
#' If `interactive=TRUE`, clicking in the left margin of the plot results
#' in switching to the previous variable and clicking in the right margin
#' results in switching to the next variable.  Clicking anywhere else on the
#' graphics device quits the interactive session.
#' 
#' @param x a vector, matrix or `data.frame`.
#' @param delimiter a character-vector to distinguish between variables and
#' imputation-indices for imputed variables (therefore, `x` needs to have
#' [colnames()]). If given, it is used to determine the corresponding
#' imputation-index for any imputed variable (a logical-vector indicating which
#' values of the variable have been imputed). If such imputation-indices are
#' found, they are used for highlighting and the colors are adjusted according
#' to the given colors for imputed variables (see `col`).
#' @param pos a numeric value giving the index of the variable of interest.
#' Additional variables in `x` are used for grouping according to
#' missingness/number of imputed missings.
#' @param selection the selection method for grouping according to
#' missingness/number of imputed missings in multiple additional variables.
#' Possible values are `"none"` (grouping according to missingness/number
#' of imputed missings in every other variable that contains missing/imputed
#' values), `"any"` (grouping according to missingness/number of imputed
#' missings in *any* of the additional variables) and `"all"`
#' (grouping according to missingness/number of imputed missings in *all*
#' of the additional variables).
#' @param col a vector of length five giving the colors to be used in the plot.
#' The first color is used for the boxplots of the available data, the
#' second/fourth are used for missing/imputed data, respectively, and the
#' third/fifth color for the frequencies of missing/imputed values in both
#' variables (see \sQuote{Details}).  If only one color is supplied, it is used
#' for the boxplots for missing/imputed data, whereas the boxplots for the
#' available data are transparent.  Else if two colors are supplied, the second
#' one is recycled.
#' @param numbers a logical indicating whether the frequencies of
#' missing/imputed values should be displayed (see \sQuote{Details}).
#' @param cex.numbers the character expansion factor to be used for the
#' frequencies of the missing/imputed values.
#' @param xlim,ylim axis limits.
#' @param main,sub main and sub title.
#' @param xlab,ylab axis labels.
#' @param axes a logical indicating whether axes should be drawn on the plot.
#' @param frame.plot a logical indicating whether a box should be drawn around
#' the plot.
#' @param labels either a logical indicating whether labels should be plotted
#' below each box, or a character vector giving the labels.
#' @param interactive a logical indicating whether variables can be switched
#' interactively (see \sQuote{Details}).
#' @param \dots for `pbox`, further arguments and graphical parameters to
#' be passed to [graphics::boxplot()] and other functions.  For
#' `TKRpbox`, further arguments to be passed to `pbox`.
#' @return a list as returned by [graphics::boxplot()].
#' @note Some of the argument names and positions have changed with version 1.3
#' due to extended functionality and for more consistency with other plot
#' functions in `VIM`.  For back compatibility, the arguments `names`
#' and `cex.text` can still be supplied to \code{\dots{}} and are handled
#' correctly.  Nevertheless, they are deprecated and no longer documented.  Use
#' `labels` and `cex.numbers` instead.
#' @author Andreas Alfons, Matthias Templ, modifications by Bernd Prantner
#' @seealso [parcoordMiss()]
#' @references M. Templ, A. Alfons, P. Filzmoser (2012) Exploring incomplete
#' data using visualization tools.  *Journal of Advances in Data Analysis
#' and Classification*, Online first. DOI: 10.1007/s11634-011-0102-y.
#' @keywords hplot
#' @family plotting functions
#' @examples
#' 
#' data(chorizonDL, package = "VIM")
#' ## for missing values
#' pbox(log(chorizonDL[, c(4,5,8,10,11,16:17,19,25,29,37,38,40)]))
#' 
#' ## for imputed values
#' pbox(kNN(log(chorizonDL[, c(4,8,10,11,17,19,25,29,37,38,40)])),
#'      delimiter = "_imp")
#' 
#' @export
pbox <- function(x, delimiter = NULL, pos = 1, selection = c("none","any","all"), 
                 col = c("skyblue","red","red4","orange","orange4"), numbers = TRUE, 
                 cex.numbers = par("cex"), xlim = NULL, ylim = NULL, 
                 main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
                 axes = TRUE, frame.plot = axes, labels = axes, 
                 interactive = TRUE, ...) {
  check_data(x)
  x <- as.data.frame(x)
    # initializations and error messages
	imputed <- FALSE # indicates if there are Variables with missing-index
    if(is.null(dim(x))) {  # vector
        n <- length(x)
        p <- 1
        if(n == 0) stop("'x' must have positive length")
    } else {  # matrix or data.frame
        if(!(inherits(x, c("data.frame","matrix")))) { 
            stop("'x' must be a data.frame or matrix")
        }
		
		##delimiterh ##
		if(!is.null(delimiter)) {
			tmp <- grep(delimiter, colnames(x)) # Position of the missing-index
			if(length(tmp) > 0) {
				imp_var <- x[, tmp, drop=FALSE]
				x <- x[, -tmp, drop=FALSE]
				
				if(ncol(x) == 0) stop("Only the missing-index is given")
				if(is.matrix(imp_var) && range(imp_var) == c(0,1)) imp_var <- apply(imp_var,2,as.logical)
				
				if(is.null(dim(imp_var))) {
					if(!is.logical(imp_var)) stop("The missing-index of imputed Variables must be of the type logical")
				} else {
					if(!any(as.logical(lapply(imp_var,is.logical)))) stop("The missing-index of imputed Variables must be of the type logical")	
				}
				imputed <- TRUE
			} else {
				warning("'delimiter' is given, but no missing-index-Variable is found", call. = FALSE)
			}
		}
		
        n <- nrow(x)
        p <- ncol(x)
        if(n == 0) stop("'x' has no rows")
        else if(p == 0) stop("'x' has no columns")
        if(is.null(colnames(x))) colnames(x) <- defaultNames(p)
    }
    if(p == 1) {
        pos <- 1
        interactive <- FALSE
    } else {
        if(!is.numeric(pos) || length(pos) != 1 ||(p < pos)) {
            stop("'pos' must be an integer specifying one column of 'x' and must be lesser than the number of colums of 'x'")
        }
        if(p == 2) selection <- "none"
        else selection <- match.arg(selection)
    }

    if(length(col) == 0) col = c("skyblue","red","red4","orange","orange4")
    else if(length(col) == 1) col <- c("transparent", rep.int(col, 4))
    else if(length(col) == 2 || length(col) == 4)  col <- c(col, rep(col[2],3))
	else if(length(col) != 5) col <- c(col[1],rep(col[2:3],2))
    
    # define local function and initialize call
    localBoxplot <- function(..., plot, log, 
        axes, frame.plot, horizontal, add) {
        boxplot(..., add=TRUE, axes=FALSE)
    }
    ca <- as.call(list(localBoxplot, ...))
    nmdots <- names(ca)[-1]
    
    # back compatibility
    if(missing(cex.numbers) && "cex.text" %in% nmdots) {
        cex.numbers <- ca$cex.text
    }
    if(missing(labels) && "names" %in% nmdots) labels <- ca$names
   