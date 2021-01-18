

#' Create a skeleton for a package extending lolog
#' @param path where to create the package
#' @seealso \code{\link{inlineLologPlugin}}
#' @details
#' lolog is a modular package, and can be extended at 
#' both the R and C++ level. This function will build a package
#' skeleton that can be used as a starting point for
#' development. To create the package in the current directory
#' run:
#' 
#' \code{lologPackageSkeleton()}
#' 
#' Build and install the package from the command line with
#' 
#' \code{R CMD build LologExtension}
#' 
#' \code{R CMD INSTALL LologExtension_1.0.tar.gz}
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#' #install package
#' lologPackageSkeleton()
#' system("R CMD build LologExtension")
#' system("R CMD INSTALL LologExtension_1.0.tar.gz")
#' 
#' library(LologExtension) #Load package
#' 
#' # Run model with new minDegree statistic
#' library(network)
#' m <- matrix(0,20,20)
#' for(i in 1:19) for(j in (i+1):20) m[i,j] <- m[j,i] <- rbinom(1,1,.1)
#' g <- network(m, directed=FALSE)
#' fit <- lologVariational(g ~ edges() + minDegree(1L))
#' summary(fit)
#' 
#' }
lologPackageSkeleton <- function(path = ".") {
  pkgPath <- find.package("lolog")
  p <- file.path(pkgPath, "examplePackage", "LologExtension")
  file.copy(p, path, recursive = TRUE)
}


#
# evaluates the formula to get all terms, which can then be used to construct a model
#
.prepModelTerms <- function(formula) {
  if (is.null(formula))
    return(NULL)
  form <- formula
  env <- environment(formula)
  
  # Parse out vertex ordering (if it exists)
  vertexOrder <- integer()
  
  if (!is.symbol(formula[[3]]) &&
      as.character(formula[[3]][[1]]) == "|") {
    tmp <- formula[[3]][[3]]
    vertexOrder <- eval(tmp, envir = env)
    if (any(is.na(vertexOrder)))
      stop("vertex order can not have any NA values")
    if (!is.numeric(vertexOrder))
      stop("vertex order must be numeric")
    form[[3]] <- formula[[3]][[2]]
  }
  
  # parse out model terms
  tmp <- form[[3]]
  lastTerm <- FALSE
  stats <- list()
  offsets <- list()
  while (!lastTerm) {
    ls <- length(stats)
    lo <- length(offsets)
    term <- if (is.symbol(tmp)) {
      lastTerm <- TRUE
      term <- tmp
    } else if (as.character(tmp[[1]]) == "+") {
      tmp[[3]]
    } else{
      lastTerm <- TRUE
      tmp
    }
    
    name <-
      if (is.symbol(term))
        as.character(term)
    else
      as.character(term[[1]])
    args <- NULL
    if (name == "offset" || name == "constraint") {
      term <- term[[2]]
      name <-
        if (is.symbol(term))
          as.character(term)
      else
        as.character(term[[1]])
      if (length(term) > 1) {
        term[[1]] <- as.name("list")
        args <- eval(term, envir = env)
      } else{
        args <- list()
      }
      offsets[[lo + 1]] <- args
      names(offsets)[lo + 1] <- name
    } else{
      if (length(term) > 1) {
        term[[1]] <- as.name("list")
        args <- eval(term, envir = env)
      } else{
        args <- list()
      }
      stats[[ls + 1]] <- args
      names(stats)[ls + 1] <- name
    }
    if (!lastTerm)
      tmp <- tmp[[2]]
  }
  
  list(stats = stats,
       offsets = offsets,
       vertexOrder = vertexOrder)
}


# Groeneveld & Meeden
.gmSkewness <- function(x){
  med <- stats::median(x)
  (mean(x) - med) / mean(abs(x - med))
}

# Add nice histogram panels to a pairs plot
.panelHist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x[-length(x)], plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey")
  abline(v=x[length(x)], col="red", lwd=3)
}
