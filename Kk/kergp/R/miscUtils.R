##*****************************************************************************
##' Make translucient colors.
##'
##' @title Make translucient colors
##'
##' @param colors A vector of colors in a format that can be
##' understood by \code{\link{col2rgb}}.
##'
##' @param alpha Level of opacity ("0" means fully transparent and
##' "max" means opaque).  After recycling to reach the required
##' length, this value or vector is used as \code{alpha} in
##' \code{\link{rgb}}.
##'
##' @return A vector of translucient (or semi-transparent) colors.
##'
translude <- function (colors, alpha = 0.6) {
    L <- pmax(length(colors), length(alpha))
    colors <- rep(colors, length.out = L)
    alpha <- rep(alpha, length.out = L)
    rgb <- as.matrix(col2rgb(colors)/255)
    colors2 <- rgb(red = rgb["red", ], green = rgb["green", ], 
                   blue = rgb["blue", ], alpha = alpha)
}

##*****************************************************************************
##' Vector of names for the general 'Symm' parameterisation.
##' 
##' @title Vector of Names for the General 'Symm' Parameterisation
##'
##' @param nlev Number of levels.
##'
##' @return Character vector of names.
##'
##' @examples
##' parNamesSymm(nlev = 4)
parNamesSymm <- function(nlev) {

    if (nlev == 1L) return(character(0))
    mat <- matrix(NA, nrow = nlev, ncol = nlev)
    a <- col(mat)[upper.tri(mat)]
    b <- row(mat)[upper.tri(mat)]
    paste("theta", a, b, sep = "_")
    
}

##*****************************************************************************
##' Check length/names for a vector of values for parameters or
##' bounds.
##'
##' @title Check Length and Names of a Vector of Values for Parameters or
##' Bounds
##'
##' @param value Numeric vector of values.
##'
##' @param parN Number of wanted values.
##'
##' @param parNames character. Names of the wanted values. 
##'
##' @param default numeric. Default value.
##'
##' @return A numeric vector.
##'
##' @examples
##' checkPar(value = c(1, 2), parN = 2L, parNames = c("theta", "sigma2"),
##'          default = 1.0)
##' checkPar(value = NULL, parN = 2L, parNames = c("theta", "sigma2"),
##'          default = 1.0)
##' checkPar(value = c("sigma2" = 100, "theta" = 1),
##'          parN = 2L, parNames = c("theta", "sigma2"),
##'          default = 1.0)
##' 
checkPar <- function(value, parN, parNames, default) {
    if (is.null(value)) {
        value <- rep(default, parN)
    }
    if (length(value) != parN) {
        stop("value must have length ", parN, " with values ",
             " for ", parNames)
    }
    if (!is.null(nm <- names(value))) {
        if (!setequal(nm, parNames)) {
            stop("bad names provided for value")
        }
        value <- value[match(nm, parNames)]
    }
    names(value) <- parNames
    value 
}

## ****************************************************************************
##' Modified Helmert contrast (or coding) matrix.
##'
##' The returned matrix is a scaled version of
##' \code{contr.helemert(A)}.
##' 
##' @title Modified Helmert Contrast Matrix
##'
##' @param n Integer.
##' 
##' @return An orthogonal matrix with \code{n} rows and \code{n - 1}
##' columns. The columns form a basis of the subspace orthogonal to
##' a vector of \code{n} ones.
##'
##' @examples
##' A <- contr.helmod(6)
##' crossprod(A)
##' 
contr.helmod <- function(n) {
    A <- stats::contr.helmert(n = n)
    n1 <- n - 1L
    norm <- sqrt((1:n1)^2 + 1L:n1)
    ## variant
    ##  scale(A, center = FALSE, scale = norm)
    sweep(A, MARGIN = 2L, STATS = norm, FUN = "/")
}
    
## *****************************************************************************
##' Vector of indices useful for symmetric or anti-symmetric matrices
##'
##' This function is intended to provide computations which are faster
##' than \code{lower.tri} and \code{upper.tri}.
##'
##' @title Vector of Indices Useful for Symmetric or Anti-Symmetric Matrices.
##' 
##' @param n Size of a square matrix. 
##'
##' @param diag Logical. When \code{FALSE} the diagonal is omitted in
##' the lower and upper triangles.
##'
##' @return A list containing the following integer vectors, each with
##' length \eqn{(n - 1) n / 2}.
##'
##' \item{i, j}{
##'
##' Row and column indices for the lower triangle to be used in a
##' two-indice style.
##'
##' }
##' \item{kL}{
##'  
##' Indices for the lower triangle, to be used in single-index
##' style. The elements are picked in column order. So if \code{X} is
##' a square matrix with size \code{n}, then \code{X[kL]} is the
##' vector containing the elements of the lower triangle of \code{X}
##' taken in colum order.
##' 
##' } \item{kU}{
##'
##' Indices for the upper triangle, to be used in a single-index
##' style.  The elements are picked in row order.  So if \code{X} is a
##' square matrix with size \code{n}, then \code{X[kU]} is the vector
##' containing the elements of the upper triangle of \code{X} taken in
##' row order.
##' 
##' }
##' 
##' @examples
##' n <- rpois(1, lambda = 10)
##' L <- symIndices(n)
##' X <- matrix(1L:(n * n), nrow = n)
##' max(abs(X[lower.tri(X, diag = FALSE)] - L$kL))
##' max(abs(t(X)[lower.tri(X, diag = FALSE)] - L$kU))
##' cbind(row = L$i, col = L$j)
symIndices <- function(n, diag = FALSE) {
    
    if (diag) stop("'diag = TRUE' not implemented yet")
    j <- rep.int(1L:(n - 1L), times = (n - 1L):1L)
    i <- sequence((n - 1L):1L) + j
    kL <- (j - 1L) * n + i
    kU <- (i - 1L) * n + j
    list(i = i, j = j, kL = kL, kU = kU)

}

## ****************************************************************************
##' Optimization methods (or algorithms) for the \code{mle} method.
##'
##' @title Optimization Methods (or Algorithms) for the \code{mle}
##' Method
##'
##' @param optimFun Value of the corresponding formal argument of the
##' \code{mle} method, or \code{"both"}. In the later case the full
##' list of algorithms will be obtained.
##'
##' @return A data frame with two columns: \code{optimFun} and
##' \code{optimMethod}.
##'
##' @references 
##' See \href{https://nlopt.readthedocs.io/en/latest/}{The NLopt website}.
##'
##' @seealso
##'
##' \code{\link{mle-methods}}.
##' 
##' @examples
##' optimMethods()
##' 
optimMethods <- function(optimFun = c("both", "nloptr::nloptr",
                             "stats::optim")) {

    name <- NULL ## to avoid NOTE at check
    optimFun <- match.arg(optimFun) 
    
    if (optimFun != "nloptr::nloptr") {
        oM <- eval(formals(optim)$method)
    } else {
        oM <- character(0)
    }
    
    if (optimFun != "stats::optim") { 
        nO <- nloptr.get.default.options()
        npV <- subset(nO, name == "algorithm")$possible_values
        npV <- strsplit(npV, ", ")[[1]]
        ## Global/Local and No derivative/Derivative Add this a
        ## columns?
        ## L <- strsplit(npV, "_")
        ## GLND <- sapply(L, function(x) x[2])
        ## GL <- substr(GLND, start = 1, stop = 1)
        ## ND <- substr(GLND, start = 2, stop = 2)
    } else {
        npV <- character(0)
    }
        
    df <- data.frame(optimFun = c(rep("stats::optim", length(oM)),
                         rep("nloptr::nloptr", length(npV))),
                     optimMethod = c(oM, npV),
                     stringsAsFactors = FALSE)

    df
}

## *****************************************************************************
##' Given a list defining groups of levels for a qualitative input,
##' two vectors of the same length are built: \code{group} gives the
##' (name of the) group and \code{nestedLevels} gives the nested
##' levels, i.e. the levels within the groups.
##' 
##' @title Translates a List Defining Groups of Levels into two
##' Character Vectors
##'
##' @param groupList A list defining the groups. This must be a list
##' containing atomic vectors, each defining a group of levels. These
##' vectors will be coerced to character. If the list is named, then
##' the names will be used as names for the groups, else default group
##' names we be given based on \code{prefix} and group numbers, see
##' \bold{Examples}.
##'
##' @param prefix A prefix to identify groups.
##'
##' @param sep Separator char used to paste groups and nested levels.
##' 
##' @return A list with the two items \code{group} and
##' \code{nestedLevels}.
##'
##' @section Caution: the levels of the wanted input must all appear
##' exactly once in \code{unlist(groupList)}. We check that the list
##' does not embed duplicated levels, but we can not tackle missing
##' levels here.
##'
##' @examples
##' 
##' gL <- list(letters[1:3], rev(letters[8:4]))
##' parseGroupList(gL)
##' parseGroupList(gL, prefix = "G", sep = "-")
##'
##' parseGroupList(list(c(1, 2, 5), c(4, 3)))
##' 
##' cities <- list("B" = c("AntWerp", "Ghent" , "Charleroi"),
##'                "F" = c("Paris", "Marseille", "Lyon"),
##'                "D" = c("Berlin", "Hamburg", "Munchen"))
##' parseGroupList(cities)
##'
##' ## duplicated levels: error
##' try(parseGroupList(list("a" = c(1, 2, 3), "b" = c(1, 4))))
##' ## not all names provided: use default names
##' try(parseGroupList(list("a" = c(1, 2, 3), c(5, 4))))
##' 
parseGroupList <- function(groupList, prefix = "gr", sep = "/") {

    flat <- unlist(groupList)
    if (any(duplicated(flat))) {
        stop("'groupList' contains duplicated elements")
    }
    
    if (is.null(names(groupList)) || any(names(groupList) == "")) {
        ng <- paste0(prefix, 1:length(groupList))
    } else {
        ng <-  names(groupList)
    }

    lg <- sapply(groupList, length)
    
    group <- character(0)
    nestedLevels <- character(0)
    for (i in seq_along(lg)) {
        group <- c(group, rep(ng[i], lg[i]))
        nestedLevels <-
            c(nestedLevels, paste(ng[i], groupList[[i]], sep = sep))
    }
    
    list(group = group,
         nestedLevels = nestedLevels,
         levels = flat)

}

