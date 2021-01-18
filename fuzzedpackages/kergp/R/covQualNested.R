setClass("covQualNested",   	
         representation(
             covLevels = "function",
             covLevMat = "matrix",
             hasGrad = "logical",           
             acceptLowerSQRT = "logical",    
             label = "character",
             d = "integer",
             inputNames = "character",     
             nlevels = "integer",                
             levels = "list",          
             parLower = "numeric", 
             parUpper = "numeric",       
             par = "numeric",             
             parN = "integer",             
             kernParNames  = "character",
             group = "integer",             ## NEW
             groupLevels = "character",     ## NEW
             between = "covQual",           ## NEW
             within = "list",               ## NEW
             parNCum = "integer",           ## NEW
             contrasts = "function"         ## New
             ),
         contains = "covQual")

## ****************************************************************************
##' Nested Qualitative Covariance
##'
##' @title Nested Qualitative Covariance
##'
##' @include q1Symm.R
##' @include q1CompSymm.R
##'
##' @param input Name of the input, i.e. name of the column in the
##' data frame when the covariance kernel is evaluated with the
##' \code{\link{covMat,covQual-method}} method.
##'
##' @param group A factor of vector giving the groups, see
##' \bold{Examples}.
##'
##' @param nestedLevels A factor or a vector giving the (nested) levels
##' within the group for each level of \code{group}. If this is
##' missing, each element of \code{group} is assumed to correspond to
##' one nested level within the group and the levels within the group
##' are taken as integers in the order of \code{group} elements. 
##' 
##' @param between Character giving the type of structure to use for
##' the \emph{between} part. For now this can be one of the three
##' choices \code{"Diag"} for the diagonal structure of
##' \code{\link{q1Diag}}, \code{"Symm"} for the general covariance of
##' \code{\link{q1Symm}}, or \code{"CompSymm"} for the Compound
##' Symmetry covariance of \code{\link{q1CompSymm}}.
##'
##' @param within Character vector giving the type of structure to use
##' for the \emph{within} part. The choices are the same as for
##' \code{between}. The character vector is recycled to have length
##' \eqn{G} so the \emph{within} covariances can differ across groups.
##'
##' @param compGrad Logical.
##'
##' @param covBet Integer indicating whether a \emph{between}
##' correlation or a covariance is wanted, as in \code{\link{q1Diag}},
##' \code{\link{q1Symm}} or \code{\link{q1CompSymm}}. See section
##' \bold{Caution}.
##'
##' @param covWith Integer indicating whether a \emph{within}
##' correlation or a covariance is wanted, as in \code{\link{q1Diag}},
##' \code{\link{q1Symm}} or \code{\link{q1CompSymm}}. See section
##' \bold{Caution}.
##' 
##' @param contrasts Object of class \code{"function"}. This function
##' is similar to the \code{\link[stats]{contr.treatment}} or
##' \code{\link[stats]{contr.treatment}} functions, but it must return
##' an \emph{orthogonal} matrix. For a given integer \code{n}, it
##' returns a matrix with \code{n} rows and \code{n - 1} columns
##' forming a basis for the supplementary of a vector of ones in the
##' \eqn{n}-dimensional Euclidean space. The
##' \code{\link{contr.helmod}} can be used to obtain an orthogonal
##' matrix hence defining an orthonormal basis.
##'
##' @return An object with class \code{"covQualNested"}.
##'
##' @section Caution: When \code{covBet} and \code{covWith} are zero,
##' the resulting matrix \emph{is not a correlation matrix}, due to
##' the mode of construction. The "between" covariance matrix is a
##' correlation but diagonal blocks are added to the extended matrix
##' obtained by resizing the "between" covariance into a \eqn{n \times
##' n}{ n * n} matrix.
##'
##' @note For now the replacement method such as \code{'coef<-'} are
##' inherited from the class \code{covQuall}. Consequently when these
##' methods are used they do not update the covariance structure in
##' the \code{between} slot nor those in the \code{within} (list)
##' slot.
##'
##' This covariance kernel involves \code{two} categorical
##' (i.e. factor) inputs, but these are nested. It could be aliased in
##' the future as \code{q1Nested} or \code{q2Nested}.
##' 
##' @examples
##' country <- c("B", "B", "B", "F", "F" ,"F", "D", "D", "D")
##' cities <- c("AntWerp", "Ghent" , "Charleroi", "Paris", "Marseille", "Lyon",
##'             "Berlin", "Hambourg", "Munchen")
##' 
##' ## create a nested covariance. Note that it will be used with ONE factor
##' ## input obtained by pasing the group and subgroup.
##'
##' nest <- covQualNested(input = "ccities",
##'                       group = country, nestedLevels = cities,
##'                       between = "Symm", within = "Diag",
##'                       compGrad = TRUE,
##'                       covBet = 0, covWith = 0)
##' 
##' ## 'show' method as automatically invocated
##' nest
##' 
##' ## 'covMat' method: if provided, 'X' must be  a data frame with a
##' ## suitable factor or integer column. The input name was given
##' ## at the creation of the covariance structure
##' 
##' Z <- sample(levels(nest)[[1]], size = 30, replace = TRUE)
##' df <- data.frame(ccities = factor(Z))
##' covMat(nest, X = df)
##' covMat(nest)
##' \dontrun{
##'     if (require(corrplot)) corrplot(cov2cor(covMat(nest)))
##' }
##'
##' ## 'simulate' method. One can give a 'X' formal, but with no
##' ## duplicated values!
##' 
##' Sim <- simulate(nest, nsim = 30)
##' levs <- levels(nest)[[1]]
##' matplot(Sim, type = "l", main = "Simulated paths", xaxt = "n")
##' axis(side = 1, at = seq_along(levs), labels = levs)
##'
##' ## another 'covMat'
##' 
##' cc <- as.factor(sample(paste(country, cities, sep = "/"),
##'                 size = 100, replace = TRUE))
##' covMat(nest, X = data.frame(ccities = cc))
##' 
covQualNested <- function(input = "x",
                          groupList = NULL,
                          group = NULL,
                          nestedLevels = NULL,
                          between = "Symm",
                          within = "Diag", 
                          covBet = c("corr", "homo", "hete"),
                          covWith = c("corr", "homo", "hete"),
                          compGrad = TRUE,
                          contrasts = contr.helmod,
                          intAsChar = TRUE) {

    ## if (intAsChar) {
    ##     warning("With this object, an input of class \"integer\" will be coerced ",
    ##             "into \"character\", not into \"factor\". Use `intAsChar = FALSE` ",
    ##             " to change this behaviour")
    ## }
    
    covBet <- match.arg(covBet)
    covWith <- match.arg(covWith, several.ok = TRUE)
  
    if ((covBet == "corr") && any(covWith != "corr")) {
        stop("If 'covBet' is 'corr', no 'covWith' element can be != 'corr'")
    }
    
    between <- match.arg(between, choices = c("Diag", "CompSymm", "Symm"))
    within <- match.arg(within, choices = c("Diag", "CompSymm", "Symm"), several.ok = TRUE)

    ## check that 'group' is given in some order (with repeted
    ## values).  Without this condition the gradient would be wrong.
    cGroup <- as.character(group)
    uGroup <- unique(cGroup)
    m <-  match(cGroup, uGroup)
    d <- diff(m)
    if (any(d < 0) || any(d > 1)) {
        stop("'group' must contain repeated values given in order")
    }
    
    contrasts <- match.fun(contrasts)

    if (!missing(group)) {
    
        group <- factor(group, levels = unique(group))
        igroup <- as.integer(group)
        
        mVec <- as.vector(table(group))
        
        m <- sum(mVec)
        G <- nlevels(group)
        
        if (!is.null(nestedLevels)) {
            if (length(nestedLevels) != length(group)) {
                stop("when 'nestedLevels' is provided it must have the ",
                     "same length as 'group'")
            }
            allLev <- paste(as.character(group),
                            as.character(nestedLevels), sep = "/")
            if (any(duplicated(allLev))) {
                stop("Comibnations of 'group' and 'nestedLevels' must not",
                     "be duplicated")
            }
            
        } else {
            nestedLevels <-  unlist(tapply(group, group, function(x) 1:length(x)))
        }
        
        nestedLevels <- paste(as.character(group), nestedLevels, sep = "/")

    } else {
        
        if (missing(groupList)) {
            stop("One of the two formals 'group' and 'groupList' must be given")
        }

        pgl <-  parseGroupList(groupList)
        group <- factor(pgl$group, levels = unique(pgl$group))
        igroup <- as.integer(group)
        mVec <- as.vector(table(group))
        
        m <- sum(mVec)
        G <- nlevels(group)

        nestedLevels <- pgl$levels
    }

    
    within <- rep(within, length.out = G)
    covWith <- rep(covWith, length.out = G)

    parN <- integer(0)
    par <- numeric(0)
    parUpper <- numeric(0)
    parLower <- numeric(0)
    parNames <- character(0)
    
    ## create the "between" covariance structure
    if (between == "Diag") {
        betweenC <- q1Diag(factor = group, input = "group", cov = covBet)
    } else if (between == "Symm") {
        betweenC <- q1Symm(factor = group, input = "group", cov = covBet)
    } else if (between == "CompSymm") {
        betweenC <- q1CompSymm(factor = group, input = "group", cov = covBet)
    }

    par <- c(par, coef(betweenC))
    parUpper <- c(parUpper, coefUpper(betweenC))
    parLower <- c(parLower, coefLower(betweenC))
    if (betweenC@parN) {
        parNames <- c(parNames, paste("bet", betweenC@kernParNames, sep = "_"))
    }
    parN <- c(parN, betweenC@parN)
    
    ## ========================================================================
    ## create the "between" covariance structure
    ## Mind the possibility that a covariance structure has zero parameters!!!
    ## ========================================================================

    withinC <- list()
    
    for (g in 1:G) {
      if (mVec[g] == 1){
        parN <- c(parN, 0L)
      } else {
        fac <- as.factor(1:(mVec[g] - 1))
        if (within[g] == "Diag") {
            withinC[[g]] <- q1Diag(factor = fac, cov = covWith[g])
        } else if (within[g] == "Symm") {
            withinC[[g]] <- q1Symm(factor = fac, cov = covWith[g])
        } else if (within[g] == "CompSymm") {
            withinC[[g]] <- q1CompSymm(factor = fac, cov = covWith[g])
        }
        par <- c(par, coef(withinC[[g]]))
        parUpper <- c(parUpper, coefUpper(withinC[[g]]))
        parLower <- c(parLower, coefLower(withinC[[g]]))
        if (withinC[[g]]@parN) {
          parNames <- c(parNames,
                        paste("with", g, withinC[[g]]@kernParNames, sep = "_"))
        }
        parN <- c(parN, withinC[[g]]@parN)
      }
    }

    parNCum <- cumsum(parN)

    thisCovLevel <- function(par, lowerSQRT = FALSE, compGrad = FALSE) {
                
        covMatBet <- betweenC@covLevels(par[1L:parNCum[1]], compGrad = compGrad)
        covMat <- covMatBet[igroup, igroup]
        
        if (compGrad) {
            grad <- array(0.0, dim = c(m, m, parNCum[G + 1L]),
                          dimnames = list(NULL, NULL, parNames))
            grad[ , , 1:parNCum[1]] <-
                attr(covMatBet, "gradient")[igroup, igroup, , drop = FALSE] 
        }
        
        for (g in 1:G) {
            ## cat("within = ", within[g], "\n")

          if (mVec[g] > 1) {
              indLev <- (igroup == g)
              nparg <- parNCum[g + 1L] - parNCum[g]
              
              if (parN[g + 1L]) {
                indPar <- (parNCum[g] + 1L):parNCum[g + 1L]
                parg <- par[indPar]
              } else parg <- numeric(0)
            
              covMatWith <- withinC[[g]]@covLevels(parg,
                                                   compGrad = compGrad)
              A <- contrasts(mVec[g])
              covMatWith2 <- A %*% covMatWith %*% t(A)
            
              if (compGrad) gradProv <- attr(covMatWith, "gradient")
              covMat[indLev, indLev] <- covMat[indLev, indLev] + covMatWith2
            
              if (compGrad && parN[g + 1L]) {
                  ## XXX il faut sweeper
                gradProv2 <- array(0.0, dim = c(mVec[g], mVec[g], parN[g + 1L]))
                for (i in 1L:parN[g + 1L]) {
                    gradProv2[ , , i] <- A %*% gradProv[ , , i] %*% t(A)
                }
                grad[indLev, indLev, indPar] <- grad[indLev, indLev, indPar,
                                                     drop = FALSE] + gradProv2
              }
            }
        }
        
        rownames(covMat) <- colnames(covMat) <- nestedLevels[[1]]
        
        if (compGrad) attr(covMat, "gradient") <- grad
        covMat
    }

    nestedLevels <- list(nestedLevels)
    names(nestedLevels) <- input
    
    covMat <- thisCovLevel(par, lowerSQRT = FALSE, compGrad = FALSE)

    new("covQualNested",            
        covLevels = thisCovLevel,
        covLevMat = covMat,
        hasGrad = TRUE,
        acceptLowerSQRT = FALSE,
        label = "Nested covariance",
        d = 1L,
        inputNames = input,
        nlevels = length(group),
        levels = nestedLevels,
        parLower = parLower,
        parUpper = parUpper,
        par = par,
        parN = length(par),
        kernParNames = parNames,
        group = igroup,                 ## NEW
        groupLevels = levels(group),    ## NEW
        between = betweenC,             ## NEW
        within = withinC,               ## NEW
        parNCum = parNCum,              ## NEW
        contrasts = contrasts,          ## NEW
        ordered = FALSE,
        intAsChar = intAsChar
        )              

}

