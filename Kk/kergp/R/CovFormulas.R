##' Parse a formula (or expression) describing a composite covariance
##' kernel.
##'
##' The formula involves existing covariance kernel objects and must
##' define a valid kernel composition rule. For instance the sum and
##' the product of kernels, the convex combination of kernels are
##' classically used. The kernels objects are used in the formula with
##' parentheses as is they where functions calls with no formal
##' arguments e.g. \code{obj( )}. Non-kernel objects used in the
##' formula must be numeric scalar parameters and are called \emph{top}
##' parameters. The covariance objects must exist in the environment
##' defined by \code{where} because their slots will be used to
##' identify the inputs and the parameters of the composite kernel
##' defined by the formula.
##' 
##' @title Parse a Formula or Expression Describing a Composite
##' Covariance Kernel
##' 
##' @param formula A formula or expression describing a covariance
##' kernel. See \bold{Examples}.
##' 
##' @param where An environment where kernel objects and top
##' parameters are searched for.
##' 
##' @param trace Integer level of verbosity.
##'
##' @return A list with the results of parsing. Although the results
##' content is easy to understand, the function is not intended to be
##' used by the final user, and the results may change in future
##' versions.
##'
##' @section Caution: Only relatively simple formulas are correctly
##' parsed. So use only formulas having a structure similar to one of
##' those given in the examples. In case of problems, error messages
##' are likely to be difficult to understand.
##'
##' @note The parsing separates covariance objects from top
##' parameters.  It retrieves information about the kernel inputs and
##' parameters from the slots. Obviously, any change in the
##' covariances objects after the parsing (e.g. change in the
##' parameters names or values) will not be reported in the results of
##' the parsing, so kernel any needed customization must be done prior
##' to the parsing.
##' 
##' @author Yves Deville
##'
##' @examples
##' ## =========================================================================
##' ## build some kernels (with their inputNames) in the global environment
##' ## =========================================================================
##' 
##' myCovExp3 <- kMatern(d = 3, nu = "1/2")
##' inputNames(myCovExp3) <- c("x", "y", "z")
##'
##' myCovGauss2 <- kGauss(d = 2)
##' inputNames(myCovGauss2) <- c("temp1", "temp2")
##'
##' k <- kMatern(d = 1)
##' inputNames(k) <- "x"
##'
##' ell <- kMatern(d = 1)
##' inputNames(ell) <- "y"
##'
##' ## =========================================================================
##' ## Parse a formula. This formula is stupid because 'myCovGauss2'
##' ## and 'myCovExp3' should be CORRELATION kernels and not
##' ## covariance kernels to produce an identifiable model.
##' ## =========================================================================
##' 
##' cov <- ~ tau2 * myCovGauss2() * myCovExp3() + sigma2 * k()
##' pf <- parseCovFormula(cov, trace = 1)
##'
##' ## =========================================================================
##' ## Parse a formula with ANOVA composition
##' ## =========================================================================
##'
##' cov1 <- ~ tau2 * myCovGauss2() * myCovExp3() + sigma2 * (1 + k()) * (1 + ell())
##' pf1 <- parseCovFormula(cov1, trace = 1) 
##' 
parseCovFormula <- function(formula, where = .GlobalEnv,
                            trace = 0) {

    ## ==========================================================================
    ## First make sur that we work with an expression
    ## ==========================================================================
    
    if (is(formula, "formula")) {
        covExpr <- TRUE
        formula <- as.expression(formula[[2]])
    } else if (is(formula, "expression")) {
        covExpr <- TRUE
    } else if (is.character(formula)) {
        covExpr <- TRUE
        formula <- parse(text = formula)
    } else {
        stop("'formula' must here be a formula, an expression or a character")
    }
    
    expr <- formula

    if (trace) {
        cat("o Initial expression\n")
        cat(deparse(expr), "\n\n")
    }
    
    ## ==========================================================================
    ## find the variables, which include kernel names.
    ## ==========================================================================

    aVar <- all.vars(expr)
    
    ## ==========================================================================
    ## find the covariance kernel names. There might be a more direct
    ## solution
    ## ==========================================================================

    kernNames <- setdiff(all.names(expr), c(aVar, "*", "+", "("))
   
    if (any(duplicated(kernNames))) {
        stop("For now a kernel name can not be repeated in the formula")
    }
    
    ## ==========================================================================
    ##  Check that all 'kernNames' are actually the names of an object
    ##  inheriting from the class "covAll", and 
    ## ==========================================================================

    inputList <- list()
    parNames <- character(0)
    parNamesList <- list()
    hasGrad <- logical(length(kernNames))
    names(hasGrad) <- kernNames
    par <- numeric(0)
    parLower <- numeric(0)
    parUpper <- numeric(0)
    kernClasses <- character(0)
    
    for (kn in kernNames) {
        obj <- get(kn, envir = where)
        if (!inherits(x = obj, what = "covAll")) {
            stop("object '", kn, "' not identified as a kernel")
        }
        kernClasses <- c(kernClasses, class(obj))
        inputList[[kn]] <- inputNames(obj)
        ## use a 'parNames' method (to be implemented)
        parNamesList[[kn]] <- obj@kernParNames
        parNames <- c(parNames, paste(kn, obj@kernParNames, sep = "."))
        hasGrad[kn] <- obj@hasGrad
        par <- c(par, obj@par)
        parLower <- c(parLower, obj@parLower)
        parUpper <- c(parUpper, obj@parUpper)
    }
      
    ## ==========================================================================
    ## Now for a new expression that will be used to compute the covMat
    ## ==========================================================================

    str <- str0 <- strD <- as.character(expr)
    for (kn in kernNames) {
        str <- gsub(pattern = paste(kn, " *\\(", sep = ""),
                    replacement = sprintf("covMat(%s, X = X", kn),
                    str)
        strD <- gsub(pattern = paste(kn, " *\\([^\\)]*\\)", sep = ""),
                    replacement = sprintf("%s", kn),
                    strD)
    }

    
    if (trace) {
        cat("o Transformed expression\n")
        cat(strD, "\n\n")
    }

    inputs <- unique(unlist(inputList))
    d <- length(inputs)

    ## XXX donner une valeur aux parametres top
    tops <- setdiff(aVar, kernNames)
    if (length(tops)) {
        parNames <- c(parNames, paste(".top", tops, sep = "."))
    }
    parNamesList[[".top"]] <- tops
    
    ## ==========================================================================
    ## Compute the gradients
    ## ==========================================================================
    
    toDer <- parse(text = strD)

    ## if (trace) {
        ## cat(sprintf("o Expression to differentiate\n%s\n\n", strD))
    ## }
    
    gradExp <- parsedGradExp <- list()
    parNamesShort <- character(0)
    
    ## differentiate w.r.t. kernel parameters
    for (kn in kernNames) {
        s <- deparse(D(toDer, kn))
        pn <- paste(kn, obj@kernParNames, sep = ".")
        for (ppn in parNamesList[[kn]]) {
            ppnLong <- paste(kn, ppn, sep = ".")
            parNamesShort <- c(parNamesShort, ppn)
            gradExp[[ppnLong]] <-
                sprintf("%s * .%s[ , , \"%s\"]", s, kn, ppn)
            parsedGradExp[[ppnLong]] <- parse(text = gradExp[[ppnLong]])
           
        }
    }
    
    ## differentiate w.r.t. top parameters
    for (pn in tops) {
        ppnLong <- paste(".top", pn, sep = ".")
        parNamesShort <- c(parNamesShort, pn)
        gradExp[[ppnLong]] <- deparse(D(toDer, pn))
        parsedGradExp[[ppnLong]] <- parse(text = gradExp[[ppnLong]])
    }
    
    names(parNamesShort) <- parNames
    
    ## allNms.inputs <- unique(c(unlist(inputs), other.inputs))
    
    list(expr = expr,
         d = d,
         inputNames = inputs,
         inputNamesList = inputList,
         hasGrad = hasGrad,
         kernNames = kernNames,
         kernClasses = kernClasses,
         parNames = parNames,
         parNamesList = parNamesList,
         parNamesShort = parNamesShort,
         par = par,
         parLower = parLower,
         parUpper = parUpper,
         covMatCall = str,
         simpleCall = strD,
         parsedSimpleCall = parse(text = strD),
         gradExp = gradExp,
         parsedGradExp = parsedGradExp)
    
}

