#' @title Generic function for conditional independence test
#' 
#' @description Generic function for conditional independence test. Specializes
#'     to specific types of data.
#'
#' @name citest-generic
#' 
#' @param x An object for which a test for conditional independence is to be
#'     made. See 'details' for valid types of \code{x}.
#' 
#' @param set A specification of the test to be made. The tests are of the form
#'     u and v are independent condionally on S where u and v are variables and
#'     S is a set of variables. See 'details' for details about specification of
#'     \code{set}.
#' 
#' @param \dots Additional arguments to be passed on to other methods.
#'
#' @return An object of class `citest` (which is a list).
#'
#' @details \code{x} can be
#'  1. a table,
#'  1. a dataframe whose columns are
#'     numerics and factors or
#'  1. a list with components \code{cov} and
#'     \code{n.obs}.
#'
#' @details
#'  \code{set} can be
#'  1. a vector,
#'  1. a right-hand sided
#'     formula in which variables are separated by '+'.
#'
#' In either case, it is tested if the first two variables in the
#' \code{set} are conditionally independent given the remaining
#' variables in \code{set}.  (Notice an abuse of the '+' operator in
#' the right-hand sided formula: The order of the variables does
#' matter.)
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{ciTest_table}},
#'     \code{\link{ciTest_df}},
#'     \code{\link{ciTest_mvn}},
#'     \code{\link{chisq.test}}
#' @keywords htest
#' @examples
#' 
#' ## contingency table:
#' data(reinis)
#' ## dataframe with only numeric variables:
#' data(carcass)
#' ## dataframe with numeric variables and factors:
#' data(milkcomp1)
#' 
#' ciTest(cov.wt(carcass, method='ML'), set=~Fat11 + Meat11 + Fat12)
#' ciTest(reinis, set=~smo + phy + sys)
#' ciTest(milkcomp1, set=~tre + fat + pro)


#' @export ciTest
ciTest <- function(x, set=NULL, ...){
  UseMethod("ciTest")
}

#' @export
ciTest.table <- function(x, set=NULL, ...){
  ciTest_table(x, set, ...)
}

#' @export
ciTest.list <- function(x, set=NULL, ...){
  ciTest_mvn(x, set, ...)
}
#' @export
ciTest.data.frame <- function(x, set=NULL, ...){
  ciTest_df(x, set, ...)
}

#' @export
print.citest <- function(x, ...){
    if (length(x$varNames) > 2){
        cat("Testing", x$varNames[1], "_|_", x$varNames[2], "|",x$varNames[-(1:2)],"\n")
    } else {
        cat("Testing", x$varNames[1], "_|_", x$varNames[2], "\n")
    }
    cat(sprintf("Statistic (%s): %8.3f df: %s p-value: %6.4f method: %s\n",
                x$statname, x$statistic, x$df, x$p.value, x$method))

    if ( !is.null(x$slice) ){
        cat("Slice information:\n")
        print( x$slice, digits=4 )
    }

    invisible( x )
}

#' @export
summary.citest <- function(object,...){
    print( object )
    if ( !is.null(object$slice) ){
        cat("Slice information:\n")
        print( object$slice, digits=4 )
    }
    invisible( object )
}



#############################################################################
#'
#' @title Test for conditional independence in a dataframe
#' @description Test for conditional independence in a dataframe.
#' @name citest-df
#'
#############################################################################
#'
#' @param x A dataframe.
#' @param set A specification of the test to be made. The tests are of
#'     the form u and v are independent condionally on S where u and v
#'     are variables and S is a set of variables. See 'details' for
#'     details about specification of \code{set}.
#' @param \dots Additional arguments.
#' @return An object of class `citest` (which is a list).
#' @details
#'
#' \code{set} can be 1) a vector or 2) a right-hand sided formula in which
#' variables are separated by '+'. In either case, it is tested if the first
#' two variables in the \code{set} are conditionally independent given the
#' remaining variables in \code{set}.  (Notice an abuse of the '+' operator in
#' the right-hand sided formula: The order of the variables does matter.)
#' 
#' If \code{set} is \code{NULL} then it is tested whether the first two
#' variables are conditionally independent given the remaining variables.
#' 
#' If \code{set} consists only of factors then \code{x[,set]} is converted to a
#' contingency table and the test is made in this table using
#' \code{ciTest_table()}.
#' 
#' If \code{set} consists only of numeric values and integers then
#' \code{x[,set]} is converted to a list with components \code{cov} and
#' \code{n.obs} by calling \code{cov.wt(x[,set], method='ML')}. This list is
#' then passed on to \code{ciTest_mvn()} which makes the test.
#' 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#'
#' @seealso \code{\link{ciTest}}, \code{\link{ciTest_table}},
#'     \code{\link{ciTest_mvn}}, \code{\link{chisq.test}}
#' @keywords htest
#' @examples
#' 
#' data(milkcomp1)
#' ciTest(milkcomp1, set=~tre + fat + pro)
#' ciTest_df(milkcomp1, set=~tre + fat + pro)

#' @export 
ciTest_df <- function(x, set=NULL,...){
    
    if (is.null(set)) set <- names(x)
    else {
        if (inherits(set, c("formula", "character"))){
            set <- unlist(rhsFormula2list(set), use.names = FALSE)
            set <- names(x)[pmatch(set, names(x))] 
        }
    }

    wdata       <- x[ , set]
    varTypes    <- unique.default(unlist(lapply(wdata, class), use.names=FALSE))

    has.factor  <- "factor" %in% varTypes
    has.numeric <- any(c("integer", "numeric") %in% varTypes)
    switch.code <- as.character(1 * has.factor + 2 * has.numeric)

    switch(switch.code,
           "0"={ ## F & F
               stop("Strange error...\n")},
           "1"={ ## T & F
               ciTest_table(xtabs(~., data=wdata),set=set, ...)},
           "2"={ ## F & T
               ciTest_mvn(cov.wt(wdata, method="ML"), set=set, ...)},
           "3"={ ## T & T
               .ciTest_df_internal(wdata, set, ...)}
           )
}

###
### If mixed data,  we test for deleting edge in model
###

### FIXME: set can be e.g. c(2,4,1) and this causes an error.
.ciTest_df_internal <- function(x, set=NULL,...){
    ##cat("CHK: ciTestmixed\n")
    if (is.numeric(set)) set <- names(x)[set]

    obj <- mmod(list(set), data=x)
    ans <- testdelete(obj, set[1:2])
    ans <- ans[c(1,3,2)] ## FIXME: This is fragile
    ans$method   <- "CHISQ"
    ans$statname <- "DEV"
    ans$varNames <- set
    class(ans)   <- "citest"
    ans
}

