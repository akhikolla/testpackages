\name{rowimpPrep}

\alias{rowimpPrep}

\title{Missing-data pattern identifier}

\description{\sQuote{rowimpPrep} identifies all missingness patterns
  within an incomplete data set. Running \code{rowimpPrep} is a
  prerequisite for \code{\link{BBPMM.row}}. 
}

\usage{rowimpPrep(data, ID=NULL, verbose=TRUE)}

\arguments{
  \item{data}{Either a data frame or matrix with missing values.}
  \item{ID}{A numeric or character string vector indicating the column
    positions or names of the ID variable (if two data sets were stacked
    that have a joint subset of variables). The first element refers to
    the 'donor ID', the second element refers to the 'recipient
    ID'. This distinction is only of relevance, if the data set is 'L-shaped',
    i.e. if the data contains only one missing-data pattern (where
    incomplete cases are 'recipients'). If \code{ID} has only one
    element, The function assumes that the identifier variables of the
    two data sets are packed into a single variable. Default=NULL is
    used, if no ID variable is specified.}
  \item{verbose}{Prints information on identified
    missing-data patterns. Default=TRUE.}
}

\value{
  \item{data}{The original data set minus the ID variable(s).}
  \item{key}{The ID variable(s) from the original data set.}
  \item{blocks}{A list containing the column positions of all identified
    missing-data patterns.}
  \item{blockNames}{A list containing the variable names corresponding
    to object \code{blocks}.}
  \item{compNames}{A character vector containing the variable names
    of the (completely observed) imputation model variables.}
  \item{ignore}{Contains positions of ignored variables.}
  \item{ignored_data}{Contains ignored variables.}
  \item{indMatrix}{A matrix with the same dimensions as the incomplete data 
   containing flags for missing values.}
}

\details{\code{rowimpPrep} identifies
  all patterns, and allows to decide, whether
  to impute all missing-data patterns with \code{\link{BBPMM.row}}or
  just some of them. This comes in handy if variables that were assumed to be
  completely observed have missing values. These variables are then
  likely to define an unexpected 'block' of their own. Of course,
  \code{BBPMM.row} can be used to impute missing data that are
  not missing-by-design as well, but  \code{\link{BBPMM}} would
  probably be the better option. Note that all variables listed in
  \code{compNames} are used for the imputation model in
  \code{BBPMM.row}, i.e. completely observed variables (ID variables
  aside) which are not to be used in the imputation model, have to be
  removed from the data set beforehand.}

\author{Florian Meinfelder, Thorsten Schnapp [ctb]}

\examples{


### sample data set with non-normal variables and a single
### missingness pattern
set.seed(1000)
n <- 50
x1 <- round(runif(n,0.5,3.5))
x2 <- as.factor(c(rep(1,10),rep(2,25),rep(3,15)))
x3 <- round(rnorm(n,0,3))
y1 <- round(x1-0.25*(x2==2)+0.5*x3+rnorm(n,0,1))
y1 <- ifelse(y1<1,1,y1)
y1 <- ifelse(y1>4,5,y1)
y2 <- y1+rnorm(n,0,0.5)
y3 <- round(x3+rnorm(n,0,2))
data1 <- as.data.frame(cbind(x1,x2,x3,y1,y2,y3))
misrow1 <- sample(n,20)
is.na(data1[misrow1, c(4:6)]) <- TRUE

### preparation step
impblock <- rowimpPrep(data1)

impblock$blockNames

}
