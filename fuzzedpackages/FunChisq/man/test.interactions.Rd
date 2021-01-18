\name{test.interactions}
\alias{test.interactions}
\title{Functional Chi-Squared Test of Functional Dependency among Many Variables in a Data Set
}
\description{
  Apply functional chi-squared tests on many-to-one combinatorial relationships for functional dependency using multivariate discrete data.
}
\usage{
test.interactions(
  x, list.ind.vars, dep.vars, var.names = rownames(x),
  index.kind = c("conditional", "unconditional")
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A numeric matrix or data frame of discrete values. Rows represent variables and columns represent samples. Thus, each row index is a variable index, used by \code{list.ind.vars} and \code{dep.vars}.
}
  \item{list.ind.vars}{
    A list of numeric or integer vectors, each vector representing independent variable indices in one interaction. Each vector (parents) forms a pair with a dependent variable (child) of the same position in \code{dep.vars} to represent a many-to-one directional interaction.
}
  \item{dep.vars}{
    A numeric vector representing indices of dependent variables (children) in multiple interactions.
}
  \item{var.names}{
    Optional. A character vector specifying names of all variables (rows). If not provided, the default is the row names of \code{x}; or \code{1:nrow(x)} if \code{x} does not have row names.
}
  \item{index.kind}{
    A character string to specify the kind of function index to return, identical to the same argument in \code{fun.chisq.test}. The value can be \code{"unconditional"} (default) or \code{"conditional"}.
}
}
\details{
  \code{test.interactions} tests functional dependencies in multiple directional interactions. Each interaction, either one-to-one or many-to-one, is a parents-child pair representing a relationship from independent variables (parents) to a dependent variable (child). The parents-child pairs are specified in two input arguments \code{list.ind.vars} (a list of parents for each interaction) and \code{dep.vars} (vector of children in each interaction).

The function automatically creates contingency tables for interactions of interest, thus convenient to use on multivariate data sets. As the function is implemented in C++ and capable of testing multiple many-to-one interactions in one call, it is much faster than calling the R function \code{fun.chisq.test} multiple times.

\code{test.interactions} implements only the \code{method="fchisq"} option in \code{fun.chisq.test}.

When a contingency table is created for each interaction, all combinations of unique values of the independent variables (parents) form the rows and the unique values of dependent variable (child) form the columns in the contingency table. The table entries are the counts of the corresponding combination of parent and child values. Either rows or columns with all zero counts are removed from the contingency table before functional chi-squared test is applied.
}
\value{
  A data frame with five columns. Each row represents the testing result of each directional interaction. The 1st column is either the indices or names (if \code{var.names} is not \code{NULL}) of independent variables (parents); The 2nd column is the indices or names of the dependent variable (child); The 3rd column named \code{p.value} are p-values; The 4th column named \code{statistic} is chi-squared values; and the 5th column named \code{estimate} is the function indices for each interaction.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Hua Zhong and Joe Song
}
\seealso{
  This function calls functional chi-squared test implemented in C++ and is thus much faster than the R version \code{\link{fun.chisq.test}}.

  For data discretization by optimal univariate \var{k}-means clustering, see \pkg{Ckmeans.1d.dp}.
}

\examples{
x <- matrix(
  c(0,0,1,0,1,
    1,0,2,1,0,
    2,2,0,0,0,
    1,2,1,1,2,
    1,0,2,1,2),
  nrow = 5, ncol = 5, byrow = TRUE)

list.ind.vars <-list(
  c(1),c(1),c(1),
  c(2),c(2),c(2),
  c(1,2), c(2,3),
  c(3,4), c(4,5))
dep.vars <- c(
  3,4,5,
  3,4,5,
  3,4,
  5,1)

# list.ind.vars and dep.vars together specify
#   the following ten interactions:
#   1 -> 3
#   1 -> 4
#   1 -> 5
#   2 -> 3
#   2 -> 4
#   2 -> 5
# 1,2 -> 3
# 2,3 -> 4
# 3,4 -> 5
# 4,5 -> 1

var.names <- paste0("var", 1:5)

test.interactions(
  x = x,
  list.ind.vars = list.ind.vars,
  dep.vars = dep.vars,
  var.names = var.names,
  index.kind = "unconditional")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\keyword{htest}
\keyword{nonparametric}

