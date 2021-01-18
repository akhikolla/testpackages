\name{covQualNested}
\alias{covQualNested}
\title{Nested Qualitative Covariance}
\usage{

covQualNested(input = "x",
              groupList = NULL,
              group = NULL,
              nestedLevels = NULL,
              between = "Symm",
              within = "Diag",
              covBet = c("corr", "homo", "hete"), 
              covWith = c("corr", "homo", "hete"),
              compGrad = TRUE,
              contrasts = contr.helmod,
              intAsChar = TRUE)

}
\arguments{

  \item{input}{
    
    Name of the input, i.e. name of the column in the
    data frame when the covariance kernel is evaluated with the
    \code{\link{covMat,covQual-method}} method.
    
  }
  \item{groupList}{

    A list giving the groups, see \bold{Examples}. Groups of
    size 1 are accepted. Note that the group values should be given in
    some order, with no gap between repeated values, see \bold{Examples}.
  }

  \item{group}{

    Inactive if \code{groupList} is used. 
    A factor or vector giving the groups, see \bold{Examples}. Groups of
    size 1 are accepted. Note that the group values should be given in
    some order, with no gap between repeated values, see
    \bold{Examples}.

  }
  \item{nestedLevels}{
    
    Inactive if \code{groupList} is used. 
    A factor or a vector giving the (nested) levels within the group for
    each level of \code{group}. If this is missing, each element of
    \code{group} is assumed to correspond to one nested level within the
    group and the levels within the group are taken as integers in the order
    of \code{group} elements.
    
  }
  \item{between}{
    
    Character giving the type of structure to use for the \emph{between}
    part. For now this can be one of the three choices \code{"Diag"},
    the diagonal structure of \code{\link{q1Diag}}, \code{"Symm"} for
    the general covariance of \code{\link{q1Symm}}, or \code{"CompSymm"}
    for the Compound Symmetry covariance of \code{\link{q1CompSymm}}.
    Default is \code{Symm}, corresponding to a specific correlation
    value for each pair of groups. On the other hand, \code{Diag}
    corresponds to a common correlation value for all pairs of groups.

  }
  \item{within}{
    
    Character vector giving the type of structure to use for the
    \emph{within} part. The choices are the same as for
    \code{between}. The character vector is recycled to have length
    \eqn{G} so the \emph{within} covariances can differ across groups.
    Default is \code{"Diag"}, corresponding to a compound symmetry matrix.
    
  }
  \item{covBet}{
  }
  \item{covWith}{
    
    Character vector indicating the type of covariance matrix to be used
    for the generator between- or within- matrices, 
    as in \code{\link{q1Diag}}, \code{\link{q1Symm}} or \code{\link{q1CompSymm}}:
    correlation ("corr"), homoscedastic ("homo") or heteroscedastic ("hete").
    Partial matching is allowed. 
    This is different from the form of the resulting covariance matrix, see section
    \bold{Caution}.
    
  }
  \item{compGrad}{
    
    Logical.
    
  }
  \item{contrasts}{
    
    Object of class \code{"function"}. This function is similar to the
    \code{\link[stats]{contr.helmert}} or
    \code{\link[stats]{contr.treatment}} functions, but it must return
    an \emph{orthogonal} matrix. For a given integer \code{n}, it
    returns a matrix with \code{n} rows and \code{n - 1} columns forming
    a basis for the supplementary of a vector of ones in the
    \eqn{n}-dimensional Euclidean space. The \code{\link{contr.helmod}}
    can be used to obtain an orthogonal matrix hence defining an
    orthonormal basis.

  } \item{intAsChar}{

    Logical. If \code{TRUE} (default), an integer-valued input will be
    coerced into a character. Otherwise, it will be coerced into a
    factor.

  }

}
\value{
  
  An object with class \code{"covQualNested"}.

}
\description{

  Nested Qualitative Covariance

}
\note{

  For now the replacement method such as \code{'coef<-'} are inherited
  from the class \code{covQuall}. Consequently when these methods are
  used they do not update the covariance structure in the \code{between}
  slot nor those in the \code{within} (list) slot.
  
  This covariance kernel involves \code{two} categorical (i.e. factor)
  inputs, but these are nested. It could be aliased in the future as
  \code{q1Nested} or \code{q2Nested}.
  
}
\section{Caution}{

  When \code{covBet} and \code{covWith} are zero, the resulting matrix
  \emph{is not a correlation matrix}, due to the mode of
  construction. The "between" covariance matrix is a correlation but
  diagonal blocks are added to the extended matrix obtained by re-sizing
  the "between" covariance into a \eqn{n \times n}{ n * n} matrix.
  
}

\examples{

### Ex 1. See the vignette "groupKernel" for an example 
###       inspired from computer experiments.

### Ex 2. Below an example in data analysis.

country <- c("B", "B", "B", "F", "F" ,"F", "D", "D", "D")
cities <- c("AntWerp", "Ghent" , "Charleroi", "Paris", "Marseille",
            "Lyon", "Berlin", "Hamburg", "Munchen")
myGroupList <- list(B = cities[1:3],
                    F = cities[4:6],
                    D = cities[7:9])

## create a nested covariance. 
# first way, with argument 'groupList':

nest1 <- covQualNested(input = "ccities",
                       groupList = myGroupList,
                       between = "Symm", within = "Diag",
                       compGrad = TRUE,
                       covBet = "corr", covWith = "corr")

# second way, with arguments 'group' and 'nestedLevels'
nest2 <- covQualNested(input = "ccities",
                       group = country, nestedLevels = cities,
                       between = "Symm", within = "Diag",
                       compGrad = TRUE,
                       covBet = "corr", covWith = "corr")


## 'show' and 'plot' method as automatically invocated
nest2
plot(nest2, type = "cor")

## check that the covariance matrices match for nest1 and nest2
max(abs(covMat(nest1) - covMat(nest2)))


## When the groups are not given in order, an error occurs!

countryBad <- c("B", "B", "F", "F", "F", "D", "D", "D", "B")
cities <- c("AntWerp", "Ghent", "Paris", "Marseille", "Lyon",
            "Berlin", "Hamburg", "Munchen", "Charleroi")

nestBad <- try(covQualNested(input = "ccities",
                             group = countryBad, nestedLevels = cities,
                             between = "Symm", within = "Diag",
                             compGrad = TRUE,
                             covBet = "corr", covWith = "corr"))


}
