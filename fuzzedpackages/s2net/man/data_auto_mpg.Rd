\name{auto_mpg}
\alias{auto_mpg}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Auto MPG Data Set
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset was taken from the UCI Machine Learning Repository \url{https://archive.ics.uci.edu/ml/datasets/Auto+MPG}, and processed for the semi-supervised setting (Ryan and Culp, 2015).
}
\usage{data("auto_mpg")}
\format{
  There are two lists that contain partitions from a data frame with 398 observations on the following 9 variables.
  \describe{
    \item{\code{mpg}}{a numeric vector}
    \item{\code{cylinders}}{an ordered factor with levels \code{3} < \code{4} < \code{5} < \code{6} < \code{8}}
    \item{\code{displacement}}{a numeric vector}
    \item{\code{horsepower}}{a numeric vector}
    \item{\code{weight}}{a numeric vector}
    \item{\code{acceleration}}{a numeric vector}
    \item{\code{year}}{a numeric vector}
    \item{\code{origin}}{a factor}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
This dataset is a slightly modified version of the dataset provided in the StatLib library. In line with the use by Ross Quinlan (1993) in predicting the attribute "mpg", 8 of the original instances were removed because they had unknown values for the "mpg" attribute. "The data concerns city-cycle fuel consumption in miles per gallon, to be predicted in terms of 3 multivalued discrete and 5 continuous attributes." (Quinlan, 1993)
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Quinlan,R. (1993). Combining Instance-Based and Model-Based Learning. In Proceedings on the Tenth International Conference of Machine Learning, 236-243, University of Massachusetts, Amherst. Morgan Kaufmann.

Dua, D. and Graff, C. (2019). UCI Machine Learning Repository [\url{http://archive.ics.uci.edu/ml}]. Irvine, CA: University of California, School of Information and Computer Science.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Ryan, K. J., & Culp, M. V. (2015). On semi-supervised linear regression in covariate shift problems. \emph{The Journal of Machine Learning Research, 16}(1), 3183-3217.
}
\examples{
data(auto_mpg)
head(auto_mpg$P1$xL)
}
\keyword{datasets}
