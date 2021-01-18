\name{ojaMedian}
\alias{ojaMedian}
\alias{ojaMedianEvo}
\alias{ojaMedianGrid}
\alias{ojaMedianEx}
\alias{ojaMedianExB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Oja Median}
\description{
  Function to compute the Oja median. Several algorithms are possible.
}
\usage{
ojaMedian(X, alg = "evolutionary", sp = 1, na.action = na.fail, 
          control = ojaMedianControl(...), ...)
          
ojaMedianEvo(X, control = ojaMedianControl(...), ...)
ojaMedianGrid(X, control = ojaMedianControl(...), ...)
ojaMedianEx(X, control = ojaMedianControl(...), ...)
ojaMedianExB(X, control = ojaMedianControl(...), ...)
}

\arguments{
  \item{X}{numeric data.frame or matrix.}
  \item{alg}{character string denoting the algorithm to be used for computing the Oja median. Options are
                \code{"exact"}, \code{"bounded_exact"}, \code{"evolutionary"} and \code{"grid"}. Default is \cr 		\code{"evolutionary"}. See Details.}
  \item{sp}{number of runs to average over.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail. } 
  \item{control}{a list specifying the control parameters of the different algorithms; use the function \code{\link{ojaMedianControl}} and see its help page.}
  \item{\dots}{can be used to specify control parameters directly instead of via \code{control}.}
}
\details{
  There are four possible algorithms to calculate the Oja median. The exact algorithm uses a gradient method. It follows intersection lines of hyperplanes until it
  reaches the minimum of an objective function. It is computationally a very intensive algorithm and it calculates the Oja median in acceptable time in the bivariate case for at least 1200
  datapoints. For a 7-dimensional dataset it is possible to calculate it for 24 datapoints.
  
  The bounded exact algorithm modifies the exact algorithm by employing bounded regions which contain the median. The regions are built using the centered rank function. The new algorithm is faster and has less complexity. 
  Parameter \code{volume} is the desired size of the bounded region, which is selected as a part of the original volume. Here the volume is calculated as the volume of a minimal multivariate circumscribed rectangle with edges parallel to the coordinate axes.
  Setting parameter \code{boundedExact} to \code{FALSE} stops the algorithm after the bounded region is found, and its center is reported as an approximation of the median.
  
  With the evolutionary algorithm it is possible to calculate an approximative solution. It starts with a random point and mutates this temporary best solution in order to gain a better one. There are several options
  to control the mutation process. If you are interested in a fast calculation of the Oja median and you tolerate a higher error rate, you should set \code{sigmaAdaption} to 1. As a second
  possibility you could limit the number of subsets used to a small number. If you use all subsets, there are in total \eqn{n} choose \eqn{k}, with \eqn{n} number of datapoints and \eqn{k} dimensions. 
  If you are interested in a precise solution, the following options have turned out to be useful:
  \code{initialSigma}: 0.5, \code{sigmaAdaptation}: 20, \code{adaptationFactor}: 0.5, \code{sigmaLog20Decrease}: 10.
  Tests have been made in the bivariate case, but these values should work for every dimension. 
  In the bivariate case it is possible to calculate the Oja median for more than \eqn{22*10^6} datapoints. In the 10-dimensional case the algorithm is still able to calculate an approximative solution
  for \eqn{10^6} datapoints. 
  Before the algorithm starts itself we transform the data with ICS in order to get a more stable version (with respect to the location of the data) and improve the quality of the approximation.
  Another reason for this was to get an affine invariant way of the approximation. 
  
  The fourth algorithm calculates the Oja median by means of a grid. The grid points are possible approximations of the Oja median. Every grid point is tested to be
  the Oja median. If the test results are not unique the algorithm will take a bigger sample of subsets into account and test it again. In comparison to the evolutionary algorithm it is
  slower and less precise. Only in special data situations it might be useful. The algorithm constitutes an earlier heuristical solution to the Oja median problem and is included mainly for historical reasons.
  
  The exact algorithm and the grid algorithm are also described in Ronkainen et al. (2002). The bounded search algorithm is described in Mosler and Pokotylo (2015).

  A lot of calculation time in the \code{ojaMedian} function might be spend for checking the input and for transforming it. So if you do time-critical calculations, e.g. with loops, you might want to take the variants \code{ojaMedianEx},  \code{ojaMedianExB}, 
  \code{ojaMedianEvo} or \code{ojaMedianGrid}. Please use this only if you know what you are doing, because there are no checks, just the \code{.Call} to the algorithm itself.
  
  If the dimension of your data is too big or if there are too many observations, it is possible that the exact algorithm will crash R. On a common PC with a 32-bit operating system the following
  combinations of dimension:amount will work fine: 2:1200, 3:300, 4:100, 5:63, 6:38, 7:24. Bigger datasets might be possible, depending on your system. 
  
  Another general restriction with this function is that there should be more data points than dimensions.
  
  There is a demo available which demonstrates graphically the Oja median in simple data situations in the bivariate case. To view the demo run \code{demo(ojaMedianDemo)}.
}
\value{
 a numeric vector containing the Oja median.
}
\references{
\cite{Oja, H. (1983), Descriptive statistics for multivariate distributions, \emph{Statistics and Probability Letters}, \bold{1}, 327--332. }

\cite{Ronkainen, T., Oja, H. and Orponen, P. (2002), Computation of the multivariate Oja median, in Dutter R.,  Filzmoser P.,Gather U. and Rousseeuw, P. J.: \emph{Developments in Robust Statistics}, Heidelberg: Springer, 344--359. }

\cite{Fischer, D. (2008), \emph{Diplomarbeit}, Statistische Eigenschaften des Oja-Medians mit einer algorithmischen Betrachtung, Dortmund: Technische \enc{Universität}{Universitat} Dortmund. In German.}

\cite{Mosler, K. and Pokotylo, O. (2015), "Computation of the Oja Median by Bounded Search." \emph{Modern Nonparametric, Robust and Multivariate Methods.} Springer International Publishing, 185--203. }

\cite{Fischer D, Mosler K, \enc{Möttönen}{Mottonen} J, Nordhausen K, Pokotylo O and Vogel D (2020). “Computing the Oja Median in R: The Package OjaNP.” \emph{Journal of Statistical Software}, \bold{92}(8), pp. 1-36. doi: 10.18637/jss.v092.i08 (URL:
http://doi.org/10.18637/jss.v092.i08).}
}

\author{Ported to R by Daniel Fischer. Original C++-code by Thorsten Bernholt, Robin Nunkesser and Tommi Ronkainen.
Bounded search algorithm by Oleksii Pokotylo.
}

\examples{

data(biochem)
X <- as.matrix(biochem[,1:2])
ojaMedian(X)
ojaMedian(X, alg = "evo")
ex <-ojaMedian(X, alg = "exact")
exb<-ojaMedian(X, alg = "bounded_exact")

ojaMedianFn(X, ex)
ojaMedianFn(X, exb)
}
\keyword{multivariate}
\keyword{nonparametric}
