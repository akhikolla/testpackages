\name{shapereg}
\alias{shapereg}
\title{Shape-Restricted Regression}
\description{The regression model \eqn{y_i = f(t_i) + x_i'\beta + \varepsilon_i, i = 1,\ldots,n} is considered, where the only assumptions about \eqn{f} concern its shape. The vector expression for the model is \eqn{y = \theta + X\beta + \varepsilon}. \eqn{X} represents a parametrically modelled covariate, which could be a categorical covariate or a linear term. The shapereg function allows eight shapes: increasing, decreasing, convex, concave, increasing-convex, increasing-concave, decreasing-convex, and decreasing-concave. This routine employs a single cone projection to find \eqn{\theta} and \eqn{\beta} simultaneously.  }

\usage{
shapereg(formula, data = NULL, weights = NULL, test = FALSE, nloop = 1e+4)
}
\arguments{
  \item{formula}{  A formula object which gives a symbolic description of the model to be fitted. It has the form "response ~ predictor". The response is a vector of length \eqn{n}.  A predictor can be a non-parametrically modelled variable with a shape restriction or a parametrically modelled unconstrained covariate. In terms of a non-parametrically modelled predictor, the user is supposed to indicate the relationship between \eqn{E(y)} and a predictor \eqn{t} in the following way:
    \itemize{
    \item{incr(t): }{\eqn{E(y)} is increasing in \eqn{t}. See \code{\link{incr}} for more details.}
    \item{decr(t): }{\eqn{E(y)} is decreasing in \eqn{t}.  See \code{\link{decr}} for more details.}
    \item{conc(t): }{\eqn{E(y)} is concave in \eqn{t}.  See \code{\link{conc}} for more details.}
    \item{conv(t): }{\eqn{E(y)} is convex in \eqn{t}.  See \code{\link{conv}} for more details.}
    \item{incr.conc(t): }{\eqn{E(y)} is increasing and concave in \eqn{t}.  See \code{\link{incr.conc}} for more details.}
    \item{decr.conc(t): }{\eqn{E(y)} is decreasing and concave in \eqn{t}. See \code{\link{decr.conc}} for more details.}
    \item{incr.conv(t): }{\eqn{E(y)} is increasing and convex in \eqn{t}. See \code{\link{incr.conv}} for more details.}
    \item{decr.conv(t): }{\eqn{E(y)} is decreasing and convex in \eqn{t}. See \code{\link{decr.conv}} for more details.}
    }
}
  \item{data}{An optional data frame, list or environment containing the variables in the model. The default is data = NULL.}
  \item{weights}{  An optional non-negative vector of "replicate weights" which has the same length as the response vector. If weights are not given, all weights are taken to equal 1. The default is weights = NULL.}
  \item{test}{The test parameter given by the user.}
  \item{nloop}{ The number of simulations used to get the p-value for the \eqn{E_{01}} test. The default is 1e+4. }
}


\details{
This routine constrains \eqn{\theta} in the equation \eqn{y = \theta + X\beta + \varepsilon} by a shape parameter. 

The constraint cone \eqn{C} has the form \eqn{\{\phi: \phi = v + \sum b_i\delta_i, i = 1,\ldots,m, b_1,\ldots, b_m \ge 0 \}}, \eqn{v} is in \eqn{V}. The column vectors of \eqn{X} are in \eqn{V}, i.e., the linear space contained in the constraint cone. 

The hypothesis test \eqn{H_0: \phi} is in \eqn{V} versus \eqn{H_1: \phi} is in \eqn{C} is an exact one-sided test, and the test statistic is \eqn{E_{01} = (SSE_0 - SSE_1)/(SSE_0)}, which has a mixture-of-betas distribution when \eqn{H_0} is true and \eqn{\varepsilon} is a vector following a standard multivariate normal distribution with mean 0. The mixing parameters are found through simulations. The number of simulations used to obtain the mixing distribution parameters for the test is 10,000. Such simulations usually take some time. For the "feet" data set used as an example in this section, whose sample size is 39, the time to get a p-value is roughly between 4 seconds. 

This routine calls coneB for the cone projection part. 
}

\value{
\item{coefs}{The estimated coefficients for \eqn{X}, i.e., the estimation for the vector \eqn{\beta}. Note that even if the user does not provide a constant vector in \eqn{X}, the coefficient for the intercept will be returned.}

\item{constr.fit}{The shape-restricted fit over the constraint cone \eqn{C} of the form \eqn{\{\phi: \phi = v + \sum b_i\delta_i, 
i = 1,\ldots,m, b_1,\ldots, b_m \ge 0 \}}, \eqn{v} is in \eqn{V}.}

\item{linear.fit}{The least-squares regression of \eqn{y} on \eqn{V}, i.e., the linear space contained in the constraint cone. If shape is 3 or shape is 4, \eqn{V} is spanned by \eqn{X} and \eqn{t}. Otherwise, it is spanned by \eqn{X}. \eqn{X} must be full column rank, and the matrix formed by combining \eqn{X} and \eqn{t} must also be full column rank.}

\item{se.beta}{The standard errors for the estimation of the vector \eqn{\beta}. The degree of freedom is returned by coneB and is multiplied by 1.5. Note that even if the user does not provide a constant vector in \eqn{X}, the standard error for the intercept will be returned.}

\item{pval}{The p-value for the hypothesis test \eqn{H_0: \phi} is in \eqn{V} versus \eqn{H_1: \phi} is in \eqn{C}. \eqn{C} is the constraint cone of the form \eqn{\{\phi: \phi = v + \sum b_i\delta_i, i = 1,\ldots,m, b_1,\ldots, b_m \ge 0 \}}, \eqn{v} is in \eqn{V}, and \eqn{V} is the linear space contained in the constraint cone. If test == TRUE, a p-value is returned. Otherwise, the test is skipped and no p-value is returned.}

\item{pvals.beta}{The approximate p-values for the estimation of the vector \eqn{\beta}. A t-distribution is used as the approximate distribution. Note that even if the user does not provide a constant vector in \eqn{X}, the approximate p-value for the intercept will be returned. }

\item{test}{The test parameter given by the user.}

\item{SSE0}{The sum of squared residuals for the linear part.}

\item{SSE1}{The sum of squared residuals for the full model.}

\item{shape}{A number showing the shape constraint given by the user in a shapereg fit.}

\item{tms}{The terms objects extracted by the generic function \emph{terms} from a shapereg fit. See the official help page (\url{http://stat.ethz.ch/R-manual/R-patched/library/stats/html/terms.html}) of the \emph{terms} function for more details.}

 \item{zid}{A vector keeping track of the position of the parametrically modelled covariate.}
  \item{vals}{A vector storing the levels of each variable used as a factor.}
  \item{zid1}{A vector keeping track of the beginning position of the levels of each variable used as a factor.}
  \item{zid2}{A vector keeping track of the end position of the levels of each variable used as a factor.}

\item{tnm}{The name of the shape-restricted predictor.}

\item{ynm}{The name of the response variable.}

\item{znms}{A vector storing the name of the parametrically modelled covariate.}

\item{is_param}{A logical scalar showing if or not a variable is a parametrically modelled covariate, which could be a factor or a linear term.}
\item{is_fac}{A logical scalar showing if or not a variable is a factor.}

\item{xmat}{A matrix whose columns represent the parametrically modelled covariate.}

\item{call}{The matched call.}
}

\references{Raubertas, R. F., C.-I. C. Lee, and E. V. Nordheim (1986) Hypothesis tests for normals
means constrained by linear inequalities. \emph{Communications in Statistics - Theory and 
Methods \bold{15 (9)}}, 2809--2833.

Robertson, T., F. Wright, and R. Dykstra (1988) \emph{Order Restricted Statistical Inference}
New York: John Wiley and Sons.

Fraser, D. A. S. and H. Massam (1989) A mixed primal-dual bases algorithm for regression
under inequality constraints application to concave regression. \emph{Scandinavian Journal of Statistics \bold{16}}, 65--74.

Meyer, M. C. (2003) A test for linear vs convex regression function using shape-restricted
regression. \emph{Biometrika \bold{90(1)}}, 223--232.

Cheng, G.(2009) Semiparametric additive isotonic regression. \emph{Journal of Statistical Planning and Inference \bold{139}}, 1980--1991.

Meyer, M.C.(2013a) Semiparametric additive constrained regression. \emph{Journal of Nonparametric Statistics \bold{25(3)}}, 715--743.  

Liao, X. and M. C. Meyer (2014) coneproj: An R package for the primal or dual cone projections with routines for constrained regression. \emph{Journal of Statistical Software \bold{61(12)}}, 1--22.
}

\author{Mary C. Meyer and Xiyue Liao}

\seealso{
\code{\link{coneB}}
}

\examples{
# load the feet data set
    data(feet)

# extract the continuous and constrained predictor
    l <- feet$length

# extract the continuous response
    w <- feet$width

# extract the categorical covariate: sex
    s <- feet$sex

# make an increasing fit with test set as FALSE
    ans <- shapereg(w ~ incr(l) + factor(s))

# check the summary table 
    summary(ans)

# make an increasing fit with test set as TRUE
    ans <- shapereg(w ~ incr(l) + factor(s), test = TRUE, nloop = 1e+3)

# check the summary table 
    summary(ans)

# make a plot comparing the unconstrained fit and the constrained fit
    par(mar = c(4, 4, 1, 1))
    ord <- order(l)
    plot(sort(l), w[ord], type = "n", xlab = "foot length (cm)", ylab = "foot width (cm)")
    title("Shapereg Example Plot")

# sort l according to sex
    ord1 <- order(l[s == "G"])
    ord2 <- order(l[s == "B"])

# make the scatterplot of l vs w for boys and girls
    points(sort(l[s == "G"]), w[s == "G"][ord1], pch = 21, col = 1)
    points(sort(l[s == "B"]), w[s == "B"][ord2], pch = 24, col = 2)

# make an unconstrained fit to boys and girls
    fit <- lm(w ~ l + factor(s))

# plot the unconstrained fit 
    lines(sort(l), (coef(fit)[1] + coef(fit)[2] * l + coef(fit)[3])[ord], lty = 2)
    lines(sort(l), (coef(fit)[1] + coef(fit)[2] * l)[ord], lty = 2, col = 2)
    legend(21.5, 9.8, c("boy","girl"), pch = c(24, 21), col = c(2, 1)) 

# plot the constrained fit
    lines(sort(l), (ans$constr.fit - ans$linear.fit + coef(ans)[1])[ord], col = 1)
    lines(sort(l), (ans$constr.fit - ans$linear.fit + coef(ans)[1] + coef(ans)[2])[ord], col = 2)
}
\keyword{cone projection routines}

