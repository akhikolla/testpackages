\name{plot}
\alias{plot.GPvam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for GPvam
}
\description{
Plot teacher effects and residuals.  The caterpillar plots use a modified version of the plotCI function from R package gplots. According to that package, "Original version [of plotCI] by Bill Venables wvenable@attunga.stats.adelaide.edu.au posted to r-help on Sep. 20, 1997. Enhanced version posted to r-help by Ben Bolker ben@zoo.ufl.edu on Apr. 16, 2001. This version was modified and extended by Gregory R. Warnes greg@warnes.net. Additional changes suggested by Martin Maechler maechler@stat.math.ethz.ch integrated on July 29, 2004."
}
\usage{
\method{plot}{GPvam}(x, ..., alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class GPvam
}
  \item{\dots}{
other arguments
}
  \item{alpha}{
the significance level for the caterpillar plots
}

}

\value{
Requires user to click window or press "enter" to progress through plots. Returns caterpillar plots (via the package gplots) and residual plots. 
}
\references{
Karl, A., Yang, Y. and Lohr, S. (2013) Efficient Maximum Likelihood Estimation of Multiple Membership Linear Mixed Models, with an Application to Educational Value-Added Assessments \emph{Computational Statistics & Data Analysis} \bold{59}, 13--27.  

Karl, A., Yang, Y. and Lohr, S. (2014) Computation of Maximum Likelihood Estimates for Multiresponse Generalized Linear Mixed Models with Non-nested, Correlated Random Effects \emph{Computational Statistics & Data Analysis} \bold{73}, 146--162.

Karl, A., Yang, Y. and Lohr, S. (2014) A Correlated Random Effects Model for Nonignorable Missing Data in Value-Added Assessment of Teacher Effects \emph{Journal of Educational and Behavioral Statistics} \bold{38}, 577--603.

Lockwood, J., McCaffrey, D., Mariano, L., Setodji, C. (2007) Bayesian Methods for Scalable Multivariate Value-Added Assesment. \emph{Journal of Educational and Behavioral Statistics} \bold{32}, 125--150. 

Mariano, L., McCaffrey, D. and Lockwood, J. (2010) A Model for Teacher Effects From Longitudinal Data Without Assuming Vertical Scaling. \emph{Journal of Educational and Behavioral Statistics} \bold{35}, 
253--279.

McCaffrey, D. and Lockwood, J. (2011) Missing Data in Value-Added Modeling of Teavher Effects, \emph{Annals of Applied Statistics} \bold{5}, 773--797
}
\author{
Andrew Karl \email{akarl@asu.edu}
Yan Yang
Sharon Lohr

Other authors as listed above for the caterpillar plots.
}


\seealso{
 \code{\link{summary.GPvam}}
}
\examples{
data(vam_data)
\donttest{
GPvam(vam_data,student.side="R",persistence="VP",
fixed_effects=formula(~as.factor(year)+cont_var+0),verbose=TRUE,max.iter.EM=1)
\donttest{
result <- GPvam(vam_data,student.side="R",persistence="VP",
fixed_effects=formula(~as.factor(year)+cont_var+0),verbose=TRUE)
 summary(result)

 plot(result)
}                
}
}
\keyword{ regression }

