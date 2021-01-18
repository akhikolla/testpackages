\name{GPvam}
\alias{GPvam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting the Generalized and Variable Persistence Models
}
\description{
An EM algorithm, Karl et al. (2013) <doi:10.1016/j.csda.2012.10.004>, is used to estimate the generalized, variable, and complete persistence models, Mariano et al. (2010) <doi:10.3102/1076998609346967>. These are multiple-membership linear mixed models with teachers modeled as "G-side" effects and students modeled with either "G-side" or "R-side" effects.
}
\usage{
GPvam(vam_data, fixed_effects = formula(~as.factor(year) + 0), 
   student.side = "R", persistence="GP", max.iter.EM = 1000, tol1 = 1e-07, 
   hessian = FALSE, hes.method = "simple", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vam_data}{
a data frame that contains at least a column "y" containing the student scores, a 
column "student" containing unique student ID's, a column "teacher" containing the 
teacher ID's, and a column "year" which contains the year (or semester, etc.) of the time period. The "y" and "year" variables needs to be numeric. If other variables are to be included 
as fixed effects, they should also be included in vam_data. See 'Note' for further discussion.
}
  \item{fixed_effects}{
an object of class \code{formula} describing the structure of the fixed effects. Categorical variables
should be wrapped in an \code{as.factor} statement.
}
  \item{student.side}{
a character. Choices are "G" or "R". See section 'Details'.
}
  \item{persistence}{
a character. Choices are "GP", "rGP", "VP", "CP", or "ZP". Only "GP" is currently 
compatible with student.side="G". See section 'Details'.
}


  \item{max.iter.EM}{
the maximum number of EM iterations
}
  \item{tol1}{
convergence tolerance for EM algorithm. The convergence criterion is specified under 'Details'.
}
  \item{hessian}{
logical indicating whether the Hessian of the variance parameters (and persistence parameters for persistence="VP") should be calculated after convergence of the EM algorithm. Standard errors for the fixed and EBLUPs are calculated by default.
}
  \item{hes.method}{
a character string indicating the method of numerical differentiation used
            to calculate the Hessian of the variance parameters. Options are "simple" or "richardson".
}
  \item{verbose}{
logical. If TRUE, model information will be printed at each iteration.
}
}
\details{
The design for the random teacher effects according to the generalized persistence model
of Mariano et al. (2010) is built into the function. The model includes correlated current- and
future-year effects for each teacher. By setting \code{student.side="R"}, the intra-student correlation is 
modeled via an unstructured, block-diagonal error covariance matrix, as specified by Mariano et al. (2010). Setting 
\code{student.side="G"} keeps the same teacher structure, but models intra-student correlation via 
random student effects. This is similar to the model used by McCaffrey and Lockwood (2011), and is appropriate when the testing 
scale is the same across years. In this case, the error covariance matrix is diagonal, although a separate variance is calculated for each year.
From a computational perspective, the model estimating the R-side student effects has better scalability properties, although the G-side function is faster (Karl et al. 2012).

The \code{persistence} option determines the type of persistence effects that are modeled. The generalized persistence model ("GP") is described above. When \code{student.side="R"}, other models for teacher persistence are available. The reduced GP model  ("rGP", Karl et al. 2012) combines each teacher's future year effects from the GP model into a single effect. The variable persistence model ("VP") assumes that teacher effects in future years are multiples of their effect in the current year (Lockwood et al. 2007). The multipliers in the VP model are called persistence parameters, and are estimated. By contrast, the complete ("CP") and zero ("ZP") persistence models fix the persistence parameters at 1 and 0, respectively (Lockwood et al. 2007). 

Convergence is declared when \eqn{(l_k-l_{k-1})/l_k < 1E-07}, where \eqn{l_k} is the log-likelihood 
at iteration \eqn{k}.


The model is estimated via an EM algorithm. For details, see Karl et al. (2012). The model was estimated through Bayesian computation in Mariano et al. (2010).

Note: When student.side="R" is selected, the first few iterations of the EM algorithm will take longer than subsequent iterations. This is a result of the hybrid gradient-ascent/Newton-Raphson method used in the M-step for the R matrix in the first two iterations (Karl et al. 2012).

\bold{Program run time and memory requirements:}
The data file GPvam.benchmark that is included with the package contains runtime and peak memory requirements for different persistence settings, using simulated data sets with different values for number of years, number of teachers per year, and number of students per teacher. These have been multiplied to show the total number of teachers in the data set, as well as the total number of students. With student.side="R", the persistence="GP" model is most sensitive to increases in the size of the data set. With student.side="G", the memory requirements increase exponentially with the number of students and teachers, and that model should not be considered scalable to extremely large data sets.

All of these benchmarks were performed with Hessian=TRUE. Calculation of the Hessian accounts for anywhere from 20\% to 75\% of those run times. Unless the standard errors of the variance components are needed, leaving Hessian=FALSE will lead to a faster run time with smaller memory requirements.

}
\value{
GPvam returns an object of class \code{GPvam}



An object of class \code{GPvam} is a list containing the following components:

\item{loglik}{the maximized log-likelihood at convergence of the EM algorithm}
\item{teach.effects}{a data frame containing the predicted teacher effects and standard errors}
\item{parameters}{a matrix of estimated model parameters and standard errors}
\item{Hessian}{if requested, the Hessian of the variance parameters}
\item{R_i}{(only when \code{student_side} is set to 'R') a matrix containing the error covariance matrix of a student}
\item{teach.cov}{a list containing the unique blocks of the covariance matrix of teacher effects}
\item{mresid}{a vector of the raw marginal residuals}
\item{cresid}{a vector of the raw conditional residuals}
\item{sresid}{a vector of the scaled conditional residuals}
\item{yhat}{a vector of the predicted values}

The function  \code{\link{summary}} provides a summary of the results. This includes the estimated model parameters and standard errors, along with the correlation matrices corresponding to the estimated correlation matrices. Summary information about scaled and raw residuals is reported.
}
\references{
Karl, A., Yang, Y. and Lohr, S. (2013) Efficient Maximum Likelihood Estimation of Multiple Membership Linear Mixed Models, with an Application to Educational Value-Added Assessments \emph{Computational Statistics & Data Analysis} \bold{59}, 13--27.  

Karl, A., Yang, Y. and Lohr, S. (2014) Computation of Maximum Likelihood Estimates for Multiresponse Generalized Linear Mixed Models with Non-nested, Correlated Random Effects \emph{Computational Statistics & Data Analysis} \bold{73}, 146--162.

Karl, A., Yang, Y. and Lohr, S. (2014) A Correlated Random Effects Model for Nonignorable Missing Data in Value-Added Assessment of Teacher Effects \emph{Journal of Educational and Behavioral Statistics} \bold{38}, 577--603.

Lockwood, J., McCaffrey, D., Mariano, L., Setodji, C. (2007) Bayesian Methods for Scalable Multivariate Value-Added Assesment. \emph{Journal of Educational and Behavioral Statistics} \bold{32}, 125--150.

Mariano, L., McCaffrey, D. and Lockwood, J. (2010) A Model for Teacher Effects From Longitudinal Data Without Assuming Vertical Scaling. \emph{Journal of Educational and Behavioral Statistics} \bold{35}, 
253--279.

McCaffrey, D. and Lockwood, J. (2011) Missing Data in Value-Added Modeling of Teacher Effects," \emph{Annals of Applied Statistics} \bold{5}, 773--797

}



\author{
Andrew Karl \email{akarl@asu.edu}, Yan Yang, Sharon Lohr
}
\note{
The model assumes that each teacher teaches only one year. If, for example, a teacher teaches in years 1 and 2, his/her first year performance is modeled independently of the second year performance.
To keep these effects separate, the progam appends "(year i)" to each teacher name, where i is the year in which the teacher taught.

The \code{fixed_effects} argument of \code{GPvam} utilizes the functionality of R's \code{formula} class. In the statement \code{fixed_effects=formula(~as.factor(year)+cont_var+0))}, \code{as.factor(year)} identifies \code{year} as a categorical variable. \code{+0} indicates that no intercept is to be fitted, and \code{+cont_var} indicates that a seperate effect is to be fitted for the continuous variable "cont_var." An interaction between "year" and "cont_var" could be specified by \code{~as.factor(year)*cont_var+0}, or equivalently, \code{~as.factor(year)+cont_var+as.factor(year):cont_var+0}. See \code{\link{formula}} for more details.

When applied to an object of class GPvam, \code{\link{plot.GPvam}} returns a caterpillar plot for each effect, as well as residual plots.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{plot.GPvam}},
 \code{\link{summary.GPvam}},
 \code{\link{vam_data}}
}
\examples{
data(vam_data)
GPvam(vam_data,student.side="R",persistence="CP",
fixed_effects=formula(~as.factor(year)+cont_var+0),verbose=TRUE,max.iter.EM=1)
\donttest{
result <- GPvam(vam_data,student.side="R",persistence="VP",
fixed_effects=formula(~as.factor(year)+cont_var+0),verbose=TRUE)
 summary(result)

 plot(result)
}
}

\keyword{ regression }

