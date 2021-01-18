\name{RealVAMS}
\alias{RealVAMS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate VAM Fitting
}
\description{
Fits a multivariate value-added model (VAM), see Broatch, Green, and Karl (2018) \cr <doi:10.32614/RJ-2018-033>, and Broatch and Lohr (2012) \cr <doi:10.3102/1076998610396900>, with normally distributed test scores and a binary outcome indicator. A pseudo-likelihood approach, Wolfinger and O'Connell (1993) \cr <doi:10.1080/00949659308811554>, is used for the estimation of this joint generalized linear mixed model.  The inner loop of the pseudo-likelihood routine (estimation of a linear mixed model) occurs in the framework of the EM algorithm presented by \cr Karl, Yang, and Lohr (2013) <DOI:10.1016/j.csda.2012.10.004>. This material is based upon work supported by the National Science Foundation under grants DRL-1336027 and DRL-1336265.  
}
\usage{
RealVAMS(score.data, outcome.data, persistence = "CP",school.effects = FALSE, 
REML = TRUE, score.fixed.effects = formula(~as.factor(year) + 0), 
outcome.fixed.effects = formula(~1), max.iter.EM = 10,
outcome.family = binomial(link = "probit"), tol1 = 1e-07, max.PQL.it = 30, 
pconv = .Machine$double.eps*1e9, var.parm.hessian = TRUE, verbose = TRUE,
independent.responses = FALSE,cpp.benchmark=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{score.data}{
a data frame that contains at least a column "y" containing the student scores, 
a column "student" containing unique student ID's, a column "teacher" containing 
the teacher ID's, and a column "year" which contains the year (or semester, etc.) 
of the time period. The "y" and "year" variables needs to be numeric. If other 
variables are to be included as fixed effects in the score model, they should also be included in 
score.data. See 'Note' for further discussion.
}
  \item{outcome.data}{
a data frame that contains at least a column "r" containing the binary student outcomes (coded 0/1),
and a column "student" containing unique student ID's. The student ID's should match those in \code{score.data}.
If other variables are to be included as fixed effects in the outcome model, they should also be included in 
outcome.data.
}
  \item{persistence}{
a character object. Choices are "CP" or "VP", for complete and variable persistence of the teacher score effects, respectively. 
The teacher outcome effects are modeled with complete persistence, regardless of the selection here.
}
\item{school.effects}{
logical. If TRUE, correlated random school-level effects are fitted in the score and outcome response models. For both responses, the
school effects are fit with zero-persistence (a student's score in each year is associated with the current school attended, and their outcome is associated with the last school the student attended). The school ID should be included
as a column \code{schoolID} in the score.data data frame. 
}
  \item{REML}{
logical. If TRUE, the pseudo-response is fit using REML. If FALSE, ML is used.
}
  \item{score.fixed.effects}{
an object of class formula describing the structure of the fixed effects for the student scores. Categorical variables should be wrapped in an as.factor statement.
}
  \item{outcome.fixed.effects}{
an object of class formula describing the structure of the fixed effects for the student outcomes. Categorical variables should be wrapped in an as.factor statement.
}
  \item{max.iter.EM}{
numeric. The maximum number of EM iterations during each pseudo-likelihood iteration
}
  \item{outcome.family}{
an object of class family describing the assumed distribution of the response. \code{binomial} is required, but any link function may be used.
}
  \item{tol1}{
numeric. Convergence tolerance for EM algorithm during each interior pseduo-likelihood iteration. The convergence criterion is specified under 'Details'.
}
  \item{max.PQL.it}{
numeric. Maximum number of outer pseudo-likelihood iterations.
}
  \item{pconv}{
numeric. Convergence criterion for outer pseudo-likelihood iterations. Compare to the PCONV option of SAS PROC GLIMMIX.
}
\item{var.parm.hessian}{
logical. If TRUE, the Hessian of the parameters in the error and random effects covariance matrices
is calculated, providing standard errors for those parameters. Setting this option to FALSE will reduce the
run time of the program: only standard errors for the fixed effects will be returned.
}
  \item{verbose}{
logical. If TRUE, model information will be printed at each iteration.
}
\item{independent.responses}{
logical. If TRUE, this option will model the responses independently by fixing the covariances in G at 0 as well as the 
covariances in the last row/column of R. The resulting estimates are the same as those that would be obtained by modelling 
the test scores in package GPvam (with REML=FALSE) and modelling the binary respones in SAS GLIMMIX (with link=probit) 
RealVAMS has been validated against these programs.
}                                                  
\item{cpp.benchmark}{
logical. If TRUE, this option will perform the calculations shown in equation (16) of 
Karl, Yang, Lohr (2013) using both R and the embedded C++ code to demonstrate the time savings of using C++. 
A summary table is printed at the end.
}
}
\details{

*The \code{persistence} option determines the type of persistence effects that are modeled. 
The variable persistence model ("VP") assumes that teacher effects in future years are multiples of their effect in the current year (Lockwood et al. 2007). 
The multipliers in the VP model are called persistence parameters, and are estimated. By contrast, the complete persistence ("CP") model fixes the persistence 
parameters at 1 and 0 (Lockwood et al. 2007). 

*Convergence is declared for each interior iteration when \eqn{(l_k-l_{k-1})/l_k < \code{tol1}}, where \eqn{l_k} is the log-likelihood 
at iteration \eqn{k}.


*The model is linearized using a pseudo-likelihood approach (Wolfinger 1993) and the resulting multiple membership linear mixed model is estimated 
via an EM algorithm (Karl et al. 2012). 


}
\value{
RealVAMS returns an object of class \code{RealVAMS}


\item{loglik}{the maximized log-likelihood at convergence of the EM algorithm. Warning: Likelihood-ratio tests are not valid with results from a PQL estimation routine.}
\item{teach.effects}{a data frame containing the predicted teacher effects and standard errors}
\item{school.effects}{if \code{school.effects=TRUE}, a data frame containing the predicted school effects and standard errors. Otherwise, NULL.}
\item{parameters}{a matrix of estimated model parameters and standard errors. Wald p-values and 95\% confidence intervals are provided for convenience, but these calculations assume infinite degrees of freedom. If working with a "small" data set, re-calculate these using the appropriate degrees of freedom.}
\item{Hessian}{the Hessian of the variance parameters. \code{var.parm.hessian} must be set to \code{TRUE}.}
\item{R_i}{a matrix containing the error covariance matrix of a student's test scores and outcome (assuming a complete response vector with no missing test scores or outcome indicator). The bottom-right component corresponds to the variance of the binary response, and is fixed at 1.}
\item{teach.cov}{a list containing the unique blocks of the covariance matrix of teacher effects (the G matrix).}
\item{mresid}{a vector of the raw marginal residuals. Can be reproduced with \cr\code{y.combined-X\%*\%fixed.effects}. }
\item{cresid}{a vector of the raw conditional residuals. Can be reproduced with \cr\code{y.combined-X\%*\%fixed.effects-Z\%*\%eblup[,2]}.}
\item{y.combined}{a vector of the pseudo-responses from the final PQL iteration, with score and outcome responses interleaved (see the notes below). The test scores will be the same as those given as an input, but the 0/1 responses for the binary distribution will be the pseudo-responses. The vector y.response.type indicates which response each component corresponds to. For components corresponding to the binary response,
the original response can be obtained from \cr\code{joined.table$y.combined.original}.}
\item{y.combined.hat}{a vector of the predicted values. This vector can be reconstructed with other variables returned in this object as \code{X\%*\%fixed.effects+Z\%*\%eblup[,2]}. The values corresponding to the binary outcome responses can be converted to probabilities (of a 1 response) via the inverse link function \code{outcome.family$linkinv()}. The covariance matrix for this vector of predictions is \code{Z\%*\%G\%*\%t(Z)+R.full}.}
\item{y.response.type}{a vector indicating the type of response in each component of \code{y.combined}}
\item{y.year}{a vector indicating the year in which each component of \code{y.combined} was recorded}
\item{num.obs}{total number of observations (test scores and binary responses)}
\item{num.student}{total number of students included in the data}
\item{num.year}{number of years over which test scores were modeled}
\item{num.teach}{a vector listing the number of teachers in each year}
\item{persistence}{a character vector indicating the persistence structure (VP or CP) used to model the teacher test-score effects}   
\item{persistence_parameters}{a matrix of the persistence parameters. The (i,j)-th component gives the persistence parameter for year-j teachers on year-i scores.}
\item{X}{the fixed effects design matrix of the interleaved score and outcome responses}
\item{Z}{the random effects design matrix of the interleaved score and outcome responses}
\item{G}{the random effects covariance matrix. This matrix is block diagonal and contains the teacher variance components and, if included, school-level variance components.}
\item{C}{the joint covariance matrix of the fixed effect estimates and the predicted random effects, that is, for c(fixed.effects,eblup[,1]). See Henderson (1975) for details. This can be used to estimate the variance of an estimable/predictable function.}
\item{R}{the error covariance matrix, with the variance of the outcome pseudo-responses restricted to 1. See the description for \code{R.full}. In most cases, \code{R.full} should be used instead of \code{R} when analyzing the results of \bold{RealVAMS}.}
\item{R.full}{the error covariance matrix, which is formed as the product \cr
\code{diag(sqrt.w)\%*\%R\%*\%diag(sqrt.w)}. The matrix R assumes a variance of 1 for all of the binomial responses, while R.full includes the variance from the binomial distribution (in Wolfinger (1993), diag(sqrt.w) is called R_{mu}).}
\item{sqrt.W}{vector of weights for the error covariance matrix. See the description for R.full above}
\item{eblup}{a matrix containing the complete random effects vector and associated standard errors. When \code{school.effects=FALSE}, this matrix is identical to the one returned by \cr \code{teach.effects}. When \code{school.effects=TRUE}, this matrix is equal to\cr \code{rbind(teach.effects,school.effects)}. Thus, \code{Z\%*\%eblup[,2]} may be used to return the subject-specific portion of the predictions.}
\item{fixed.effects}{a vector containing the fixed effect estimates. This is a subset of \code{parameters}, and provided for compatibility with \code{X}. That is, \code{X\%*\%fixed.effects} will yield the marginal means.}
\item{joined.table}{a data frame containing the interleaved score and outcome data sets. See the notes below.}
\item{outcome.family}{returns information about the distribution and link function used for the outcomes.}
}
\references{
Broatch, J. and Lohr, S. (2012) <DOI:10.3102/1076998610396900> Multidimensional Assessment of Value Added by Teachers to Real-World Outcomes. \emph{Journal of Educational and Behavioral Statistics} \bold{37}, 256--277.

Broatch, J., Green, J., Karl, A. (2018) <DOI:10.32614/RJ-2018-033> RealVAMS: An R Package for Fitting a Multivariate Value-added Model (VAM). \emph{The R Journal} 10/1, 22--30.

Karl, A., Yang, Y. and Lohr, S. (2013) <DOI:10.1016/j.csda.2012.10.004> Efficient Maximum Likelihood Estimation of Multiple Membership Linear Mixed Models, with an Application to Educational Value-Added Assessments. \emph{Computational Statistics & Data Analysis} \bold{59}, 13--27.  

Karl, A., Yang, Y. and Lohr, S. (2013) <DOI:10.3102/1076998613494819> A Correlated Random Effects Model for Nonignorable Missing Data in Value-Added Assessment of Teacher Effects. \emph{Journal of Educational and Behavioral Statistics} \bold{38}, 577--603.

Karl, A., Yang, Y. and Lohr, S. (2014) <DOI:10.1016/j.csda.2013.11.019> Computation of Maximum Likelihood Estimates for Multiresponse Generalized Linear Mixed Models with Non-nested, Correlated Random Effects. \emph{Computational Statistics & Data Analysis} \bold{73}, 146--162.

Henderson, C.R. (1975) Best linear unbiased estimation and prediction under a selection model. \emph{Biometrics} \bold{31(2)}, 423-447.

Lockwood, J., McCaffrey, D., Mariano, L., Setodji, C. (2007) <DOI:10.3102/1076998606298039> Bayesian Methods for Scalable Multivariate Value-Added Assessment. \emph{Journal of Educational and Behavioral Statistics} \bold{32}, 125--150.

Wolfinger, R. (1993) <DOI:10.1080/00949659308811554> Generalized linear mixed models a pseudo-likelihood approach. \emph{Journal of Statistical Computation and Simulation} \bold{48} 233--243.
}
\author{
Andrew Karl \email{akarl@asu.edu}, Jennifer Broatch, Jennifer Green
}
\note{
The first few iterations of the EM algorithm will take longer than subsequent iterations. This is a result of the hybrid gradient-ascent/Newton-Raphson method used in the M-step for the R matrix in the first two iterations (Karl et al. 2012).


The model assumes that each teacher teaches only one year. If, for example, a teacher teaches in years 1 and 2, his/her first year performance is modeled independently of the second year performance.
To keep these effects separate, the program appends "(year i)" to each teacher name, where i is the year in which the teacher taught.

To fit the model and allow correlation between test scores and outcomes (at both the student and teacher levels), the score and outcome response vectors are interleaved into a single response vector. For example, if there are three years of test scores modeled with a binary outcome indicator, the binary indicator for a student is inserted immediately after that student's test scores. The \code{joined.table} that is returned by RealVAMS shows how this was done for a particular data set. Row \emph{i} of \code{joined.table} corresponds to row \emph{i} of \code{X}, \code{Z}, \code{R.full}, \code{y.combined}, \code{y.combined.hat}, and \code{eblup}.

The \code{fixed.effects} arguments of \code{RealVAMS} utilizes the functionality of R's \code{formula} class. In the statement \cr\code{score.fixed.effects=formula(~as.factor(year)+cont_var+0))}, \code{as.factor(year)}\cr identifies \code{year} as a categorical variable. \code{+0} indicates that no intercept is to be fitted, and \code{+cont_var} indicates that a separate effect is to be fitted for the continuous variable "cont_var." An interaction between "year" and "cont_var" could be specified by \code{~as.factor(year)*cont_var+0}, or equivalently, \code{~as.factor(year)+cont_var+as.factor(year):cont_var+0}. See \code{\link{formula}} for more details.


}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(example.score.data)
data(example.outcome.data)
#The next line exists to show that the function can run and that the package
#installed correctly
RealVAMS(example.score.data,example.outcome.data,max.PQL.it=1,max.iter.EM=2,
var.parm.hessian=FALSE)

\donttest{
res<-RealVAMS(example.score.data,example.outcome.data)}
}
\keyword{regression}

