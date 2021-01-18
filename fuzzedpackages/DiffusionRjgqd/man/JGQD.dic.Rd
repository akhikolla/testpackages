\name{JGQD.dic}
\alias{JGQD.dic}
\title{
Summarize MCMC Selection Output for a List of JGQD.mcmc or BiJGQD.mcmc Objects.
}
\description{\code{JGQD.dic()}  summarizes the MCMC output from a list of \code{JGQD.mcmc()} objects. This may be used to neatly summarize the MCMC output of various models fitted to a given dataset.
}
\usage{
JGQD.dic(model.list, type = "col")
}

\arguments{
  \item{model.list}{
A list of \code{JGQD.mcmc()} objects.
}
\item{type}{Shoould output be of row (\code{'row'}) or column (\code{'col'}) format.
}
}
\details{
\code{JGQD.dic()} summarizes the output from various models fitted via \code{JGQD.mcmc()}. By ranking them according to DIC. \code{[=]} indicates which model has the minimal DIC.
\tabular{lrrrrrr}{
        \tab \code{Elapsed_Time}\tab \code{Time_Homogeneous}\tab           \code{p}\tab    \code{DIC}\tab     \code{pD}\tab   \code{N}\cr
\code{Model 1} \tab     \code{00:00:47}\tab               \code{No}\tab   \code{5.00     }\tab \code{    389.30}\tab   \code{3.92}\tab 201\cr
\code{Model 2} \tab     \code{00:01:00}\tab               \code{No}\tab   \code{5.00     }\tab \code{[=]386.45}\tab   \code{4.13}\tab 201\cr
\code{Model 3} \tab    \code{ 00:02:50}\tab               \code{No}\tab   \code{5.00     }\tab \code{    391.37}\tab   \code{3.94}\tab 201\cr
}}


\value{
  \item{}{A data frame with summary of model output. See Details.}
}
\references{
Updates available on GitHub at \url{https://github.com/eta21}.
}
\author{
Etienne A.D. Pienaar: \email{etiannead@gmail.com}
}


\seealso{
  \code{\link{JGQD.mcmc}}
}

\examples{
\donttest{

#===============================================================================

 # For detailed notes and examples on how to use the BiJGQD.dic() function, see
 # the following vignette:

 RShowDoc('Part_4_Likelihood_Inference',type='html','DiffusionRjgqd')

#===============================================================================
}
}
\keyword{deviance information criterion (DIC)}
