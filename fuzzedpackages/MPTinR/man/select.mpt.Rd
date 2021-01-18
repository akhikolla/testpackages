\name{select.mpt}
\alias{select.mpt}
\title{
Model Selection with MPTinR
}
\description{
This function performs model selection for results produced by MPTinR's \code{\link{fit.mpt}}. It takes multiple results from \code{\link{fit.mpt}} as a list and returns a \code{data.frame} comparing the models using various model selection criteria (e.g., FIA) and AIC and BIC weights. For model selection of multi-dataset fits \code{select.mpt} will additionally count how often each model provided the best fit.
}
\usage{
select.mpt(mpt.results, output = c("standard", "full"), round.digit = 6, dataset)
}

\arguments{
  \item{mpt.results}{
A \code{list} containing results from \code{\link{fit.mpt}}.
}
  \item{output}{
\code{"standard"} or \code{"full"}. If \code{"full"} additionally returns original FIA, AIC, and BIC values, and, for multi-individual fits, compares the model-selection criteria for the aggregated data.
}
  \item{round.digit}{
Integer specifying to which decimal place the results should be rounded. Default is 6. Is also used for rounding FIA, AIC, and BIC values before counting the best fitting values per individual datasets.
}
  \item{dataset}{
Integer vector specifying whether or not to restrict the individual comparison top certain dataset(s). Aggregated results will not be displayed if this argument is present. 
}
}
\details{
\code{select.mpt} is the second major function of MPTinR, next to \code{\link{fit.mpt}}. It takes a list of results produced by \code{fit.mpt} and returns a \code{data.frame} comparing the models using the information criteria obtained by \code{fit.mpt}. That is, if FIA was not obtained for the models, \code{select.mpt} only uses AIC and BIC. We strongly recommend using FIA for model selection (see e.g., Gruenwald, 2000).

The outputs follows the same principle for all information criteria. The lowest value is taken as the reference value and the differences to this value (i.e., the \code{delta}) are reported for all models (e.g., \code{delta.FIA}). If one additionally wants the original values, \code{output} needs to be set to \code{"full"}.

For AIC and BIC, AIC and BIC weights are reported as \code{wAIC} and \code{wBIC} (Wagenmakers & Farrell, 2004).

For multi-individual fit, \code{select.mpt} will additionally return how often each model provided the best fit (e.g., \code{FIA.best}). Values are rounded before determining which is the best fitting model. Note that there can be ties so that two models provide the best fit. Furthermore, if \code{output} is \code{"standard"}, only results for the summed information criteria are returned (indicated by the postfix \code{.sum}). To obtain model selection results for the aggregated data (indicated by postfix \code{.aggregated}), \code{output} needs to be set to \code{"full"}.

select.mpt will check if the data of the results returned from \code{fit.mpt} are equal. (If they are not equal model selection can not be done.)

Note that the values in the returned \code{data.frame} are rounded to the \code{round.digit}th decimal place.

}
\value{
A \code{data.frame} containing the model selection values:\cr
\code{model}: Name or number of model (names are either taken from \code{mpt.results} or obtained via \code{match.call}).\cr
\code{n.parameters}: Number of parameters for each model.\cr
\code{G.Squared}: G.Squared values of the model (from summed fits for multiple datasets).\cr
\code{df}: df values of the model (from summed fits for multiple datasets).\cr
\code{p.value}: p values of the model (from summed fits for multiple datasets).\cr 
\code{p.smaller.05}: How many of the individual data sets have p < .05 (for multiple datasets only).\cr 
For the information criteria (i.e., FIA, AIC, BIC) \code{X}, \code{delta.X}, \code{X.best}, \code{X}, \code{wX} represent: The difference from the reference model, how often each model provided the best fit (only for multi-individual fit), the absolute value, the weights (only AIC and BIC).\cr
For multi-indivudal fit the postfix indicates whether the results refer to the summed information criteria from individual fit \code{.sum} or the information criteria from the aggregated data \code{.aggregated}.
}

\author{
Henrik Singmann
}

\note{
As of March 2015 BIC and FIA are calculated anew if the results are displayed for multiple data sets as BIC and FIA cannot directly be summed across participants due to the \eqn{log(n)} terms in their formula (while AIC can be summed). Instead one first needs to sum the \eqn{G^2} values, \eqn{n}, and the number of parameters, and only then can BIC and FIA be calculated for those summed values.

If any of the models is fitted with \code{fit.aggregated = FALSE} no aggregated results are presented.
}

\references{
Gruenwald, P.D. (2000). Model selection based on minimum description length. \emph{Journal of Mathematical Psychology}, 44, 133-152.

Wagenmakers, E.J. & Farrell, S. (2004). AIC model selection using Akaike weights. \emph{Psychonomic Bulletin & Review}, 11, 192-196.
}

\seealso{
\code{\link{fit.mpt}} for obtaining the results needed here and an example using multi-individual fit and FIA.
}
\examples{

# This example compares the three versions of the model in 
# Riefer and Batchelder (1988, Figure 2)

data(rb.fig2.data)
model2 <- system.file("extdata", "rb.fig2.model", package = "MPTinR")
model2r.r.eq <- system.file("extdata", "rb.fig2.r.equal", package = "MPTinR")
model2r.c.eq <- system.file("extdata", "rb.fig2.c.equal", package = "MPTinR")

# The full (i.e., unconstrained) model
ref.model <- fit.mpt(rb.fig2.data, model2)
# All r equal
r.equal <- fit.mpt(rb.fig2.data, model2, model2r.r.eq)
# All c equal
c.equal <- fit.mpt(rb.fig2.data, model2, model2r.c.eq)

select.mpt(list(ref.model, r.equal, c.equal))



\dontrun{

# Example from Broder & Schutz (2009)

data(d.broeder, package = "MPTinR")
m.2htm <- system.file("extdata", "5points.2htm.model", package = "MPTinR")
r.2htm <- system.file("extdata", "broeder.2htm.restr", package = "MPTinR")
r.1htm <- system.file("extdata", "broeder.1htm.restr", package = "MPTinR")

br.2htm.fia <- fit.mpt(d.broeder, m.2htm, fia = 50000, fit.aggregated = FALSE)
br.2htm.res.fia <- fit.mpt(d.broeder, m.2htm, r.2htm, fia = 50000, fit.aggregated = FALSE)
br.1htm.fia <- fit.mpt(d.broeder, m.2htm, r.1htm, fia = 50000, fit.aggregated = FALSE)

select.mpt(list(br.2htm.fia, br.2htm.res.fia, br.1htm.fia))
# This table shows that the n (number of trials) is too small to correctly compute 
# FIA for the 1HT model (as the penalty for the 1HTM is larger than for the 2HTM, 
# although the former is nested in the latter).
# This problem with FIA can only be overcome by collecting more trials per participant,
# but NOT by collecting more participants (as the penalties are simply summed).

# using the dataset argument we see the same
select.mpt(list(br.2htm.fia, br.2htm.res.fia, br.1htm.fia), dataset = 4, output = "full")

select.mpt(list(br.2htm.fia, br.2htm.res.fia, br.1htm.fia),	dataset = 1:10)
}


}
