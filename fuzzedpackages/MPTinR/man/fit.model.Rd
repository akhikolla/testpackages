\name{fit.model}
\alias{fit.model}
\title{
Fit cognitive models for categorical data using model files
}
\description{
\code{fit.model} fits MPT and other cognitive models for categorical data (e.g., SDT models) that can be specified in a model file.
}
\usage{

fit.model(
	data,
	model.filename, 
	restrictions.filename = NULL, 
	n.optim = 5,
	fia = NULL,
	ci = 95, 
	starting.values = NULL,
	lower.bound = 0,
	upper.bound = 1,
	output = c("standard", "fia", "full"),
	reparam.ineq = TRUE,
	fit.aggregated = TRUE,
	sort.param = TRUE,
	show.messages = TRUE,
	model.type = c("easy", "eqn", "eqn2"),
	multicore = c("none", "individual", "n.optim", "fia"), sfInit = FALSE, nCPU = 2,
	control = list(),
	use.gradient = TRUE, use.hessian = FALSE, check.model = TRUE, 
    args.fia = list(), numDeriv = TRUE
)
}

\arguments{
  \item{data}{
Either a \emph{numeric} \code{vector} for individual fit or a \emph{numeric} \code{matrix} or \code{data.frame} for multi-individual fit. The data on each position (column for multi-individual fit) must correspond to the respective line in the model file. Fitting for multiple individuals can be parallelized via \code{multicore}.
}
  \item{model.filename}{
A character \code{vector} specifying the location and name of the model file. 
}
  \item{restrictions.filename}{
\code{NULL} or a character \code{vector} or a \code{list} of characters. The default is \code{NULL} which corresponds to no restrictions. A character \code{vector} specifies the location or name of the restrictions file. A \code{list} of characters contains the restrictions directly. See Details and Examples.
}
  \item{n.optim}{
Number of optimization runs. Can be parallelized via \code{multicore}. Default is 5. If the number is high, fitting can take long time for \emph{large} models.
}
\item{fia}{
Number of random samples to be drawn in the Monte Carlo algorithm to estimate the Fisher Information Approximation (FIA) for MPTs only. See Details at \code{\link{fit.mpt}}
}
  \item{ci}{
A scalar corresponding to the size of the confidence intervals for the parameter estimates. Default is 95 which corresponds to 95\% confidence intervals.
}
  \item{starting.values}{
A \code{vector}, a \code{list}, or \code{NULL} (the default). If \code{NULL} starting values for parameters are randomly drawn from a uniform distribution with the interval \code{(0.1 - 0.9)}. See Details for the other options.
}
  \item{output}{
If "full" \code{fit.mpt} will additionally return the return values of \code{\link{nlminb}} and the Hessian matrices. (If "fia", \code{fit.mpt} will additionally return the results from \code{\link{get.mpt.fia}} (if \code{fia} not equal \code{NULL}).)
}
\item{reparam.ineq}{
  Logical. Indicates whether or not inequality restrictions (when present in the model file) should be enforced while fitting the model. If \code{TRUE} (the default) inequality restricted parameters will be reparameterized, if \code{FALSE} not. Probably irrelevant for none MPTs.
}
\item{fit.aggregated}{
Logical. Only relevant for multiple datasets (i.e., \code{matrix} or \code{data.frame}). Should the aggregated dataset (i.e., data summed over rows) be fitted? Default (\code{TRUE}) fits the aggregated data. 
}
\item{sort.param}{
Logical. If TRUE, parameters are alphabetically sorted in the parameter table. If FALSE, the first parameters in the parameter table are the non-restricted ones, followed by the restricted parameters. Default is TRUE.
}
  \item{show.messages}{
Logical. If TRUE the time the fitting algorithms takes is printed to the console.
}
  \item{model.type}{
Character vector specifying whether the model file is formatted in the easy way (\code{"easy"}; i.e., each line represents all branches corresponding to a response category) or the traditional EQN syntax (\code{"eqn"} or \code{"eqn2"}; see Details and e.g., Stahl & Klauer, 2007). If \code{model.filename} ends with .eqn or .EQN, \code{model.type} is automatically set to \code{"eqn"}. Default is \code{"easy"}.
}
  \item{multicore}{
Character vector. If not \code{"none"}, uses \code{snowfall} for parallelization (which needs to be installed separately via \code{install.packages(snowfall)}). If \code{"individual"}, parallelizes the optimization for each individual (i.e., data needs to be a \code{matrix} or \code{data.frame}). If \code{"n.optim"}, parallelizes the \code{n.optim} optimization runs. If not \code{"none"} (e.g., \code{"fia"}) calculation of FIA is parallelized (if FIA is requested). Default is \code{"none"} which corresponds to no parallelization. Note that you need to initialize \code{snowfall} in default settings. See \code{sfInit} and Details.
}
  \item{sfInit}{
Logical. Relevant if \code{multicore} is not \code{"none"}. If \code{TRUE}, \code{fit.mpt} will initialize and close the multicore support. If \code{FALSE}, (the default) assumes that \code{sfInit()} was initialized before. See Details.
}
  \item{nCPU}{
Scalar. Only relevant if \code{multicore} is not \code{"none"} and \code{sfInit} is TRUE. Number of CPUs used by \code{snowfall}. Default is 2.
}
  \item{lower.bound}{
numeric scalar or vector. Can be used in \code{fit.model} to set the lower bounds of the parameter space. See Details.
}
  \item{upper.bound}{
numeric scalar or vector. Can be used in \code{fit.model} to set the upper bounds of the parameter space. See Details.
}
  \item{control}{
list containing control arguments passed on to \code{\link{nlminb}}. See there.
}
  \item{use.gradient}{
logical. Whether or not the symbolically derived function returning the gradient should be used for fitting. Default is \code{TRUE} meaning gradient function is used.
}
  \item{use.hessian}{
logical. Whether or not the symbolically derived function returning the Hessian matrix should be used for fitting. Default is \code{FALSE} meaning hessian function is not used.
}
  \item{check.model}{
logical. Should model be checked with random values whether or not the expected values sum to one per tree? Default is \code{TRUE}. (This also controls whether other model checks during optimization are performed. If \code{FALSE} the most permissive fitting is performed.)
}
  \item{args.fia}{
named list of further arguments passed to \code{\link{get.mpt.fia}}, such as \code{mConst} to avoid numerical problems in the FIA function.
}
  \item{numDeriv}{
logical. Should the Hessian matrix of the maximum likelihood estimates be estimated numerically using \code{numDeriv::hessian} in case it cannot be estimated analytically? This can be extremely time and memory consuming for larger models. Default is TRUE.
  }
}
\details{
This functions should be used when fitting a model that is not an MPT model or when fitting using \code{fit.mpt} fails. For fitting MPT models and information on fitting MPT models see \code{\link{fit.mpt}}.

The model file for non-MPT models should be of the \code{easy} format. That is the ordinal number or rank of each line should correspond to this column/position in the data object. Model files can contain any visible function (i.e., including self-defined functions). However, note that the derivation that is needed for the gradient and Hessian function can only be done for those functions that \code{\link{D}} can handle. If derivation fails a warning will be given and fitting will be done without gradient and/or Hessian function.

Equations that correspond to one item type/category must be not be separated by an empty line. Equations that do not correspond to the same item type/category must be separated by at least one empty line.

Note that names of parameters in the model file should NOT start with \code{hank}. Variables with these names can lead to unforeseen problems as variables starting with these letters are internally used.

The restrictions file may contain (sequential) equality (i.e., =) and inequality (i.e., <) restrictions (see \code{\link{fit.mpt}} for more general info on the restrictions files). Note that inequality restrictions usually will lead to catastrophic results when used for non-MPT models. Our recommendation: Do never use inequality restrictions for non-MPT models. Equality restrictions or fixing parameters should be no problem though.

For equality restrictions, the equality restricted parameters are simply exchanged with their restrictions (i.e., another parameter or a number) before the fitting.

Restrictions or model files can contain comments (i.e., everything to the right of a # will be ignored; new behavior since version 0.9.2)

Both models and restrictions can be specified as \code{\link{textConnection}}s instead of as external files (see examples). Note that textConnections get "consumed" so you may need to specify them each time you fit a model using a connection (see \code{Examples} for how to avoid this).

Confidence intervals (CI) are based on the Hessian matrix produced by the symbolically derived function for the Hessian (i.e., the second derivative of the likelihood function). If it is based on a numerically estimated Hessian, a warning will be given. \cr

To set the starting values for the fitting process (e.g., to avoid local minima) one can set \code{starting.values} to a vector of length 2 and \code{n.optim > 1}. Then, starting values are randomly drawn from a uniform distribution from \code{starting.values[1]} to \code{starting.values[2]}.

Alternatively, one can supply a list with two elements to \code{starting.values}. Both elements need to be either of length 1 or of length equal to the number of parameters (if both are of length 1, it is the same as if you supply a vector of length 2). For each parameter n (in alphabetical order), a starting value is randomly drawn from a uniform distribution \code{starting.values[[1]][n]} to \code{starting.values[[2]][n]} (if length is 1, this is the border for all parameters).

The least interesting option is to specify the starting values individually by supplying a vector with the same length as the number of parameters. Starting values must be ordered according to the alphabetical order of the parameter names. Use \code{\link{check.mpt}} for a function that returns the alphabetical order of the parameters. If one specifies the starting values like that, \code{n.optim} will be set to 1 as all other values would not make any sense (the optimization routine will produce identical results with identical starting values).

The \code{lower.bound} and \code{upper.bound} needs to be of length 1 or equal to the number of \code{free} parameters. If length > 1, parameters are mapped to the bounds in alphabetic order of the parameters. Use \code{\link{check.mpt}} to obtain the alphabetical order of parameters for your model.

This function is basically a comfortable wrapper for \code{\link{fit.mptinr}} producing the appropriate objective, gradient, hessian, and prediction function from the model equations (passed via \code{model.filename}) whilst allowing for custom lower or upper bounds on the parameters. You can specify whether or not gradient or hessian function should be used for fitting with \code{use.gradient} or \code{use.hessian}, respectively.

Multicore fitting is achieved via the \code{snowfall} package and needs to be initialized via \code{sfInit}. As initialization needs some time, you can either initialize multicore facilities yourself using \code{sfInit()} and setting the \code{sfInit} argument to \code{FALSE} (the default) or let MPTinR initialize multicore facilities by setting the \code{sfInit} argument to \code{TRUE}. The former is recommended as initializing \code{snowfall} takes some time and only needs to be done once if you run \code{fit.mpt} multiple times. If there are any problems with multicore fitting, first try to initialize \code{snowfall} outside MPTinR (e.g., \code{sfInit( parallel=TRUE, cpus=2 )}). If this does not work, the problem is not related to MPTinR but to snowfall (for support and references visit: \url{http://www.imbi.uni-freiburg.de/parallel/}).\cr
Note that you should \emph{close} snowfall via \code{sfStop()} after using MPTinR.

}
\value{
For individual fits (i.e., \code{data} is a \code{vector}) a \code{list} containing one or more of the following components from the best fitting model:
\item{goodness.of.fit}{A \code{data.frame} containing the goodness of fit values for the model. \code{Log.Likelihood} is the Log-Likelihood value. \code{G.Squared}, \code{df}, and \code{p.value} are the \eqn{G^2}{G^2} goodness of fit statistic.}
\item{information.criteria}{A \code{data.frame} containing model information criteria based on the \eqn{G^2}{G^2} value. The FIA values(s) are presented if \code{fia} is not \code{NULL}.}
\item{model.info}{A \code{data.frame} containing other information about the model. If the rank of the Fisher matrix (\code{rank.fisher}) \emph{does not} correspond to the number of parameters in the model (\code{n.parameters}) this indicates a serious issue with the identifiability of the model. A common reason is that one of the parameter estimates lies on the bound of the parameter space (i.e., 0 or 1).}
\item{parameters}{A data.frame containing the parameter estimates and corresponding confidence intervals. If a restriction file was present, the restricted parameters are marked.}
\item{data}{A \code{list} of two matrices; the first one (\code{observed}) contains the entered data, the second one (\code{predicted}) contains the predicted values.}

For multi-dataset fits (i.e., \code{data} is a \code{matrix} or \code{data.frame}) a \code{list} with similar elements, but the following differences:\cr
The first elements, \code{goodness.of.fit}, \code{information.criteria}, and \code{model.info}, contain the same information as for individual fits, but each are \code{lists} with three elements containing the respective values for: each individual in the list element \code{individual}, the sum of the individual values in the list element \code{sum}, and the values corresponding to the fit for the aggregated data in the list element \code{aggregated}.\cr
\code{parameters} is a list containing:
\item{individual}{A 3-dimensional array containing the parameter estimates ([,1,]), confidence intervals [,2:3,], and, if restrictions not \code{NULL}, column 4 [,4,] is 0 for non-restricted parameters, 1 for equality restricted parameters, and 2 for inequality restricted parameters. The first dimension refers to the parameters, the second to the information on each parameter, and the third to the individual/dataset.}
\item{mean}{A \code{data.frame} with the mean parameter estimates from the individual estimates. No confidence intervals can be provided for these values.}
\item{aggregated}{A data.frame containing the parameter estimates and corresponding confidence intervals for the aggregated data. If a restriction file was present, the restricted parameters are marked.}
The element \code{data} contains two matrices, one with the \code{observed}, and one with the \code{predicted} data (or is a list containing lists with \code{individual} and \code{aggregated} \code{observed} and \code{predicted} data).

If \code{n.optim} > 1, the \code{\link{summary}} of the vector (matrix for multi-individual fit) containing the Log-Likelihood values returned by each run of \code{optim} is added to the output: \code{fitting.runs}

When \code{output == "full"} the list contains the additional items:\cr
\item{optim.runs}{A list (or list of lists for multiple datasets) containing the outputs from all runs by \code{nlminb} (including those runs produced when fitting did not converge)}
\item{best.fits}{A list (or list of lists for multiple datasets) containing the outputs from the runs by \code{nlminb} that had the lowest likelihood (i.e., the successful runs)}
\item{hessian}{A list  containing the Hessian matrix or matrices of the final parameter estimates.}
}
\note{
Warnings may relate to the optimization routine (e.g., \code{Optimization routine [...] did not converge successfully}).\cr
In these cases it is recommended to rerun the model fitting to check if the results are stable.

The likelihood returned does not include the factorial constants of the multinomial probability-mass functions.

All (model or restriction) files should end with an empty line, otherwise a warning will be shown.
}
\references{
Broeder, A., & Schuetz, J. (2009). Recognition ROCs are curvilinear-or are they? On premature arguments against the two-high-threshold model of recognition. \emph{Journal of Experimental Psychology: Learning, Memory, and Cognition}, 35(3), 587. doi:10.1037/a0015279

Wickens, T. D. (2002). \emph{Elementary Signal Detection Theory}. Oxford; New York: Oxford University Press.
}
\author{
Henrik Singmann and David Kellen.
}

\examples{

\dontrun{

#####################################
## Fit response-bias or payoff ROC ##
#####################################
  
# Example from Broder & Schutz (2009)
# We fit the data from the 40 individuals from their Experiment 3
# We fit three different models:
# 1. Their SDT Model: br.sdt
# 2. Their 2HTM model: br.2htm
# 3. A restricted 2HTM model with Dn = Do: br.2htm.res
# 4. A 1HTM model (i.e., Dn = 0): br.1htm

data(d.broeder, package = "MPTinR")
m.2htm <- system.file("extdata", "5points.2htm.model", package = "MPTinR")


# We specify the SDT model in the code using a textConnection.
# However, textConnection is only called in the function call on the string.

m.sdt <- "
1-pnorm((cr1-mu)/ss)
pnorm((cr1-mu)/ss)

1-pnorm(cr1)
pnorm(cr1)

1-pnorm((cr2-mu)/ss)
pnorm((cr2-mu)/ss)

1-pnorm(cr2)
pnorm(cr2)

1-pnorm((cr3-mu)/ss)
pnorm((cr3-mu)/ss)

1-pnorm(cr3)
pnorm(cr3)

1-pnorm((cr4-mu)/ss)
pnorm((cr4-mu)/ss)

1-pnorm(cr4)
pnorm(cr4)

1-pnorm((cr5-mu)/ss)
pnorm((cr5-mu)/ss)

1-pnorm(cr5)
pnorm(cr5)
"

# How does the model look like?
check.mpt(textConnection(m.sdt))

# fit the SDT (unequal variance version)
br.uvsdt <- fit.model(d.broeder, textConnection(m.sdt), 
  		lower.bound = c(rep(-Inf, 5), 0, 1), upper.bound = Inf)

# Is there any effect of studying the items?
br.uvsdt.2 <- fit.model(d.broeder, textConnection(m.sdt), 
			restrictions.filename = list("mu = 0", "ss = 1"), 
			lower.bound = -Inf, upper.bound = Inf)

(diff.g2 <- br.uvsdt.2[["goodness.of.fit"]][["sum"]][["G.Squared"]] -
			br.uvsdt[["goodness.of.fit"]][["sum"]][["G.Squared"]])
(diff.df <- br.uvsdt.2[["goodness.of.fit"]][["sum"]][["df"]] - 
			br.uvsdt[["goodness.of.fit"]][["sum"]][["df"]])
1 - pchisq(diff.g2, diff.df)

# fit the equal variance SDT model:
br.evsdt <- fit.model(d.broeder, textConnection(m.sdt), 
			lower.bound = c(rep(-Inf, 5), 0), upper.bound = Inf, 
			restrictions.filename = list("ss = 1"))

# fit the MPTs (see also ?fit.mpt).
# In contrast to ?fit.mpt we specify the restrictions using a textConnection or a list!
br.2htm <- fit.mpt(d.broeder, m.2htm)
br.2htm.res <- fit.mpt(d.broeder, m.2htm, textConnection("Do = Dn"))
br.1htm <- fit.mpt(d.broeder, m.2htm, list("Dn = 0"))

select.mpt(list(uvsdt = br.uvsdt, evsdt = br.evsdt, two.htm = br.2htm, 
			two.htm.res = br.2htm.res, one.htm = br.1htm), output = "full")

# the restricted 2HTM "wins" for individual data (although evsdt does not perform too bad),
# but the 2htm and restricted 2htm restricted "win" for aggregated data.


###################################
## Fit confidence rating ROC SDT ##
###################################
#(see ?roc6 for more examples)

# We fit example data from Wickens (2002, Chapter 5)
# The example data is from Table 5.1, p. 84
# (data is entered in somewhat different order).

# Note that criteria are defined as increments to 
# the first (i.e., leftmost) criterion!
# This is the only way to do it in MPTinR.

# Data
dat <- c(47, 65, 66, 92, 136, 294, 166, 161, 138, 128, 63, 43)

# UVSDT
m.uvsdt <- "
pnorm(cr1, mu, sigma)
pnorm(cr1+cr2, mu, sigma) - pnorm(cr1, mu, sigma)
pnorm(cr3+cr2+cr1, mu, sigma) - pnorm(cr2+cr1, mu, sigma)
pnorm(cr4+cr3+cr2+cr1, mu, sigma) - pnorm(cr3+cr2+cr1, mu, sigma)
pnorm(cr5+cr4+cr3+cr2+cr1, mu, sigma) - pnorm(cr4+cr3+cr2+cr1, mu, sigma)
1 - pnorm(cr5+cr4+cr3+cr2+cr1, mu, sigma)
  
pnorm(cr1)
pnorm(cr2+cr1) - pnorm(cr1)
pnorm(cr3+cr2+cr1) - pnorm(cr2+cr1)
pnorm(cr4+cr3+cr2+cr1) - pnorm(cr3+cr2+cr1)
pnorm(cr5+cr4+cr3+cr2+cr1) - pnorm(cr4+cr3+cr2+cr1)
1 - pnorm(cr5+cr4+cr3+cr2+cr1)
"
check.mpt(textConnection(m.uvsdt))

# Model fitting
(cr_sdt <- fit.model(dat, textConnection(m.uvsdt),
            lower.bound=c(-Inf, rep(0, 5), 0.1), upper.bound=Inf))

# To obtain the criteria (which match those in Wickens (2002, p. 90)
# obtain the cumulative sum:

cumsum(cr_sdt$parameters[paste0("cr",1:5), 1, drop = FALSE])

}

}


\seealso{
\code{\link{check.mpt}} for a function that can help in constructing models.

\code{\link{fit.mptinr}} for a function that can fit arbitrary objective functions.

\code{\link{fit.mpt}} for the function to fit MPTs (it should be slightly faster for MPTs).

\link{roc6} for more examples fitting different SDT models.
}
  
\keyword{models}
\keyword{tree}

