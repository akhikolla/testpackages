\name{fit.mpt.old}
\alias{fit.mpt.old}
\title{
Function to fit MPT models (old)
}
\description{
\code{fit.mpt.old} function fits \emph{binary} multinomial processing tree models (MPT models; e.g., Riefer & Batchelder, 1988). However, this function is an old version  using the L-BFGS-B optimization routine. See \code{\link{fit.mpt}} for the new version. 
}
\usage{
fit.mpt.old(
	data,
	model.filename, 
	restrictions.filename = NULL, 
	n.optim = 5,
	fia = NULL,
	ci = 95, 
	starting.values = NULL,
	output = c("standard", "fia", "full"),
	reparam.ineq = TRUE,
	sort.param = TRUE,
	model.type = c("easy", "eqn", "eqn2"),
	multicore = c("none", "individual", "n.optim"), sfInit = FALSE, nCPU = 2
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
\code{NULL} or a character \code{vector} or a \code{list} of characters. The default is \code{NULL} which corresponds to no restrictions. A character \code{vector} specifies the location or name of the restrictions file. A \code{list} of characters contains the restrictions directly.
}
  \item{n.optim}{
Number of optimization runs. Can be parallelized via \code{multicore}. Default is 5. If the number is high, fitting can take long for \emph{large} models.
}
\item{fia}{
Number of random samples to be drawn in the Monte Carlo algorithm to estimate the Fisher Information Approximation (FIA), a minimum description length based measure of model complexity (see Wu, Myung & Batchelder, 2010). The default is \code{NULL} which corresponds to no computation of the FIA. Reasonable values (e.g., > 200000) can lead to long computation times (minutes to hours) depending on the size of the model. See Details.
}
  \item{ci}{
A scalar corresponding to the size of the confidence intervals for the parameter estimates. Default is 95 which corresponds to 95\% confidence intervals.
}
  \item{starting.values}{
A \code{vector} or \code{NULL}. If \code{NULL} (the default), starting values for parameters are randomly drawn from a uniform distribution with the interval \code{(0.05 - 0.95)}. If \code{length(starting.values)==2}, starting values  are randomly drawn from a uniform distribution with the interval \code{starting.values[1] - starting.values[2])}. If \code{length(starting.values)} matches the number of parameters in the model, \code{starting.values} will be used as the starting values for fitting and \code{n.optim} will be set to 1. See Details.
}
  \item{output}{
If "fia", \code{fit.mpt} will additionally return the results from \code{\link{get.mpt.fia}} (if calculated \code{fia} not equal \code{NULL}). If "full" \code{fit.mpt} will additionally return the results from \code{\link{get.mpt.fia}} and the output of \code{\link{optim}}.
}
  \item{reparam.ineq}{
  Logical. Indicates whether or not inequality restrictions (when present in the model file) should be enforced while fitting the model. If \code{TRUE} (the default) inequality restricted parameters will be reparameterized, if \code{FALSE} not. See Details.
}
  \item{sort.param}{
Logical. If TRUE, parameters are alphabetically sorted in the parameter table. If FALSE, the first parameters in the parameter table are the non-restricted ones, followed by the restricted parameters. Default is TRUE.
}
  \item{model.type}{
Character vector specifying whether the model file is formatted in the easy way (\code{"easy"}; i.e., each line represents all branches corresponding to a response category) or the traditional EQN syntax (\code{"eqn"} or \code{"eqn2"}; see Details and e.g., Stahl & Klauer, 2007). If \code{model.filename} ends with .eqn or .EQN, \code{model.type} is automatically set to \code{"eqn"}. Default is \code{"easy"}.
}
  \item{multicore}{
Character vector. If not \code{"none"}, uses \code{snowfall} for parallelization (which needs to be installed separately via \code{install.packages(snowfall)}). If \code{"individual"}, parallelizes the optimization for each individual (i.e., data needs to be a \code{matrix} or \code{data.frame}). If \code{"n.optim"}, parallelizes the \code{n.optim} optimization runs. Default is \code{"none"} which corresponds to no parallelization. Note that you need to initialize \code{snowfall} in default settings. See \code{sfInit} and Details.
}
  \item{sfInit}{
Logical. Relevant if \code{multicore} is not \code{"none"}. If \code{TRUE}, \code{fit.mpt} will initialize and close the multicore support. If \code{FALSE}, (the default) assumes that \code{sfInit()} was initialized before. See Details.
}
  \item{nCPU}{
Scalar. Only relevant if \code{multicore} is not \code{"none"} and \code{sfInit} is TRUE. Number of CPUs used by \code{snowfall}. Default is 2.
}
}
\details{
There is a new version of this function using \code{\link{nlminb}} and the analytically derived gradient and hessian. See \code{\link{fit.mpt}}. We recommend using the new version \code{\link{fit.mpt}}, only use this version if you are sure on what to do.

The model file is either of the easy format or the "classical" EQN format (see below).\cr
In the easy format (the default) the model file contains all trees of the model. Trees are separated by at least one empty line. Everything to the right of a hash (#) is ignored (this behavior is new since version 0.9.2). Lines starting with a # are treated as empty. Each line in each tree corresponds to all branches of this tree (concatenated by a +) that correspond to one of the possible response categories. The position of each line must correspond to the position of this response category in the data object (for multi-individual fit to the respective column).

The difference between both types of EQN format (\code{"eqn"} or\code{"eqn2"}) is the way the first line of the model file is treated. If \code{model.file} is set to \code{"eqn"}, \code{MPTinR} will ignore the first line of the model file and will read the rest of the file (as does multiTree; Moshagen, 2010). If \code{model.file} is set to \code{"eqn2"} MPTinR will only read as many lines as indicated in the first line of the EQN model file (as does e.g., HMMTree; Stahl & Klauer, 2007). As default \code{fit.mpt} expects the easy format, but if the filename ends with .eqn or .EQN and \code{model.type} is \code{"easy"}, \code{model.type} is set to \code{"eqn"}\cr
For the EQN format consult one of the corresponding papers (see e.g., Moshagen, 2010; Stahl & Klauer, 2007). The positions in the data object (number of column for multi-individual fit) must correspond to the category number in the EQN file.\cr

Note that names of parameters in the model file should not start with \code{hank.}. Variables with these names can lead to unforeseen problems as variables starting with these letters are internally used.

The restrictions file may contain (sequential) equality (i.e., =) and inequality (i.e., <) restrictions and must adhere to the following rules:\cr
	1. Inequalities first.\cr
	2. If a variable appears in an inequality restriction, it can not be on the LHS of any further restriction.\cr
	3. If a variable appears on RHS of an equality restriction, it can not appear on LHS of an equality restriction.\cr
Note that only "<" is supported as inequality operator but not ">"!\cr
Examples of restrictions are (the following could all appear in one restrictions file):\cr
\code{D1 < D2 < D3}\cr
\code{D4 = D3}\cr
\code{B1 = B3 = 0.3333}\cr
\code{X4 = X5 = D3}\cr
Restrictions file may contain comments (i.e., everything to the right of a # will be ignored; new behavior since version 0.9.2)

For equality restrictions, the equality restricted parameters are simply exchanged with their restrictions before the fitting.\cr
For inequality restricted parameters, the model is reparameterized so that only the rightmost parameter of an inequality restriction remains the original parameter. Each instance of the other parameters in this restriction is replaced by the product of the rightmost parameter and dummy parameters (see Knapp & Batchelder, 2004). This procedure (which is equivalent to method A described in Knapp & Batchelder, 2004) leads to an equivalent model (although the binary MPT structure is not apparent in the resulting equations).\cr
To prohibit this reparameterization (i.e., if the inequality restrictions hold without reparameterization), you can set \code{reparam.ineq} to \code{FALSE}. This can be useful for obtaining the FIA (see examples in Wu, Myung, & Batchelder, 2010).

The fitting/optimization is achieved via \code{\link{optim}}'s L-BFGS-B method by Byrd et al. (1995) with random starting values. To avoid local minima it is useful to set \code{n.optim} > 1. If \code{n.optim} > 1, the \code{\link{summary}} of the vector containing the Log-Likelihood values returned by each run of \code{optim} is added to the output (to check whether local minima were present). If the model is rather big, \code{n.optim} > 1 can be slow.

To obtain a measure of the model's complexity beyond the number of parameters (and taking inequality restrictions into account), set \code{fia} to a (reasonably high) scalar integer (i.e., a number). Then, \code{fit.mpt} will obtain the Fisher information approximation (FIA), a minimum description based measure of model complexity, using the algorithm provided by Wu, Myung, & Batchelder (2010a, 2010b) ported from Matlab to R. When performing model-selection, this measure is superior to other methods such as the Akaike information criterion (AIC) or Bayesian information criterion (BIC) which basically only take the number of parameters into account.\cr
To get the FIA, \code{fit.mpt.old} performs the following steps:\cr
1. The representation of the model as equations is transformed into the string representation of the model in the context-free language of MPT models (L-BMPT; Purdy & Batchelder, 2009). For this step to be successful it is \emph{absolutely necessary} that the equations representing the model perfectly map the tree structure of the MPT. That is, the model file is only allowed to contain parameters, their negations (e.g., \code{Dn} and \code{(1 - Dn)}) and the operators + and *, but nothing else. Simplifications of the equations will seriously distort this step. This step is achieved by \code{\link{make.mpt.cf}}.\cr
2. The context free representation of the model is then fed into the MCMC function computing the FIA (the port of BMPTFIA provided by Wu, Myung & Batchelder (2010a), see \code{\link{bmpt.fia}}).\cr
(Actually, both steps are achieved by a call to \code{\link{get.mpt.fia}})

Once again: If one wants to compute the FIA, it is \emph{absolutely necessary}, that the representation of the model via equations in the model file exactly maps on the structure of the binary MPT (see \code{\link{make.mpt.cf}} for more details).

Confidence intervals (CI) are based on the observed Hessian matrix returned by the minimization algorithm (\code{\link{optim}}).\cr
For inequality restricted parameters, the CIs are computed using the parameter estimates' variance bounds (see Baldi & Batchelder, 2003; especially equation 19). Note that these bounds represent the "worst case scenario" variances, and can lead to CIs outside parameter boundaries if the set of inequalities is large and/or the variances for the reparameterized model are large (Note that CIs for non-restricted parameters can be outside the parameter boundaries as well due to large variances).

To set the starting values for the fitting process (e.g., to avoid local minima) one can set \code{starting.values} to a vector of length 2. Then, starting values are randomly drawn from a uniform distribution from \code{starting.values[1]} to \code{starting.values[2]}.\cr
Furthermore, one can specify the starting values individually by supplying a vector with the same length as the number of parameters. Starting values must be ordered according to the alphabetical order of the parameters. Use \code{\link{check.mpt}} for a function that returns the alphabetical order of the parameters. If one specifies the starting values like that, \code{n.optim} will be set to 1 as all other values would not make any sense (the optimization routine will produce identical results with identical starting values).

Multicore fitting is achieved via the \code{snowfall} package and needs to be initialized via \code{sfInit}. As initialization needs some time, you can either initialize multicore facilities yourself using \code{sfInit()} and setting the \code{sfInit} argument to \code{FALSE} (the default) or let MPTinR initialize multicore facilities by setting the \code{sfInit} argument to \code{TRUE}. The former is recommended as initializing \code{snowfall} takes some time and only needs to be done once if you run \code{fit.mpt.old} multiple times. If there are any problems with multicore fitting, first try to initialize \code{snowfall} outside MPTinR (e.g., \code{sfInit( parallel=TRUE, cpus=2 )}). If this does not work, the problem is not related to MPTinR but to snowfall (for support and references visit: \url{http://www.imbi.uni-freiburg.de/parallel/}).\cr
Note that you need to \emph{close} snowfall via \code{sfStop()} after using MPTinR.

\code{fit.model()} is essentially a copy of \code{fit.mpt.old} that allows the user to specify the upper and lower bounds of the parameters. This function can be used to fit other models than MPT models that can be described in a model file. That is, the model file can contain any type of valid R expressions including R functions (potentially self-written) visible in the global environment (i.e., not only +, *, and - as operators).  Currently \code{fit.model} should be viewed as experimental.

Note that \code{fit.model()} is usually slower than \code{fit.mpt.old} as there are some more checks in the critical function calculating the likelihood of the model. \cr


The \code{lower.bound} and \code{upper.bound} needs to be of length 1 or equal to the number of \code{free} parameters. If length > 1, parameters are mapped to the bounds in alphabetic order of the parameters. Use \code{\link{check.mpt}} to obtain the alphabetical order of parameters for your model.

While it should be possible to specify equality or fixed restrictions it will probably lead to unforeseen consequences to specify inequality restrictions for non-MPT models.

}
\value{
For individual fits (i.e., \code{data} is a \code{vector}) a \code{list} containing one or more of the following components from the best fitting model:
\item{goodness.of.fit}{A \code{data.frame} containing the goodness of fit values for the model. \code{Log.Likelihood} is the Log-Likelihood value. \code{G.Squared}, \code{df}, and \code{p.value} are the \eqn{G^2}{G^2} goodness of fit statistic.}
\item{information.criteria}{A \code{data.frame} containing model information criteria based on the G^2 value. The FIA values(s) are presented if \code{fia} is not \code{NULL}.}
\item{model.info}{A \code{data.frame} containing other information about the model. If the rank of the Hessian matrix (\code{rank.hessian}) \emph{does not} correspond to the number of parameters in the model (\code{n.parameters}) this indicates a serious issue with the identifiability of the model.}
\item{parameters}{A data.frame containing the parameter estimates and corresponding confidence intervals. If a restriction file was present, the restricted parameters are marked.}
\item{data}{A \code{list} of two matrices; the first one (\code{observed}) contains the entered data, the second one (\code{predicted}) contains the predicted values.}

For multi-individual fits (i.e., \code{data} is a \code{matrix} or \code{data.frame}) a \code{list} with similar elements, but the following differences.\cr
The first elements, \code{goodness.of.fit}, \code{information.criteria}, and \code{model.info}, contain the same information as for individual fits, but each are \code{lists} with three elements containing the respective values for: each individual in the list element \code{individual}, the sum of the individual values in the list element \code{sum}, and the values corresponding to the fit for the aggregated data in the list element \code{aggregated}.\cr
\code{parameters} is a list containing:
\item{individual}{A 3-dimensional array containing the parameter estimates ([,1,]), confidence intervals [,2:3,], and, if restrictions not \code{NULL}, column 4 [,4,] is 0 for non-restricted parameters, 1 for equality restricted parameters, and 2 for inequality restricted parameters. The first dimension refers to the parameters, the second to the information on each parameter, and the third to the individuals.}
\item{mean}{A \code{data.frame} with the mean parameter estimates from the individual estimates. No confidence intervals can be provided for these values.}
\item{aggregated}{A data.frame containing the parameter estimates and corresponding confidence intervals for the aggregated data. If a restriction file was present, the restricted parameters are marked.}
The element \code{data} contains two matrices, one with the \code{observed}, and one with the \code{predicted} data.

If \code{n.optim} > 1, the \code{\link{summary}} of the vector (matrix for multi-individual fit) containing the Log-Likelihood values returned by each run of \code{optim} is added to the output.

When using R (>= 2.13.0) compiling \code{fit.mpt.old} using \code{compilers} \code{cmpfun} can significantly improve fitting time.
}
\note{
There may be several warnings \code{fit.mpt.old} throws while fitting MPT models. Most of them are not problematic and related to matrix operations needed for confidence intervals. Examples:\cr
\code{In sqrt(var.params) : NaNs produced}\cr
\code{In sqrt(min(var.bound.tmp)) : NaNs produced}\cr
These warnings are not critical.

Other warnings may relate to the optimization routine (e.g., \code{Optimization routine [...] did not converge successfully}).\cr
In these cases it is recommended to rerun \code{fit.mpt.old} to check if the results are stable.
}
\references{
Baldi, P. & Batchelder, W. H. (2003). Bounds on variances of estimators for multinomial processing tree models. \emph{Journal of Mathematical Psychology}, 47, 467-470. 

Broeder, A., & Schuetz, J. (2009). Recognition ROCs are curvilinear-or are they? On premature arguments against the two-high-threshold model of recognition. \emph{Journal of Experimental Psychology: Learning, Memory, and Cognition}, 35(3), 587. doi:10.1037/a0015279


Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. \emph{SIAM J. Scientific Computing}, 16, 1190-1208.

Knapp, B. R., & Batchelder, W. H. (2004). Representing parametric order constraints in multi-trial applications of multinomial processing tree models. \emph{Journal of Mathematical Psychology}, 48, 215-229.

Moshagen, M. (2010). multiTree: A computer program for the analysis of multinomial processing tree models. \emph{Behavior Research Methods}, 42, 42-54.

Purdy, B. P., & Batchelder, W. H. (2009). A context-free language for binary multinomial processing tree models. \emph{Journal of Mathematical Psychology}, 53, 547-561.

Riefer, D. M., & Batchelder, W. H. (1988). Multinomial modeling and the measurement of cognitive processes. \emph{Psychological Review}, 95, 318-339.

Stahl, C. & Klauer, K. C. (2007). HMMTree: A computer program for latent-class hierarchical multinomial processing tree models. \emph{Behavior Research Methods}, 39, 267- 273.

Wu, H., Myung, J.I., & Batchelder, W.H. (2010a). Minimum description length model selection of multinomial processing tree models. \emph{Psychonomic Bulletin & Review}, 17, 275-286.

Wu, H., Myung, J.I., & Batchelder, W.H. (2010b). On the minimum description length complexity of multinomial processing trees. \emph{Journal of Mathematical Psychology}, 54, 291-303.
}
\author{
Henrik Singmann and David Kellen with help from Karl Christoph Klauer and Fabian Hoelzenbein.
}
\note{
All (model or restriction) files should end with an empty line, otherwise a warning will be shown.
}

\examples{

\dontrun{

# The first example fits the MPT model presented in Riefer and Batchelder (1988, Figure 1)
# to the data presented in Riefer and Batchelder (1988, Table 1)
# Note that Riefer and Batchelder (1988, pp. 328) did some hypotheses tests, that are not done here.
# Rather, we use each condition (i.e., row in Table 1) as a different individual.
# We try to use n.optim = 1 here, but this can lead to local minima
# In general we recommend to set n.optim >= 5

# load the data
data(rb.fig1.data)

#get the character string with the position of the model:
model1 <- system.file("extdata", "rb.fig1.model", package = "MPTinR")
model1.eqn <- system.file("extdata", "rb.fig1.model.eqn", package = "MPTinR")

# just fit the first "individual":
fit.mpt.old(rb.fig1.data[1,], model1, n.optim = 1)

#fit all "individuals":
fit.mpt.old(rb.fig1.data, model1, n.optim = 1)

#fit all "individuals" using the .EQN model file:
fit.mpt.old(rb.fig1.data, model1.eqn, n.optim = 1)



# The second example fits the MPT model presented in Riefer and Batchelder (1988, Figure 2)
# to the data presented in Riefer and Batchelder (1988, Table 3)
# First, the model without restrictions is fitted: ref.model
# Next, the model with all r set equal is fitted: r.equal
# Then, the model with all c set equal is fitted: c.equal
# Finally, the inferential tests reported by Riefer & Batchelder, (1988, p. 332) are executed.
# Note, that n.optim = 10, because of frequent local minima.

# get the data
data(rb.fig2.data)

# positions of model and restriction files:
model2 <- system.file("extdata", "rb.fig2.model", package = "MPTinR")
model2r.r.eq <- system.file("extdata", "rb.fig2.r.equal", package = "MPTinR")
model2r.c.eq <- system.file("extdata", "rb.fig2.c.equal", package = "MPTinR")

# The full (i.e., unconstrained) model
(ref.model <- fit.mpt.old(rb.fig2.data, model2, n.optim = 10))

# All r equal
(r.equal <- fit.mpt.old(rb.fig2.data, model2, model2r.r.eq, n.optim = 10))

# All c equal
(c.equal <- fit.mpt.old(rb.fig2.data, model2, model2r.c.eq, n.optim = 10))

# is setting all r equal a good idea?
(g.sq.r.equal <- r.equal[["goodness.of.fit"]][["G.Squared"]] - 
			ref.model[["goodness.of.fit"]][["G.Squared"]])
(df.r.equal <- r.equal[["goodness.of.fit"]][["df"]] - 
			ref.model[["goodness.of.fit"]][["df"]])
(p.value.r.equal <- pchisq(g.sq.r.equal, df.r.equal , lower.tail = FALSE))

# is setting all c equal a good idea?
(g.sq.c.equal <- c.equal[["goodness.of.fit"]][["G.Squared"]] - 
		ref.model[["goodness.of.fit"]][["G.Squared"]])
(df.c.equal <- c.equal[["goodness.of.fit"]][["df"]] - 
		ref.model[["goodness.of.fit"]][["df"]])
(p.value.c.equal <- pchisq(g.sq.c.equal, df.c.equal , lower.tail = FALSE))


# Example from Broeder & Schuetz (2009)
# We fit the data from the 40 individuals from their Experiment 3
# We fit three different models:
# 1. Their 2HTM model: br.2htm
# 2. A restricted 2HTM model with Dn = Do: br.2htm.res
# 3. A 1HTM model (i.e., Dn = 0): br.1htm
# We fit the models with, as well as without, applied inequality restrictions (see Details)
# that is, for some models (.ineq) we impose: G1 < G2 < G3 < G4 < G5 
# As will be apparent, the inequality restrictions do not hold for all individuals.
# Finally, we compute the FIA for all models, taking inequalities into account.
# Note: The following examples will take some time (> 1 hour).

data(d.broeder)
m.2htm <- system.file("extdata", "5points.2htm.model", package = "MPTinR")
r.2htm <- system.file("extdata", "broeder.2htm.restr", package = "MPTinR")
r.1htm <- system.file("extdata", "broeder.1htm.restr", package = "MPTinR")
i.2htm <- system.file("extdata", "broeder.2htm.ineq", package = "MPTinR")
ir.2htm <- system.file("extdata", "broeder.2htm.restr.ineq", package = "MPTinR")
ir.1htm <- system.file("extdata", "broeder.1htm.restr.ineq", package = "MPTinR")

# fit the original 2HTM
br.2htm <- fit.mpt.old(d.broeder, m.2htm)
br.2htm.ineq <- fit.mpt.old(d.broeder, m.2htm, i.2htm)

# do the inequalities hold for all participants?
br.2htm.ineq[["parameters"]][["individual"]][,"estimates",]
br.2htm[["parameters"]][["individual"]][,"estimates",]
# See the difference between forced and non-forced inequality restrictions:
round(br.2htm[["parameters"]][["individual"]][,"estimates",] -
		br.2htm.ineq[["parameters"]][["individual"]][,"estimates",],2)

# The same for the other two models
# The restricted 2HTM
br.2htm.res <- fit.mpt(d.broeder, m.2htm, r.2htm)
br.2htm.res.ineq <- fit.mpt(d.broeder, m.2htm, ir.2htm)
round(br.2htm.res[["parameters"]][["individual"]][,"estimates",] -
		br.2htm.res.ineq[["parameters"]][["individual"]][,"estimates",],2)
# The 1HTM
br.1htm <- fit.mpt(d.broeder, m.2htm, r.1htm)
br.1htm.ineq <- fit.mpt(d.broeder, m.2htm, ir.1htm)
round(br.2htm.res[["parameters"]][["individual"]][,"estimates",] -
		br.2htm.res.ineq[["parameters"]][["individual"]][,"estimates",],2)

# These results show that we cannot compute inequality constraints for the non inequality
#  imposed models (It would look differently if we excluded critical cases, 
# i.e., 2, 6, 7, 10, 18, 21, 25, 29, 32, 34, 35, 37, 38)
# Therefore, we get the FIA for the models as computed above 
# WARNING: The following part will take a long time!

br.2htm.fia <- fit.mpt.old(d.broeder, m.2htm, fia = 200000)
br.2htm.ineq.fia <- fit.mpt.old(d.broeder, m.2htm, i.2htm, fia = 200000)
br.2htm.res.fia <- fit.mpt.old(d.broeder, m.2htm, r.2htm, fia = 200000 )
br.2htm.res.ineq.fia <- fit.mpt.old(d.broeder, m.2htm, ir.2htm, fia = 200000)
br.1htm.fia <- fit.mpt.old(d.broeder, m.2htm, r.1htm, fia = 200000)
br.1htm.ineq.fia <- fit.mpt.old(d.broeder, m.2htm, ir.1htm, fia = 200000)

# Model selection using the FIA
(br.select <- select.mpt(list(orig.2htm = br.2htm.fia, orig.2htm.ineq = br.2htm.ineq.fia,
		res.2htm = br.2htm.res.fia, res.2htm.ineq = br.2htm.res.ineq.fia, 
		orig.1htm = br.1htm.fia, orig.1htm.ineq = br.1htm.ineq.fia)))
# The same results, ordered by FIA
br.select[order(br.select[,"delta.FIA.sum"]),]

# Compare this with the model selection not using FIA:
select.mpt(list(orig.2htm = br.2htm, orig.2htm.ineq = br.2htm.ineq, 
		res.2htm = br.2htm.res, res.2htm.ineq = br.2htm.res.ineq, 
		orig.1htm = br.1htm, orig.1htm.ineq = br.1htm.ineq))




# compare speed of no multicore versus multicore for multiple optimization runs:

require(snowfall)
# change number of CPUs if more are available
nCPU = 2
sfInit( parallel=TRUE, cpus=nCPU, type = "SOCK" )

# NO multicore
system.time(fit.mpt.old(d.broeder, m.2htm))

# multicore:
system.time(fit.mpt.old(d.broeder, m.2htm, multicore = "n.optim"))

sfStop()
}

  }
\seealso{
\code{\link{fit.mpt}} for the current and recommended function fot fitting MPTs
}
  
\keyword{models}
\keyword{tree}

