\name{fit.mpt}
\alias{fit.mpt}
\title{
Function to fit MPT models
}
\description{
\code{fit.mpt} fits \emph{binary} multinomial processing tree models (MPT models; e.g., Riefer & Batchelder, 1988) from an external model file and (optional) external restrictions using the general-purpose quasi-Newton box-constraint optimization routine provided by Byrd et al. (1995). Additionally, measures for model selection (AIC, BIC, FIA) can be computed. 

}
\usage{
fit.mpt(
	data,
	model.filename, 
	restrictions.filename = NULL, 
	n.optim = 5,
	fia = NULL,
	ci = 95, 
	starting.values = NULL,
	output = c("standard", "fia", "full"),
	reparam.ineq = TRUE,
	fit.aggregated = TRUE,
	sort.param = TRUE,
	show.messages = TRUE,
	model.type = c("easy", "eqn", "eqn2"),
	multicore = c("none", "individual", "n.optim", "fia"), sfInit = FALSE, nCPU = 2,
	control = list(), args.fia = list(), numDeriv = TRUE
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
Number of optimization runs. Can be parallelized via \code{multicore}. Default is 5. If the number is high, fitting can take long for \emph{large} models.
}
\item{fia}{
Number of random samples to be drawn in the Monte Carlo algorithm to estimate the Fisher Information Approximation (FIA), a minimum description length based measure of model complexity (see Wu, Myung & Batchelder, 2010). The default is \code{NULL} which corresponds to no computation of the FIA. Reasonable values (e.g., > 200000) can lead to long computation times (minutes to hours) depending on the size of the model. See Details.
}
  \item{ci}{
A scalar corresponding to the size of the confidence intervals for the parameter estimates. Default is 95 which corresponds to 95\% confidence intervals.
}
  \item{starting.values}{
A \code{vector}, a \code{list}, or \code{NULL} (the default). If \code{NULL} starting values for parameters are randomly drawn from a uniform distribution with the interval \code{(0.1 - 0.9)}. See Details of \code{\link{fit.mptinr}} for the other options.
}
  \item{output}{
If "fia", \code{fit.mpt} will additionally return the results from \code{\link{get.mpt.fia}} (if \code{fia} not equal \code{NULL}). If "full" \code{fit.mpt} will additionally return the results from \code{\link{get.mpt.fia}} and the output of \code{\link{nlminb}} and the Hessian matrix/matrices.
}
\item{reparam.ineq}{
  Logical. Indicates whether or not inequality restrictions (when present in the model file) should be enforced while fitting the model. If \code{TRUE} (the default) inequality restricted parameters will be reparameterized, if \code{FALSE} not. See Details.
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
  \item{control}{
list containing control arguments passed on to \code{\link{nlminb}}. See there.
}
  \item{args.fia}{
named list of further arguments passed to \code{\link{get.mpt.fia}}, such as \code{mConst} to avoid numerical problems in the FIA function.
}
  \item{numDeriv}{
logical. Should the Hessian matrix of the maximum likelihood estimates be estimated numerically using \code{numDeriv::hessian} in case it cannot be estimated analytically? This can be extremely time and memory consuming for larger models. Default is TRUE.
  }
}
\details{
The model file is either of the easy format or the "classical" EQN format (see below).\cr
In the easy format (the default) the model file contains all trees of the model. Trees are separated by at least one empty line. Everything to the right of a hash (#) is ignored (this behavior is new since version 0.9.2). Lines starting with a # are treated as empty. Each line in each tree corresponds to all branches of this tree (concatenated by a +) that correspond to one of the possible response categories. The position of each line must correspond to the position of this response category in the data object (for multi-individual fit to the respective column).

The difference between both types of EQN format (\code{"eqn"} or\code{"eqn2"}) is the way the first line of the model file is treated. If \code{model.file} is set to \code{"eqn"}, \code{MPTinR} will ignore the first line of the model file and will read the rest of the file (as does multiTree; Moshagen, 2010). If \code{model.file} is set to \code{"eqn2"} MPTinR will only read as many lines as indicated in the first line of the EQN model file (as does e.g., HMMTree; Stahl & Klauer, 2007). As default \code{fit.mpt} expects the easy format, but if the filename ends with .eqn or .EQN and \code{model.type} is \code{"easy"}, \code{model.type} is set to \code{"eqn"}\cr
For the EQN format consult one of the corresponding papers (see e.g., Moshagen, 2010; Stahl & Klauer, 2007). The positions in the data object (number of column for multi-individual fit) must correspond to the category number in the EQN file.\cr

Note that names of parameters in the model file should not start with \code{hank.}. Variables with these names can lead to unforeseen problems as variables starting with these letters are internally used. Furthermore, any \code{\link{reserved}} names (e.g., \code{NA}) are not allowed in model files of any types (i.e., also not as category labels in .eqn files). All names in models need to be valid R variable names (see \code{\link{make.names}}).

The restrictions file may contain (sequential) equality (i.e., =) and inequality (i.e., <) restrictions and must adhere to the following rules:\cr
	1. Inequalities first.\cr
	2. If a variable appears in an inequality restriction, it can not be on the left hand side (LHS) of any further restriction.\cr
	3. If a variable appears on the right hand side (RHS) of an equality restriction, it can not appear on LHS of an equality restriction.\cr
Note that only "<" is supported as inequality operator but not ">"!\cr
Examples of restrictions are (the following could all appear in one restrictions file):\cr
\code{D1 < D2 < D3}\cr
\code{D4 = D3}\cr
\code{B1 = B3 = 0.3333}\cr
\code{X4 = X5 = D3}\cr
Restrictions file may contain comments (i.e., everything to the right of a # will be ignored; new behavior since version 0.9.2)

Restrictions can also be specified in line as a list. The same restrictions as the one above as a \code{list} would be \code{list("D1 < D2 < D3", "D4 = D3", "B1 = B3 = 0.3333", "X4 = X5 = D3")} (simply use this \code{list} as the \code{restrictions.filename} argument).

For equality restrictions, the equality restricted parameters are simply exchanged with their restrictions before the fitting.\cr
For inequality restricted parameters, the model is reparameterized so that only the rightmost parameter of an inequality restriction remains the original parameter. Each instance of the other parameters in this restriction is replaced by the product of the rightmost parameter and dummy parameters (see Knapp & Batchelder, 2004). This procedure (which is equivalent to method A described in Knapp & Batchelder, 2004) leads to an equivalent model (although the binary MPT structure is not apparent in the resulting equations).\cr
To prohibit this reparameterization (i.e., if the inequality restrictions hold without reparameterization), you can set \code{reparam.ineq} to \code{FALSE}. This can be useful for obtaining the FIA (see examples in Wu, Myung, & Batchelder, 2010).

Both models and restrictions can be specified as \code{\link{textConnection}}s instead of as external files.\cr
Furthermore, restrictions can be specified directly as a \code{list} containing the restrictions (quoted, i.e. as characters).\cr
\code{\link{fit.model}} contains additional examples showing model and restrictions specification within the code.

Note that when setting some parameters equal and also restricting their order, the parameters set equal which are not the rightmost element in the order (i.e., inequality) restriction, are computed correctly, but are marked as inequality restricted instead of equality restricted in the output (this did not work at all before v1.0.1). An example: For the restrictions \code{list("G2 < G3 < G5", "G1 = G2", "G4 = G5")}, \code{G1} would be computed correctly, but marked as inequality restricted. In contrast, \code{G4} would be marked as equal to \code{G5} (and also computed correctly).

To obtain a measure of the model's complexity beyond the number of parameters (and taking inequality restrictions into account), set \code{fia} to a (reasonably high) scalar integer (i.e., a number). Then, \code{fit.mpt} will obtain the Fisher Information Approximation (FIA), a Minimum Description Length (MDL) based measure of model complexity, using the algorithm provided by Wu, Myung, & Batchelder (2010a, 2010b) ported from Matlab to R. When performing model-selection, this measure is superior to other methods such as the Akaike Information Criterion (AIC) or Bayesian Information Criterion (BIC) which basically only take the number of parameters into account.\cr
To get the FIA, \code{fit.mpt} performs the following steps:\cr
1. The representation of the model as equations is transformed into the string representation of the model in the context-free language of MPT models (L-BMPT; Purdy & Batchelder, 2009). For this step to be successful it is \emph{absolutely necessary} that the equations representing the model perfectly map the tree structure of the MPT. That is, the model file is only allowed to contain parameters, their inverse (e.g., \code{Dn} and \code{(1 - Dn)}) and the operators + and *, but nothing else. Simplifications of the equations will seriously distort this step. Similarly, unnecessary brackets will distort the results. Brackets must only be used to indicate the inverse of a parameter (i.e. (1 - parameter)). This step is achieved by \code{\link{make.mpt.cf}}.\cr
2. The context free representation of the model is then fed into the MCMC function computing the FIA (the port of BMPTFIA provided by Wu, Myung & Batchelder (2010a), see \code{\link{bmpt.fia}}).\cr
(Actually, both steps are achieved by a call to \code{\link{get.mpt.fia}})

Note that FIA can sometimes be non-consistent (i.e., larger FIA penalty values for restricted versions of a model than for the superordinate model; see Navarro, 2004). This may specifically happens for small ns and is for example the case for the Broder & Schutz example below. In these cases FIA cannot be used! Therefore, always check for consistency of the FIA penalty terms.

Once again: If one wants to compute the FIA, it is \emph{absolutely necessary}, that the representation of the model via equations in the model file exactly maps on the structure of the binary MPT (see \code{\link{make.mpt.cf}} for more details).

Confidence intervals (CI) are based on the observed Hessian matrix produced by the symbolically derived function for the Hessian (i.e., the second derivative of the likelihood function). If it is based on a numerically estimated Hessian, a warning will be given. \cr
For inequality restricted parameters, the CIs are computed using the parameter estimates' variance bounds (see Baldi & Batchelder, 2003; especially Equation 19). Note that these bounds represent the "worst case scenario" variances, and can lead to CIs outside parameter boundaries if the set of inequalities is large and/or the variances for the reparameterized model are large (Note that CIs for non-restricted parameters can be outside the parameter boundaries as well due to large variances).

To avoid local minima and instead find the maximum likelihood estimates it is useful to set \code{n.optim} > 1 with random starting values (see below). If \code{n.optim} > 1, the \code{\link{summary}} of the vector containing the Log-Likelihood values returned by each run of \code{nlminb} is added to the output (to check whether local minima were present). If the model is rather big, \code{n.optim} > 1 can be slow.


Multicore fitting is achieved via the \code{snowfall} package and needs to be initialized via \code{sfInit}. As initialization needs some time, you can either initialize multicore facilities yourself using \code{sfInit()} and setting the \code{sfInit} argument to \code{FALSE} (the default) or let MPTinR initialize multicore facilities by setting the \code{sfInit} argument to \code{TRUE}. The former is recommended as initializing \code{snowfall} takes some time and only needs to be done once if you run \code{fit.mpt} multiple times. If there are any problems with multicore fitting, first try to initialize \code{snowfall} outside MPTinR (e.g., \code{sfInit( parallel=TRUE, cpus=2 )}). If this does not work, the problem is not related to MPTinR but to snowfall (for support and references visit: \url{http://www.imbi.uni-freiburg.de/parallel/}).\cr
Note that you should  \emph{close} snowfall via \code{sfStop()} after using MPTinR.

The fitting/optimization is achieved via \code{\link{nlminb}} (Fox, Hall, & Schryer, 1978) a Newton based algorithm using the analytically derived gradient. In some cases (e.g., in case of empty cells) \code{nlminb} will not converge successfully in which \code{fit.mpt} will retry fitting using a numerically estimated gradient (with warning). 

\code{fit.mpt} is just a comfortable wrapper around the workhorse \code{\link{fit.mptinr}}. \code{fit.mpt} produces the appropriate objective function, gradient function, hessian function, and prediction function that are handed over to \code{fit.mptinr} (functions are produced by symbolical derivation, see \code{\link{D}}). A function similar to \code{fit.mpt} is \code{\link{fit.model}} which has the additional arguments \code{lower.bound} and \code{upper.bound} allowing to fit other models than just MPTs and the possibility to indicate whether or not to use the analytically derived gradient or hessian for fitting (here this is automatically handled). Note that for MPTs (where upper and lower bounds of parameters are set to 0 and 1, respectively) \code{fit.mpt} is probably faster as the objective function is slightly faster (i.e., more optimized). However, for datasets with many empty cells trying \code{\link{fit.model}} with or without gradient or hessian can be worth a try.

Note that \code{\link{fit.mptinr}} can fit models with arbitrary (i.e., custom) objective functions.

The old version of this function using \code{\link{optim}}'s L-BFGS-B algorithm is \code{\link{fit.mpt.old}}.

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
In these cases it is recommended to rerun fit.mpt to check if the results are stable.
}
\references{
Baldi, P. & Batchelder, W. H. (2003). Bounds on variances of estimators for multinomial processing tree models. \emph{Journal of Mathematical Psychology}, 47, 467-470. 

Broeder, A., & Schuetz, J. (2009). Recognition ROCs are curvilinear-or are they? On premature arguments against the two-high-threshold model of recognition. \emph{Journal of Experimental Psychology: Learning, Memory, and Cognition}, 35(3), 587. doi:10.1037/a0015279

Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. \emph{SIAM J. Scientific Computing}, 16, 1190-1208.

Fox, P. A., Hall, A. P., & Schryer, N. L. (1978). The PORT Mathematical Subroutine Library. \emph{CM Trans. Math. Softw.}, 4, 104-126. \url{http://doi.acm.org/10.1145/355780.355783}

Knapp, B. R., & Batchelder, W. H. (2004). Representing parametric order constraints in multi-trial applications of multinomial processing tree models. \emph{Journal of Mathematical Psychology}, 48, 215-229.

Moshagen, M. (2010). multiTree: A computer program for the analysis of multinomial processing tree models. \emph{Behavior Research Methods}, 42, 42-54.

Navarro, D. J. (2004). A Note on the Applied Use of MDL Approximations. \emph{Neural Computation}, 16(9), 1763-1768.

Purdy, B. P., & Batchelder, W. H. (2009). A context-free language for binary multinomial processing tree models. \emph{Journal of Mathematical Psychology}, 53, 547-561.

Riefer, D. M., & Batchelder, W. H. (1988). Multinomial modeling and the measurement of cognitive processes. \emph{Psychological Review}, 95, 318-339.

Stahl, C. & Klauer, K. C. (2007). HMMTree: A computer program for latent-class hierarchical multinomial processing tree models. \emph{Behavior Research Methods}, 39, 267- 273.

Wu, H., Myung, J.I., & Batchelder, W.H. (2010a). Minimum description length model selection of multinomial processing tree models. \emph{Psychonomic Bulletin & Review}, 17, 275-286.

Wu, H., Myung, J.I., & Batchelder, W.H. (2010b). On the minimum description length complexity of multinomial processing trees. \emph{Journal of Mathematical Psychology}, 54, 291-303.
}
\author{
Henrik Singmann and David Kellen with help from Karl Christoph Klauer.
}
\note{
All (model or restriction) files should end with an empty line, otherwise a warning will be shown.

The likelihood returned does not include the factorial constants of the multinomial probability-mass functions.
}

\examples{
# The first example fits the MPT model presented in Riefer and Batchelder (1988, Figure 1)
# to the data presented in Riefer and Batchelder (1988, Table 1)
# Note that Riefer and Batchelder (1988, pp. 328) did some hypotheses tests not replicated here.
# Instead, we use each condition (i.e., row in Table 1) as a different dataset.

# load the data
data(rb.fig1.data, package = "MPTinR")

#get the character string with the position of the model:
model1 <- system.file("extdata", "rb.fig1.model", package = "MPTinR")
model1.eqn <- system.file("extdata", "rb.fig1.model.eqn", package = "MPTinR")

# just fit the first dataset:
fit.mpt(rb.fig1.data[1,], model1, n.optim = 1)
fit.model(rb.fig1.data[1,], model1, n.optim = 1)

#fit all datasets:
fit.mpt(rb.fig1.data, model1, n.optim = 1)
fit.model(rb.fig1.data, model1, n.optim = 1)

#fit all datasets using the .EQN model file:
fit.mpt(rb.fig1.data, model1.eqn, n.optim = 1)

#fit using a textConnection (i.e., you can specify the model in your script/code):
model1.txt <- "p * q * r
p * q * (1-r)
p * (1-q) * r
p * (1-q) * (1-r) + (1-p)"
fit.mpt(rb.fig1.data, textConnection(model1.txt), n.optim = 1)



# The second example fits the MPT model presented in Riefer and Batchelder (1988, Figure 2)
# to the data presented in Riefer and Batchelder (1988, Table 3)
# First, the model without restrictions is fitted: ref.model
# Next, the model with all r set equal is fitted: r.equal
# Then, the model with all c set equal is fitted: c.equal
# Finally, the inferential tests reported by Riefer & Batchelder, (1988, p. 332) are executed.

# get the data
data(rb.fig2.data, package = "MPTinR")

# positions of model and restriction files:
model2 <- system.file("extdata", "rb.fig2.model", package = "MPTinR")
model2r.r.eq <- system.file("extdata", "rb.fig2.r.equal", package = "MPTinR")
model2r.c.eq <- system.file("extdata", "rb.fig2.c.equal", package = "MPTinR")

# The full (i.e., unconstrained) model
(ref.model <- fit.mpt(rb.fig2.data, model2))

# All r equal
(r.equal <- fit.mpt(rb.fig2.data, model2, model2r.r.eq))

# All c equal
(c.equal <- fit.mpt(rb.fig2.data, model2, model2r.c.eq))

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

# You can specify restrictions also via a list instead of an external file:
# All r equal
r.equal.2 <- fit.mpt(rb.fig2.data, model2, list("r0 = r1 = r2= r3 = r4"), n.optim = 5)
all.equal(r.equal, r.equal.2)

# All c equal
c.equal.2 <- fit.mpt(rb.fig2.data, model2, list("c0 = c1 = c2 = c3= c4"))
all.equal(c.equal, c.equal.2)


\dontrun{

# Example from Broder & Schutz (2009)
# We fit the data from the 40 individuals from their Experiment 3
# We fit three different models:
# 1. Their 2HTM model: br.2htm
# 2. A restricted 2HTM model with Dn = Do: br.2htm.res
# 3. A 1HTM model (i.e., Dn = 0): br.1htm
# We fit the models with, as well as without, applied inequality restrictions (see Details)
# that is, for some models (.ineq) we impose: G1 < G2 < G3 < G4 < G5 
# As will be apparent, the inequality restrictions do not hold for all individuals.
# Finally, we compute the FIA for all models, taking inequalities into account.

data(d.broeder, package = "MPTinR")
m.2htm <- system.file("extdata", "5points.2htm.model", package = "MPTinR")
r.2htm <- system.file("extdata", "broeder.2htm.restr", package = "MPTinR")
r.1htm <- system.file("extdata", "broeder.1htm.restr", package = "MPTinR")
i.2htm <- system.file("extdata", "broeder.2htm.ineq", package = "MPTinR")
ir.2htm <- system.file("extdata", "broeder.2htm.restr.ineq", package = "MPTinR")
ir.1htm <- system.file("extdata", "broeder.1htm.restr.ineq", package = "MPTinR")

# fit the original 2HTM
br.2htm <- fit.mpt(d.broeder, m.2htm)
br.2htm.ineq <- fit.mpt(d.broeder, m.2htm, i.2htm)

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

# identical to the last fit of the 1HTM (using a list as restriction):
br.1htm.ineq.list <- fit.mpt(d.broeder, m.2htm, list("G1 < G2 < G3 < G4 < G5", "Dn = 0"))
all.equal(br.1htm.ineq, br.1htm.ineq.list)  # TRUE

# These results show that inequality restrictions do not hold for all datasets.
# (It would look differently if we excluded critical cases, 
# i.e., 2, 6, 7, 10, 18, 21, 25, 29, 32, 34, 35, 37, 38)
# Therefore, we get the FIA for the models as computed above 

br.2htm.fia <- fit.mpt(d.broeder, m.2htm, fia = 200000)
br.2htm.ineq.fia <- fit.mpt(d.broeder, m.2htm, i.2htm, fia = 200000)
br.2htm.res.fia <- fit.mpt(d.broeder, m.2htm, r.2htm, fia = 200000 )
br.2htm.res.ineq.fia <- fit.mpt(d.broeder, m.2htm, ir.2htm, fia = 200000)
br.1htm.fia <- fit.mpt(d.broeder, m.2htm, r.1htm, fia = 200000)
br.1htm.ineq.fia <- fit.mpt(d.broeder, m.2htm, ir.1htm, fia = 200000)

# Model selection using the FIA
(br.select <- select.mpt(list(br.2htm.fia, br.2htm.ineq.fia, br.2htm.res.fia, 
                              br.2htm.res.ineq.fia, br.1htm.fia, br.1htm.ineq.fia)))
                              
# The same results, ordered by FIA
br.select[order(br.select[,"delta.FIA.sum"]),]

# Note that FIA for individual data (.sum) is not consistent (i.e., the penalty
# for the nested model br.1htm.ineq.fia is not really smaller than the penalty
# for the superordinate model br.2htm.ineq.fia).
# Hence, one should use the aggregated data for this analysis (not shown here)! 

# Compare this with the model selection not using FIA:
select.mpt(list(br.2htm, br.2htm.ineq, br.2htm.res, br.2htm.res.ineq, br.1htm, br.1htm.ineq))

# Only use the aggregated data:
d.broeder.agg <- colSums(d.broeder)
br.2htm.agg <- fit.mpt(d.broeder.agg, m.2htm)
br.2htm.res.agg <- fit.mpt(d.broeder.agg, m.2htm, r.2htm)
br.1htm.agg <- fit.mpt(d.broeder.agg, m.2htm, r.1htm)

select.mpt(list(br.2htm.agg, br.2htm.res.agg, br.1htm.agg), output = "full")


# compare speed of no multicore versus multicore for multiple datasets:

require(snowfall)
# change number of CPUs if more are available
nCPU = 2
sfInit( parallel=TRUE, cpus=nCPU, type = "SOCK" )

# NO multicore
system.time(fit.mpt(d.broeder, m.2htm))

# multicore:
system.time(fit.mpt(d.broeder, m.2htm, multicore = "individual"))

sfStop()
}

  }
\seealso{
\code{\link{check.mpt}} for a function that can help in constructing models.

\code{\link{select.mpt}} for the function that performs model selection on results from \code{fit.mpt}.

\code{\link{fit.model}} for a similar wrapper for which you can specify upper and lower bounds of parameters (and whether or not \code{nlminb} uses the symbolically derived gradient and hessian)

\code{\link{fit.mptinr}} is the workhorse with which you can also fit your own objective functions.
}
  
\keyword{models}
\keyword{tree}

