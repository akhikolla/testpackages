\name{fit.mptinr}
\alias{fit.mptinr}
\title{
Fit cognitive models for categorical data using an objective function
}
\description{
Fitting function for package MPTinR. Can fit any model for categorical data specified in an objective function. Fitting can be enhanced with gradient and or Hessian. Predicted values will be added when a prediction function is present.
}
\usage{

fit.mptinr(
	data,
	objective, 
	param.names,
	categories.per.type,
	gradient = NULL, use.gradient = TRUE,
	hessian = NULL, use.hessian = FALSE,
	prediction = NULL,
	n.optim = 5,
	fia.df = NULL,
	ci = 95, 
	starting.values = NULL,
	lower.bound = 0,
	upper.bound = 1,
	output = c("standard", "fia", "full"),
	fit.aggregated = TRUE,
	sort.param = TRUE,
	show.messages = TRUE,
	use.restrictions = FALSE,
	orig.params = NULL,
	restrictions = NULL,	
	multicore = c("none", "individual", "n.optim"), sfInit = FALSE, nCPU = 2,
	control = list(),
    numDeriv = TRUE,
	...
)
}

\arguments{
  \item{data}{
Either a \emph{numeric} \code{vector} for individual fit or a \emph{numeric} \code{matrix} or \code{data.frame} for multi-dataset fit. The data on each position (column for multi-dataset fit) must correspond to the respective line in the model file. Fitting for multiple datasets can be parallelized via \code{multicore}.
}
  \item{objective}{
the objective function used for fitting. Needs to return a scalar likelihood value.
}
  \item{param.names}{
character vector giving the parameters present in the model. The order of this vector determines the order with which the output from the fitting routine is interpreted.
}
\item{categories.per.type}{
numeric vector indicating how many response categories each item type has.
}
  \item{gradient}{
the gradient function used for fitting. Needs to return a vector of same length as \code{param.names}.
}
  \item{use.gradient}{
logical. indicating whether or not \code{gradient} should be used for fitting. Default is \code{TRUE}
}
  \item{hessian}{
the Hessian function used for fitting. Needs to return a matrix with \code{dim = c(length(param.names), length(param.names)}.
}
  \item{use.hessian}{
logical. indicating whether or not \code{hessian} should be used for fitting. Default is \code{FALSE}
}
\item{prediction}{
the prediction function. Needs to return a vector of equal length as the response categories or data. Needs to return probabilities!
}
  \item{n.optim}{
Number of optimization runs. Can be parallelized via \code{multicore}. Default is 5. 
}
\item{fia.df}{
needed for handling MPTs with computation of FIA coming from \code{\link{fit.mpt}} or \code{fit.model}. Do not use.
}
  \item{ci}{
A scalar corresponding to the size of the confidence intervals for the parameter estimates. Default is 95 which corresponds to 95\% confidence intervals.
}
  \item{starting.values}{
A \code{vector}, a \code{list}, or \code{NULL} (the default). If \code{NULL} starting values for parameters are randomly drawn from a uniform distribution with the interval \code{(0.1 - 0.9)}. See Details for the other options.
} 
  \item{lower.bound}{
numeric scalar or vector. Can be used in \code{fit.model} to set the lower bounds of the parameter space. See Details.
}
  \item{upper.bound}{
numeric scalar or vector. Can be used in \code{fit.model} to set the upper bounds of the parameter space. See Details.
}
  \item{output}{
If "full" \code{fit.mpt} will additionally return the return values of \code{\link{nlminb}} and the Hessian matrices. (If "fia", \code{fit.mpt} will additionally return the results from \code{\link{get.mpt.fia}} (if \code{fia} not equal \code{NULL}).)
}
\item{fit.aggregated}{
logical. Only relevant for multiple datasets (i.e., \code{matrix} or \code{data.frame}). Should the aggregated dataset (i.e., data summed over rows) be additionally fitted? Default (\code{TRUE}) fits the aggregated data. 
}
\item{sort.param}{
Logical. If TRUE, parameters are alphabetically sorted in the parameter table. If FALSE, the first parameters in the parameter table are the non-restricted ones, followed by the restricted parameters. Default is TRUE.
}
    \item{show.messages}{
Logical. If TRUE the time the fitting algorithms takes is printed to the console.
}
\item{use.restrictions}{
needed for handling MPTs coming from \code{\link{fit.mpt}}. Do not use, unless you are sure what you are doing.
}
\item{orig.params}{
needed for handling models coming from \code{\link{fit.mpt}} or \code{\link{fit.model}}. Do not use, unless you are sure what you are doing.
}
\item{restrictions}{
needed for handling models coming from \code{\link{fit.mpt}} or \code{\link{fit.model}}. Do not use, unless you are sure what you are doing.
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

  \item{control}{
list containing control arguments passed on to \code{\link{nlminb}}. See there.
}
  \item{numDeriv}{
logical. Should the Hessian matrix of the maximum likelihood estimates be estimated numerically using \code{numDeriv::hessian} in case it cannot be estimated analytically? This can be extremely time and memory consuming for larger models. Default is TRUE.
  }
  \item{...}{
arguments passed on to the objective function, the gradient function, the hessian function and the prediction function.
}

}
\details{
This functions can be used to fit any model for categorical data that can be specified via a (objective) function (i.e., especially models that are not MPTs). For fitting MPTs or other similar models such as SDTs see \code{\link{fit.mpt}} or \code{\link{fit.model}}.

The only mandatory arguments are: \code{data}, \code{objective}, \code{param.names}, and \code{categories.per.type}. Adding a function calculating the \code{gradient} will usually significantly speed up the fitting. However, in extreme cases (i.e., many empty cells) using the gradient can interfere with finding the global minima. Adding the function computing the \code{hessian} matrix is usually only useful for obtaining the accurate confidence intervals (usually the numerically estimated Hessian matrix is equivalent unless there are many empty cells or parameters at the boundary).

The \code{objective} (and \code{gradient} and \code{hessian}) function need to take as the first argument a numerical vector of \code{length(param.names)} representing the parameters. The other mandatory arguments for these functions are:\cr
\code{data}: A vector containing the data for the dataset being fitted.\cr
\code{param.names}: The character vector containing the parameter names is handed over to the objective.\cr
\code{n.params}: = \code{length(param.names)}. To speed up computation the number of parameters is also handed over to the objective on each iteration.\cr
\code{tmp.env}: A \code{\link{environment}} (created with \code{new.env}). The objective function produced by \code{fit.mpt} assign the parameter values into this environment using the following statement:\cr
\code{for (i in seq_len(n.params))  assign(param.names[i],Q[i], envir = tmp.env)}\cr
Furthermore, \code{fit.mptinr} assigns the data points before fitting each dataset into \code{tmp.env} with the variables names \code{hank.data.x} where \code{x} is the ordinal number of that data point (i.e., position or column). In other words, you can use \code{tmp.env} to \code{\link{eval}} you model within this environment and access both parameters and data in it. \cr
\code{lower.bound} and \code{upper.bound}: both \code{lower.bound} and \code{upper.bound} will be passed on to the user-supplied functions as when nlminb fits without \code{gradient} it can try to use parameter values outside the bounds. This can be controlled with these arguments isnide the objective function.

Furthermore, note that all arguments passed via \code{...} will be passed to \code{objective}, \code{gradient}, and \code{hessian}. And that these three functions need to take the same arguments. Furthermore \code{gradient} must return a vector as long as \code{param.names} and \code{hessian} must return a square matrix of order \code{length(param.names}). See \code{\link{nlminb}} for (slightly) more info.

Usage of \code{gradient} and/or \code{hessian} can be controlled with \code{use.gradient} and \code{use.hessian}.

\code{prediction} is a function similar to \code{objective} with the difference that it should return a vector of length \code{sum(categories.per.type} giving the probabilities for each item type. This function needs to take the same arguments as \code{objective} with the only exception that it does not take \code{lower.bound} and \code{upper.bound} (but \code{...} is passed on to it).

Note that parameters names should not start with \code{hank.}. 

To set the starting values for the fitting process (e.g., to avoid local minima) one can set \code{starting.values} to a vector of length 2 and \code{n.optim > 1}. Then, starting values are randomly drawn from a uniform distribution from \code{starting.values[1]} to \code{starting.values[2]}.

Alternatively, one can supply a list with two elements to \code{starting.values}. Both elements need to be either of length 1 or of length equal to the number of parameters (if both are of length 1, it is the same as if you supply a vector of length 2). For each parameter n (in alphabetical order), a starting value is randomly drawn from a uniform distribution \code{starting.values[[1]][n]} to \code{starting.values[[2]][n]} (if length is 1, this is the border for all parameters).

The least interesting option is to specify the starting values individually by supplying a vector with the same length as the number of parameters. Starting values must be ordered according to the alphabetical order of the parameter names. Use \code{\link{check.mpt}} for a function that returns the alphabetical order of the parameters. If one specifies the starting values like that, \code{n.optim} will be set to 1 as all other values would not make any sense (the optimization routine will produce identical results with identical starting values).

Multicore fitting is achieved via the \code{snowfall} package and needs to be initialized via \code{sfInit}. As initialization needs some time, you can either initialize multicore facilities yourself using \code{sfInit()} and setting the \code{sfInit} argument to \code{FALSE} (the default) or let MPTinR initialize multicore facilities by setting the \code{sfInit} argument to \code{TRUE}. The former is recommended as initializing \code{snowfall} takes some time and only needs to be done once if you run \code{fit.mpt} multiple times. If there are any problems with multicore fitting, first try to initialize \code{snowfall} outside MPTinR (e.g., \code{sfInit( parallel=TRUE, cpus=2 )}). If this does not work, the problem is not related to MPTinR but to snowfall (for support and references visit: \url{http://www.imbi.uni-freiburg.de/parallel/}).\cr
Note that you need to \emph{close} snowfall via \code{sfStop()} after using MPTinR.

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
In these cases it is recommended to rerun the model to check if the results are stable.
}
\references{
Kellen, D., Klauer, K. C., & Singmann, H. (2012). On the Measurement of Criterion Noise in Signal Detection Theory: The Case of Recognition Memory. \emph{Psychological Review}. doi:10.1037/a0027727

}
\author{
Henrik Singmann and David Kellen.
}
\note{
All (model or restriction) files should end with an empty line, otherwise a warning will be shown.
}

\examples{
\dontrun{
# the example may occasionally fail due to a starting values - integration mismatch.

# Fit an SDT for a 4 alternative ranking task (Kellen, Klauer, & Singmann, 2012).

ranking.data <- structure(c(39, 80, 75, 35, 61, 54, 73, 52, 44, 63, 40, 48, 80,
49, 43, 80, 68, 53, 81, 60, 60, 65, 49, 58, 69, 75, 71, 47, 44,
85, 23, 9, 11, 21, 12, 21, 14, 20, 19, 15, 29, 13, 14, 15, 22,
11, 12, 16, 13, 20, 20, 9, 26, 19, 13, 9, 14, 15, 24, 9, 19,
7, 9, 26, 16, 14, 6, 17, 21, 14, 20, 18, 5, 19, 17, 5, 11, 21,
4, 9, 15, 17, 7, 17, 11, 11, 9, 19, 20, 3, 19, 4, 5, 18, 11,
11, 7, 11, 16, 8, 11, 21, 1, 17, 18, 4, 9, 10, 2, 11, 5, 9, 18,
6, 7, 5, 6, 19, 12, 3), .Dim = c(30L, 4L))

expSDTrank <- function(Q, param.names, n.params, tmp.env){
   
    e <- vector("numeric",4)

    mu <- Q[1]
    ss <- Q[2]
       
    G1<-function(x){
        ((pnorm(x)^3)*dnorm(x,mean=mu,sd=ss))
    }

    G2<-function(x){
        ((pnorm(x)^2)*dnorm(x,mean=mu,sd=ss)*(1-pnorm(x)))*3
    }

     G3<-function(x){
        (pnorm(x)*dnorm(x,mean=mu,sd=ss)*(1-pnorm(x))^2)*3
    }
 

    e[1] <- integrate(G1,-Inf,Inf,rel.tol = .Machine$double.eps^0.5)$value    
    e[2] <- integrate(G2,-Inf,Inf,rel.tol = .Machine$double.eps^0.5)$value
    e[3] <- integrate(G3,-Inf,Inf,rel.tol = .Machine$double.eps^0.5)$value
    e[4] <- 1-e[1]-e[2]-e[3]  
   
    return(e)
}



SDTrank <- function(Q, data, param.names, n.params, tmp.env, lower.bound, upper.bound){
   
    e<-vector("numeric",4)

    mu <- Q[1]
    ss <- Q[2]
       
    G1<-function(x){
        ((pnorm(x)^3)*dnorm(x,mean=mu,sd=ss))
    }

    G2<-function(x){
        ((pnorm(x)^2)*dnorm(x,mean=mu,sd=ss)*(1-pnorm(x)))*3
    }

     G3<-function(x){
        (pnorm(x)*dnorm(x,mean=mu,sd=ss)*(1-pnorm(x))^2)*3
    }
 

    e[1] <- integrate(G1,-Inf,Inf,rel.tol = .Machine$double.eps^0.5)$value    
    e[2] <- integrate(G2,-Inf,Inf,rel.tol = .Machine$double.eps^0.5)$value
    e[3] <- integrate(G3,-Inf,Inf,rel.tol = .Machine$double.eps^0.5)$value
    e[4] <- 1-e[1]-e[2]-e[3]  
   
    LL <- -sum(data[data!=0]*log(e[data!=0]))
    return(LL)
}

fit.mptinr(ranking.data, SDTrank, c("mu", "sigma"), 4, prediction = expSDTrank, 
		lower.bound = c(0,0.1), upper.bound = Inf)
 }
}

\seealso{
\code{\link{fit.model}} or \code{\link{fit.mpt}} for a function that can fit model represented in a model file.
}
  
\keyword{models}
\keyword{tree}

