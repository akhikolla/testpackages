\name{fDMA}
\alias{fDMA}

\title{Computes Dynamic Model Averaging.}

\description{
The function estimates Dynamic Model Averaging (and some of its variations). The method is described in Raftery et al. (2010).
}

\usage{
fDMA(y,x,alpha,lambda,initvar,W=NULL,initial.period=NULL,V.meth=NULL,kappa=NULL,
gprob=NULL,omega=NULL,model=NULL,parallel=NULL,m.prior=NULL,mods.incl=NULL,
DOW=NULL,DOW.nmods=NULL,DOW.type=NULL,DOW.limit.nmods=NULL,progress.info=NULL,
forced.models=NULL,forbidden.models=NULL,forced.variables=NULL,bm=NULL,
small.c=NULL,fcores=NULL,mods.check=NULL,red.size=NULL,av=NULL)
}

\arguments{
\item{y}{\code{\link[base]{numeric}} or a column \code{\link[base]{matrix}} of a dependent variable, if \code{y} is \code{\link[xts]{xts}} object, then plots will have time index on the x axis}
\item{x}{\code{\link[base]{matrix}} of independent variables, different columns correspond to different variables}
\item{alpha}{\code{\link[base]{numeric}}, a forgetting factor \eqn{\alpha} between 0 and 1 used in probabilities estimations}
\item{lambda}{\code{\link[base]{numeric}}, a forgetting factor \eqn{\lambda} between 0 and 1 used in variance approximations}
\item{initvar}{\code{\link[base]{numeric}}, initial variance in the state space equation, i.e., the number by which the unit matrix is multiplied}
\item{W}{optional, a method for setting the initial values of variance for the models equations, \code{W="reg"} corresponds to the method based on the linear regression as in the paper by Raftery et al. (2010), alternatively an arbitrary positive number (\code{\link[base]{numeric}}) can be specified, by default the method of Raftery et al. (2010) is used}
\item{initial.period}{optional, \code{\link[base]{numeric}}, a number of observation since which MSE (Mean Squared Error) and MAE (Mean Absolute Error) are computed, by default the whole sample is used, i.e., \code{initial.period=1}}
\item{V.meth}{optional, a method for the state space equation variance updating, \code{V.meth="rec"} corresponds to the recursive moment estimator, as in the paper by Raftery et al. (2010), \code{V.meth = "ewma"} corresponds to the exponentially weighted moving average as in, for example, Koop and Korobilis (2012), by default \code{V.meth = "rec"} is used}
\item{kappa}{optional, \code{\link[base]{numeric}}, a parameter in the exponentially weighted moving average, between 0 and 1, used if \code{V.meth = "ewma"}}
\item{gprob}{optional, \code{\link[base]{matrix}} of Google probabilities as in Koop and Onorante (2014), columns should correspond to columns of x}
\item{omega}{optional, \code{\link[base]{numeric}}, a parameter between 0 and 1 used in probabilities estimations, used if \code{gprob} is specified}
\item{model}{optional, \code{model="dma"} for Dynamic Model Averaging, \code{model="dms"} for Dynamic Model Selection, or \code{model="med"} for Median Probability Model as in Barbieri and Berger (2004), by default \code{model="dma"} is used}
\item{parallel}{optional, \code{\link[base]{logical}}, indicate whether parallel computations should be used, by default \code{parallel=FALSE}}
\item{m.prior}{optional, \code{\link[base]{numeric}}, a parameter for general model prior (Mitchell and Beauchamp, 1988), by default \code{m.prior=0.5}, which corresponds to the uniform distribution, i.e., non-informative priors, see also Eicher et al. (2011)}
\item{mods.incl}{optional, \code{\link[base]{matrix}} indicating which models should be used for estimation, the first column indicates inclusion of a constant, by default all possible models with a constant are used, inclusion of a variable is indicated by 1, omitting by 0}
\item{DOW}{optional, \code{\link[base]{numeric}}, a threshold for Dynamic Occam's Window (Onorante and Raftery, 2016), should be a number between 0 and 1, if \code{DOW=0}, then no Dynamic Occam's Window is applied, by default \code{DOW=0}, Dynamic Occam's Window can be applied only to Dynamic Model Averaging, i.e., when \code{model="dma"}}
\item{DOW.nmods}{optional, \code{\link[base]{numeric}}, initial number of models for Dynamic Occam's Window, should be less than the number of all possible models and larger than or equal to 2, they are randomly chosen, if \code{DOW.nmods=0}, then initially models with exactly one variable are taken, by default \code{DOW.nmods=0}}
\item{DOW.type}{optional, \code{DOW.type="r"} corresponds to DMA-R from Onorante and Raftery (2016), \code{DOW.type="e"} to DMA-E, by default \code{DOW.type="r"}}
\item{DOW.limit.nmods}{optional, \code{\link[base]{numeric}}, maximum number of models selected by Dynamic Occam's Window, an additional limitation to the threshold given by \code{DOW}, by default no limit is set}
\item{progress.info}{optional, \code{\link[base]{logical}}, applicable only if Dynamic Occam's Window is used, otherwise ignored, if \code{progress.info=TRUE} number of the current recursive DMA computation round and number of models selected for this round are printed, by default \code{progress.info=FALSE}}
\item{forced.models}{optional, \code{\link[base]{matrix}}, applicable only if Dynamic Occam's Window is used, otherwise ignored, indicates models that have to be always included in the set of expanded models, similar as \code{mods.incl}, by default \code{forced.models=NULL}}
\item{forbidden.models}{optional, \code{\link[base]{matrix}}, applicable only if Dynamic Occam's Window is used, otherwise ignored, indicates models that cannot be used in the set of expanded models, similar as \code{mods.incl}, by default \code{forbidden.models=NULL}}
\item{forced.variables}{optional, \code{\link[base]{vector}}, applicable only if Dynamic Occam's Window is used, otherwise ignored, indicates variables that have to be always included in models constituting the set of expanded models, similar as \code{mods.incl}, first slot indicates inclusion of constant, by default \code{forced.variables=NULL}}
\item{bm}{optional, \code{\link[base]{logical}}, indicate whether benchmark forecast should be computed, these benchmarks are naive forecast (all forecasts are set to be the value of the last observation) and Auto Arima \code{\link[forecast]{auto.arima}}, by default \code{bm=FALSE}}
\item{small.c}{optional, \code{\link[base]{numeric}}, small constant added to posterior model proabilities as in Raftery et al. (2010) to prevent potential reduction them to 0 due to the computational issues, if not specified the value computed as in Raftery et al. (2010) is taken}
\item{fcores}{optional, \code{\link[base]{numeric}}, used only if \code{parallel=TRUE}, otherwise ignored, indicates the number of cores that should not be used, by default \code{fcores=1}}
\item{mods.check}{optional, \code{\link[base]{logical}}, indicates if \code{mods.incl} should be checked for duplicated entries, etc., by default \code{mods.check=FALSE}}
\item{red.size}{optional, \code{\link[base]{logical}}, indicates if outcomes should be reduced to save memory, by default \code{red.size=FALSE}} 
\item{av}{optional, \code{av="dma"} corresponds to the original DMA averaging scheme, \code{av="mse"} corresponds to averaging based on Mean Squared Error, \code{av="hr1"} corresponds to averaging based on Hit Ratio, assuming time-series are in levels, \code{av="hr2"} corresponds to averaging based on Hit Ratio, assuming time-series represent changes, by default \code{av="dma"}}
}

\details{It is possible to use \code{\link[base]{numeric}} \code{\link[base]{vector}} for \code{lambda}. Its values are automatically ordered in descending order and if numbers are not unique they are reduced to become unique. If more than one value is given for \code{lambda}, then model state space, i.e., \code{mods.incl}, is expanded by considering all these models with given values of \code{lambda}. The outcomes are then ordered by columns in a way that first outcomes from models with first value of \code{lambda} are presented, then from models with second value of \code{lambda}, etc. (Raftery et al., 2010).

If \code{nrow(gprob)<length(y)}, then the method by Koop and Onorante (2014) is used for the last \code{nrow(gprob)} observations. For the preceding ones the original method by Raftery et al. (2010) is used. In such case a \code{warning} is generated.
}

\value{class \code{dma} object, \code{\link[base]{list}} of
\item{$y.hat}{forecasted values}
\item{$post.incl}{posterior inclusion probabilities for independent variables}
\item{$MSE}{Mean Squared Error of forecast}
\item{$MAE}{Mean Absolute Error of forecast}
\item{$models}{models included in estimations, or models used in the last step of Dynamic Occam's Window method (if this method has been selected)}
\item{$post.mod}{posterior probabilities of all used models, or \link{NA} if Dynamic Occam's Window method has been selected}
\item{$exp.var}{expected number of variables (incl. constant)}
\item{$exp.coef.}{expected values of regression coefficients}
\item{$parameters}{parameters of the estimated model}
\item{$yhat.all.mods}{predictions from every sub-model used in estimations}
\item{$y}{\code{y}, dependent variable}
\item{$benchmarks}{Root Mean Squared Error and Mean Absolute Error of naive and auto ARIMA forecast}
\item{$DOW.init.mods}{models initially selected to Dynamic Occam's Window, if this method has been selected}
\item{$DOW.n.mods.t}{number of models used in Dynamic Model Averaging at time \eqn{t}, if Dynamic Occam's Window method has been selected}
\item{$p.dens.}{predicitive densities from the last period of all sub-models used in estimations}
\item{$exp.lambda}{expected values of \code{lambda} parameter}
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10)
m2 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,model="dms")
m3 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,V.meth="ewma",kappa=0.9)
m4 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,DOW=0.7)
}
\dontshow{
s <- ld.drivers[1:20,1:4]
mod1 <- fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10)
mod2 <- fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,model="dms")
mod3 <- fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,V.meth="ewma",kappa=0.9)
mod4 <- fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,DOW=0.7)
}
# Google trends are available since 2004
\dontrun{
gp <- trends/100
s <- ld.drivers['2004-01-01/']
m5 <- fDMA(y=ld.wti['2004-01-01/'],x=s,alpha=0.99,lambda=0.90,initvar=10,gprob=gp,omega=0.5)
}
\dontshow{
s <- ld.drivers['2004-01-01/'][1:20,1:4]
gp <- trends/100
gp <- gp[1:20,1:4]
mod5 <- fDMA(y=ld.wti['2004-01-01/'][1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,gprob=gp,omega=0.5)
}
# models just with one independent variable and a constant will be averaged
mds <- diag(1,ncol(ld.drivers),ncol(ld.drivers))
mds <- cbind(rep(1,ncol(ld.drivers)),mds)
\dontrun{
m6 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,mods.incl=mds)
}
\dontshow{
s <- ld.drivers[1:20,1:4]
mod6 <- fDMA(y=ld.wti[1:20,1],x=s,alpha=0.99,lambda=0.90,initvar=10,mods.incl=mds[1:5,1:5])
}
# models just with one independent variable (without a constant) will be averaged
mds.nc <- diag(1,ncol(ld.drivers),ncol(ld.drivers))
mds.nc <- cbind(rep(0,ncol(ld.drivers)),mds.nc)
\dontrun{
m7 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,mods.incl=mds.nc)
}
# model with multiple lambda
\dontrun{
m8 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=c(0.99,0.95,0.90),initvar=10)
}
}

\source{
Raftery, A. E., Karny, M., Ettler, P., 2010. Online prediction under model uncertainty via Dynamic Model Averaging: Application to a cold rolling mill. \emph{Technometrics} \bold{52}, 52--66.
}

\references{
Barbieri, M. M., Berger, J. O., 2004. Optimal predictive model selection. \emph{The Annals of Statistics} \bold{32}, 870--897.

Eicher, T. S., Papageorgiou, C., Raftery, A. E., 2011. Default priors and predictive performance in Bayesian Model Averaging, with application to growth determinants. \emph{Journal of Applied Econometrics} \bold{26}, 30--55.

Koop, G., Korobilis, D., 2012. Forecasting inflation using Dynamic Model Averaging. \emph{International Economic Review} \bold{53}, 867--886.

Koop, G., Korobilis, D., 2018. Variational Bayes inference in high-dimensional time-varying parameter models. \url{https://arxiv.org/pdf/1809.03031}

Koop, G., Onorante, L., 2014. Macroeconomic nowcasting using Google probabilities. \url{https://goo.gl/ATsBN9}

Mitchell, T. J., Beauchamp, J. J., 1988. Bayesian variable selection in linear regression (with discussion). \emph{Journal of the American Statistical Association} \bold{83}, 1023--1036.

Onorante, L., Raftery, A. E., 2016. Dynamic model averaging in large model spaces using dynamic Occam's window. \emph{European Economic Review} \bold{81}, 2--14.

Yin, X., Peng, J., Tang, T., 2018. Improving the forecasting accuracy of crude oil prices. \emph{Sustainability} \bold{10}, 454. \url{https://doi.org/10.3390/su10020454}
}

\seealso{\code{\link{grid.DMA}}, \code{\link{print.dma}}, \code{\link{summary.dma}}, \code{\link{plot.dma}}, \code{\link{hit.ratio}}.}
