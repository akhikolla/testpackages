\name{OptimPriceStorage}
\alias{ComputeMarketPrices}
\alias{OptimPriceMarket_l}
\alias{OptimPriceMarket_q}
\alias{OptimPriceStorage}
\alias{OptimPriceStorage_}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimisation of storage operation with market prices taking into acount storage efficiency and network taxes. 
}
\description{
Optimisation of storage operation with market prices taking into acount storage efficiency and network taxes. 
}
\usage{
OptimPriceStorage(Prices,Pplus,Pmoins,Cplus,Cmoins=0,
                            efficiencyS=0,efficiencyP=efficiencyS,networkTax=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Prices}{A vector of prices}
  \item{Pplus}{A value for the upper power constraint or a vector of values with the same size as Prices}
\item{Pmoins}{A value for the lower power constraint or a vector of values with the same size as Prices}
\item{Cplus}{A value for the upper capacity constraint or a vector of values with the same size as Prices}
\item{Cmoins}{A value for the lower capacity constraint or a vector of values with the same size as Prices}

\item{efficiencyS}{storage efficiency when storing electricity}
\item{efficiencyP}{storage efficiency when producing electricity}
\item{networkTax}{networkTax}

}
\details{
 function OptimPriceStorage solves 
#             min_x sum_i=1^n Y_i*efficiencyP x_i*(x_i<0)+(Y_i*efficiencyS +networkTax)*x_i*(x_i>0)
#                   Pmoins_i<=  x_i  			<=Pplus_i 		i=1,...,n
#					Cmoins_i<=	sum_j=1^i x_j	<=Cplus_i 	i=1,...,n

 when efficiency=1 and networkTax=0 this gives
#             min_x sum_i=1^n Y_i x_i
#                   Pmoins_i<=  x_i				<=Pplus_i 		i=1,...,n
#					Cmoins_i<=	sum_j=1^i x_j	<=Cplus_i 	i=1,...,n

}
\value{
A list with

\item{Operation }{the optimal operation for each time step}
\item{Revenue }{the revenue for each time step}
%% ...
}
\references{
TODO%% ~put references to the literature/web site here ~
}
\author{
Robin Girard
}
\note{
TODO %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 to See Also \code{\link{cplfunction}} (method OptimMargInt that is more general)
}
\examples{



n=8760
Prices=runif(n,1,100) ##uniform random prices in [1;100] in Euro/MWh
Pmax=1; Pmin=-1; Cmax=5; ## 1MW maximum  during 5 hours.
res=OptimPriceStorage(Prices,Pmax,Pmin,Cmax) # solving the optimization problem
sum(res$Revenue)## Revenue
res=OptimPriceStorage(Prices,Pmax,Pmin,Cmax,efficiencyS=0.8) # solving the optimization problem
sum(res$Revenue)## Revenue

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
