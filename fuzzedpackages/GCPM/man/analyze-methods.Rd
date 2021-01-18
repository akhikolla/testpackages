\name{analyze-methods}
\docType{methods}
\alias{analyze-methods}
\alias{analyze,GCPM-method}
\alias{analyze,GCPM,data.frame,missing,missing-method}
\alias{analyze,GCPM,data.frame,numeric,missing-method}
\alias{analyze,GCPM,data.frame,missing,numeric-method}
\alias{analyze,GCPM,data.frame,numeric,numeric-method}
\alias{analyze}

\title{Analyze a Credit Portfolio}
\description{
The method analyzes a given portfolio with a predefined portfolio model
(i.e. a \code{GCPM} object). Portfolio key numbers such as the number of portfolio
positions, sum of EAD and PL or the expected loss are calculated. Afterwards the
loss distribution is estimated according to \code{model.type}.
}
\usage{
analyze(this,portfolio,alpha,Ncores)
}
\arguments{
\item{this}{object of class \code{GCPM}}
\item{portfolio}{data frame containing portfolio data. The following columns
have to be defined (please be aware of the correct spelling of the column names):\cr
  \code{Number}: identification number for each portfolio position (numeric)\cr
  \code{Name}: counterparty name (character)\cr
  \code{Business}: business information (character/factor)\cr
  \code{Country}: country information (character/factor)\cr
  \code{EAD}: exposure at default (numeric)\cr
  \code{LGD}: loss given default (numeric in [0,1])\cr
  \code{PD}: probability of default (numeric in [0,1])\cr
  \code{Default}: default distribution either \dQuote{Bernoulli} or
  \dQuote{Poisson} (employable for pools)\cr
  sectors: starting with the 9th column, the sector weights have to be defined..
  }
  \item{alpha}{loss levels for risk measures economic capital, value at risk
  and expected shortfall (optional)}
  \item{Ncores}{number of (virtual) cores used to perfom Monte Carlo simulation
  (requires package \pkg{parallel}, default=1)}
}
\value{
object of class \code{GCPM}.
}
\details{
In case of an analytical CreditRisk+ model, a modified version of the
algorithm described in Gundlach & Lehrbass (2003) is used. For a simulative
model, the loss distribution is estimated based on \code{N} simulations with
sector drawings specified by \code{random.numbers} (see \code{\link{init}}).
The sector names (column names) should not include any white spaces. In case of a
CreditMetrics type model, the values of R (not R^2) have to be provided as sector
weights. In the standard CreditMetrics or CreditRisk+ framework a counterparty can
be assigned to more than one sector. Within a analytical CreditRisk+ model, the
sector names have to match the names of \code{sec.var} or in a simulative model
the column names of \code{random.numbers} (see \code{\link{init}})
}
\section{Methods}{
\describe{

\item{\code{signature(this = "GCPM", portfolio = "data.frame",
alpha = "missing")}}{
If loss levels \code{alpha} are not provided, risk measures such as economic
capital, value at risk and expected shortfall are not calculated by default.
However, they can be calculated afterwards by calling the corresponding methods
(see \code{\link{VaR}}, \code{\link{ES}}, \code{\link{EC}})
}

\item{\code{signature(this = "GCPM", portfolio = "data.frame",
alpha = "numeric")}}{
If loss levels \code{alpha} are provided, risk measures such as economic
capital, value at risk and expected shortfall are calculated and printed.
To extract these risk measures into a separate variable you can use the
corresponding methods.
}
}}
\references{
Jakob, K. & Fischer, M. "GCPM: A flexible package to explore credit portfolio risk" Austrian Journal of Statistics 45.1 (2016): 25:44\cr
Morgan, J. P. "CreditMetrics-technical document." JP Morgan, New York, 1997\cr 
First Boston Financial Products, "CreditRisk+", 1997 \cr
Gundlach & Lehrbass, "CreditRisk+ in the Banking Industry", Springer, 2003\cr
}
\seealso{\code{\link{init}}}
\keyword{methods}
\keyword{GCPM}
\examples{
#create a random portfolio with NC counterparties
NC=100
#assign business lines and countries randomly
business.lines=c("A","B","C")
CP.business=business.lines[ceiling(runif(NC,0,length(business.lines)))]
countries=c("A","B","C","D","E")
CP.country=countries[ceiling(runif(NC,0,length(countries)))]

#create matrix with sector weights (CreditRisk+ setting)
#according to business lines
NS=length(business.lines)
W=matrix(0,nrow = NC,ncol = length(business.lines),
dimnames = list(1:NC,business.lines))
for(i in 1:NC){W[i,CP.business[i]]=1}

#create portfolio data frame
portfolio=data.frame(Number=1:NC,Name=paste("Name ",1:NC),Business=CP.business,
                     Country=CP.country,EAD=runif(NC,1e3,1e6),LGD=runif(NC),
                     PD=runif(NC,0,0.3),Default=rep("Bernoulli",NC),W)

#draw sector variances randomly
sec.var=runif(NS,0.5,1.5)
names(sec.var)=business.lines

#draw N sector realizations (independent gamma distributed sectors)
N=5e4
random.numbers=matrix(NA,ncol=NS,nrow=N,dimnames=list(1:N,business.lines))
for(i in 1:NS){
random.numbers[,i]=rgamma(N,shape = 1/sec.var[i],scale=sec.var[i])}

#create a portfolio model and analyze the portfolio
TestModel=init(model.type = "simulative",link.function = "CRP",N = N,
loss.unit = 1e3, random.numbers = random.numbers,LHR=rep(1,N),loss.thr=5e6,
max.entries=2e4)
TestModel=analyze(TestModel,portfolio)

#plot of pdf of portfolio loss (in million) with indicators for EL, VaR and ES
alpha=c(0.995,0.999)
plot(TestModel,1e6,alpha=alpha)

#calculate portfolio VaR and ES
VaR=VaR(TestModel,alpha)
ES=ES(TestModel,alpha)

#Calculate risk contributions to VaR and ES
risk.cont=cbind(VaR.cont(TestModel,alpha = alpha),
ES.cont(TestModel,alpha = alpha))

#Use parallel computing for Monte Carlo simulation
TestModel=analyze(TestModel,portfolio,Ncores=2)
}
