\name{init}
\alias{init}
\title{
Initialize an Object of Class \code{GCPM} 
}
\description{
The function helps to create a new object of class \code{GCPM}. The arguments 
of the function are passed to the object after performing some plausibility checks.
}
\usage{
init(model.type = "CRP", link.function = "CRP", N, seed,
loss.unit, alpha.max = 0.9999, loss.thr = Inf, sec.var,
random.numbers = matrix(), LHR, max.entries=1e3)
}
\arguments{
  \item{model.type}{Character value, specifying the model type. One can choose
  between \dQuote{simulative} and \dQuote{CRP} which corresponds to the
  analytical version of the CreditRisk+ model (see First Boston Financial
  Products, 1997)}
  \item{link.function}{
character value, specifying the type of the link function. One can choose
between \dQuote{CRP}, which corresponds to \eqn{\overline{PD}=PD\cdot (w^Tx)}
and \dQuote{CM} which corresponds to
\eqn{\overline{PD}=\Phi\left(\frac{\Phi^{-1}PD-w^Tx}{\sqrt{1-w^T\Sigma w}}\right)},
where PD is the original PD from portfolio data, x is the vector of sector
drawings, \eqn{\Phi} is the CDF of the standard normal distribution, w is the
vector of sector weights given in the portfolio data and \eqn{\Sigma} is the
correlation matrix of the sector variables estimated from \code{random.numbers}.
\dQuote{CRP} will be used automatically if \code{model.type} == "CRP".
}
  \item{N}{
numeric value, defining the number of simulations if
\code{model.type} == "simulative". If \code{N} is greater than the number of
scenarios provided via \code{random.numbers}, scenarios are reused.
This parameter is used only if \code{model.type} == "simulative".
}
  \item{seed}{
numeric value used to initialize the random number generator. If \code{seed} is 
not provided a value based on current system time will be used. This parameter 
is used only if \code{model.type} == "simulative".
}
  \item{loss.unit}{
numeric positive value used to discretize potential losses.
}
  \item{alpha.max}{
numeric value between 0 and 1 defining the maximum CDF-level which will be
computed in case of an analytical CreditRisk+ type model.
}
  \item{loss.thr}{
numeric value specifying a lower bound for portfolio losses to be stored in
order to derive risk contributions on counterparty level. Using a lower value
needs a lot of memory but will be necessary in order to calculate risk
contributions on lower CDF levels. This parameter is used only if
\code{model.type} == "simulative".}
  \item{sec.var}{
named numeric vector defining the sector variances in case of a CreditRisk+
type model. The names have to correspond to the sector names given in the
portfolio. This parameter is used only if \code{model.type} == "CRP".
}
  \item{random.numbers}{
matrix with sector drawings. The columns represent the sectors,
whereas the rows represent the scenarios (number of different simulations).
The column names must correspond to the names used in the portfolio data
(see \code{\link{analyze}}) and to the names of \code{sec.var} if
\code{model.type} == "CRP". This parameter is used only if
\code{model.type} == "simulative".
}
  \item{LHR}{
numeric vector of length equal to \code{nrow(random.numbers)} defining the
likelihood ratio of each scenario. If not provided, all scenarios are assumed
to be equally likely.  This parameter is used only if
\code{model.type} == "simulative".
}
\item{max.entries}{
numeric value defining the maximum number of loss scenarios stored to
calculate risk contributions.
}
}
\value{
object of class \code{GCPM}
}
\references{
Jakob, K. & Fischer, M. "GCPM: A flexible package to explore credit portfolio risk" Austrian Journal of Statistics 45.1 (2016): 25:44\cr
Morgan, J. P. "CreditMetrics-technical document." JP Morgan, New York, 1997\cr 
First Boston Financial Products, "CreditRisk+", 1997 \cr
Gundlach & Lehrbass, "CreditRisk+ in the Banking Industry", Springer, 2003\cr
}
\author{
Kevin Jakob
}
\seealso{
\code{\link{GCPM}}, \code{\link{GCPM-class}}, \code{\link{analyze}}
}
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
}

