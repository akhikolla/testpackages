\name{GCPM-package}
\alias{GCPM-package}
\alias{GCPM}
\docType{package}
\title{
Generalized Credit Portfolio Model
}
\description{
The package helps to analyze the default risk of credit portfolios. Commonly known models, 
like CreditRisk+ or the CreditMetrics model are implemented in their very basic settings.
The portfolio loss distribution can be achieved either by simulation or analytically 
in case of the classic CreditRisk+ model. Models are only implemented to respect losses
caused by defaults, i.e. migration risk is not included. The package structure is kept
flexible especially with respect to distributional assumptions in order to quantify the
sensitivity of risk figures with respect to several assumptions. Therefore the package
can be used to determine the credit risk of a given portfolio as well as to quantify
model sensitivities.
}
\details{
\tabular{ll}{
Package: \tab GCPM\cr
Type: \tab Package\cr
Version: \tab 1.2.2\cr
Date: \tab 2016-12-29\cr
License: \tab GPL-2\cr
}
}
\author{
Kevin Jakob

Maintainer: Kevin Jakob <Kevin.Jakob.Research@gmail.com>
}
\references{
Jakob, K. & Fischer, M. "GCPM: A flexible package to explore credit portfolio risk" Austrian Journal of Statistics 45.1 (2016): 25:44\cr
Morgan, J. P. "CreditMetrics-technical document." JP Morgan, New York, 1997\cr 
First Boston Financial Products, "CreditRisk+", 1997 \cr
Gundlach & Lehrbass, "CreditRisk+ in the Banking Industry", Springer, 2003\cr
}
\keyword{ package }
\seealso{
\code{\link{GCPM-class}}, \code{\link{init}}, \code{\link{analyze}}
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
