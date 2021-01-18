\name{survSNP.power.table}
\alias{survSNP.power.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Table for Asymptotic and Empirical Power
}
\description{
  This function produces a table with the resulting empirical and
  asymptotic power over a given range of Genotype Hazard Ratios,
  relative allelic frequencies, sample sizes, and event rates.
}
\usage{
survSNP.power.table(GHRs, ns, rafs, erates, pilm, lm, model, test, alpha, 
                    exactvar = FALSE, B = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GHRs}{
    A vector of Genotype Hazard Ratios
}
\item{ns}{
  A vector of sample sizes
}
\item{rafs}{
  A vector of relative allelic frequencies (for the risk allele) 
}
\item{erates}{
  A vector of event rates
}
\item{pilm}{
  Probability that the time-to-event is greater than lm
}
\item{lm}{
  Landmark time used for powering the study
}
\item{model}{
    True genetic risk model (choices are "additive", "recessive" or
    "dominant". For the asymptotic calculations only the "additive"
    should be used until further notice)
}
  \item{test}{
    Hypothesized genetic risk model (choices are "additive", "recessive"
    or "dominant". For the asymptotic calculations only the "additive" should
    be used until further notice)

}
\item{alpha}{
  Nominal two-sided type I error rate
}
  \item{exactvar}{
    Indicator for using the exact variance formula
}
\item{B}{
   Number of simulation replicates (set to 0 if no empirical
    calculations are desired)

}
}
\details{
  This version only supports additive models.
}
\value{
  See output of \link{sim.snp.expsurv.power}.
}
\references{

  Kouros Owzar, Zhiguo Li, Nancy Cox and Sin-Ho Jung.
  Power and Sample Size Calculations for SNP Association
  Studies with Censored Time-to-Event Outcomes.
 \url{http://onlinelibrary.wiley.com/doi/10.1002/gepi.21645/abstract}

}
\author{
  Kouros Owzar, Zhiguo Li, Nancy Cox, Sin-Ho Jung and Chanhee Yi
}

\examples{
GRRs<-seq(1.5,2,by=0.25)
ns<-c(100,500,1000)
rafs<-c(0.3,0.5,0.7)
erates=c(0.5,0.7,0.9)
res<-survSNP.power.table(GRRs,ns,rafs,erates,pilm=0.5,lm=1,model="additive",
                         test="additive",alpha=0.05)

# Create key for illustration
KEY=paste("q=",levels(factor(res$raf)),sep="")
KEY<-list(lines=list(col=1:length(KEY),lty=1:length(KEY)),
          text=list(labels=paste("q=",levels(factor(res$raf)),sep="")),
          column=3)


# Illustrate Power
print(xyplot(pow0~GHR|factor(erate)*factor(n),group=factor(raf),
             data=res,type="l",lty=KEY$lines$lty,col=KEY$lines$col,
             key=KEY,
             xlab="Genotype Hazard Ratio",ylab="Power"))

# Illustrate Power (restricted to n=100)

print(xyplot(pow0~GHR|factor(erate),group=factor(raf),
             data=subset(res,n==ns[1]),
             type="l",lty=KEY$lines$lty,col=KEY$lines$col,
             key=KEY,
             xlab="Genotype Hazard Ratio",ylab="Power",
             sub=paste("n=",ns[1],", alpha=",round(unique(res$alpha),2))))


}

