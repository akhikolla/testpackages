\name{Graphs.Bayesian}
\alias{Graphs.Bayesian}
\title{ Graphs }
\description{ Plot graphs to visualize the results of \code{ASP.Bayesian} }
\usage{ Graphs.Bayesian(M, burn=0, xbins=200, ORlim=c(1,5),
         conf.int=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95), print=TRUE) }
\arguments{
  \item{M}{ Object given by the function ASP.Bayesian }
  \item{burn}{ The first \code{burn} values of the sampling are removed }
  \item{xbins}{ The number of bins which partition the range of graph variables }
  \item{ORlim}{ OR limits in graphs }
  \item{conf.int}{ Chosen credibility intervals }
  \item{print}{ Logical, if TRUE the plots are printed }
}
\details{ Plot two graphs and give associated hexbinplot objects. This two graphs summarize the results of the Bayesian method.
          The first graph shows the linkage disequilibrium between observed and causal SNPs in abscissae and the OR of causal SNP in ordinates.
          The second graph displays the alternative allele frequency of causal SNP in abscissae and the alternative allele frequency of observed SNP in ordinates.
		  Before plotting the graphs, the causal odds ratio is transformed. The value of OR is kept if it is superior to 1, otherwise it is inversed.
		  The alternative causal allele frequency is transformes accordingly: if the OR is inferior to 1, the frequency is replaced by its complement to 1.
		  With this transformations, we avoid to obtain two peaks corresponding to equivalent parameter values.}
\value{ List of 2 hexbin objects:
 \item{hex_r2_OR}{ Hexbinplot object with the linkage disequilibrium between observed and causal SNPs in abscissae and the OR of causal SNP in ordinates. }
 \item{hex_fa_fb}{ Hexbinplot object with the alternative allele frequency of causal SNP in abscissae and the alternative allele frequency of observed SNP in ordinates. }
}
\references{ Dandine-Roulland, Claire and Perdry, Herve. \emph{Where is the causal variant? On the advantage of the family design over the case-control design in genetic association studies.} 
Submitted to Eur J Hum Genet }
\author{Claire Dandine-Roulland}

\seealso{ \code{\link{ASP.Bayesian}} }
\examples{
data(ASPData)
B <- ASP.Bayesian(1e5, ASPData$Control, ASPData$Index,
                  ASPData$IBD, 15)
G <- Graphs.Bayesian(B, burn = 5000, xbins=100)
}
\keyword{ Graphs }
