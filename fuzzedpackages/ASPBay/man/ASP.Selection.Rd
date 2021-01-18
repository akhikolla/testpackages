\name{ASP.Selection}
\alias{ASP.Selection}
\title{ Select a subset of SNPs }
\description{ Select a subset of SNPs using discrimination method for affected sib pairs }
\usage{ ASP.Selection(Tem_Gen, Index_Gen, IBD, k = log(10000)) }
\arguments{
  \item{Tem_Gen}{ Genotypes of controls (denoted by the number of alternative allele) }
  \item{Index_Gen}{ Genotypes of index cases }
  \item{IBD}{ IBD states for each affected sib pair }
  \item{k}{ Selection threshold (by default log(1e4)) }
}
\details{
Take the genotypes of controls and index cases and the IBD states. Give the score statistics, discrimination statistics and the subset of selected SNPs with the chosen threshold.
}
\value{ List of 4 vectors with components:
  \item{score}{ The values of the score statistic for each SNPs }
  \item{stat}{ The values of discrimination statistic comparing each SNPs with the most associated SNP }
  \item{SNP_subset}{ The indexes (numbers of columns) of selected SNPs }
  \item{SNPname_subset}{ The names (names of columns) of selected SNPs }
}
\references{ Dandine-Roulland, Claire and Perdry, Herve. \emph{Where is the causal variant? On the advantage of the family design over the case-control design in genetic association studies.} 
Submitted to Eur J Hum Genet }
\author{Claire Dandine-Roulland}
\seealso{\code{\link{ASP.Score}}, \code{\link{ASP.Bayesian}}}
\examples{
data(ASPData)
ASP.Selection(ASPData$Control, ASPData$Index, ASPData$IBD )
}
\keyword{ Selection }
\keyword{ Affected sib pairs }
