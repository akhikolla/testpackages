\name{ ASP.Score }
\alias{ ASP.Score }
\title{ Score test of association }
\description{ Calculate score statistics and the associated P-value for each SNPs }
\usage{ ASP.Score(Tem_Gen, Index_Gen, IBD) }
\arguments{
  \item{Tem_Gen}{ Genotypes of controls (denoted by the number of alternative allele) }
  \item{Index_Gen}{ Genotypes of Index cases }
  \item{IBD}{ IBD states for each affected sib pair }
}
\details{ Give the values of statistic and p-value of the association score test. }
\value{List of 2 vectors of length the number of SNPs:
  \item{Value}{ Statistic values for each SNPs }
  \item{Pvalue}{ P-values of the score test for each SNPs }
}
\references{ Perdry, Herve, Muller-Myhsok, Bertram, et Clerget-Darpoux, Francoise. \emph{Using affected sib-pairs to uncover rare disease variants.} Hum Hered, 2013. }
\author{Claire Dandine-Roulland}
\seealso{ \code{\link{ASP.Selection}} }
\examples{
data(ASPData)
ASP.Score(ASPData$Control, ASPData$Index, ASPData$IBD )
}
\keyword{ Affected sib pairs }
\keyword{ Score test }
