\name{Jaccard}
\alias{Jaccard}
\title{
Jaccard index
}
\description{
Calculates the Jaccard index for each pair of individuals using a bed.matrix
}
\usage{
Jaccard(x, maf.threshold = 0.01)
}
\arguments{
  \item{x}{A bed.matrix}
  \item{maf.threshold}{The MAF used for the definition of a rare variant, set at 0.01 by default}
}
\details{
The individuals carrying no rare variants will have a null Jaccard index with all the individuals including themselves.
}
\value{
A squared matrix giving the Jaccard index for each pair of individuals
}
\references{ Jaccard, P. (1908) \emph{Nouvelles researches sur la distribution florale}, Bulletin de la Société vaudoise des sciences naturelles, \bold{44, 223-270} } 
\examples{
#Simulation of genetic data with GRR values according to the SKAT formula
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, 
                      n.case.groups = 2, select.gene = "R1",
                      GRR.multiplicative.factor=2)

#Simulation of one group of 1,000 controls and two groups of 500 cases, 
#50% of causal variants, 5 genomic regions are simulated.
x <- random.bed.matrix(genes.maf=Kryukov, size = c(1000, 500, 500),
                       prev = c(0.001, 0.001), select.gene = "R1", 
                       GRR.matrix.del = GRR.del, p.causal = 0.5,
                       genetic.model = "multiplicative", replicates = 5)
                           
#Calculate the Jaccard matrix
J <- Jaccard(x, maf.threshold = 0.01)
}
