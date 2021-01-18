\name{GRR.matrix}
\alias{GRR.matrix}
\title{
  GRR matrix for genetic data simulation
}
\description{
  Computes a GRR matrix based on a simulation model
}
\usage{
GRR.matrix(genes.maf = Kryukov, n.case.groups = 2, 
           GRR = c("SKAT", "constant", "variable"), 
           GRR.value, GRR.function, GRR.multiplicative.factor = 2, select.gene)
}
\arguments{
  \item{genes.maf}{ A dataframe containing at least the MAF in the general population (column \code{maf}) with their associated gene (column \code{gene}). By default, \code{maf} from the file Kryukov are used}
  \item{n.case.groups}{ The number of cases groups (set at 2 by default), i.e. the number of groups where variants will have a GRR greater than 1}
  \item{GRR}{ How to calculate the GRR }
  \item{GRR.value}{ GRR value if \code{GRR="constant"}}
  \item{GRR.function}{ A function indicating how to calculate the GRR depending on MAF in the general population, only needed if \code{GRR="variable"}}
  \item{GRR.multiplicative.factor}{ A vector of size (\code{n.case.groups}-1) containing the multiplicative factor for the GRR for each group of cases compared to the first group of cases}
  \item{select.gene}{ The gene(s) to be selected from the file \code{genes.maf} if multiple genes are present. If missing, the first level of \code{genes.maf$gene} is kept.}
}
\value{
  A matrix containing the GRR values with one column per variant and one line per cases group
}
\details{
  The GRR can be computed in three ways using the argument \code{GRR}.

  If \code{GRR="constant"}, the same GRR is given to all the variants, its value being specified to \code{GRR.value}.  
  If \code{GRR="SKAT"}, the GRR are calculating using the formula from the paper presenting the SKAT method and thus depend on MAF.  
  If \code{GRR="variable"}, the GRR are calculating using a function given by the user to \code{GRR.function} depending only on the MAF in the general population.  

  The argument \code{multiplicative.factor} contains \code{n.case.groups}-1 values; if \code{multiplicative.factor=1}, GRR will be the same between the different groups of cases.  

  The two dataframes \code{Kryukov} (used by default) and \code{GnomADgenes} (containing MAF in the NFE population) can be used as \code{genes.maf}.

  \code{GRR.matrix} returns a matrix that can be used in other simulation functions such as \code{random.bed.matrix}.
}
\seealso{ \code{\link{random.bed.matrix}}, \code{\link{GnomADgenes}}, \code{\link{Kryukov}} }
\examples{
#GRR calculated on the MAF from the first unit of the file Kryukov
#using the formula from the SKAT paper, with the second group of cases 
#having GRR values twice as high as the first one
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov,
                      GRR.multiplicative.factor=2, select.gene = "R1")
}
