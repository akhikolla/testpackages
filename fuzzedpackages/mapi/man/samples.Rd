\name{samples}
\alias{samples}
\title{'samples' test dataset}
\usage{
data("samples")
}
\value{
A data.table object with 200 simulated individuals (one per row) and 4 columns ("ind", "x", "y", "errRad") 
including the sample names, coordinates x and y and an error circle radius on coordinates.
}
\description{
Test dataset provided with the MAPI package. We used generation-by-generation coalescent algorithms 
(Hudson et al., 1990) to simulate 10 microsatellite genotypes at migration-mutation-drift equilibrium 
for 200 diploid individuals, distributed on the nodes of a 20x10 lattice.

Mutations for each locus followed a symmetric generalized stepwise model with a variance equal to 0.36
(Estoup et al., 2001) and a maximum range of allelic states of 40. The mutation rate was fixed so that
heterozygosity ranged from 0.6 to 0.8 as frequently observed at microsatellite
markers (Chapuis et al., 2012).

Two panmictic populations are separated by a barrier. We used
Simcoal2 (Laval and Excoffier, 2004) to generate two panmictic populations of equal effective size
\eqn{N_e = 100} and exchanging \eqn{N_{e}m = 0.1} migrants at each generation. The barrier to gene
flow bisected the lattice from north to south in its center.

The differentiation values computed between these samples is described in the \code{\link{metric}} section.
}
\examples{
data("samples")
}
\references{
Chapuis, M.-P., Streiff, R., and Sword, G. (2012). Long microsatellites and unusually high levels 
of genetic diversity in the orthoptera. Insect Molecular Biology, 21(2):181--186.

Estoup, A., Wilson, I. J., Sullivan, C., Cornuet, J.-M., and Moritz, C. (2001). Inferring 
population history from microsatellite and enzyme data in serially introduced cane toads, 
Bufo marinus. Genetics, 159(4):1671--1687.

Hudson, R. R. et al. (1990). Gene genealogies and the coalescent process. 
Oxford Surveys in Evolutionary Biology, 7(1):44.

Laval, G. and Excoffier, L. (2004). SIMCOAL 2.0: a program to simulate genomic diversity 
over large recombining regions in a subdivided population with a complex history. 
Bioinformatics, 20(15):2485--2487.
}
