\name{genepop-package}
\alias{genepop-package}
\alias{genepop}
\docType{package}
\title{
  Population genetic analyses using the Genepop software
}
\description{
  A distribution of the Genepop software as an R package. The included C++ sources are suitable for compilation as a stand-alone executable. A shiny interface is also included.    Genepop performs three main tasks: it computes \bold{exact tests} for Hardy-Weinberg equilibrium (\code{\link{test_HW}}), for population differentiation  (\code{\link{test_diff}}) and for genotypic disequilibrium among pairs of loci  (\code{\link{test_LD}}); it computes \bold{estimates} of F-statistics  (\code{\link{Fst}}), null allele frequencies (\code{\link{nulls}}), allele size-based statistics for microsatellites, etc., and of number of immigrants by Barton & Slatkin's (1986) private allele method  (\code{\link{Nm_private}}); It performs analyses of \bold{isolation by distance} from pairwise comparisons of individuals or groups (\code{\link{ibd}}), including confidence intervals for "neighborhood size". It also provides various data \code{\link{conversion}} and \code{\link{manipulation}} utilities.  
}
\author{
R package originally developed by Jimmy Lopez and Khalid Belkhir from the C++ sources of the Genepop executable version 4.6 (2016; Rousset 2008). 
}
\references{
\bold{Main reference for current maintained version of Genepop}:

Rousset, F. (2008). Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. Mol. Ecol. Res. 8: 103-106.

Original Genepop publication: 

Raymond, M. & Rousset, F., 1995b. GENEPOP Version 1.2: population genetics software for exact tests and ecumenicism. J. Hered. 86: 248-249.

Methods implemented in Genepop:

Barton, N. H. & Slatkin, M., 1986. A quasi-equilibrium theory of the distribution of rare alleles in a subdivided population. Heredity 56: 409-415.

Brookfield, J. F. Y., 1996. A simple new method for estimating null allele frequency from heterozygote deficiency. Mol. Ecol. 5: 453-455.

Goudet, J., Raymond, M., de Meeus, T. & Rousset, F., 1996. Testing differentiation in diploid populations. Genetics 144: 1931-1938.

Guo, S. W. & Thompson, E. A., 1992. Performing the exact test of Hardy-Weinberg proportion for multiple alleles. Biometrics 48: 361-372.

Kalinowski, S. T. & Taper, M. L., 2006. Maximum likelihood estimation of the frequency of null alleles at microsatellite loci. Conserv. Genetics 7:991-995.

Louis, E. J. & Dempster, E. R., 1987. An exact test for Hardy-Weinberg and multiple alleles. Biometrics 43: 805-811.

Mantel, N., 1967. The detection of disease clustering and a generalized regression approach. Cancer Research 27: 209-220.

Michalakis, Y. & Excoffier, L., 1996. A generic estimation of population subdivision using distances between alleles with special interest to microsatellite loci. Genetics 142: 1061-1064.

Raymond, M. & Rousset, F., 1995a. An exact test for population differentiation.Evolution 49: 1283-1286.

Robertson, A. & Hill, W. G., 1984. Deviations from Hardy-Weinberg proportions: sampling variances and use in estimation of inbreeding coefficients. Genetics 107: 703-718.

Rousset, F., 1996. Equilibrium values of measures of population subdivision for stepwise mutation processes. Genetics 142: 1357-1362.

Rousset, F., 1997. Genetic differentiation and estimation of gene flow from F-statistics under isolation by distance. Genetics 145: 1219-1228.

Rousset, F., 2000. Genetic differentiation between individuals. J. Evol. Biol. 13:58-62.

Rousset, F. & Raymond, M., 1995. Testing heterozygote excess and deficiency. Genetics 140: 1413-1419.

Watts, P. C., Rousset, F., Saccheri, I. J., Leblois, R., Kemp, S. J. & Thompson, D. J., 2007. Compatible genetic and ecological estimates of dispersal rates in insect (Coenagrion mercuriale: Odonata: Zygoptera) populations: analysis of 'neighbourhood size' using a more precise estimator. Mol. Ecol. 16: 737-751.

Weir, B. S., 1996. Genetic Data Analysis II. Sinauer, Sunderland, Mass.

Weir, B. S. & Cockerham, C. C., 1984. Estimating F-statistics for the analysis of population structure. Evolution 38: 1358-1370.
}
\keyword{ package }
