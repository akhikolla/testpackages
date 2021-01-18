\name{plot.eqtl}
\docType{methods}
\alias{plot.eqtl}
\alias{plot,eqtl-method}

\title{Plot an eqtl Object}
\description{
 The function offers informative plots for an \code{eqtl} object.
}
\usage{
 \method{plot}{eqtl}(x, file = NULL, which = NULL, sig = 0.01, verbose = TRUE, 
                     centered = TRUE, log = FALSE, mc.cores = 1, genome = NULL, ...)
}

\arguments{
  \item{x}{Object of class \code{eqtl}.}
  \item{file}{Store set of graphics under that file name.}
  \item{which}{Specifies for which genes should the plot be created.}
  \item{sig}{Chosen significance level.}
  \item{verbose}{Logical, extended feedback of the function.}
  \item{centered}{Logical, plot should be centered around center gene.}
  \item{log}{Logical, y-axis scale is log(base=10)-scaled.}
  \item{mc.cores}{Amount of cores for parallel computing.}
  \item{genome}{The name of an existing genome.}
  \item{...}{Additional plotting parameters.}
}

\details{
 This function plots the test results of an \code{eqtl} object. Typically is the tested gene in the center and the p-values of associated SNPs are visualized. 
 Monomorphic SNPs and those that were missing are separately plotted. Test results that are smaller than the value given to \code{sig} are marked in red. 
 The y.axis can be switched to log10 scale by setting the logical parameter \code{log=TRUE} in that case are bars instead of dots plotted. If the y-axis is on log-scale
 it is also possible to give a second \code{eqtl} object to the function and plot the test results for both.
 
 The annotation feature is currently under development and only available in limited form. 
}

\author{Daniel Fischer}

\examples{
# Load the data
  data("annotTrack")   
  data("geneEXP")      
  data("genotData") 
\dontrun{
# A cis eQTL for 10 different genes:
  EQTL1 <- eQTL(gex=geneEXP[,1:10] , xAnnot = annotTrack,
                geno= genotData)

# A trans-eQTl for two different genes:
  EQTL2   <- eQTL(gex=geneEXP[,c(1,1000)] , xAnnot = annot.bed, 
                  geno= genotData, windowSize = NULL)

# Visualize the second cis-eQTL:
  plot(EQTL1, which=2)
  
# Visualize the trans-eQTL
  plot(EQTL2)
}
}

\keyword{methods}
\keyword{hplot}
