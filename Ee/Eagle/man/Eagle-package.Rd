\name{Eagle-package}
\alias{Eagle-package}
\alias{Eagle}
\docType{package}
\title{
  Eagle for Genome-wide Association Mapping 
}
\description{
An implementation of multiple-locus association mapping on a genome-wide scale. 'Eagle' can handle inbred and  outbred study populations, populations of arbitrary unknown complexity, and data larger than the memory capacity of the computer. Since 'Eagle' is based on linear mixed models, it is best suited to the analysis of data on continuous traits. However, it can tolerate non-normal data. 'Eagle' reports, as its findings, the best set of snp in strongest association with a trait. For users unfamiliar with R, to perform an analysis, run 'OpenGUI()'. This  opens a web browser to the menu-driven user interface for the input of data, and for performing genome-wide analysis.
}
\details{


\emph{Motivation:}  Data from genome-wide association studies are analyzed, commonly, with single-locus 
 models. That is,  analyzes are performed on a locus-by-locus basis. Multiple-locus approaches 
 that model the association between a trait and multiple loci simultaneously are more powerful. However, 
 these methods do not scale well with study size and many of the packages that implement these methods are not easy 
 to use. 
 Eagle was specifically designed to make 
 genome-wide association mapping with multiple-locus models simple and practical.  

\emph{Assumptions}
 \enumerate{
 \item Individuals are diploid but they can be inbred or outbred.
 \item The marker and phenotype data are in separate files.
 \item  Marker loci are snps. Dominant and multi-allelic loci will need to be converted into biallelic (snp-like) loci. 
 \item The trait is continuous and normally distributed. Eagle can handle non-normally distributed trait data but 
  there may be a loss of power
 to detect marker-trait associations. 
 }

\emph{Important Functions:}
\enumerate{
 \item \code{\link{ReadMarker}} for reading in the snp data.
 \item \code{\link{ReadPheno}} for reading in the phenotypic data (traits and features/covariates)
 \item \code{\link{ReadMap}} for reading in the marker map.
 \item \code{\link{FPR4AM}} for calculating the value of the \code{lambda} parameter to be used by \code{\link{AM}} that will give a desired false positive rate for detecting SNP-trait associations. 
 \item \code{\link{AM}} for performing association mapping on the data.
 \item \code{\link{OpenGUI}} which opens the GUI. 
}

\emph{Output:} 
The key output from \code{\link{AM}} is a list of snp. Each snp identifies a separate genomic region of interest, housing 
genes that are affecting the trait. Additional summary information such as the size of the snp effects, their statistical significance, and 
how much phenotypic variation they explain can be obtained by running \code{\link{SummaryAM}}.


\emph{Where to get help: }
 A variety of different help options are available. 
 \itemize{
  \item At the R prompt, type  \preformatted{library(, "Eagle")} for an overview of the package and its functions. 
 \item For detailed help on a function called "foo" say, type  \preformatted{help("foo")} 
 \item Visit the Eagle website at \url{http://eagle.r-forge.r-project.org/} where 
 you can find a quick start guide, instructions on getting the most out of Eagle, 
 video tutorials, and other useful information. 
}
}
\author{
Andrew W. George (Data61, CSIRO)  with a lot of support from Joshua Bowden (IM&T, CSIRO)

Maintainer: Andrew W. George  <andrew.george@csiro.au>
}


