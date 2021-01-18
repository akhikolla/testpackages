\name{PINSPlus-package}
\alias{PINSPlus-package}
\alias{PINSPlus}
\docType{package}
\title{
Perturbation Clustering for data INtegration and disease Subtyping
}

\description{
This package implements clustering algorithms proposed by Nguyen et al. (2017, 2019). Perturbation Clustering for data INtegration and disease Subtyping (PINS) is an approach for integraton of data and classification of diseases into various subtypes. PINS+ provides algorithms supporting both single data type clustering and multi-omics data type. PINSPlus is an improved version of PINS by allowing users to customize the based clustering algorithm and perturbation methods. Furthermore, PINSPlus is fast and supports the analysis or large datasets with millions of samples and features.
}

\details{
PINS+ provides \code{\link{PerturbationClustering}} and \code{\link{SubtypingOmicsData}} functions for single data type clustering and multi-omics data type clustering. PINS makes use of different clustering algorithms such as \code{kmeans} and \code{pam} to perform clustering actions. The principle of PINS is to find the optimum number of clusters and location of each sample in the clusters based on perturbation methods such as \code{noise} or \code{subsampling}. PINS+ allows users to pass their own clustering algorithm and perturbation method.
}

\references{

H Nguyen, S Shrestha, S Draghici, & T Nguyen. PINSPlus: a tool for tumor subtype discovery in integrated genomic data. Bioinformatics, 35(16), 2843-2846, (2019).

T Nguyen, R Tagett, D Diaz, S Draghici. A novel method for data integration and disease subtyping. 
Genome Research, 27(12):2025-2039, 2017.

Nguyen, H., Shrestha, S., Draghici, S., & Nguyen, T. (2019). PINSPlus: a tool for tumor subtype discovery in integrated genomic data. Bioinformatics, 35(16), 2843-2846.
}

\keyword{ package }
\seealso{
\code{\link{PerturbationClustering}}, \code{\link{SubtypingOmicsData}}
}
