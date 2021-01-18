\name{TSDFGS-package}
\alias{TSDFGS-package}
\alias{TSDFGS}
\docType{package}
\title{
\packageTitle{TSDFGS}
}
\description{
\packageDescription{TSDFGS}
}
\details{
The package is used to determine the optimal training set in a highly structured, mild structured and diverse population. The function "optTrain" use a genetic algorithm or simple exchange algorithm to evaluate an optimal solution using one of the criteria (r-score (Ou J.H., Liao C.T. (2018) <DOI:10.6342/NTU201802290>), PEV-score (Akdemir D. et al. (2015) <DOI:10.1186/s12711-015-0116-6>), CD-score(Laloe D. (1993) <DOI:10.1186/1297-9686-25-6-557>)).
%The DESCRIPTION file:
%\packageDESCRIPTION{TSDFGS}
%\packageIndices{TSDFGS}
%optTrain(geno, cand, n.train, subpop=NULL, test=NULL, method="rScore", min.iter=NULL)
}
\author{
\packageAuthor{TSDFGS}

Maintainer: \packageMaintainer{TSDFGS}
}
\references{
Akdemir D., Sanchez JI., Jannink JL. (2015), Optimization of genomic selection training populations with a genetic algorithm. Genetic Selection Evolution 47:38.\\

Laloe D. (1993), Precision and information in linear models of genetic evolution. Genetics Selection Evolution 25:557.\\

Ou J.H., Liao C.T. (2018), Training set determination for genomic selection. National Taiwan University Master Thesis.
}
\keyword{ package }
\seealso{
\code{\link[STPGA:STPGA-package]{STPGA}}
}
%\examples{}
