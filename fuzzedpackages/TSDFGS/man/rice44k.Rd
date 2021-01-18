\name{rice44kPCA}
\alias{rice44kPCA}
\docType{data}
\title{
44k genome rice data
}
\description{
This data set was provided by Zhao et al. (2011) <DOI:10.1038/ncomms1467> which genotyping 44,100 SNP variants across 413 diverse accessions of O. sativa from 82 countries. We converted the genomic information into a PC matrix.
}
\usage{data("rice44kPCA")}
\format{
  \describe{
    \item{\code{geno}}{A numeric matrix of principal components}
    \item{\code{subpop}}{A character vector of subpopulation's group name.}
  }
}
%\details{}
\source{
http://www.ricediversity.org/data/index.cfm
}
\references{
Keyan Zhao, Chih-Wei Tung, Georgia C. Eizenga, Mark H. Wright, M. Liakat Ali, Adam H. Price, Gareth J. Norton, M. Rafiqul Islam, Andy Reynolds, Jason Mezey, Anna M. McClung, Carlos D. Bustamante & Susan R. McCouch (2011). Genome-wide association mapping reveals a rich genetic architecture of complex traits in Oryza sativa. Nat Comm 2:467 | DOI: 10.1038/ncomms1467, Published Online 13 Sep 2011.
}
\examples{
data("rice44kPCA")
}
%\keyword{datasets}
