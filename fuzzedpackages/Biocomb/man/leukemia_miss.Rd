\encoding{UTF-8}
\name{leukemia_miss}
\docType{data}
\alias{leukemia_miss}
\title{
 desease data
}
\description{
  Leukemia dataset includes the bone marrow samples obtained from acute leukemia
patients at the time of diagnosis: 25 acute myeloid leukemia (AML) samples; 9
T-lineage acute lymphoblastic leukemia (ALL) samples; and 38 B-lineage ALL
samples. After preprocessing, the 100 genes with the largest variation across samples are selected. This dataset is the same as leukemia72 with the 5 percent of missing values.
}
\usage{leukemia_miss}

\format{A data.frame containing 72 observations of 101 variables: 100 features and class with three diagnosis: 38 B-lineage ALL, 9 T-lineage ALL and 25 AML. It has 5 percent missing values.}

\source{
Handl J, Knowles J, Kell DB, Computational cluster validation in post-genomic data analysis, Bioinformatics 21:3201-3212, 2005.
}

\seealso{
  \code{\link{data_test}}, \code{\link{datasetF6}}, \code{\link{leukemia72}}, \code{\link{leukemia72_2}}
}

\examples{
# load the dataset
data(leukemia_miss)
# X95735_at
with(leukemia_miss, by(X95735_at,Class,mean,na.rm=TRUE))

# M27891_at
with(leukemia_miss,tapply(M27891_at, Class, FUN = mean,na.rm=TRUE))
}

\keyword{datasets}
