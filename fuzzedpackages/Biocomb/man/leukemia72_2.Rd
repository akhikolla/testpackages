\encoding{UTF-8}
\name{leukemia72_2}
\docType{data}
\alias{leukemia72_2}
\title{
 desease data
}
\description{
  Leukemia dataset includes the bone marrow samples obtained from acute leukemia
patients at the time of diagnosis: 25 acute myeloid leukemia (AML) samples; 9
T-lineage acute lymphoblastic leukemia (ALL) samples; and 38 B-lineage ALL
samples. After preprocessing, the 100 genes with the largest variation across samples are selected.
}
\usage{leukemia72_2}

\format{A data.frame containing 72 observations of 101 variables: 100 features and class with two diagnosis: 47 ALL and 25 AML).}

\source{
Handl J, Knowles J, Kell DB, Computational cluster validation in post-genomic data analysis, Bioinformatics 21:3201-3212, 2005.
}

\seealso{
  \code{\link{data_test}}, \code{\link{datasetF6}}, \code{\link{leukemia72}}, \code{\link{leukemia_miss}}
}

\examples{
# load the dataset
data(leukemia72_2)
# X95735_at
with(leukemia72_2, by(X95735_at,Class,mean))

# M27891_at
with(leukemia72_2,tapply(M27891_at, Class, FUN = mean))
with(leukemia72_2, table(M27891_at=ifelse(M27891_at<=mean(M27891_at), "1", "2"), Class))
}

\keyword{datasets}
