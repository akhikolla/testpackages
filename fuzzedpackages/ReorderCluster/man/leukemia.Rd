\encoding{UTF-8}
\name{leukemia}
\alias{leukemia}
\docType{data}
\title{
Real biological dataset to perform the analysis.
}
\description{
The Leukemia dataset includes the bone marrow samples obtained from acute leukemia patients at the time of diagnosis: 25 acute myeloid leukemia (AML) samples and 47 acute lymphoblastic leukemia (ALL) (24 acute 9 T-lineage acute lymphoblastic leukemia samples; and 38 B-lineage ALL samples), The most informative 100 genes were selected to proceed with reordering algorithm. The dataset can be downloaded from the website http://www.broadinstitute.org/cgi-bin/cancer/datasets.cgi.The last column stores the class labels for dataset objects.
}
\usage{data(leukemia)}
\format{
  A data frame with 72 observations on the 102 variables.
}

\source{
Molecular classification of Cancer:class discovery and class prediction by gene expression monitoring / Golub T.R., Slonim D.K., Tamayo P., Huard C., Gaasenbeek M., Mesirov J.P., Coller H., Loh M.L., Downing J.R., Caligiuri M.A., et al. // Nature. – Vol. 286, 1999, P.531-537.
}

\examples{
data(leukemia)
}

\keyword{datasets}