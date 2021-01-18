\name{TCPAprad}
\alias{TCPAprad}
\docType{data}

\title{
  Protein expression data.
}

\description{Level 3 normalized proteomic data (v3.0) from The Cancer Proteome Alas (http://tcpaportal.org/tcpa). The data comprise the measurements of 189 antibodies obtained from 164 tumor tissue samples (prostate adenocarcinoma) using reverse phase protein arrays (RPPA).}

\usage{
TCPAprad
}

\format{A 164 by 189 matrix}

\source{
The Cancer Proteome Alas (http://tcpaportal.org/tcpa)
}

\references{
Li J, Lu Y, Akbani R, Ju Z, Roebuck PL, Liu W, Yang J-Y, Broom BM, Verhaak RGW, Kane DW, Wakefield C, Weinstein JN, Mills GB, Liang H. (2013). TCPA: A Resource for Cancer Functional Proteomics Data. \emph{Nature Methods} 10(11), 1046-1047.
}

\examples{
data(TCPAprad)
dim(TCPAprad)
TCPAprad[1:5, 1:5]
}
