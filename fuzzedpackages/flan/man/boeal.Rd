\docType{data}
\name{boeal}
\alias{boeal}

\title{Number of nalidixic acid-resistant mutants from Boe et al. (1994)}

\description{
  Samples from Table 4 p. 2784-2785 of Boe and al (1994).
  
  Plates with more than 512 mutants could not be counted with precision, hence the value 512 must be understood as ``512 or more''.
}

\usage{
  data(boeal)
}

\format{
  A list of 23 samples of mutants counts, each named "B<index of the sample>".
  
  The i-th sample of the list includes the i-th column of the table.
}

\examples{
b <- unlist(boeal)              # concatenate all samples
ml <- mutestim(b)               # maximum likelihood
gf <- mutestim(b,method="GF")   # generating function
p0 <- mutestim(b,method="P0")   # P0 method
cbind(ml,gf,p0)                 # compare 3 methods

                                # test values of mutations and fitness
flan.test(b,alternative=c("greater","less"),mutations0=0.6,fitness0=1)

b1 <- unlist(boeal[1:10])       # first 10 samples
b2 <- unlist(boeal[11:20])      # next 10 samples
flan.test(list(b1,b2))          # test equality
}

\source{L. Boe, T. Tolker-Nielsen, K. M. Eegholm, H. Spliid, and A. Vrang:  Fluctuation analysis of mutations to nalidixic acid resistance in \emph{Escherichia Coli}, \emph{J. Bacteriol.}, 176(10):2781-2787 (1994)
}
\keyword{datasets}
