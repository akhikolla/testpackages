\name{AAToCodon}
\alias{AAToCodon}
\title{Amino Acid to codon set}
\usage{
 AAToCodon(aa, focal = FALSE)
}
\arguments{
\item{aa}{Amino acid in single character notation}

\item{focal}{logical, Include the alphabetically last (focal) codon}
}
\value{
Returns the names of the codon encoding the give amino acid
}
\description{
Converts one character amino acid code to the set of codon encoding that amino acid
}
\seealso{
\code{\link{codonToAA}}
}

