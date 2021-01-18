\name{prePed}
\Rdversion{1.1}
\alias{prePed}
\title{\bold{Pre}pares a \bold{Ped}igree 
}
\description{
Prepares a pedigree by sorting and adding founders and pruning the pedigree.
}
\usage{
prePed(Pedig, keep=NULL, thisBreed=NA, lastNative=NA, addNum=FALSE)}

\arguments{
\item{Pedig}{Data frame containing the pedigree where the first 3 columns correspond to: Individual ID,  Sire, and Dam. More columns can be passed in the \code{Pedig} argument including columns named \code{Sex}, \code{Breed} (with breed names), and \code{Born} (with years of birth or generation numbers). Missing parents are coded as NA, 0, or "0".}
\item{keep}{Vector with IDs of individuals, or \code{NULL}, or a logical vector indicating the individuals to be kept. If this parameter is not \code{NULL}, then only these individuals and their ancestors will be kept in the pedigree.}
\item{thisBreed}{Name of the breed.}
\item{lastNative}{Last year of birth for which individuals with unknown pedigree are considered native.}
\item{addNum}{If \code{TRUE}, then columns with IDs of individuals, sires, and dams in integer form will be added.
}
}

\details{

This function takes a pedigree, adds missing founders, and sorts the pedigree. If parameter \code{keep} contains IDs of individuals then only these individuals and their ancestors will be kept in the pedigree.

If the pedigree contains loops, then the loops will be broken by setting the parents of one animal in each loop to \code{NA}. 

If the pedigree contains columnn \code{Sex} then the sexes will be recoded as  \code{'male'} and \code{'female'}. Missing sexes will be determined from pedigree structure if possible. 

If the pedigree contains column \code{Breed} then for ancestors with missing breed the breed name is estimated. If parameter \code{lastNative} is not \code{NA} then for each animal with one missing parent an imaginary founder is added to the pedigree in order to enable classifying the breed names of all founders as follows: 
In general animals with missing breed are assumed to have the same breed as most of their offspring. But there is one exception: For founders belonging to \code{thisBreed} who are born after \code{lastNative} the breed name will be set to \code{"unknown"}. Moreover for founders from \code{thisBreed} with unknown year of birth the breed name will be set to \code{"unknown"} if all their descendants are born after \code{lastNative+I}. 
}

\value{
Data frame containing the pedigree with columns:

\item{Indiv}{Character column with IDs of the individuals}
\item{Sire}{Character column with IDs of the sires}
\item{Dam}{Character column with IDs of the dams}
\item{Sex}{Character column with sexes of the individuals denoted as \code{"male"} and \code{"female"}.}
\item{Breed}{Character column with adjusted breed names of the individuals (only if \code{Pedig} has column \code{Breed}.).}
\item{Born}{Numeric column with Year-of-Birth of the individuals (only if \code{Pedig} has column \code{Born}.).}
\item{I}{Numeric column with the average age of the parents when the respective individual was born (only if \code{Pedig} has column \code{Born}.).}
\item{numIndiv}{Numeric IDs of the individuals, which are equal to the row numbers (only if \code{addNum=TRUE}). }
\item{numSire}{Numeric IDs of the sires (only if \code{addNum=TRUE}).}
\item{numDam}{Numeric IDs of the dams (only if \code{addNum=TRUE}).}
\item{Offspring}{Logical column indicating the individuals with offspring.}
}


\examples{
data(PedigWithErrors)
Pedig <- prePed(PedigWithErrors)

tail(Pedig)
hist(Pedig$I, freq=FALSE, ylim=c(0,0.2))
}



\author{Robin Wellmann}
