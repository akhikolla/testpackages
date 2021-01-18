\name{bed.matrix.split.genomic.region}
\alias{bed.matrix.split.genomic.region}
\title{ Bed matrix for variants associated to multiple genomic regions}
\description{
Creates a new bed matrix with variants associated to multiple genomic regions being duplicated
}
\usage{
bed.matrix.split.genomic.region(x, changeID=TRUE, genomic.region=NULL, 
                                split.pattern=",")
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{changeID}{ TRUE/FALSE: whether to change the variants ID by including the gene name}
  \item{genomic.region}{ A vector containing the genomic region of each variant }
  \item{split.pattern}{ The character separating the genomic regions}
}
\value{
A bed matrix with variants assigned to multiple genomic regions being duplicated and the corresponding genomic regions separated and transformed into factors.
}
\details{
  If \code{changeID=TRUE}, variants will have new IDs being CHR:POS:A1:A2:genomic.region.  

  The genomic region(s) associated to each varaint should be in \code{x@snps$genomic.region} or given as a vector to \code{genomic.region}. If both are present, \code{genomic.region} is used.
}
\examples{
#Example bed matrix with 4 variants
  x.ex <- as.bed.matrix(x=matrix(0, ncol=4, nrow=10), 
                       bim=data.frame(chr=1:4, id=paste("rs", 1:4, sep=""), dist = rep(0,4), 
                                      pos=c(150,150,200,250), A1=rep("A", 4), A2=rep("T", 4)))

#Example genes dataframe
  genes.ex <- data.frame(Chr=c(1,1,3,4), Start=c(10,110,190,220), End=c(170,180,250,260), 
                         Gene_Name=factor(letters[1:4]))

#Attribute genomic regions
  x.ex <- set.genomic.region(x.ex, regions = genes.ex)

#Split genomic regions
  x.ex.split <- bed.matrix.split.genomic.region(x.ex, split.pattern = ",")
}

