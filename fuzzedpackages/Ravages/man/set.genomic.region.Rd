\name{set.genomic.region}
\alias{set.genomic.region}
\title{ Variants annotation based on gene positions }
\description{
Attributes regions to variants based on given region positions
}
\usage{
set.genomic.region(x, regions = genes.b37, flank.width = 0L)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{regions}{
    A dataframe containing the fields : \code{Chr} (the chromosome of the gene),
    \code{Start} (the start position of the gene), \code{End} (the end position of the gene), and
    \code{Gene_Name} (the name of the gene - a factor),
  }
  \item{flank.width}{ An integer: width of the flanking regions in base pairs downstream and upstream the regions. }
}
\details{ 
Warnings: \code{regions$Gene_Name} should be a factor containing the names of the regions, ordered in the genome order. 

We provide two data sets of autosomal humain genes, \code{genes.b37} and \code{genes.b38}.  

If \code{x@snps$chr} is not a vector of integers, it should be a factor with same levels as \code{regions$Chr}.

If \code{flank.width} is null, only the variants having their position between the \code{regions$Start} and the \code{regions$End} of a gene will be attributed to the corresponding gene.  
When two regions overlap, variants in the overlapping zone will be assigned to those two regions, separated by a comma.

If \code{flank.width} is a positive number, variants \code{flank.width} downstream or upstream a gene will be annotated annotated to this gene. You can use \code{flank.width = Inf}
to have each variant attributed to the nearest gene.  
}
\value{
The same bed matrix as x with an additional column \code{x@snps$genomic.region} containing the annotation of each variant (character).
}
\seealso{
 \code{\link{genes.b37}}, \code{\link{genes.b38}}
}
\examples{
#Import 1000Genome data from region around LCT gene
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

#Group variants within known genes
x <- set.genomic.region(x)

#Group variants within know genes +/- 500bp
x <- set.genomic.region(x, flank.width=500)
}
