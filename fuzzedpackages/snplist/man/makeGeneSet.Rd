\name{makeGeneSet}
\alias{makeGeneSet}
\title{
  makeGeneSet
}
\description{
    This function uses existing SQLite tables (from \code{setGeneTable} and \code{setSNPTable}) to make SNP sets. The SNP Set for each gene is the collection of SNPs located either between the start and end locations of the gene, or within a specified neighborhood around the gene. The SNP Sets are stored in the SQLite database, and returned as a list object.
}
\usage{
    makeGeneSet(annoInfo=NULL,margin=0,annoTable='anno',geneTable='gene',
                        allTable='allchrpos',db='snplistdb',dbCleanUp=FALSE)
}
\arguments{
  \item{annoInfo}{
    A \code{vector} of rsIDs, a \code{data.frame} with an 'rsid' column, or a file with one rsID per line. The SNP sets will be restricted to contain only the SNPs listed here. Default is NULL, in which case all SNPs present in the SNP table in the SQLite database will be used.
}
  \item{margin}{
    A number, indicating the size of the neighborhood (in base pairs) surrounding a genes start and end positions in which a SNP will be included in that genes SNP set. Default is 0.
}
  \item{annoTable}{
    A string indicating the name of the SQLite table for the rsIDs from \code{annoInfo}. Also used in naming the resulting table of SNP sets ('<name>ToGene').  Default is 'anno'.
}
  \item{geneTable}{
    Name of the SQLite table containg chromosome, start and end positions for each gene, as previously created by \code{setGeneTable}. Default is 'gene'.
}
  \item{allTable}{
   Name of the SQLite table containg chromosome and position for each SNP, as previously created by \code{setSNPTable}. Default is 'allchrpos'.
}
  \item{db}{
   Name of the SQLite database in which to find the gene and SNP tables and create the SNP set table. Default is 'snplistdb'.
}
  \item{dbCleanUp}{
    Boolean indicating if the tables and views created by the function should be dropped after the SNP set is returned. Default is FALSE.
}
}
\details{
  Note: This function relies on the prior execution of the \code{setGeneTable} and \code{setSNPTable} functions and the SQLite databes and tables they create. If the \code{table} or \code{db} argument in either of those functions is changed from the default value, it must also be changed here.
}
\value{
  Returns a \code{list} of SNP sets of the form:
    \item{<gene name>}{Vector of rsIDs of SNPs within <gene> (or the neighborhood around it)}
}

\seealso{
    \code{\link{setGeneTable}}, \code{\link{setSNPTable}}, \code{\link{snplist-package}}
}
\examples{
    # Please see the vignette or the package description 
    # for an example of using this function.
}
\keyword{ database }
