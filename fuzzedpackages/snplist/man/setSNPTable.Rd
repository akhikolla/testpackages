\name{setSNPTable}
\alias{setSNPTable}
\title{
    setSNPTable
}
\description{
    Takes a file or \code{data.frame} object with columns 'chr','pos', and 'rsid', and creates an SQLite table of the information. Returns a count of the number of SNPs in the table.
}
\usage{
    setSNPTable(snpInfo,table='allchrpos',db='snplistdb')
}
\arguments{
  \item{snpInfo}{
    A \code{data.frame} object of SNP location info with columns 'chr','pos', and 'rsid', or a tab-delimited file with those columns and one record per row.
}
  \item{table}{
    Name of the SQLite table to be created. Default is 'allchrpos'.
}
  \item{db}{
    Name of the SQLite database in which to create \code{table}. Default is 'snplistdb'.
}
}

\value{
    Count of genes included in \code{table}.
}

\examples{
    snpInfo <- cbind(c(17,17,13,13),
		     c(41211653, 41213996, 32890026,32890572),
		     c("rs8176273","rs8176265","rs9562605","rs1799943") )
    colnames(snpInfo) <- c('chr','pos','rsid')
    setSNPTable(as.data.frame(snpInfo))
}
\keyword{ database }
