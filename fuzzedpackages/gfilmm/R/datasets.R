#' @title Krishnamoorthy & Mathew's example 4.1
#' @description The dataset used in Krishnammorthy & Mathew's example 4.1.
#' @name KM41
#' @docType data
#' @references Krishnamoorthy and Mathew, Statistical Tolerance Regions, Wiley 2009.
#' @keywords data
#' @usage data(KM41)
#' @format A data frame with 25 rows and 2 columns.
#' @examples
#' data(KM41)
#' str(KM41)
#' table(KM41$Batch)
NULL

#' @title pH dataset
#' @description A dataset from ?? (I don't remember).
#' @name pHdata
#' @docType data
#' @keywords data
#' @usage data(pHdata)
#' @format A data frame with 160 rows and 4 columns. 
#' Column \code{SIRE} is a factor nested in column \code{DAM}.
#' @examples
#' data(pHdata)
#' str(pHdata)
#' table(droplevels(pHdata[pHdata$DAM=="D1","SIRE"]))
#' table(droplevels(pHdata[pHdata$DAM=="D2","SIRE"]))
#' table(droplevels(pHdata[pHdata$DAM=="D3","SIRE"]))
NULL