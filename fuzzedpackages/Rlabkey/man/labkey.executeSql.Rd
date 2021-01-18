\name{labkey.executeSql}
\alias{labkey.executeSql}
\title{Retrieve data from a LabKey Server using SQL commands}
\description{
Use Sql commands to specify data to be imported into R. Prior to import, data can
be manipulated through standard SQL commands supported in LabKey SQL.
}
\usage{
labkey.executeSql(baseUrl, folderPath, schemaName, sql,
    maxRows = NULL, rowOffset = NULL, colSort=NULL,
    showHidden = FALSE, colNameOpt='caption',
    containerFilter=NULL, parameters=NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the  \code{schemaName} for the query}
  \item{sql}{a string containing the \code{sql} commands to be executed}
  \item{maxRows}{(optional) an integer specifying the maximum number of rows to return. If no value is specified, all rows are returned.}
  \item{rowOffset}{(optional) an integer specifying which row of data should be the first row in the retrieval.
If no value is specified, rows will begin at the start of the result set.}
  \item{colSort}{(optional) a string including the name of the column to sort preceeded by a \dQuote{+} or
                 \dQuote{-} to indicate sort direction}
  \item{showHidden}{(optional) a logical value indicating whether or not to return data columns that would
normally be hidden from user view. Defaults to FALSE if no value provided.}
  \item{colNameOpt}{(optional) controls the name source for the columns of the output dataframe, with valid values of 'caption', 'fieldname', and 'rname'
  See labkey.selectRows for more details.}
  \item{containerFilter}{(optional) Specifies the containers to include in the scope of selectRows request. A value of NULL is equivalent to "Current".   Valid values are
     \itemize{
  	      \item{ "Current": Include the current folder only}
  	      \item{ "CurrentAndSubfolders": Include the current folder and all subfolders }
  	      \item{ "CurrentPlusProject": Include the current folder and the project that contains it }
  	      \item{ "CurrentAndParents": Include the current folder and its parent folders }
  	      \item{ "CurrentPlusProjectAndShared": Include the current folder plus its project plus any shared folders }
  	      \item{ "AllFolders": Include all folders for which the user has read permission }
	}
  }
  \item{parameters}{(optional) List of name/value pairs for the parameters if the SQL references underlying queries that are parameterized. For example, parameters=c("X=1", "Y=2").}
}

\details{
A full dataset or any portion of a dataset can be imported into an R
data frame using the \code{labkey.executeSql}
function. Function arguments are components of the url that identify the location of the
data and the SQL actions that should be taken on the data prior to import.

See \code{\link{labkey.selectRows}} for a discussion of the valid options and defaults for \code{colNameOpt}.
}
\value{The requested data are returned in a data frame with stringsAsFactors set to FALSE.  Column names are set as determined by the colNameOpt parameter.
}
\author{Valerie Obenchain}
\seealso{
\code{\link{labkey.selectRows}}, \code{\link{makeFilter}}, \code{\link{labkey.insertRows}}, 
\code{\link{labkey.importRows}}, \code{\link{labkey.updateRows}}, \cr
\code{\link{labkey.deleteRows}}, \code{\link{getRows}}
}
\examples{
\donttest{

## Example of an expicit join and use of group by and aggregates
# library(Rlabkey)

sql<- "SELECT AllTypesCategories.Category AS Category, 
    SUM(AllTypes.IntFld) AS SumOfIntFld,
    AVG(AllTypes.DoubleFld) AS AvgOfDoubleFld
    FROM AllTypes LEFT JOIN AllTypesCategories
    ON (AllTypes.Category = AllTypesCategories.TextKey)
    WHERE AllTypes.Category IS NOT NULL
    GROUP BY AllTypesCategories.Category"

sqlResults <- labkey.executeSql(
    baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples",
    schemaName="lists",
    sql = sql)

sqlResults

}
}
\keyword{IO}
