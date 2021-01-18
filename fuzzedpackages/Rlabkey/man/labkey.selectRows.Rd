\name{labkey.selectRows}
\alias{labkey.selectRows}
\title{Retrieve data from a labkey database}
\description{
Import full datasets or selected rows into R. The data can be sorted and filtered prior to import.
}
\usage{
labkey.selectRows(baseUrl = NULL, folderPath, schemaName, queryName,
    viewName = NULL, colSelect = NULL, maxRows = NULL,
    rowOffset = NULL, colSort = NULL, colFilter = NULL,
    showHidden = FALSE, colNameOpt="caption",
    containerFilter = NULL, parameters = NULL,
    includeDisplayValues = FALSE, method = "POST")
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl}for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{schemaName}{a string specifying the  \code{schemaName} for the query}
  \item{queryName}{a string specifying the \code{queryName}}
  \item{viewName}{(optional) a string specifying the \code{viewName} associated with the query. If not specified, the default view determines the rowset returned.  }
  \item{colSelect}{(optional) a vector of strings specifying which columns of a dataset or view to import.
       \itemize{
    	      \item{ The wildcard character ("*") may also be used here to get all columns including those not in the default view. }
    	      \item{ If you include a column that is a lookup (foreign key) the value of the primary key for that target will be returned. }
    	      \item{ Use a backslash character ("/") to include non-primary key columns from a lookup target (foreign key), e.g "LookupColumnName/targetColumn". }
    	      \item{ When using a string to specify the colSelect set, the column names must be separated by a comma and not include spaces between the names. }
        }
    }
  \item{maxRows}{(optional) an integer specifying how many rows of data to return. If no value is specified, all rows are returned.}
  \item{rowOffset}{(optional) an integer specifying which row of data should be the first row in the retrieval. If no
value is specified, the retrieval starts with the first row.}
  \item{colSort}{(optional) a string including the name of the column to sort preceded by a \dQuote{+} or
\dQuote{-} to indicate sort direction}
  \item{colFilter}{(optional) a vector or array object created by the \code{makeFilter} function which contains the
column name, operator and value of the filter(s) to be applied to the retrieved data.}
  \item{showHidden}{(optional) a logical value indicating whether or not to return data columns that would normally be hidden from user view. Defaults to FALSE if no value provided.}
  \item{colNameOpt}{(optional) controls the name source for the columns of the output dataframe, with valid values of 'caption', 'fieldname', and 'rname'}
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
  \item{includeDisplayValues}{(optional) a logical value indicating if display values should be included in the response object for lookup fields.}
  \item{method}{(optional) HTTP method to use for the request, defaults to POST.}
}
\details{
A full dataset or any portion of a dataset can be downloaded into an R data frame using the
\code{labkey.selectRows} function. Function arguments are the components of the url that identify
the location of the data and what actions should be taken on the data prior to import
(ie, sorting, selecting particular columns or maximum number of rows, etc.).

Stored queries in LabKey Server have an associated default view and may have one or more named views.  Views determine the column set of
the return data frame.  View columns can be a subset or superset of the columns of the underlying query (a subset if columns from the
query are left out of the view, and a superset if lookup columns in the underlying query are used to include columns from related queries).
Views can also include filter and sort properties that will make their result set different from the underlying query.
If no view is specified, the columns and rows returned are determined by the default view, which may not be the same as the 
result rows of the underlying query.  Please see the topic on Saving Views in the LabKey online documentation.

In the returned data frame, there are three different ways to have the columns named:  
\code{colNameOpt='caption'} uses the caption value, and is the default option for backward compatibility.  
It may be the best option for displaying to another user, but may make scripting more difficult.
\code{colNameOpt='fieldname'} uses the field name value, so that the data frame colnames are the same names
that are used as arguments to labkey function calls.  It is the default for the new \code{\link{getRows}} session-based
function.
\code{colNameOpt='rname'} transforms the field name value into valid R names by substituting an underscore for both
spaces and forward slash (/) characters and lower casing the entire name.  This option is the way a data frame is passed to a
script running in a LabKey server in the R View feature of the data grid.  If you are writing scripts for running in an
R view on the server, or if you prefer to work with legal r names in the returned grid, this option may be useful.

For backward compatibility, column names returned by \code{labkey.executeSql} and \code{labkey.selectRows} 
are field captions by default.  The \code{getRows} function has the same
colNameOpt parameter but defaults to field names instead of captions.

}
\value{
The requested data are returned in a data frame with stringsAsFactors set to FALSE.  Column names are set as determined by the colNameOpt parameter.
}
\references{https://www.labkey.org/Documentation/wiki-page.view?name=savingViews}

\author{Valerie Obenchain}
\seealso{
{Retrieve data:  } \code{\link{makeFilter}}, \code{\link{labkey.executeSql}} \cr
{Modify data:  } \code{\link{labkey.updateRows}}, \code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \code{\link{labkey.deleteRows}}\cr
{List available data: } \code{\link{labkey.getSchemas}}, \code{\link{labkey.getQueries}}, \code{\link{labkey.getQueryViews}},
\code{\link{labkey.getQueryDetails}},\code{\link{labkey.getDefaultViewDetails}},\code{\link{labkey.getLookupDetails}},

}
\examples{
\donttest{

## select from a list named AllTypes
# library(Rlabkey)

rows <- labkey.selectRows(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples",
	schemaName="lists", 
	queryName="AllTypes")
	
## select from a view on that list
viewrows <- labkey.selectRows(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="Lists", queryName="AllTypes",
    viewName="rowbyrow")

## select a subset of columns
colSelect=c("TextFld", "IntFld")
subsetcols <- labkey.selectRows(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    colSelect=colSelect)

## including columns from a lookup (foreign key) field
lookupcols <- labkey.selectRows(baseUrl="http://localhost:8080/labkey",
    folderPath="/apisamples", schemaName="lists", queryName="AllTypes",
    colSelect="TextFld,IntFld,IntFld/LookupValue")

}
}
\keyword{IO}
