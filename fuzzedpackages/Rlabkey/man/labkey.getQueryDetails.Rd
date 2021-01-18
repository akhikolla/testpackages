\name{labkey.getQueryDetails}
\alias{labkey.getQueryDetails}
\title{Retrieve detailed information on a LabKey query}
\description{
Fetch a list of output columns and their attributes that are avaialble from a given query}
\usage{
labkey.getQueryDetails(baseUrl, folderPath, schemaName, queryName)
}
\arguments{
  \item{baseUrl}{a string specifying the address of the LabKey Server, including the context root}
  \item{folderPath}{a string specifying the hierarchy of folders to the current folder (container) for the operation, starting with the project folder }
  \item{schemaName}{a string specifying the schema name in which the query object is defined}
  \item{queryName}{a string specifying the name of the query}
}
\details{
Queries have a default output list of fields defined by the "default view" of the query.  
To retrieve that set of fields with their detailed properties such as type and nullability, use
\code{labkey.getQueryDetails} function. Function arguments are the components of the url that identify
the location of the server, the folder path, the schema, and the name of the query.
\cr \cr
The results from getQueryDetails describe the \dQuote{field names} that are used to build the colSelect, 
colFilter and colSort parameters to selectRows. Each column in the data frame returned from selectRows
corresponds to a field in the colSelect list.
\cr
There are two types of fieldNames that will be reported by the server in the output of this function. For fields that are
directly defined in the query corresponding the queryName parameter for this function, the fieldName 
is simply the name assigned by the query. Because selectRows returns the results specified by the default view, however,
there may be cases where this default view incorporates data from other queries that have a defined 1-M relationship
with the table designated by the queryName. Such fields in related tables are referred to as \dQuote{lookup} fields.
Lookup fields have multi-part names using a forward slash as the delimiter. For example, in a samples data set, if the ParticipantId
identifies the source of the sample, \code{ParticipantId/CohortId/CohortName} could be a reference to a \code{CohortName} field 
in a \code{Cohorts} data set.

These lookup fieldNames can appear in the default view and show up in the selectRows result. If a field from a lookup table
is not in the default view, it can still be added to the output column list of \code{labkey.selectRows}. Use the
labkey.getLookups to discover what additional fields are available via lookups, and then put their multipart fieldName values
into the colSelect list. Lookup fields have the semantics of a LEFT JOIN in SQL, such that
every record from the target queryName appears in the output whether or not there is a matching lookup field value.

}
\value{
The available schemas are returned as a data frame,

\code{queryName} {the name of the query, repeated n times, where n is the number of output fields from the query} \cr
\code{fieldName} {the fully qualified name of the field, relative to the specified queryName.  } \cr 
\code{caption}  {a more readable label for the data field, appears as a column header in grids} \cr 
\code{fieldKey} {the name part that identifies this field within its containing table, independent of its use as a lookup target.}\cr 
\code{type} {a string specifying the field type, e.g. Text, Number, Date, Integer}\cr 
\code{isNullable} {TRUE if the field can be left empty (null)} \cr
\code{isKeyField} {TRUE if the field is part of the primary key} \cr
\code{isAutoIncrement} {TRUE if the system will automatically assign a sequential integer in this on inserting a record}\cr 
\code{isVersionField} {TRUE if the field issued to detect changes since last read} \cr
\code{isHidden} {TRUE if the field is not displayed by default} \cr
\code{isSelectable} {reserved for future use. }\cr
\code{isUserEditable} {reserved for future use.}\cr
\code{isReadOnly} {reserved for future use}\cr
\code{isMvEnabled} {reserved for future use}\cr
\code{lookupKeyField}  {for a field defined as a lookup the primary key column of the query referenced by the lookup field;  NA for non-lookup fields} \cr
\code{lookupSchemaName} {the schema of the query referenced by the lookup field;  NA for non-lookup fields} \cr
\code{lookupDisplayField} {the field from the query referenced by the lookup field that is shown by default in place of the lookup field;  NA for non-lookup fields} \cr 
\code{lookupQueryName} {the query referenced by the lookup field;  NA for non-lookup fields.  A non-NA value indicates that you can use this field in a call to getLookups} \cr 
\code{lookupIsPublic} {reserved for future use}\cr
}
\author{Peter Hussey, peter@labkey.com}
\seealso{
{Retrieve data:  } \code{\link{labkey.selectRows}}, \code{\link{makeFilter}}, \code{\link{labkey.executeSql}} \cr
{Modify data:  } \code{\link{labkey.updateRows}}, \code{\link{labkey.insertRows}}, \code{\link{labkey.importRows}}, \code{\link{labkey.deleteRows}} \cr
{List available data: } \code{\link{labkey.getSchemas}}, \code{\link{labkey.getQueries}}, \code{\link{labkey.getQueryViews}}, 
\code{\link{labkey.getDefaultViewDetails}}, \code{\link{labkey.getLookupDetails}}
}
\examples{
\donttest{

## Details  of fields of a query
# library(Rlabkey)

queryDF<-labkey.getQueryDetails(
	baseUrl="http://localhost:8080/labkey",
	folderPath="/apisamples",
	schemaName="lists",
	queryName="AllTypes")

}
}
\keyword{IO}
