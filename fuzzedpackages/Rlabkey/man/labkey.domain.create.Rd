\name{labkey.domain.create}
\alias{labkey.domain.create}
\title{Create a new LabKey domain}
\description{
Create a domain of the type specified by the domainKind and the domainDesign. A LabKey domain represents a table in a specific schema.
}
\usage{
labkey.domain.create(baseUrl=NULL, folderPath,
    domainKind=NULL, domainDesign=NULL, options=NULL,
    module=NULL, domainGroup=NULL, domainTemplate=NULL,
    createDomain=TRUE, importData=TRUE)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{domainKind}{(optional) a string specifying the type of domain to create }
  \item{domainDesign}{(optional) a list containing the domain design to create }
  \item{options}{(optional) a list containing options specific to the domain kind }
  \item{module}{(optional) the name of the module that contains the domain template group }
  \item{domainGroup}{(optional) the name of a domain template group }
  \item{domainTemplate}{(optional) the name of a domain template within the domain group }
  \item{createDomain}{(optional) when using a domain template, create the domain. Defaults to TRUE }
  \item{importData}{(optional) when using a domain template, import initial data asociated in the template. Defaults to TRUE }
}
\details{
When creating a domain using a domainKind parameter, the domainDesign parameter will be required. If a domain template is being
used, then module, domainGroup, and domainTemplate are required.

Will create a domain of the specified domain type, valid types are
     \itemize{
  	      \item{ "IntList": A list with an integer key field }
  	      \item{ "VarList": A list with a string key field }
  	      \item{ "StudyDatasetVisit": A dataset in a visit based study }
  	      \item{ "StudyDatasetDate": A dataset in a date based study }
  	      \item{ "IssueDefinition": An issue list domain }
  	      \item{ "SampleSet": Sample set }
  	      \item{ "DataClass": Data class }
	}
The domain design parameter describes the set of fields in the domain, see \code{labkey.domain.createDesign} for the helper
function that can be used to construct this data structure. The options parameter should contain a list of attributes that are specific
to the domain kind specified. The list of valid options for each domain kind are:
     \itemize{
  	      \item{ IntList and VarList }
  	      \itemize{
  	        \item{ \code{keyName (required)} : The name of the field in the domain design which identifies the key field}
  	      }
  	      \item{ StudyDatasetVisit and StudyDatasetDate }
  	      \itemize{
  	        \item{ \code{datasetId} : Specifies a dataset ID to use, the default is to auto generate an ID}
  	        \item{ \code{categoryId} : Specifies an existing category ID}
  	        \item{ \code{categoryName} : Specifies an existing category name}
  	        \item{ \code{demographics} : (TRUE | FALSE) Determines whether the dataset is created as demographic}
  	        \item{ \code{keyPropertyName} : The name of an additional key field to be used in conjunction with participantId and (visitId or date) to create unique records}
  	        \item{ \code{useTimeKeyField} : (TRUE | FALSE) Specifies to use the time portion of the date field as an additional key }
  	        \item{ \code{isManagedField} : (TRUE | FALSE) Specifies whether the field from keyPropertyName should be managed by LabKey. }
  	      }
  	      \item{ IssueDefinition }
  	      \itemize{
  	        \item{ \code{providerName} : The type of issue list to create (IssueDefinition (default) or AssayRequestDefinition)}
  	        \item{ \code{singularNoun} : The singular name to use for items in the issue definition (defaults to issue)}
  	        \item{ \code{pluralNoun} : The plural name (defaults to issues)}
  	      }
  	      \item{ SampleSet }
  	      \itemize{
  	        \item{ \code{idCols} : The columns to use when constructing the concatenated unique ID. Can be up to 3 numeric IDs which represent the zero-based position of the fields in the domain.}
  	        \item{ \code{parentCol} : The column to represent the parent identifier in the sample set. This is a numeric value representing the zero-based position of the field in the domain.}
  	        \item{ \code{nameExpression} : The name expression to use for creating unique IDs}
  	      }
  	      \item{ DataClass }
  	      \itemize{
  	        \item{ \code{sampleSet} : The ID of the sample set if this data class is associated with a sample set.}
  	        \item{ \code{nameExpression} : The name expression to use for creating unique IDs}
  	      }
	}
}
\value{
A list containing elements describing the newly created domain.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.domain.get}},
\code{\link{labkey.domain.inferFields}},
\code{\link{labkey.domain.createDesign}},
\code{\link{labkey.domain.createIndices}},
\code{\link{labkey.domain.save}},
\code{\link{labkey.domain.drop}},
\code{\link{labkey.domain.createConditionalFormat}},
\code{\link{labkey.domain.createConditionalFormatQueryFilter}},
\code{\link{labkey.domain.FILTER_TYPES}}
}
\examples{
\donttest{

## create a data frame and infer it's fields, then create a domain design from it
library(Rlabkey)

df <- data.frame(ptid=c(1:3), age = c(10,20,30), sex = c("f", "m", "f"))
fields <- labkey.domain.inferFields(baseUrl="http://labkey/", folderPath="home", df=df)
dd <- labkey.domain.createDesign(name="test list", fields=fields)

## create a new list with an integer key field
labkey.domain.create(baseUrl="http://labkey/", folderPath="home",
    domainKind="IntList", domainDesign=dd, options=list(keyName = "ptid"))

## create a domain using a domain template
labkey.domain.create(baseUrl="http://labkey/", folderPath="home",
    domainTemplate="Priority", module="simpletest", domainGroup="todolist")
}
}
\keyword{IO}
