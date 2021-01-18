\name{labkey.domain.createAndLoad}
\alias{labkey.domain.createAndLoad}
\title{Create a new LabKey domain and load data}
\description{
Create a domain of the type specified by the domainKind. A LabKey domain represents a table in a specific schema. Once the domain is created
the data from the data frame will be imported.
}
\usage{
labkey.domain.createAndLoad(baseUrl=NULL, folderPath,
    name, description="", df, domainKind, options=NULL, schemaName=NULL)
}
\arguments{
  \item{baseUrl}{a string specifying the \code{baseUrl} for the labkey server}
  \item{folderPath}{a string specifying the \code{folderPath} }
  \item{name}{a string specifying the name of the domain to create }
  \item{description}{(optional) a string specifying the domain description }
  \item{df}{a data frame specifying fields to infer. The data frame must have column names as well as row data to infer the type of the field from.}
  \item{domainKind}{a string specifying the type of domain to create }
  \item{options}{(optional) a list containing options specific to the domain kind }
  \item{schemaName}{(optional) a string specifying the schema name to import the data into }
}
\details{
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

The options parameter should contain a list of attributes that are specific
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
A list containing the newly uploaded data frame rows.
}

\author{Karl Lum}
\seealso{
\code{\link{labkey.domain.get}},
\code{\link{labkey.domain.inferFields}},
\code{\link{labkey.domain.createDesign}},
\code{\link{labkey.domain.createIndices}},
\code{\link{labkey.domain.save}},
\code{\link{labkey.domain.drop}}
}
\examples{
\donttest{

library(Rlabkey)

## Prepare a data.frame
participants = c("0001","0001","0002","0002","0007","0008")
Visit = c("V1", "V2", "V2", "V1", "V2", "V1")
IntValue = c(256:261)
dataset = data.frame("ParticipantID" = participants, Visit,
    "IntegerValue" = IntValue, check.names = FALSE)

## Create the dataset and import
labkey.domain.createAndLoad(baseUrl="http://labkey", folderPath="home",
    name="demo dataset", df=dataset, domainKind="StudyDatasetVisit")
}
}
\keyword{IO}
