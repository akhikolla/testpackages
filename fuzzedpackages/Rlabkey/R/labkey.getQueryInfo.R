##
#  Copyright (c) 2010-2018 LabKey Corporation
# 
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
##

##labkey.getLookupDetails
labkey.getLookupDetails <- function(baseUrl=NULL, folderPath, schemaName, queryName, lookupKey)
{
	if(missing("lookupKey") )
		{stop ("You must supply the key (name) value of a query field defined as a lookup type field.")}

	lookupFields <- getQueryInfo(baseUrl=baseUrl, folderPath=folderPath, schemaName=schemaName, queryName=queryName,showDefaultView=FALSE, lookupKey=lookupKey)
	return(lookupFields)
}

##Public getQueryDetails
labkey.getQueryDetails <- function(baseUrl=NULL, folderPath, schemaName, queryName)
{
	queryDetails <- getQueryInfo(baseUrl=baseUrl, folderPath=folderPath, schemaName=schemaName, queryName=queryName,showDefaultView=FALSE)
	return(queryDetails)
}

## Public getDefaultViewDetails
labkey.getDefaultViewDetails <- function(baseUrl=NULL, folderPath, schemaName, queryName)
{
	viewDetails <- getQueryInfo(baseUrl=baseUrl, folderPath=folderPath, schemaName=schemaName, queryName=queryName,showDefaultView=TRUE)
	return(viewDetails)
}

## internal reoutine that handles all of these
getQueryInfo <- function(baseUrl=NULL, folderPath, schemaName, queryName, showDefaultView=FALSE, lookupKey=NULL)
{
	baseUrl=labkey.getBaseUrl(baseUrl)

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(schemaName)) stop (paste("A value must be specified for schemaName."))
    if (missing(queryName)) stop (paste("A value must be specified for queryName."))

	if(is.null(lookupKey)==FALSE) {char <- nchar(lookupKey); if(char<1) {lookupKey<-NULL} }

	## URL encoding (if not already encoded)
	if(schemaName==URLdecode(schemaName)) {schemaName <- URLencode(schemaName)}
	if(queryName==URLdecode(queryName)) {queryName <- URLencode(queryName)}
	if(is.null(lookupKey)==FALSE) {if(lookupKey==URLdecode(lookupKey)) lookupKey <- URLencode(lookupKey)}

	## normalize the folder path
	folderPath <- encodeFolderPath(folderPath)

	## Construct url
	myurl <- paste(baseUrl,"query",folderPath,"getQueryDetails.api?schemaName=", schemaName, "&queryName=", queryName, "&apiVersion=8.3", sep="")
	if(is.null(lookupKey)==FALSE) {myurl <- paste(myurl,"&fk=",lookupKey,sep="")}

	## Execute via our standard GET function
	mydata <- labkey.get(myurl)

	decode <- fromJSON(mydata, simplifyVector=FALSE, simplifyDataFrame=FALSE)

	## If querying the default view, the metadata is in a differnt object in the json stream
	if (showDefaultView==TRUE) {qcs<-decode$defaultView$columns}
	else {qcs <- decode$columns}

	## parsed JSON stream has two types of problems related to nulls:
	## the value NULL as as named element of the parent node
	## the absence of either a value or a name for some columns on some records
	## etiher one can be detected by checking for is.null on a row-by row basis against the total set of column names

	baseColumns <- c("name", "caption", "fieldKey", "type", "isNullable","isKeyField",
				"isAutoIncrement", "isVersionField","isHidden","isSelectable",
				"isUserEditable", "isReadOnly", "isMvEnabled","description")
	lookupColumns <- c("keyColumn", "schemaName", "displayColumn", "queryName", "isPublic")

	dmall <- matrix(nrow=0, ncol=20, byrow=TRUE)
	if(length(qcs)>0)
	{
		for (j in 1:length(qcs))
		{
			dmqrow<- matrix(data=decode$name[[1]], nrow=1, ncol=1, byrow=FALSE)
			for (nm in baseColumns) {
				if (is.null(qcs[[j]][[nm]])) {qcs[[j]][[nm]] <- NA}

				dmqrow<- matrix(data=cbind(dmqrow, qcs[[j]][[nm]]), nrow=1, byrow=FALSE)
			}


			if (is.null(qcs[[j]]$lookup))
			{
				lookupinfo <- matrix(data=cbind(NA,NA,NA,NA,NA), ncol=5, byrow=FALSE)
			}
			else
			{
				for (nm in lookupColumns) {
					if (is.null(qcs[[j]]$lookup[[nm]])) {qcs[[j]]$lookup[[nm]] <- NA}
				}
				nm <- lookupColumns[1]
				lookupinfo<- as.matrix(qcs[[j]]$lookup[[nm]], nrow=1, byrow=FALSE)
				for (nm in lookupColumns[-1]) {
					lookupinfo<- matrix(data=cbind(lookupinfo, qcs[[j]]$lookup[[nm]]), nrow=1, byrow=FALSE)
				}
			}

			dmqrow<-cbind(dmqrow, lookupinfo)
			dmall <- rbind(dmall,dmqrow)
		}
	}
	dfall <- as.data.frame(dmall, stringsAsFactors=FALSE)
	colnames(dfall)<-c("queryName", "fieldName", "caption", "fieldKey", "type", "isNullable","isKeyField",
				"isAutoIncrement", "isVersionField","isHidden","isSelectable",
				"isUserEditable", "isReadOnly", "isMvEnabled", "description",
				"lookupKeyField","lookupSchemaName","lookupDisplayField", "lookupQueryName", "lookupIsPublic")

	return(dfall)
}
