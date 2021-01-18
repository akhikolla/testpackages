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

labkey.selectRows <- function(baseUrl=NULL, folderPath, schemaName, queryName, viewName=NULL, colSelect=NULL,
        maxRows=NULL, rowOffset=NULL, colSort=NULL, colFilter=NULL, showHidden=FALSE, colNameOpt='caption',
        containerFilter=NULL, parameters=NULL, includeDisplayValues=FALSE, method='POST')
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## Empty string/NULL checking
    if(is.null(viewName)==FALSE) {char <- nchar(viewName); if(char<1){viewName<-NULL}}
    if(is.null(colSelect)==FALSE) {char <- nchar(colSelect[1]); if(char<1){colSelect<-NULL}}
    if(is.null(maxRows)==FALSE) {char <- nchar(maxRows); if(char<1){maxRows<-NULL}}
    if(is.null(rowOffset)==FALSE) {char <- nchar(rowOffset); if(char<1){rowOffset<-NULL}}
    if(is.null(colSort)==FALSE) {char <- nchar(colSort); if(char<1){colSort<-NULL}}
    if(is.null(colFilter)==FALSE) {char <- nchar(colFilter[1]); if(char<1){colFilter<-NULL}}
    if(is.null(showHidden)==FALSE) {char <- nchar(showHidden); if(char<1){showHidden<-FALSE}}
    if(is.null(containerFilter)==FALSE) {char <- nchar(containerFilter[1]); if(char<1){containerFilter<-NULL}}
    if(is.null(parameters)==FALSE) {char <- nchar(parameters[1]); if(char<1){parameters<-NULL}}
    if(is.null(includeDisplayValues)==FALSE) {char <- nchar(includeDisplayValues); if(char<1){includeDisplayValues<-FALSE}}

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(schemaName)) stop (paste("A value must be specified for schemaName."))
    if (missing(queryName)) stop (paste("A value must be specified for queryName."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    apiVersion = "8.3"

    ## Format colSelect
    colSelect2=NULL
    if(is.null(colSelect)==FALSE) {
        lencolSel <- length(colSelect)
        holder <- NULL
        for(i in 1:length(colSelect)) {
            holder <-paste(holder,URLencode(colSelect[i]),",",sep="")
        }
        colSelect2 <- substr(holder, 1, nchar(holder)-1)
        colSelect <- paste(colSelect, collapse=",")

        # when using colSelect, always set showHidden to TRUE
        showHidden = TRUE
    }

    if(is.null(method) == FALSE && method == "GET")
    {
        ## URL encoding of schema, query, view, etc. (if not already encoded)
        if(schemaName==URLdecode(schemaName)) {schemaName <- URLencode(schemaName)}
        if(queryName==URLdecode(queryName)) {queryName <- URLencode(queryName)}
        if(is.null(viewName)==FALSE) {if(viewName==URLdecode(viewName)) viewName <- URLencode(viewName)}
        if(is.null(containerFilter)==FALSE) {if(containerFilter==URLdecode(containerFilter)) containerFilter<- URLencode(containerFilter)}
        if(is.null(colSort)==FALSE) {if(colSort==URLdecode(colSort)) colSort <- URLencode(colSort)}

        ## Construct url
        myurl <- paste(baseUrl,"query",folderPath,"selectRows.api?schemaName=",schemaName,"&query.queryName=",queryName,"&apiVersion=",apiVersion,sep="")
        if (!is.null(includeDisplayValues) && includeDisplayValues == TRUE) {myurl <- paste(myurl,"&includeDisplayValues=true",sep="")}
        if(is.null(viewName)==FALSE) {myurl <- paste(myurl,"&query.viewName=",viewName,sep="")}
        if(is.null(colSelect2)==FALSE) {myurl <- paste(myurl,"&query.columns=",colSelect2,sep="")}
        if(is.null(maxRows)==FALSE) {myurl <- paste(myurl,"&query.maxRows=",maxRows,sep="")}
        if(is.null(maxRows)==TRUE) {myurl <- paste(myurl,"&query.maxRows=-1",sep="")}
        if(is.null(rowOffset)==FALSE) {myurl <- paste(myurl,"&query.offset=",rowOffset,sep="")}
        if(is.null(colSort)==FALSE) {myurl <- paste(myurl,"&query.sort=",colSort,sep="")}
        if(is.null(colFilter)==FALSE) {for(j in 1:length(colFilter)) myurl <- paste(myurl,"&query.",colFilter[j],sep="")}
        if(is.null(parameters)==FALSE) {for(k in 1:length(parameters)) myurl <- paste(myurl,"&query.param.",parameters[k],sep="")}
        if(is.null(containerFilter)==FALSE) {myurl <- paste(myurl,"&containerFilter=",containerFilter,sep="")}

        ## Execute via our standard GET function
        mydata <- labkey.get(myurl);
    }
    else
    {
        ## Construct url and parameters
        myurl <- paste(baseUrl, "query", folderPath, "selectRows.api", sep="")
        params <- list(schemaName=schemaName, queryName=queryName, apiVersion=apiVersion)
        if (!is.null(includeDisplayValues) && includeDisplayValues == TRUE) {params <- c(params, list(includeDisplayValues="true"))}
        if(is.null(containerFilter)==FALSE) {params <- c(params, list(containerFilter=containerFilter))}
        if(is.null(viewName)==FALSE) {params <- c(params, list(viewName=viewName))}
        if(is.null(colSelect)==FALSE) {params <- c(params, list("query.columns"=colSelect))}
        if(is.null(maxRows)==FALSE) {params <- c(params, list("query.maxRows"=maxRows))}
        if(is.null(maxRows)==TRUE) {params <- c(params, list("query.maxRows"=-1))}
        if(is.null(rowOffset)==FALSE) {params <- c(params, list("query.offset"=rowOffset))}
        if(is.null(colSort)==FALSE) {params <- c(params, list("query.sort"=colSort))}
        if(is.null(colFilter)==FALSE) {for(j in 1:length(colFilter)) {
            # note that the makFilter call uses URLencode() so we need to unescape here
            key = paste("query.",URLdecode(strsplit(colFilter[j],"=")[[1]][1]),sep="")
            value = URLdecode(strsplit(colFilter[j],"=")[[1]][2])
            params[key] = value
        }}
        if(is.null(parameters)==FALSE) {for(k in 1:length(parameters)) {
            key = paste("query.param.",strsplit(parameters[k],"=")[[1]][1],sep="")
            value = strsplit(parameters[k],"=")[[1]][2]
            params[key] = value
        }}

        ## Execute via our standard POST function
        mydata <- labkey.post(myurl, toJSON(params, auto_unbox=TRUE))
    }

    newdata <- makeDF(mydata, colSelect, showHidden, colNameOpt)

    ## Check for less columns returned than requested
    if(is.null(colSelect)==FALSE){if(ncol(newdata)<lencolSel)warning("Fewer columns are returned than were requested in the colSelect variable. The column names may be invalid. Be sure to use the column name and not the column caption. See the documentation for further explaination.")}

    return(newdata)
}

