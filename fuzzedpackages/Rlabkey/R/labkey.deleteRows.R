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

labkey.deleteRows <- function(baseUrl=NULL, folderPath, schemaName, queryName, toDelete, provenanceParams=NULL)
{  
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## Default showAllRows=TRUE
    showAllRows=TRUE

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(schemaName)) stop (paste("A value must be specified for schemaName."))
    if (missing(queryName)) stop (paste("A value must be specified for queryName."))
    if (missing(toDelete)) stop (paste("A value must be specified for toDelete."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    ## URL encode folder path, JSON encode post body (if not already encoded)
    toDelete <- convertFactorsToStrings(toDelete);
    nrows <- nrow(toDelete)
    ncols <- ncol(toDelete)

    params <- list(schemaName=schemaName, queryName=queryName, apiVersion=8.3)
    if (!missing(provenanceParams))
        params$provenance = provenanceParams

    p1 <- toJSON(params, auto_unbox=TRUE)
    cnames <- colnames(toDelete)
    p3 <- NULL
    for(j in 1:nrows)
    {
        cvalues <- as.list(toDelete[j,])
        names(cvalues) <- cnames
        p2 <- toJSON(cvalues, auto_unbox=TRUE)
        p3 <- c(p3, p2)
    }
    p3 <- paste(p3, collapse=",")
    pbody <- paste(substr(p1, 1, nchar(p1) - 1),', \"rows\":[',p3,"] }",sep="")

    myurl <- paste(baseUrl, "query", folderPath, "deleteRows.api", sep="")

    ## Execute via our standard POST function
    mydata <- labkey.post(myurl, pbody)
    newdata <- fromJSON(mydata, simplifyVector=FALSE, simplifyDataFrame=FALSE)

    return(newdata)
}
                                                              
labkey.truncateTable <- function(baseUrl=NULL, folderPath, schemaName, queryName)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## check required parameters
    if (missing(baseUrl) || is.null(baseUrl) || missing(folderPath))
        stop (paste("A value must be specified for each of baseUrl and folderPath."))

    ## Validate required parameters
    if (missing(schemaName)) stop (paste("A value must be specified for schemaName."))
    if (missing(queryName)) stop (paste("A value must be specified for queryName."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    url <- paste(baseUrl, "query", folderPath, "truncateTable.api", sep="")

    params <- list(schemaName=schemaName, queryName=queryName)
    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

    return (fromJSON(response))
}
