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

labkey.insertRows <- function(baseUrl=NULL, folderPath, schemaName, queryName, toInsert, na=NULL, provenanceParams=NULL)
{  
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(schemaName)) stop (paste("A value must be specified for schemaName."))
    if (missing(queryName)) stop (paste("A value must be specified for queryName."))
    if (missing(toInsert)) stop (paste("A value must be specified for toInsert."))

    ## Default showAllRows=TRUE
    showAllRows=TRUE

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    ## URL encode folder path, JSON encode post body (if not already encoded)
    toInsert <- convertFactorsToStrings(toInsert);
    nrows <- nrow(toInsert)
    ncols <- ncol(toInsert)

    params <- list(schemaName=schemaName, queryName=queryName, apiVersion=8.3)
    if (!missing(provenanceParams))
        params$provenance = provenanceParams

    p1 <- toJSON(params, auto_unbox=TRUE)
    cnames <- colnames(toInsert)
    p3 <- NULL
    for(j in 1:nrows)
    {
        cvalues <- as.list(toInsert[j,])
        names(cvalues) <- cnames
        cvalues[is.na(cvalues)] = na
        p2 <- toJSON(cvalues, auto_unbox=TRUE)
        p3 <- c(p3, p2)
    }
    p3 <- paste(p3, collapse=",")
    pbody <- paste(substr(p1, 1, nchar(p1)-1), ', \"rows\":[' ,p3, "] }", sep="")

    myurl <- paste(baseUrl, "query", folderPath, "insertRows.api", sep="")

    ## Execute via our standard POST function
    mydata <- labkey.post(myurl, pbody)
    newdata <- fromJSON(mydata, simplifyVector=FALSE, simplifyDataFrame=FALSE)

    return(newdata)
}
                                                              
