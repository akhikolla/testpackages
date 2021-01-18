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

labkey.executeSql <- function(baseUrl=NULL, folderPath, schemaName, sql, maxRows=NULL,
        rowOffset=NULL, colSort=NULL, showHidden=FALSE, colNameOpt='caption',
        containerFilter=NULL, parameters=NULL)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(schemaName)) stop (paste("A value must be specified for schemaName."))
    if (missing(sql)) stop (paste("A value must be specified for sql."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    ## Construct url
    myurl <- paste(baseUrl, "query", folderPath, "executeSql.api", sep="")

    ## Construct parameters
    params <- list(schemaName=schemaName, apiVersion=8.3, sql=sql)
    if(is.null(maxRows)==FALSE) {params <- c(params, list(maxRows=maxRows))}
    if(is.null(maxRows)==TRUE) {params <- c(params, list(maxRows="-1"))}
    if(is.null(rowOffset)==FALSE) {params <- c(params, list(offset=rowOffset))}
    if(is.null(colSort)==FALSE) {params <- c(params, list(query.sort=colSort))}
    if(is.null(parameters)==FALSE) {for(k in 1:length(parameters)) params <- c(params, list("query.param."=parameters[k]))}
    if(is.null(containerFilter)==FALSE) {params <- c(params, list("containerFilter"=containerFilter))}

    ## Execute via our standard POST function
    mydata <- labkey.post(myurl, toJSON(params, auto_unbox=TRUE))

    newdata <- makeDF(rawdata=mydata, showHidden=showHidden, colNameOpt=colNameOpt)

    return(newdata)
}

