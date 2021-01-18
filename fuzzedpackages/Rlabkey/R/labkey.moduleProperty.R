##
#  Copyright (c) 2018 LabKey Corporation
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


## Get module property for folder
##
labkey.getModuleProperty <- function(baseUrl=NULL, folderPath, moduleName, propName)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(moduleName)) stop (paste("A value must be specified for moduleName."))
    if (missing(propName)) stop (paste("A value must be specified for propName."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    url <- paste(baseUrl, "project", folderPath, "getContainers.api", sep="")
    params <- list(moduleProperties=c(moduleName))
    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

    if (exists("response"))
    {
        result <- (fromJSON(response))
        if (is.null(result$moduleProperties))
        {
            return (paste("Module property does not exist for", moduleName, "module in folder", folderPath))
        }

        moduleProperties <- result$moduleProperties

        for(i in 1:nrow(moduleProperties))
        {
            row <- moduleProperties[i,]
            if (tolower(row$module) == tolower(moduleName) && row$name == propName)
            {
                return (row$effectiveValue)
            }
        }

        return (paste("Module property", propName, "does not exist for", moduleName, "module in folder", folderPath))
    }

    return (paste("Failed to get", moduleName, "module property for folder", folderPath))
}

## Set module property for folder
##
labkey.setModuleProperty <- function(baseUrl=NULL, folderPath, moduleName, propName, propValue)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(moduleName)) stop (paste("A value must be specified for moduleName."))
    if (missing(propName)) stop (paste("A value must be specified for propName."))
    if (missing(propValue)) stop (paste("A value must be specified for propValue."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    property <- list()
    property$moduleName = moduleName
    property$userId = 0 ##currently do not support individualized properties
    property$propName = propName
    property$value = propValue
    property$currentContainer = TRUE

    params <- list(properties=list(property))

    url <- paste(baseUrl, "core", folderPath, "saveModuleProperties.api", sep="")
    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

	return (fromJSON(response, simplifyVector=FALSE, simplifyDataFrame=FALSE))
}