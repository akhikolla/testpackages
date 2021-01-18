##
#  Copyright (c) 2020 LabKey Corporation
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

## Helper function to create a provenance param object that can be used in subsequent
## provenance functions
##
labkey.provenance.createProvenanceParams <- function(recordingId=NULL, name=NULL, description=NULL,
        materialInputs=NULL, materialOutputs=NULL, dataInputs=NULL, dataOutputs=NULL,
        inputObjectUriProperty=NULL, outputObjectUriProperty=NULL, objectInputs=NULL, objectOutputs=NULL,
        provenanceMap=NULL)
{
    param <- list()

    if (!missing(recordingId))
        param$recordingId = recordingId
    if (!missing(name))
        param$name = name
    if (!missing(description))
        param$description = description

    if (!missing(materialInputs))
    {
        if (!is.list(materialInputs))
            stop (paste("The 'materialInputs' parameter must be a list of material inputs."))

        param$materialInputs = materialInputs
    }

    if (!missing(materialOutputs))
    {
        if (!is.list(materialOutputs))
            stop (paste("The 'materialOutputs' parameter must be a list of material outputs."))

        param$materialOutputs = materialOutputs
    }

    if (!missing(dataInputs))
    {
        if (!is.list(dataInputs))
            stop (paste("The 'dataInputs' parameter must be a list of data inputs."))

        param$dataInputs = dataInputs
    }

    if (!missing(dataOutputs))
    {
        if (!is.list(dataOutputs))
            stop (paste("The 'dataOutputs' parameter must be a list of data outputs."))

        param$dataOutputs = dataOutputs
    }

    if (!missing(inputObjectUriProperty))
        param$inputObjectUriProperty = inputObjectUriProperty
    if (!missing(outputObjectUriProperty))
        param$outputObjectUriProperty = outputObjectUriProperty

    if (!missing(objectInputs))
        param$objectInputs = objectInputs
    if (!missing(objectOutputs))
        param$objectOutputs = objectOutputs

    if (!missing(provenanceMap))
    {
        if (!is.list(provenanceMap))
            stop (paste("The 'provenanceMap' parameter must be a list of provenance map entries."))

        param$provenanceMap = provenanceMap
    }

    return (param)
}

labkey.provenance.startRecording <- function(baseUrl=NULL, folderPath, provenanceParams = NULL)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## check required parameters
    if (missing(baseUrl) || is.null(baseUrl) || missing(folderPath))
        stop (paste("A value must be specified for each of baseUrl and folderPath."))

    if (is.null(provenanceParams))
        stop (paste("Provenance start recording must include the provenanceParams."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    url <- paste(baseUrl, "provenance", folderPath, "startRecording.api", sep="")
    response <- labkey.post(url, toJSON(provenanceParams, auto_unbox=TRUE))

    return (fromJSON(response))
}

labkey.provenance.addRecordingStep <- function(baseUrl=NULL, folderPath, provenanceParams = NULL)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## check required parameters
    if (missing(baseUrl) || is.null(baseUrl) || missing(folderPath))
        stop (paste("A value must be specified for each of baseUrl and folderPath."))

    if (is.null(provenanceParams))
        stop (paste("Provenance start recording must include the provenanceParams."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    url <- paste(baseUrl, "provenance", folderPath, "addRecordingStep.api", sep="")
    response <- labkey.post(url, toJSON(provenanceParams, auto_unbox=TRUE))

    return (fromJSON(response))
}

labkey.provenance.stopRecording <- function(baseUrl=NULL, folderPath, provenanceParams = NULL)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## check required parameters
    if (missing(baseUrl) || is.null(baseUrl) || missing(folderPath))
        stop (paste("A value must be specified for each of baseUrl and folderPath."))

    if (is.null(provenanceParams))
        stop (paste("Provenance start recording must include the provenanceParams."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    url <- paste(baseUrl, "provenance", folderPath, "stopRecording.api", sep="")
    response <- labkey.post(url, toJSON(provenanceParams, auto_unbox=TRUE))

    return (fromJSON(response))
}