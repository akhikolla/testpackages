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
labkey.experiment.SAMPLE_DERIVATION_PROTOCOL <- "Sample Derivation Protocol"

## Create an ExpData object
##
labkey.experiment.createData <- function(config, dataClassId = NULL, dataClassName = NULL, dataFileUrl = NULL)
{
    ## check required parameters
    if (missing(config))
        stop (paste("A list of ExpObject config params must be specified in the config parameter."))

    data <- config

    if (!is.null(dataFileUrl))
        data$dataFileURL = dataFileUrl

    if (!is.null(dataClassId) || !is.null(dataClassName))
    {
        dataClass <- list()
        if (!is.null(dataClassId))
            dataClass$id = dataClassId
        if (!is.null(dataClassName))
            dataClass$name = dataClassName

        data$dataClass = dataClass
    }
    return (data)
}

## Create an ExpMaterial object
##
labkey.experiment.createMaterial <- function(config, sampleSetId = NULL, sampleSetName = NULL)
{
    ## check required parameters
    if (missing(config))
        stop (paste("A list of ExpObject config params must be specified in the config parameter."))

    material <- config

    if (!is.null(sampleSetId) || !is.null(sampleSetName))
    {
        sampleSet <- list()
        if (!is.null(sampleSetId))
            sampleSet$id = sampleSetId
        if (!is.null(sampleSetName))
            sampleSet$name = sampleSetName

        material$sampleSet = sampleSet
    }
    return (material)
}

## Create an ExpRun object
##
labkey.experiment.createRun <- function(config, dataInputs = NULL, dataOutputs = NULL, dataRows = NULL, materialInputs = NULL, materialOutputs = NULL, plateMetadata = NULL)
{
    ## check required parameters
    if (missing(config))
        stop (paste("A list of ExpObject config params must be specified in the config parameter."))

    run <- config

    if (!is.null(dataInputs))
    {
        if (!is.list(dataInputs))
            stop (paste("dataInputs must be a list of data objects, see labkey.experiment.createData."))

        ## ensure dataInputs is serialized as an array of objects
        run$dataInputs = ensureNestedList(dataInputs)
    }

    if (!is.null(dataOutputs))
    {
        if (!is.list(dataOutputs))
            stop (paste("dataOutputs must be a list of data objects, see labkey.experiment.createData."))

        run$dataOutputs = ensureNestedList(dataOutputs)
    }

    if (!is.null(materialInputs))
    {
        if (!is.list(materialInputs))
            stop (paste("materialInputs must be a list of material objects, see labkey.experiment.createMaterial."))

        run$materialInputs = ensureNestedList(materialInputs)
    }

    if (!is.null(materialOutputs))
    {
        if (!is.list(materialOutputs))
            stop (paste("materialOutputs must be a list of material objects, see labkey.experiment.createMaterial."))

        run$materialOutputs = ensureNestedList(materialOutputs)
    }

    if (!is.null(dataRows))
    {
        if (!is.data.frame(dataRows))
            stop (paste("dataRows must be a data frame."))

    	## build Assay object tree based on R lists
    	nrows <- nrow(dataRows)
    	ncols <- ncol(dataRows)
    	cnames <- colnames(dataRows)
    	rowsVector <- vector(mode="list", length=nrows)
    	for (j in 1:nrows)
    	{
    		cvalues <- as.list(dataRows[j,])
    		names(cvalues) <- cnames
    		rowsVector[[j]] <- cvalues
    	}

        run$dataRows <- rowsVector
    }

    if (!is.null(plateMetadata))
    {
        run$plateMetadata = plateMetadata;
    }
    return (run)
}

## Helper to ensure the passed object is serialized to JSON as an array of objects
##
ensureNestedList <- function(data)
{
    if (is.null(names(data)))
        data
    else
        data <- list(data)

    return (data)
}

labkey.experiment.saveBatch <- function(baseUrl=NULL, folderPath, assayConfig = NULL, protocolName = NULL, batchPropertyList = NULL, runList)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(runList)) stop (paste("A value must be specified for runList."))

    if (is.null(assayConfig) && is.null(protocolName))
        stop (paste("Either an assay config list or protocolName must be specified. The assay configuration must contain either an assayId or both assayName and providerName"))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

	## Now post form with batch object filled out
    url <- paste(baseUrl, "assay", folderPath, "saveAssayBatch.api", sep="")

    if (!is.null(assayConfig))
        params = assayConfig
    else if (!is.null(protocolName))
    {
        params = list()
        params$protocolName = protocolName
    }
    params$batch = c(batchPropertyList, list(runs = ensureNestedList(runList)))

    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

	return (fromJSON(response, simplifyVector=FALSE, simplifyDataFrame=FALSE))
}

labkey.experiment.saveRuns <- function(baseUrl=NULL, folderPath, protocolName, runList)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(protocolName)) stop (paste("A value must be specified for protocolName."))
    if (missing(runList)) stop (paste("A value must be specified for runList."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    ## Now post form with runs object filled out
    url <- paste(baseUrl, "assay", folderPath, "saveAssayRuns.api", sep="")

    if (!is.null(runList))
    {
        params = list()
        params$protocolName = protocolName
    }
    params$runs = ensureNestedList(runList)

    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

    return (fromJSON(response, simplifyVector=FALSE, simplifyDataFrame=FALSE))
}
