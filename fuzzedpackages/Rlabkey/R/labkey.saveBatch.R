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

labkey.saveBatch <- function(baseUrl=NULL, folderPath, assayName, resultDataFrame, batchPropertyList=NULL, runPropertyList=NULL)
{
    .Deprecated(msg = "This function is deprecated and will be removed in a future release. Please use the newer replacement function : labkey.experiment.saveBatch, as it supports additional features")

	baseUrl = labkey.getBaseUrl(baseUrl)

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(assayName)) stop (paste("A value must be specified for assayName."))
    if (missing(resultDataFrame)) stop (paste("A value must be specified for resultDataFrame."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

	## Translate assay name to an ID
	myurl <- paste(baseUrl,"assay",folderPath,"assayList.api", sep="")
    params <- list(name=assayName)
    assayInfoJSON <- labkey.post(myurl, toJSON(params, auto_unbox=TRUE))
	assayDef <- NULL
	assayInfo<- fromJSON(assayInfoJSON, simplifyVector=FALSE, simplifyDataFrame=FALSE)
	if (length(assayInfo) == 1 && length(assayInfo[[1]]) == 1)
	{
		assayDef <- assayInfo[[1]][[1]]
		if (assayDef$name != assayName)
			{assayDef <- NULL}
		## TODO:  check assay domain def against dataframe
	}
	if (is.null(assayDef))
		{stop(paste("Could not find an assay matching that name." ,sep=""))}

	# build Assay object tree based on R lists
	nrows <- nrow(resultDataFrame)
	ncols <- ncol(resultDataFrame)
	cnames <- colnames(resultDataFrame)
	rowsVector <- vector(mode="list", length=nrows)
	for(j in 1:nrows) {
		cvalues <- as.list(resultDataFrame[j,])
		names(cvalues) <- cnames
		rowsVector[[j]] <- cvalues
	}

	dataInputsArray <- vector(mode="list", length=0)

	runsArray <- vector(mode="list", length=1)
	runPropertyList <- c(runPropertyList, list("dataInputs" = dataInputsArray))

	runsArray[[1]] <- c(runPropertyList, list("dataRows" = rowsVector))

	batchPropertyList <- c(batchPropertyList, list("runs" = runsArray))

	baseAssayList <- list(assayId=assayDef$id)
	baseAssayList <- c(baseAssayList, list(batch=batchPropertyList))

	## Now post form with batch object filled out
	myurl <- paste(baseUrl, "assay", folderPath, "saveAssayBatch.api", sep="")
	pbody <- toJSON(baseAssayList, auto_unbox=TRUE)

	## Execute via our standard POST function
	mydata <- labkey.post(myurl, pbody)
	newAssayInfo <- fromJSON(mydata, simplifyVector=FALSE, simplifyDataFrame=FALSE)

	return(newAssayInfo)
}
                                                              
