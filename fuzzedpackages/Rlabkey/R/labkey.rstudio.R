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


## initialize a RStudio session for LabKey user
##
labkey.rstudio.initSession <- function(requestId, baseUrl)
{
    ## check required parameters
    if(missing(requestId) || missing(baseUrl))
        stop (paste("A value must be specified for each of requestId and baseUrl."))

    url <- paste(baseUrl, "rstudio-fetchCmd.api?id=", requestId, sep="")
    response <- labkey.get(url)
    lkResult <- (fromJSON(response))
    if (lkResult$success == TRUE)
    {
        # do not print initialization cmd that contain apikey
        eval(parse(text=lkResult$initializationCmd))

        additionalCmds <- lkResult$additionalCmds
        if (length(additionalCmds) > 0)
        {
            for (i in 1:length(additionalCmds))
            {
                get(".rs.api.sendToConsole")(additionalCmds[i]);
            }
        }
    }
    else
    {
        warning("Unable to initialize integration with Labkey. Invalid requestId or baseUrl.")
    }
}

## initialize a RStudio session for LabKey user
##
labkey.rstudio.initRStudio <- function(apiKey="", baseUrl="", folderPath, skipViewer=FALSE)
{
    labkey.setDefaults(apiKey, baseUrl);

    if (missing(skipViewer) || skipViewer == FALSE)
    {
        if (!is.null(.lkdefaults[["baseUrl"]]))
        {
            folderPath <- encodeFolderPath(folderPath)
            get(".rs.api.viewer")(paste(.lkdefaults[["baseUrl"]], folderPath, "rstudio-viewer.view", sep=""))
        }
    }

    if (exists("labkey.rstudio.extend", mode="function")) get("labkey.rstudio.extend")()
}


## initialize a RStudio session for LabKey R report source editing
##
labkey.rstudio.initReport <- function(apiKey="", baseUrl="", folderPath, reportEntityId, skipViewer=FALSE, skipEdit=FALSE)
{
    labkey.rstudio.initRStudio(apiKey, baseUrl, folderPath, skipViewer);

    ## check required parameters
    if(missing(folderPath) || missing(reportEntityId))
        stop (paste("A value must be specified for each of folderPath and reportEntityId."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    url <- paste(baseUrl, "rstudio", folderPath, "getRReportContent.api", sep="")

    params <- list(entityId=reportEntityId)
    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

    result <- (fromJSON(response))

    if (result$success == TRUE)
    {
        ## reset working directory to home directory
        setwd('~/')

        ## create dir for report
        dir.create(file.path("LabKeyReports"), showWarnings = FALSE)
        dir.create(file.path("LabKeyReports", reportEntityId), showWarnings = FALSE)

        ## change working directory to report directory
        setwd(paste("LabKeyReports", reportEntityId, sep="/"))

        ## create props.JSON for folderPath, filename and timestamp
        labkey.rstudio.updateProp("folderPath", folderPath)
        labkey.rstudio.updateProp("reportFilename", result$filename)
        labkey.rstudio.updateProp("lastModified", result$lastModified)

        ## create prolog script and update its content
        prologFileConn <- file("prolog.R", open="w")
        writeLines(result$prolog, prologFileConn)
        close(prologFileConn)

        ## create report file and update its content
        fileConn <- file(result$filename, open="w")
        writeLines(result$reportSource, con=fileConn, sep="")
        close(fileConn)

        ## create input data file
        if (!is.null(result$queryName))
        {
            inputData <- labkey.selectRows(folderPath=folderPath, schemaName=result$schemaName, queryName=result$queryName, viewName=result$viewName, colNameOpt="rname", showHidden = TRUE, includeDisplayValues = TRUE)
            write.table(inputData, file="input_data.tsv", append=FALSE, sep="\t", quote=TRUE)
        }

        ## open report for editing
        if (missing(skipEdit) || skipEdit == FALSE)
            get("file.edit")(result$filename)
    }
    else
    {
        warning(result$errorMsg)
    }
}

## Update RStudio report source back to LabKey
##
labkey.rstudio.saveReport <- function(folderPath, reportEntityId, reportFilename, useWarning=FALSE)
{
    ## check required parameters
    if(missing(reportEntityId) || missing(reportFilename))
        stop (paste("A value must be specified for each of reportEntityId and reportFilename."))

    ## check working directory
    if (!grepl(reportEntityId, getwd()))
    {
        return("Working directory is currently not set to report's directory. Skip saving source to LabKey.")
    }

    if (!file.exists(reportFilename))
    {
        stop (paste("File doesn't exist: ", reportFilename))
    }

    targetFilename <- labkey.rstudio.getSavedProp("reportFilename")
    if (reportFilename != targetFilename)
    {
        return("Skip saving non LabKey report file.")
    }

    if(missing(folderPath))
    {
        folderPath <- labkey.rstudio.getSavedProp("folderPath")
        if (is.null(folderPath) || folderPath == "NULL")
        {
            return("Unable to determine report folderPath. Skip saving source to LabKey.")
        }
    }

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    baseUrl=labkey.getBaseUrl(NULL)

    ## check valid report
    url <- paste(baseUrl, "rstudio", folderPath, "ValidateRStudioReport.api", sep="")
    params <- list(entityId=reportEntityId)
    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))
    lkResult <- (fromJSON(response))

    if (lkResult$isValid == TRUE)
    {
        localLastModified <- labkey.rstudio.getSavedProp("lastModified")
        doSave <- TRUE
        if (localLastModified == lkResult$lastModified)
        {
            if (useWarning == TRUE)
                doSave <- get(".rs.api.showQuestion")("Save to LabKey?", "Do you want to update report content to LabKey Server?")
        }
        else
        {
            doSave <- get(".rs.api.showQuestion")("Save to LabKey? (Potential conflicting edit)", "The report source was modified in LabKey Server and the content might have diverged from local copy. Do you want to save local changes to LabKey Server?")
        }
        if (!doSave)
        {
            return("Skipped saving updated source to LabKey Server");
        }
        url <- paste(baseUrl, "rstudio", folderPath, "SaveRReportContent.api", sep="")

        script <- readChar(reportFilename, file.info(reportFilename)$size)

        params <- list(entityId=reportEntityId, runScript=script)
        response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

        result <- (fromJSON(response))

        if (result$success == TRUE)
        {
            labkey.rstudio.updateProp("lastModified", result$lastModified)
            return("Successfully updated report source to LabKey Server.")
        }
        else
            warning(result$errorMsg)
    }
    else
    {
        warning("Failed to update report source to LabKey Server. Report doesn't exist.")
    }
}

## Read property value form props.JSON
##
labkey.rstudio.getSavedProp <- function(propName)
{
    propsFilepath <- "props.JSON"
    if (!file.exists(propsFilepath))
        return (NULL)
    props <- fromJSON(txt = propsFilepath)
    return (props[[propName]])
}

## Update property value to props.JSON
##
labkey.rstudio.updateProp <- function(propName, propValue)
{
    propsFilepath <- "props.JSON"
    if (!file.exists(propsFilepath))
        props <- list()
    else
        props <- fromJSON(txt = propsFilepath)
    props[propName] = propValue
    write(toJSON(props), file=propsFilepath)
}

## check valid rlabkey session
##
labkey.rstudio.isInitialized <- function()
{
    return (!is.null(.lkdefaults[["baseUrl"]]))
}
