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

## Returns the domain design (as a dataframe) for the specified domain
##
labkey.domain.get <- function(baseUrl=NULL, folderPath, schemaName, queryName)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## check required parameters
    if(missing(baseUrl) || is.null(baseUrl) || missing(folderPath) || missing(schemaName) || missing(queryName))
        stop (paste("A value must be specified for each of baseUrl, folderPath, schemaName and queryName."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    url <- paste(baseUrl, "property", folderPath, "getDomain.api", sep="")

    params <- list(schemaName=schemaName, queryName=queryName)
    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

    result <- fromJSON(response)

    # Issue 36894: translate NULL values for string props to NA
    if (is.null(result$description)) result$description = NA
    if (is.null(result$schemaName)) result$schemaName = NA
    if (is.null(result$queryName)) result$queryName = NA
    if (is.null(result$templateDescription)) result$templateDescription = NA
    if (is.null(result$instructions)) result$instructions = NA
    if (is.null(result$domainKindName)) result$domainKindName = NA

    return (result)
}

## Update an existing domain
##
labkey.domain.save <- function(baseUrl=NULL, folderPath, schemaName, queryName, domainDesign)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## check required parameters
    if (missing(baseUrl) || is.null(baseUrl) || missing(folderPath) || missing(schemaName) || missing(queryName) || missing(domainDesign))
        stop (paste("A value must be specified for each of baseUrl, folderPath, schemaName, queryName and domainDesign."))

    if (!is.list(domainDesign))
        stop (paste("domainDesign must be a list data structure."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)
    params <- list(schemaName = schemaName, queryName = queryName, domainDesign = domainDesign)

    url <- paste(baseUrl, "property", folderPath, "saveDomain.api", sep="")
    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

    return (fromJSON(response))
}

## Helper function to create the domain design list
##
labkey.domain.createDesign <- function(name, description = NULL, fields, indices = NULL)
{
    ## check required parameters
    if (missing(name) || missing(fields))
        stop (paste("A value must be specified for each of name and fields."))

    if (!is.list(fields))
        stop (paste("The 'fields' parameter must be a list of field definitions."))

    dd <- list(name = name, fields = fields$fields)

    if (!missing(description))
        dd$description = description

    if (!missing(indices)) {
        if (!is.list(indices))
            stop (paste("The 'indices' parameter must be a list of index definitions (including a list of 'columnNames' and a 'unique' boolean)."))

        dd$indices = indices
    }

    return (dd)
}

## Helper function to create the domain design indices list
##
labkey.domain.createIndices <- function(colNames, asUnique, existingIndices = NULL)
{
    if (!is.list(colNames))
        stop (paste("The 'colNames' parameter must be a list of column names."))

    if (!is.logical(asUnique))
        stop (paste("The 'asUnique' parameter must be a logical value of either TRUE or FALSE."))

    columnNames <- list(colNames)
    unique <- list(tolower(toString(asUnique)))
    indices = as.data.frame(cbind(columnNames, unique))

    if (!missing(existingIndices))
        indices = rbind(existingIndices, indices)

    return (indices)
}

labkey.domain.create <- function(baseUrl=NULL, folderPath, domainKind=NULL, domainDesign=NULL, options=NULL,
        module=NULL, domainGroup=NULL, domainTemplate=NULL, createDomain=TRUE, importData=TRUE)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## check required parameters
    if (missing(baseUrl) || is.null(baseUrl) || missing(folderPath))
        stop (paste("A value must be specified for each of baseUrl and folderPath."))

    if (is.null(domainKind) && is.null(domainTemplate))
        stop (paste("Domain creation must use either a domain kind or a domain template."))

    if (!is.null(domainKind))
    {
        if (is.null(domainDesign))
            stop (paste("If domainKind is specified, then domainDesign must also be included."))

        if (!is.list(domainDesign))
            stop (paste("domainDesign must be a list data structure."))

        params <- list(kind = domainKind, domainDesign = domainDesign)
        if (!missing(options))
        {
            if (!is.list(options))
                stop (paste("options must be a list data structure."))
            params$options = options
        }
    }

    if (!is.null(domainTemplate))
    {
        if (is.null(module) || is.null(domainGroup))
            stop (paste("If domainTemplate is specified, module and domainGroup are required."))

        params <- list(domainTemplate = domainTemplate, module = module, domainGroup = domainGroup,
                createDomain = createDomain, importData = importData)
    }

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    url <- paste(baseUrl, "property", folderPath, "createDomain.api", sep="")
    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

    return (fromJSON(response))
}

labkey.domain.drop <- function(baseUrl=NULL, folderPath, schemaName, queryName)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## Validate required parameters
    if (missing(folderPath)) stop (paste("A value must be specified for folderPath."))
    if (missing(schemaName)) stop (paste("A value must be specified for schemaName."))
    if (missing(queryName)) stop (paste("A value must be specified for queryName."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    url <- paste(baseUrl, "property", folderPath, "deleteDomain.api", sep="")

    params <- list(schemaName=schemaName, queryName=queryName)
    response <- labkey.post(url, toJSON(params, auto_unbox=TRUE))

    return (fromJSON(response))
}

labkey.domain.inferFields <- function(baseUrl=NULL, folderPath, df)
{
    baseUrl=labkey.getBaseUrl(baseUrl)

    ## check required parameters
    if (missing(baseUrl) || is.null(baseUrl) || missing(folderPath) || missing(df))
        stop (paste("A value must be specified for each of baseUrl, folderPath and df."))

    ## normalize the folder path
    folderPath <- encodeFolderPath(folderPath)

    ## write the dataframe to a tempfile to post to the server
    tf <- tempfile(fileext=".tsv")
    write.table(df, file=tf, sep="\t", quote=FALSE, row.names=FALSE)

    ## Execute via our standard POST function
    url <- paste(baseUrl, "property", folderPath, "inferDomain.api", sep="")

    rawdata <- labkey.post(url, list(file=upload_file(tf)), encoding="multipart")
    ## delete the temp file
    file.remove(tf)
    response <- fromJSON(rawdata)

    return (response)
}

labkey.domain.createAndLoad <- function(baseUrl=NULL, folderPath, name, description="", df, domainKind, options=NULL, schemaName=NULL)
{
    ## check required parameters
    if (missing(baseUrl) || is.null(baseUrl) || missing(folderPath) || missing(name) || missing(df) || missing(domainKind))
        stop (paste("A value must be specified for each of baseUrl, folderPath, name, df or domainKind."))

    if (is.null(options))
        options <- list(strictFieldValidation = FALSE)
    else
        options <- c(options, list(strictFieldValidation = FALSE))

    if (is.null(schemaName))
    {
        if (domainKind == "StudyDatasetVisit" || domainKind == "StudyDatatsetDate")
            schemaName <- "study"
        else if (domainKind == "IntList" || domainKind == "VarList")
            schemaName <- "lists"
        else if (domainKind == "IssueDefinition")
            schemaName <- "issues"
        else if (domainKind == "SampleSet")
            schemaName <- "samples"
        else if (domainKind == "DataClass")
            schemaName <- "exp.data"
    }

    if (is.null(schemaName))
        stop (paste("A value must be specified for schemaName."))

    fields = labkey.domain.inferFields(baseUrl = baseUrl, folderPath = folderPath, df = df[,colnames(df)])

    design = labkey.domain.createDesign( fields = fields, name = name, description = description)
    labkey.domain.create(baseUrl = baseUrl, folderPath = folderPath, domainKind = domainKind,
        domainDesign = design, options = options)

    labkey.insertRows(baseUrl = baseUrl, folderPath = folderPath, schemaName = schemaName, queryName= name, df)
}

labkey.domain.createConditionalFormat <- function(queryFilter, bold=FALSE, italic=FALSE, strikeThrough=FALSE, textColor="", backgroundColor="")
{
    data.frame(filter = queryFilter, bold = bold, strikethrough = strikeThrough, italic = italic, textColor = textColor, backgroundColor = backgroundColor)
}

labkey.domain.createConditionalFormatQueryFilter <- function(filterType, value, additionalFilter=NULL, additionalValue=NULL)
{
    qf1 <- makeFilter(c("format.column", filterType, value))[1]
    qf2 <- NULL

    if (!is.null(additionalValue) || !is.null(additionalFilter))
        qf2 <- makeFilter(c("format.column", additionalFilter, additionalValue))[1]

    qf <- if(is.null(qf2)) qf1 else paste0(qf1, "&", qf2)
    return(qf)
}

labkey.domain.FILTER_TYPES <-
  list(
        HAS_ANY_VALUE = '',

        EQUAL = 'eq',
        DATE_EQUAL = 'dateeq',

        NEQ = 'neq',
        NOT_EQUAL = 'neq',
        DATE_NOT_EQUAL = 'dateneq',

        NEQ_OR_NULL = 'neqornull',
        NOT_EQUAL_OR_MISSING = 'neqornull',

        GT = 'gt',
        GREATER_THAN = 'gt',
        DATE_GREATER_THAN = 'dategt',

        LT = 'lt',
        LESS_THAN = 'lt',
        DATE_LESS_THAN = 'datelt',

        GTE = 'gte',
        GREATER_THAN_OR_EQUAL = 'gte',
        DATE_GREATER_THAN_OR_EQUAL = 'dategte',

        LTE = 'lte',
        LESS_THAN_OR_EQUAL = 'lte',
        DATE_LESS_THAN_OR_EQUAL = 'datelte',

        STARTS_WITH = 'startswith',
        DOES_NOT_START_WITH = 'doesnotstartwith',

        CONTAINS = 'contains',
        DOES_NOT_CONTAIN = 'doesnotcontain',

        CONTAINS_ONE_OF = 'containsoneof',
        CONTAINS_NONE_OF = 'containsnoneof',

        IN = 'in',

        EQUALS_ONE_OF = 'in',

        NOT_IN = 'notin',
        EQUALS_NONE_OF = 'notin',

        BETWEEN = 'between',
        NOT_BETWEEN = 'notbetween',

        IS_BLANK = 'isblank',
        IS_NOT_BLANK = 'isnonblank',

        MEMBER_OF = 'memberof'
    )
