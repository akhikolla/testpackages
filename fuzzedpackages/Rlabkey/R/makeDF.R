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

makeDF <- function(rawdata, colSelect=NULL, showHidden, colNameOpt)
{
    decode <- fromJSON(rawdata, simplifyVector=FALSE, simplifyDataFrame=FALSE)

	## Check to make sure at least one column exists in the result set
	if(length(decode$columnModel)==0){
		stop('No columns exist in the result set. Be sure you are using the column name for colNameOpt="fieldname" and the column label for colNameOpt="caption" in the colSelect vector. See the documentation for more details.')
	}

	colModelNames = c()
	colModelLabels = c()
	colModelRNames = c()
	for(i in 1:length(decode$columnModel)){
		if(!is.null(decode$columnModel[[i]]$dataIndex)){
			colModelNames = c(colModelNames, decode$columnModel[[i]]$dataIndex)
			colModelRNames = c(colModelRNames, .getRNameFromName(decode$columnModel[[i]]$dataIndex, existing=colModelRNames))
		}
		if(!is.null(decode$columnModel[[i]]$header)){
			colModelLabels = c(colModelLabels, decode$columnModel[[i]]$header)
		}
	}

	## Check for invalid colSelect name
	colSelectVector <- NULL
	if(!is.null(colSelect)){
		colSelectVector <- strsplit(colSelect, ",")[[1]]
	}
	if(!is.null(colSelectVector) & length(colSelectVector)>0){
		for(i in 1:length(colSelectVector)){
			if(colSelectVector[[i]] != "*" & !(colSelectVector[[i]] %in% colModelNames) & !(colSelectVector[[i]] %in% colModelLabels) & !(colSelectVector[[i]] %in% colModelRNames)){
				stop(paste('The column "',colSelectVector[[i]],'" specified in the colSelect variable does not exist in the result set. Be sure you are using the column name for colNameOpt="fieldname" and the column label for colNameOpt="caption". See the documentation for more details.',sep=''))
			}
		}
    }

  	## Get column names in proper order, associated header index, hidden tag, and data type
  	cnames <- NULL
  	hindex <- NULL
  	hide <- NULL
  	for(j in 1:length(decode$columnModel))
  	{   	
		## three different ways to refer to columns exist today
		## selectRows and executeSQL by default return the field caption (also called the  "label")
		## When sepcifying colSelect, colFilter, or ExecuteSql, you must use field names.
		## when running R "views"  at the server, the field_name is modified to use underscores and lower cased.  
		## This also makes it a legal name in R, which can be useful.

  		if (colNameOpt == "caption") {
			cname <- decode$columnModel[[j]]$header
		} else if (colNameOpt == "fieldname") {
			cname <- decode$columnModel[[j]]$dataIndex
		} else if (colNameOpt == "rname" ) {
			cname <- .getRNameFromName(decode$columnModel[[j]]$dataIndex, existing=cnames)
		} else {
			stop("Invalid colNameOpt option.  Valid values are caption, fieldname, and rname.")
		}
  			
		cnames <- c(cnames, cname)
        hindex <- c(hindex, decode$columnModel[[j]]$dataIndex)

        colHidden <- decode$columnModel[[j]]$hidden
        # issue 26561: include column if it is specifically included in the colSelect vector
        if(!is.null(colHidden)) {
            if(colHidden & !is.null(colSelectVector)) {
               colHidden <- !(decode$columnModel[[j]]$header %in% colSelectVector | decode$columnModel[[j]]$dataIndex %in% colSelectVector)
            }
        }
        hide <- c(hide, colHidden)
    }
  	refdf <- data.frame(cnames,hindex,hide)

	## Check for no rows returned, put data in data frame 
  	if(length(decode$rows)<1){
		tohide <- length(which(refdf$hide==TRUE))
		totalcol <- length(refdf$cnames)
		if(showHidden==FALSE){
			emptydf <- as.data.frame(rep(list(num=double(0)), each=(totalcol-tohide)))
			colnames(emptydf) <- refdf$cnames[refdf$hide==FALSE]
			warning("Empty data frame was returned. Query may be too restrictive.", call.=FALSE)
			return(emptydf)
		} else {
			emptydf <- as.data.frame(rep(list(num=double(0)), each=(totalcol)))
			colnames(emptydf) <- refdf$cnames
			warning("Empty data frame was returned. Query may be too restrictive.", call.=FALSE)
		}
		return(emptydf)
	}
		
	if(length(decode$rows)>0) {
		hold.dat <- NULL
		tmprow <- filterrow(decode$rows[[1]])
		hold.dat<-listToMatrix(decode$rows, names(tmprow))
        hold.dat <- as.data.frame(hold.dat,stringsAsFactors=FALSE)
		names(hold.dat) <- names(tmprow)   			
	}

	## Order data
	oindex <- NULL
	## number of cols selected may be more or less than described in metadata
  	for(k in 1:length(cnames)){oindex <- rbind(oindex, which(names(hold.dat)==refdf$hindex[k]))}

  	refdf$oindex <- oindex
	newdat <- as.data.frame(hold.dat[,refdf$oindex],stringsAsFactors=FALSE)

  	refdf$type <- NULL
  	for(p in 1:dim(refdf)[1]) {   
  	    ind <- which(refdf$hindex==decode$metaData$fields[[p]]$name)
  	    refdf$type[ind] <- decode$metaData$fields[[p]]$type
  	}

  	## Delete hidden column(s) unless showHidden=TRUE
	if (showHidden==TRUE || is.null(decode$metaData$id)) {}  else {
		hide.ind <- which(refdf$hide==TRUE)
		if(length(hide.ind)>0 ) {
			if(length(hide.ind) == length(newdat)) {
				stop("No visible columns selected.  Use the showHidden=TRUE to see these columns.") 
			} else {
				newdat <- newdat[,-hide.ind]
				refdf <- refdf[-hide.ind,]
				cnames <- cnames[-hide.ind] 
			}
		}
	}

	## Set mode for multiple columns of data (this also removes list factor)
	if(is.null(dim(newdat))==FALSE)
  	{
  	    for(j in 1:ncol(newdat))
  	    {
  	        mod <- refdf$type[j]
  	        try(
                if(mod=="date") { newdat[,j] <- .parseDate(newdat[,j])} else
                if(mod=="string"){	suppressWarnings(mode(newdat[,j]) <- "character")} else
                if(mod=="int"){ suppressWarnings(mode(newdat[,j]) <- "integer")} else
                if(mod=="boolean"){suppressWarnings(mode(newdat[,j]) <- "logical")} else
                if(mod=="float"){suppressWarnings(mode(newdat[,j]) <- "numeric")} else
                {print("MetaData field type not recognized.")}
            , silent=TRUE)
        }
	    newdat <- as.data.frame(newdat, stringsAsFactors=FALSE); colnames(newdat)<-cnames
    }

	## Set mode for single column of data
	if(is.null(dim(newdat))==TRUE & length(newdat)>1)
	{
	    mod <- refdf$type
	    try(
            if(mod=="date"){ newdat <- .parseDate(newdat)}else
            if(mod=="string"){suppressWarnings(mode(newdat) <- "character")} else
            if(mod=="int"){ suppressWarnings(mode(newdat) <- "integer")} else
            if(mod=="boolean"){suppressWarnings(mode(newdat) <- "logical")} else
            if(mod=="float"){suppressWarnings(mode(newdat) <- "numeric")} else
            {print("MetaData field type not recognized.")}
        , silent=TRUE)
	    newdat <- as.data.frame(newdat, stringsAsFactors=FALSE); colnames(newdat)<-cnames[1]
    }

return(newdat)
}

## need to get rid of hidden hrefs within the row, R doesn't use them and their presence causes problems
## also consolidate null handling here
filterrow<-function(row)
{
	filtered <- NULL
	for (x in 1:length(row)) {		
		valname <- names(row[x])
		if ((nchar(valname)>11) && (substr(valname,1,11) == as.character("_labkeyurl_"))) {
			next
		}
		if (is.null(row[x][[valname]])) { row[x][[valname]]<-NA }
		filtered <- c(filtered, row[x])
	}
return(filtered)

}

.getRNameFromName <- function(lkname, existing=NULL)
{
	rname <- gsub("::", "_", lkname)
	rname <- tolower(chartr(" /", "__", rname))

	if (length(existing)>0)
	{ 
		for (i in 1:99)
		{
			if(length(existing[rname == existing]) ==0)
				{break;}
			else 
				{rname<- c(rname + as.character(i))}
        } 
  	}    	
  	return (rname)
}

.parseDate <- function(s)
{
    s <- as.character(s);

    ## format from DateUtil.getJsonDateTimeFormatString ("yyyy/MM/dd HH:mm:ss")
    d <- tryCatch(as.POSIXct(s),error = function(e) NA);

    ## if none of the values have time part, convert the variable to a date
    t <- format(d, "%H-%M-%S")
    if (all(t == "00-00-00", na.rm = TRUE)) {
        d <- as.Date(d, tz = attr (d, which = "tzone"));
    }

    return(d);
}

# converts the factors of a data frame to characters
convertFactorsToStrings <- function(df)
{
    factors <- sapply(df, is.factor)
    df[factors] <- lapply(df[factors], as.character)
    return(df);
}