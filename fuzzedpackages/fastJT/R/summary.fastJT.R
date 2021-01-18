summary.fastJT <- function(object, Y2Print = 1:10, X2Print =1:10, printP = TRUE, outTopN=NA, subObj=FALSE,...)
{
	#-----------------------------------------------------------------------#
	# Check if the user set X2Print and Y2Print 
	SNP.is.default.for.top.hit <- .checkSNPDefault(X2Print)
	marker.is.default <- .checkMarkerDefault(Y2Print)
	# get the colnames 
	markerNames <- colnames(object$J)
	
	#-----------------------------------------------------------------------#
	# check the validaty of marker range
	Y2Print <- .markerRangeCheck(object, Y2Print)
	# check the validaty of SNP range
	X2Print <-.SNPRangeCheck(object,X2Print)
	
	#-----------------------------------------------------------------------#
	# compute the p-values
	if(attr(object,'standardized'))
		p.values <- 2*pnorm(-abs(object$J))
	else
		printP = FALSE
	#-----------------------------------------------------------------------#
	# print out note
	.note.title.print(object, outTopN,printP, marker.is.default, SNP.is.default.for.top.hit)
	
	#-----------------------------------------------------------------------#
	# check print out the summary
	#-----------------------------------------------------------------------#
	if(is.na(attr(object,'outTopN'))&&is.na(outTopN))
	{
		# Determine the SNPs that are required for out print 
		# for special case that the user provide the SNP names.
		if(is.character(X2Print))
		{
			# SNP name not in the range warning
			if(length(X2Print[!(X2Print %in% object$XIDs)]))
			{	
				cat("\nWARNING: specified X IDs/range: ")	
				cat(X2Print[!(X2Print %in% object$XIDs)]) 
				cat("\nare not found, results for them are omitted.\n")
			}
			# keep SNP that appears in the input data.
			X2Print <- intersect(X2Print, object$XIDs)
		}
		if(!attr(object,'standardized'))
			 print(object$J[X2Print,Y2Print])
		else{
			if(printP)
           		print(p.values[X2Print,Y2Print])
       		else
           		print(object$J[X2Print,Y2Print])
		}
		object$J <-object$J[X2Print,Y2Print]
		object$XIDs <- object$XIDs[X2Print]
		colnames(object$J) <-colnames(object$J)[Y2Print]
		if(subObj)	
			return(object)
	}
	
	#-----------------------------------------------------------------------#	
	if(is.na(attr(object,'outTopN'))&& !is.na(outTopN))		
	{
		if(is.numeric(Y2Print))
			marker.name.to.print <- markerNames[Y2Print]
		else
			marker.name.to.print <- Y2Print
		
		summary.mat <- .sortTopN(object,Y2Print, outTopN, printP)
		if(!attr(object,'standardized'))
			.printJStatisticsFormat(summary.mat, marker.name.to.print)	
		else{
			if(printP)
       	    	.printPvaluesFormat(summary.mat, marker.name.to.print)
       		else
           		.printJstarFormat(summary.mat, marker.name.to.print)
		}
		subobject <- .sort2saveObj(object, Y2Print,outTopN)
		colnames(subobject$J) <- marker.name.to.print
		colnames(subobject$XIDs) <- marker.name.to.print
		if(subObj)	
			return(subobject)
	}		
	#-----------------------------------------------------------------------#
    if(!is.na(attr(object,'outTopN')))	
	{
		if(SNP.is.default.for.top.hit)
			X2Print = seq(1:nrow(object$XIDs))
		
		if(is.numeric(Y2Print))
			marker.name.to.print <- markerNames[Y2Print]
		else
			marker.name.to.print <- Y2Print
		
		summary.mat <- .summary.mat(object, Y2Print, X2Print, printP)	
		
		if(!attr(object,'standardized'))
			.printJStatisticsFormat(summary.mat, marker.name.to.print)	
		else{
			if(printP)
       	    	.printPvaluesFormat(summary.mat, marker.name.to.print)
       		else
           		.printJstarFormat(summary.mat, marker.name.to.print)
		}
		object$J <- object$J[X2Print,Y2Print]
		object$XIDs <- object$XIDs[X2Print,Y2Print]	
		colnames(object$J) <-colnames(object$J)[Y2Print]
		colnames(object$XIDs) <-colnames(object$XIDs)[Y2Print] 
		if(subObj)
			return(object)
	}
}
#------------------------------------------------------------------------#
.checkSNPDefault <- function(X2Print)
{
	SNP.is.default.for.top.hit =T
	if(is.numeric(X2Print) && !(X2Print[length(X2Print)]==10&&length(X2Print)==10))
		SNP.is.default.for.top.hit = F
	if(length(X2Print)==1)
	{
		if(is.na(X2Print))
			SNP.is.default.for.top.hit =F
	}
  if(is.character(X2Print))
		SNP.is.default.for.top.hit = F
	return(SNP.is.default.for.top.hit)
}
#------------------------------------------------------------------------#
.checkMarkerDefault <- function(Y2Print)
{
	marker.is.default = T
	if(is.numeric(Y2Print) && !(Y2Print[length(Y2Print)]==10&&length(Y2Print)==10))
        marker.is.default = F
	if(is.character(Y2Print))
		marker.is.default = F
	return(marker.is.default)
}

#------------------------------------------------------------------------#
.sortTopN <- function(object, Y2Print, outTopN, printP)
{
 	markerNames <- colnames(object$J)
	col.name <- NULL
	summary.mat <-NULL
	nSNP <- nrow(object$J)
    SNPIDs <-rownames(object$J)
	pvalues <- 2*pnorm(-abs(object$J))
	
	for(i in Y2Print )
    {
        # generate colname names for the dataframe
        col.name <- cbind(col.name, "SNPID")
        if(printP)
            col.name <- cbind(col.name, "P-value")
        else{
            if(attr(object,'standardized'))
				col.name <- cbind(col.name, "J*")
			else
				col.name <- cbind(col.name, "J")
		}
        # create summary data frames    
        temp <- object$J[,i]
        temp <- cbind(temp,pvalues[,i])

        rownames(temp) <- SNPIDs
        colnames(temp) <- colnames(temp) <- c("J","pvalues")
        temp <-as.data.frame(temp)
        temp <- temp[order(-abs(temp$J)),]

        #rownames(order(object$J[,i]))
        summary.mat <- cbind(summary.mat, rownames(temp)[1:outTopN])
        if(printP)
            summary.mat <- cbind(summary.mat,temp$pvalues[1:outTopN])
        else
           summary.mat <- cbind(summary.mat,temp$J[1:outTopN])

    }
	colnames(summary.mat) <- col.name
    summary.mat = as.data.frame(summary.mat, stringsAsFactors = FALSE)
	return(summary.mat)
}

#-------------------------------------------------------------------------#
.sort2saveObj <- function(object,Y2Print,outTopN)
{
    markerNames <- colnames(object$J)
    col.name <- markerNames[Y2Print]
    XIDs <-NULL
    Js <-NULL

    for(i in Y2Print )
    {
        # create summary data frames    
        temp <- object$J[,i]
        temp <- temp[order(-abs(temp))]

        XIDs <- cbind(XIDs, names(temp)[1:outTopN])
        Js <- cbind(Js,temp[1:outTopN])

    }
    colnames(Js) <- col.name
    colnames(XIDs) <- col.name

    res <- NULL
    res$J <- as.matrix(Js)
    rownames(res$J) <- NULL
    res$XIDs <- as.matrix(XIDs)
    # set class for the result
    class(res) <- "fastJT"

    # add attribute for the result
    attr(res, 'outTopN') <- outTopN
    attr(res, 'standardized') <- attr(object, 'standardized')
    return(res)
}

#-------------------------------------------------------------------------#
# build summary matrix
.summary.mat <- function(object, Y2Print, X2Print, printP)
{
	markerNames <- colnames(object$J)
	p.values <- 2*pnorm(-abs(object$J))
	summary.mat <-NULL
	col.name <- NULL
	for(i in Y2Print )
	{
		# generate colname names for the dataframe
		col.name <- cbind(col.name, "SNPID")
		if(printP)
			col.name <- cbind(col.name, "P-value")
		else{
			if(attr(object,'standardized'))
				col.name <- cbind(col.name, "J*")
			else
				col.name <- cbind(col.name, "J")
		}

		# create summary data frames	
		summary.mat <- cbind(summary.mat,object$XIDs[X2Print,i])
		if(printP)
			summary.mat <- cbind(summary.mat,p.values[X2Print,i])
		else
			summary.mat <- cbind(summary.mat,object$J[X2Print,i])
	} 
   
	colnames(summary.mat) <- col.name
	summary.mat = as.data.frame(summary.mat, stringsAsFactors = FALSE)
	return(summary.mat)
}

## format printing function for p-values
.printPvaluesFormat <- function(object, markerNames)
{
	
	# process the input
    summary.mat <-object
	col.name <- colnames(summary.mat)

	# printing parameters
	# width: maximum number of charactor allow to print in each line default to be 100
	# n.print.repeat: number of table to be cut due to the width 
	
	n = floor(getOption("width")/20)
	width = n*20
    n.marker = ncol(summary.mat)/2
    n.print.repeat  <- ceiling(n.marker/n)
		
	# some basic element printint functions 
    marker.name.print <- function(x) cat(sprintf("%19s|",x,quote=F))
    space.print <- function(x) cat(" ")
    dash.print <- function(x) cat("-")
    double.dash.print <- function(x) cat("=")
    col.name.print <- function(i, col.name)
                    {
                        cat(sprintf("%11s", col.name[2*i-1]))
                        cat(sprintf("%8s", col.name[2*i]))
                        cat("|")
                    }
    SNPid.jstat.print <- function(j, i, summary.mat)
                    {
                        cat(sprintf("%11s", summary.mat[i, 2*j-1]))
                        cat(sprintf("%8.1e", as.numeric(summary.mat[i, 2*j])))
                        cat("|")
                    }

	# start of printing
    cat("\n\n")
    # Print out " ***Johckheere-Terpstra Test for Large Matrices****"
    title.string <- "Johckheere-Terpstra Test for Large Matrices\n"
    dummy <- sapply(1:(floor((width - nchar(title.string))/2)),space.print)
    cat(title.string)

    # Print out " *** Top Standardized Statistics and The Corresponding Variables IDs ***"
    subtitle.string <- "P-values for Top Standardized Statistics\n"
    dummy <- sapply(1:(ceiling((width - nchar(subtitle.string))/2)),space.print)
    cat(subtitle.string)

    # Print out "=======...="    
    dummy <- sapply(1:width,double.dash.print)
    cat("\n\n")
    for(n.repeat in 1:n.print.repeat)
    {
        start.marker = (n.repeat-1)*n + 1
        end.marker = min((n.repeat)*n, n.marker)

        # Print out "|  marker: x | marker: x| ...|"
    #    cat("|")
        dummy <- sapply(markerNames[start.marker:end.marker],marker.name.print)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|  ID   J* | SNPID   J*| ....|"
     #   cat("|")
        dummy <- sapply(start.marker:end.marker, col.name.print, col.name = col.name)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|"rnID"     "numeric value of J*" | ... |" 
        for(i in 1: nrow(summary.mat)){
     #       cat("|")
            dummy <- sapply(start.marker:end.marker, SNPid.jstat.print, i=i,summary.mat = summary.mat)
            cat("\n")
        }
        cat("\n\n")
    }
}

# print out statistics.

.printJStatisticsFormat <- function(summary.mat, markerNames)
{
	
	col.name <- colnames(summary.mat)
	
	# printing parameters
	# width: maximum number of charactor allow to print in each line default to be 100
	# n.print.repeat: number of table to be cutted due to the width 
	n = floor(getOption("width")/20)
	width = n*20
    n.marker = ncol(summary.mat)/2
    n.print.repeat  <- ceiling(n.marker/n)
   
	# some basic element printint functions  
	marker.name.print <- function(x) cat(sprintf("%19s|",x,quote=F))
    space.print <- function(x) cat(" ")
    dash.print <- function(x) cat("-")
    double.dash.print <- function(x) cat("=")
    col.name.print <- function(i, col.name)
                    {
                        cat(sprintf("%11s", col.name[2*i-1]))
                        cat(sprintf("%8s", col.name[2*i]))
                        cat("|")
                    }
    SNPid.jstat.print <- function(j, i, summary.mat)
                    {
                        cat(sprintf("%11s", summary.mat[i, 2*j-1]))
			            cat(sprintf("%8.1e", as.numeric(summary.mat[i, 2*j])))
						cat("|")
                    }

    cat("\n\n")
    # Print out " ***Johckheere-Terpstra Test for Large Matrices****"
    title.string <- "Johckheere-Terpstra Test for Large Matrices\n"
    dummy <- sapply(1:(floor((width - nchar(title.string))/2)),space.print)
    cat(title.string)

    # Print out " *** Top Standardized Statistics and The Corresponding Variables IDs ***"
    subtitle.string <- "Top Statistics\n"
    dummy <- sapply(1:(ceiling((width - nchar(subtitle.string))/2)),space.print)
    cat(subtitle.string)

    # Print out "=======...="    
    dummy <- sapply(1:width,double.dash.print)
    cat("\n\n")
    for(n.repeat in 1:n.print.repeat)
    {
        start.marker = (n.repeat-1)*n + 1
        end.marker = min((n.repeat)*n, n.marker)

        # Print out "|  marker: x | marker: x| ...|"
 #       cat("|")
        dummy <- sapply(markerNames[start.marker:end.marker],marker.name.print)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|  ID   J* | SNPID   J*| ....|"
  #      cat("|")
        dummy <- sapply(start.marker:end.marker, col.name.print, col.name = col.name)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|"rnID"     "numeric value of J*" | ... |" 
        for(i in 1:nrow(summary.mat)){
 #           cat("|")
            dummy <- sapply(start.marker:end.marker, SNPid.jstat.print, i=i,summary.mat = summary.mat)
            cat("\n")
        }
        cat("\n\n")
    }
}
# format printing function for standardized statistics 

.printJstarFormat <- function(object, markerNames)
{
	
	# process input
    summary.mat <-object
	col.name <- colnames(summary.mat)
	
	# printing parameters
	# width: maximum number of charactor allow to print in each line default to be 100
	# n.print.repeat: number of table to be cutted due to the width 
	n = floor(getOption("width")/20)
	width = n*20
    n.marker = ncol(summary.mat)/2
    n.print.repeat  <- ceiling(n.marker/n)
   
	# some basic element printint functions  
	marker.name.print <- function(x) cat(sprintf("%19s|",x,quote=F))
    space.print <- function(x) cat(" ")
    dash.print <- function(x) cat("-")
    double.dash.print <- function(x) cat("=")
    col.name.print <- function(i, col.name)
                    {
                        cat(sprintf("%11s", col.name[2*i-1]))
                        cat(sprintf("%8s", col.name[2*i]))
                        cat("|")
                    }
    SNPid.jstat.print <- function(j, i, summary.mat)
                    {
                        cat(sprintf("%11s", summary.mat[i, 2*j-1]))
			            cat(sprintf("%8.3f", as.numeric(summary.mat[i, 2*j])))
						cat("|")
                    }

    cat("\n\n")
    # Print out " ***Johckheere-Terpstra Test for Large Matrices****"
    title.string <- "Johckheere-Terpstra Test for Large Matrices\n"
    dummy <- sapply(1:(floor((width - nchar(title.string))/2)),space.print)
    cat(title.string)

    # Print out " *** Top Standardized Statistics and The Corresponding Variables IDs ***"
    subtitle.string <- "Top Standardized Statistics\n"
    dummy <- sapply(1:(ceiling((width - nchar(subtitle.string))/2)),space.print)
    cat(subtitle.string)

    # Print out "=======...="    
    dummy <- sapply(1:width,double.dash.print)
    cat("\n\n")
    for(n.repeat in 1:n.print.repeat)
    {
        start.marker = (n.repeat-1)*n + 1
        end.marker = min((n.repeat)*n, n.marker)

        # Print out "|  marker: x | marker: x| ...|"
 #       cat("|")
        dummy <- sapply(markerNames[start.marker:end.marker],marker.name.print)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|  ID   J* | SNPID   J*| ....|"
  #      cat("|")
        dummy <- sapply(start.marker:end.marker, col.name.print, col.name = col.name)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|"rnID"     "numeric value of J*" | ... |" 
        for(i in 1:nrow(summary.mat)){
   #         cat("|")
            dummy <- sapply(start.marker:end.marker, SNPid.jstat.print, i=i,summary.mat = summary.mat)
            cat("\n")
        }
        cat("\n\n")
    }
}


# function to check valid the requested Y2Print
.markerRangeCheck <- function(object,Y2Print)	
{
	# get the colnames and rownames 
	markerNames <- colnames(object$J)	
	
	# for the NA cases, which means full print out the results.
	if(is.na(Y2Print[1]))
		Y2Print = seq(1:ncol(object$J))

	# out of range
	if(is.character(Y2Print))
	{
		# if the user specify an wired names for the marker, we will tell
		# the user we ignore it.
		if(length(Y2Print[!(Y2Print %in% markerNames)]))
		{	
			cat("\nWARNING: specified Y ID/range: ")
			cat(Y2Print[!(Y2Print %in% markerNames)]) 
			cat(" are not found, \nresults for them are omitted.\n")
		}
		Y2Print <- intersect(Y2Print, markerNames)
	}
	
	if(is.numeric(Y2Print))
	{
		# if the user specify an wired numbers for the marker, we will tell
        # the user we ignore it.
		
        if(length(Y2Print[!(Y2Print %in% seq( 1: ncol(object$J)) )]))
        {   
			if( length(Y2Print) == 10&& Y2Print[length(Y2Print)] == 10)
			{}
			else
			{	
				cat("\nWARNING: specified marker ID/range: ")
            	cat(Y2Print[!(Y2Print %in% seq( 1: ncol(object$J)))])
            	cat(" are not found, \nresults for them are omitted.\n")
			}
        }
		Y2Print <- intersect(Y2Print, seq( 1: ncol(object$J)) )
	}
	return(Y2Print)	
}

# check the validaty of the requested SNP range
.SNPRangeCheck <- function(object, X2Print)
{	
    if(is.na(X2Print[1]))
        X2Print <- seq(1:nrow(object$J))

	if(is.numeric(X2Print))
    {
        # if the user specify an wired range for the SNP, we will tell
        # the user we ignore it.
        if(length(X2Print[!(X2Print %in% seq( 1: nrow(object$J)) )]))
        {
            cat("WARNING: specified X ID/range: ")
            cat(X2Print[!(X2Print %in% seq( 1: nrow(object$J)))])
            cat(" are not found, results for them are omitted.\n")
        }
        X2Print <- intersect(X2Print, seq(1: nrow(object$J)) )
    }
	return(X2Print)
}


# warning and tile print function

.note.title.print <- function(object, outTopN, printP, marker.is.default, SNP.is.default.for.top.hit)
{
	# print out note
	if(marker.is.default && SNP.is.default.for.top.hit)
	{
		if(is.na(attr(object,'outTopN'))&&is.na(outTopN))
		{
			cat("\n\nNote: By default, only part of the results may be printed!\n")
			cat("Please specify names/ranges using 'Y2Print' and 'X2Print',\n")
			cat("or set to 'NA' to print the full results.\n\n")
		}
		else
		{
			cat("\n\nNote: By default, all requested top hits are printed for up to 10 markers.\n")
			cat("Please specify names/ranges using 'Y2Print' and 'X2Print',\n")
			cat("or set to 'NA' to print the full results.\n\n")
			
		}
	}
	
	# print out titles
	if(is.na(attr(object,'outTopN'))&&is.na(outTopN))
	{           
		if(printP){
         	title.String <- "\n      Johckheere-Terpstra Test \n    P-values Based on Standardized Statistics\n\n"
     	}else{
			if(attr(object,'standardized'))
				title.String <- "\n    Johckheere-Terpstra Test Standardized Statistics\n\n"
			else
				title.String <- "\n    Johckheere-Terpstra Test Statistics\n\n"
     	}
     	cat(title.String)
	}

}




















