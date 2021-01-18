
#' @title Read phenotype file
#' @description Read in the phenotype data. 
#' @param filename contains the name of the phenotype  file. The file name needs to be in quotes. If the file is not in the working directory, then the full path 
#' to the file is required.
#' @param header a logical value. When \code{TRUE}, the first row of the file contains the names of the columns.  Default is \code{TRUE}.
#' @param csv   a logical value. When \code{TRUE}, a csv file format is assumed. When \code{FALSE}, a space separated format is assumed. Default
#'              is \code{FALSE}.
#' @param missing the number or character for a missing phenotype value.
#' @param ...   arguments to be passed to read.table such as \code{skip}, \code{sep}. See \code{\link{read.table}} so the list 
#'             of arguments. 




#' @details  
#' 
#' \code{ReadPheno} reads in the phenotype data 
#' which are data measured on traits and any fixed effects (or predictors/features/explanatory variables). 
#' A space separated plain text file is assumed. Each row in this file 
#' corresponds to an individual. The number of rows in the phenotype file must be the same as the number of rows in 
#' the marker data file. Also, the ordering of the individuals must be the same in the two files.  A space separated file with 
#' column headings is the default but can be changed with the \code{header} and \code{csv} options. 
#'
#' The phenotype file may contain multiple traits and fixed effects variables. 
#'
#' Missing values are allowed. Eagle is told which value should be treated as missing by setting the \code{missing} 
#' parameter to the value. 
#'
#' For example, suppose we have three individuals for which we have collected data on two quantitative traits (y1 and y2), and 
#' four explanatory variables (age, weight, height, and sex). The data looks like  
#' \tabular{cccccc}{
#'     y1      \tab y2      \tab  age  \tab weight \tab   height  \tab sex \cr
#'     112.02  \tab -3.123  \tab  26    \tab 75   \tab   168.5     \tab M \cr
#'     156.44  \tab 1.2     \tab  45    \tab 102   \tab   NA     \tab NA  \cr
#'     10.3    \tab NA   \tab  28     \tab 98   \tab   189.4     \tab F 
#'}
#' where the first row has the column headings and the next three rows contain the observed data on three 
#' individuals. 
#'
#' To load these data, we would use the command 
#'
#' \preformatted{pheno_obj <- ReadPheno(filename='pheno.dat', missing='NA')}
#'
#' where pheno.dat is the name of the phenotype file, and \code{pheno_obj} is the R object that contains the 
#' results from reading in the phenotype data.    The file is located in the working directory so there is no need to specify the full path, just the file name is suffice. 
#'
#'
#' \subsection{Dealing with missing trait data}{
#'
#'  \code{AM} deals automatically with individuals with missing trait data. 
#' These individuals are removed  from the analysis and a warning message is generated.
#' }
#' 
#' \subsection{Dealing with missing fixed effects values}{
#'
#' \code{AM} deals automatically with individuals with missing fixed effects values. 
#' These individuals are removed from the analysis and a warning message is generated
#' }
#'
#'
#' @seealso \code{\link{ReadMarker}} for reading in marker data, \code{\link{AM}} for performing association mapping.
#' @return 
#' a data frame is returned of the phenotype data. If \code{header} is true, the 
#' names of the columns will be as specified by the first row of the phenotype file. If \code{header} is \code{FALSE}, 
#' generic names are supplied by R in the form of V1, V2, etc.  If no column headings are given, these 
#' generic names will need to be used in the \code{trait} and \code{fformula} parameters in 
#' \code{\link{AM}}.  You can print out the column names of the data frame by using
#'
#' \preformatted{names(pheno_obj)}
#'
#' The column names are also printed along with other summary information when \code{ReadPheno} is run. 
#'
#'
#'
#' @examples
#' # Read in  phenotype data from ./extdata/
#' 
#' # find the full location of the phenotype data 
#' complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
#'
#' pheno_obj <- ReadPheno(filename=complete.name)
#'   
#'  ## print a couple of lines of the data file
#'  head(pheno_obj)
#'
ReadPheno <- function(filename = NULL, header=TRUE, csv=FALSE, missing = "NA" , ... ){

  phenofile <- fullpath(filename)



  error.code <- check.inputs(file_phenotype=phenofile)
  if(error.code){
     message(" ReadPheno has terminated with errors.")
     return(FALSE)
  }
  sep <- ""
  if(csv) sep=","
  message("\n\n Loading Phenotype file ... \n\n")

  # argument list - guarding against multiple uses of the same argument by mistake
  ar <- list(...)
  exist.missing  <- "na.strings" %in% names(ar) 
  exist.sep      <- "sep" %in% names(ar)

  if (exist.missing & exist.sep )
    args <- list(file=filename, header=header, ...)
  if (exist.missing & !exist.sep )
    args <- list(file=filename, header=header, sep=sep, ...)
  if (!exist.missing & exist.sep )
    args <- list(file=filename, header=header, na.strings=missing, ...)
  if (!exist.missing & !exist.sep )
    args <- list(file=filename, header=header, sep=sep, na.strings=missing, ...)

    phenos <- try(do.call(read.table, args))
if(class(phenos) == "try-error") {
        message(" ReadPheno as encountered errors in the reading in of the phenotype file")
        return(FALSE)
    }

  ## check for factors with only one level which will cause contrast code to crash
  for(ii in names(phenos)){
    if(is.factor(phenos[[ii]])){
       if(length(levels(phenos[[ii]]))==1){
          message(" The phenotype file contains factors that only have a single value. \n")
          message(" Please remove this factor. \n") 
          message(" ReadPheno has terminated with errors.")
          return(FALSE)

       }
    }
  }


message("               Summary of Phenotype File  \n")
message("              ~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
message(" File name:                   ",  phenofile, "\n")
message(" Number of individuals:       ", nrow(phenos), "\n")
message(" Number of columns:           ", ncol(phenos), "\n\n")
message(" First 5 rows of the phenotype file are \n")
if(nrow(phenos) > 5){
  mat <- as.matrix(phenos[1:5,])
  for(ii in 1:5){
  message(paste(as.character(mat[ii,]), collapse="   "))
  }
} else {
  mat <- as.matrix(phenos)
  for(ii in 1:nrow(phenos) ){
  message(paste(as.character(mat[ii,]), collapse="   "))
  }
}




message("\n Column classes are  \n")
for(ii in 1:ncol(phenos))
  message(c( sprintf("%20s   %15s", names(phenos)[ii], class(phenos[[ii]]) ), "\n"))



message("The phenotype file has been loaded.")


return(phenos)


}



