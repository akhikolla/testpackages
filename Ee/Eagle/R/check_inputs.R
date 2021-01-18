check.inputs <- function(ncpu=NULL, availmemGb=NULL,
                         file_genotype=NULL,
                         file_phenotype=NULL )
{
 ## internal function to AM

if(!is.null(ncpu)){
 if(!is.numeric(ncpu)){
   message("Error:  ncpu is not a numeric value. It is of class ", class(ncpu), "It should be the number of cpu.\n")
   return(TRUE)
 }
 if(ncpu < 1){
    message("Error: ncpu cannot be a zero or a negative number. It should be the number of cpu. \n")
    return(TRUE)
 }

}

if(!is.null(availmemGb))
{
 if(!is.numeric(availmemGb)){
   message("Error: availmemGb is not a numeric value. It is of class ", class(availmemGb), "It should be the number of gigabytes of RAM available. \n")
   return(TRUE)
 }
 if(availmemGb <= 0){
    message("Error: availmemGb cannot be zero or a negative number.  It should be the number of gigabytes of RAM available. \n")
    return(TRUE)
  }
}


if(!is.null(file_genotype))
{
  genofile <- fullpath(file_genotype)

  if(!file.exists(genofile)){
    message("Error: Cannot find marker file ", genofile, "\n")
    message("       This could be a problem with the name of the file and/or the location of the file. \n")
    message("       Perhaps specify the full name of the file (i.e. absolute directory path and file name) \n")
    message("       Type help(ReadMarker) and go to the examples section for an example of this. \n")
    return(TRUE)
  }
}


if(!is.null(file_phenotype))
{
  phenofile <- fullpath(file_phenotype)

  if(!file.exists(phenofile)){
    message("Error: Cannot find phenotype file ", phenofile, "\n")
    message("       This could be a problem with the name of the file and/or the location of the file. \n")
    message("       Perhaps specify the full name of the file (i.e. absolute directory path and file name) \n")
    message("       Type help(ReadPheno) and go to the examples section for an example of this. \n")
    return(TRUE)
  }
}

  return(FALSE)


}



