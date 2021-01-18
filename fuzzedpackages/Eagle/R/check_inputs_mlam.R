
check.inputs.mlam <- function (ncpu,  colname.trait, map, pheno,
                  geno, Zmat, lambda, falseposrate=NULL )
{
  ## internal function for AM

if(!is.data.frame(pheno)){
   message("Error: the phenotypic data needs to be stored as a data.frame. The command as.data.frame may help.     \n")
   return(TRUE)
}


if(!is.null(falseposrate)){
  if(!is.numeric(falseposrate)){
    message("Error: the falseposrate parameter must be a value between 0 and 1. \n")
    return(TRUE)
  } else {
   if(falseposrate < 0 | falseposrate > 1){
    message("Error: the falseposrate parameter must be a value between 0 and 1. \n")
    return(TRUE)
   }
 }
}
if(is.factor(pheno[[colname.trait]])){
   message("Error: the trait data in ", colname.trait, " is of class factor. Eagle can only handle quantitative traits.   \n")
   message("       Factors are assumed by R if the trait data contains alphanumeric values. \n")
   message("       To proceed, change the data in ", colname.trait, " in the input file so that the trait has only numeric values. \n")
   return(TRUE)
}

if(is.null(colname.trait)){
   message("Error: the name of the column containing the trait data must be given. \n")
   return(TRUE)
}

if(is.null(pheno)){
   message("Error: the pheno parameter has not been set. ")
   message("       Set this parameter to the object that contains ")
   message("       the phenotype data. This object is the result of running  ")
   message("       ReadPheno. \n")
   return(TRUE)
}

if(is.null(geno)){
   message("Error: the marker data has not been specified.")
   message("       If you are using the GUI, then go back and read in your marker data. ")
   message("       If you are using 'Eagle' functions, then set the geno parameter to the ")
   message("       output from running ReadMarker(). ")
   return(TRUE)
}


if(class(try(class(geno), silent=TRUE)) == "try-error"){
   message("Error: the object supplied to the geno parameter does not exist. ")
   message("       This object is set by running ReadMarker. Type help(ReadMarker) for help ")
   message("       on running this command. \n")
   return(TRUE)
}


if(class(try(class(pheno), silent=TRUE)) == "try-error"){
   message("Error: the object supplied to the pheno parameter does not exist. ")
   message("       This object is set by running ReadPheno. Type help(ReadPheno) for help. ")
   message("       on running this command. \n")
   return(TRUE)
}


## checking list structure of geno
if(!is.list(geno)){
  message("Error: the geno object is not a list object. ")
  message("     The geno object is obtained from running ReadMarker.Type help(ReadMarker) for help. \n")
  return(TRUE)
}

## checking if pheno is a data frame 
if(!is.data.frame(pheno)){
  message("Error: the pheno object is not a data frame. ")
  message("      It is a ", class(pheno), "\n")
  return(TRUE)
}




nms <- names(geno)
indx <- match(nms, c("tmpM", "tmpMt", "dim_of_M" , "dim_of_Mt", "availmemGb"))
if(any(is.na(indx))){
  message("Error: there is a problem with the list structure of the geno object. ")
  message("       It should contain the elements tmpM, tmpMt, and dim_of_M. ")
  message("       The object supplied contains the elements ", names(geno) )
  return(TRUE)
}

if(is.null(map)){
    message("WARNING: no map object has been specified. A generic map ")
    message("         will be assumed.                                ")
    map <- data.frame(Mrk= paste("M", 1:geno[["dim_of_M"]][2]), Chrm=1, Pos=1:geno[["dim_of_M"]][2])
}






 ## checks for colname.trait
 if(is.null(colname.trait)){
    message("Error: the column name for the trait/response has not been specified.")
    message("       Please set trait to the column name of the trait data in ")
    message("       the phenotype file. The allowable column names are ", names(pheno) )
    return(TRUE)
 }

 if(length(colname.trait)>1){
    message("Error: multiple column names for the trait have been specified. ")
    message("       Only a single column name should be  assigned to trait. ")
    return(TRUE)
 }

 indx <- match(colname.trait, names(pheno))
 if(any(is.na(indx))){
   message("Error: the trait column name does not match any of the column names in the phenotype file. ")
   message("       The name that has been supplied is ", colname.trait)
   message("       The column names of the phenotype file are ", names(pheno))
   return(TRUE)
 }





 ## when Zmat is null
 if(is.null(Zmat)){
    ## check that geno and pheno contain the same number of individuals
    if(nrow(pheno) !=  geno[["dim_of_M"]][1])
    {
      message("Error: the number of individuals specified in the phenotype file is ", nrow(pheno))
      message("       the number of individuals specified in the genotype file is ",  geno[["dim_of_M"]][1])
      message("       The number of individuals should be the same in the two files.")
      return(TRUE)
    }
 }

 ## make sure Zmat is of the dimension n > n_g
 if(!is.null(Zmat)){
   if(dim(Zmat)[1] <= dim(Zmat)[2]){
      message("Error: the dimensions of the Z matrix are incorrect. ")
      message("       the number of columns specified in the Z matrix file is ", ncol(Zmat))
      message("       the number of rows specified in the Z matrix file is ", nrow(Zmat))
      message("       the number of rows must be greater than the number of columns. ") 
      return(TRUE)
   }
 }


 ## when Zmat is not null
 if(!is.null(Zmat)){
   ## check that the number of columns in Zmat match the number of rows in geno
   if(ncol(Zmat) != geno[["dim_of_M"]][1])
   {
      message("Error: the number of columns specified in the Z matrix file is ", ncol(Zmat))
      message("       the number of rows specified in the genotype file is ",  geno[["dim_of_M"]][1])
      message("       The number of columns in the Z matrix should be the same as the number of rows in the genotype file.")
      return(TRUE)
   }
 }




 ## when Zmat is not null
 if(!is.null(Zmat)){
   ## check that the number of rows in Zmat match the number of rows in pheno
   if(nrow(Zmat) != nrow(pheno))
   {
      message("Error: the number of rows specified in the Z matrix file is ", nrow(Zmat))
      message("       the number of rows specified in the phenotype file is ",  nrow(pheno) )
      message("       The number of rows in the Z matrix file and phenotype file must be the same.")
      return(TRUE)
   }
 }




 ## check that map and geno contain the same number of snp
 if(nrow(map) != geno[["dim_of_M"]][2])
 {
   message("Error: the number of marker loci in the map file is ", nrow(map))
   message("       The number of marker loci in the genotype file is ", geno[["dim_of_M"]][2])
   message("       The number of marker loci in the two files should be the same." )
   return(TRUE)
 }


## check lambda is between 0 and 1 if specified. 
if(!is.null(lambda)){
  if(!is.numeric(lambda))
  { 
    message("Error: a non-numeric value has been set for the lambda parameter. ")
    message("       Valid values are between 0 and 1. ")
    return(TRUE)
  } else {

     if(lambda < 0)
     {
        message("Error: the lambda parameter cannot be negative. ")
        message("       Valid values are between 0 and 1. ")
        return(TRUE)
     } 
     if(lambda > 1 )
     {
        message("Error: the lambda parameter cannot be greater than 1. ")
        message("       Valid values are between 0 and 1. ")
        return(TRUE)
    }

 }  ## end if else

} ## end outer if



  return(FALSE)

}




