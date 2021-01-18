ReadVCF <- function( filename=NULL, availmemGb=16, quiet=TRUE ){


 if (nargs() == 0){
    ## checking that function has arguments
    message(" Please supply arguments to function \n")
    return(NULL)
 }

       
 if (is.null(filename)){
            message(" The name of the vcf file is missing.")
            message(" ReadVCF has terminated with errors.")
            return(NULL)
 }
 if (!file.exists(fullpath(filename) )){
            message(" The vcf file ", filename, " could not be found. ")
            message(" ReadVCF has terminated with errors ")
            return(NULL)
 }

 ## Rcpp function to create binary packed M and Mt file 
 #dim_of_M <- create.vcf.bin(file_genotype=fullpath(filename), availmemGb=availmemGb,  quiet=quiet  )
 liststr <- create.vcf.bin(file_genotype=fullpath(filename), availmemGb=availmemGb,  quiet=quiet  )


    if(.Platform$OS.type == "unix") {
       binfileM <- paste(tempdir(), "/", "M.bin", sep="")
       binfileMt <- paste(tempdir(), "/", "Mt.bin", sep="")
     } else {
       binfileM <- paste(tempdir()  , "\\", "M.bin", sep="")
       binfileMt <- paste(tempdir() , "\\", "Mt.bin", sep="")
     }


 geno <- list("tmpM"=binfileM, "tmpMt"=binfileMt,
               "dim_of_M" = liststr[["dim_of_M"]],
               "dim_of_Mt" = c( liststr[["dim_of_M"]][2], liststr[["dim_of_M"]][1]),
               "availmemGb" = availmemGb, 
               "map" = liststr[["map"]] )



  if(.Platform$OS.type == "unix") {
       RDatafile <- paste(tempdir() , "/", "M.RData", sep="")
  } else {
       RDatafile <- paste( tempdir() , "\\", "M.RData", sep="")
  }




  save(geno, file=RDatafile)
  ## create M.Rdata file in current directory
  return(geno)


}  ## end function call ReadVCF



