create.vcf.bin  <- function(file_genotype=NULL,  availmemGb=8,  quiet=TRUE){
 ## an Rcpp function to create the no-space file of the genotype data M and Mt
 ## from vcf marker data. 
 ## Args
 ## file_genotype    absolute path and file name of genotype file
 ## availmemGb     available memory for conversion to reformatted file

 if(.Platform$OS.type == "unix") {
       binMfile <- paste(tempdir() , "/", "M.bin", sep="")
       binMtfile <- paste(tempdir()  , "/", "Mt.bin", sep="")
 } else {
       binMfile <- paste(tempdir() , "\\", "M.bin", sep="")
       binMtfile <- paste(tempdir() , "\\", "Mt.bin", sep="")
 }


    # need to create M file  
      liststr <- create_vcf_BIN_rcpp(f_name = file_genotype,  f_name_bin_M = binMfile,   f_name_bin_Mt = binMtfile,
               max_memory_in_Gbytes=availmemGb,  quiet = quiet, message=message )

 message("\n ReadVCF is complete ... \n")



 return(liststr)

}



