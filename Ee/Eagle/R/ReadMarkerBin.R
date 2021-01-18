#' @title Read  marker data.
#' 
#' @description
#' A function for reading in different types of snp marker data. 
#' @param filename contains the name of the marker  file. The file name needs to be in quotes. If the file is not in the working directory, then the full path 
#' to the file is required.
#' @param type  specify the type of file. Choices are 'text' (the default) , PLINK, and vcf.
#' @param missing the number or character for a missing genotype in the text file. There is no need to specify this for a vcf or PLINK ped file. Missing 
#' allele values in a vcf file are coded as "." and missing 
#' allele values in a PLINK file must be coded as '0' or '-'.  
#' @param AA     the character or number corresponding to the 'AA' snp genotype in the marker genotype file. 
#' This need only be specified if the file type is 'text'.  If a character then it must be in quotes.
#' @param AB     the  character or number  corresponding to the 'AB' snp genotype in the marker genotype file. 
#' This need only be specified if the file type is 'text'.
#'This can be left unspecified 
#'               if there are no heterozygous genotypes (i.e. the individuals are inbred). Only a single 
#'  heterozygous genotype is allowed ('Eagle' does not distinguish between 'AB' and 'BA').  
#' If specified and a character, it must be in quotes. 
#' @param BB        the character or number corresponding to the 'BB' snp genotype in the marker genotype file. 
#' This need only be  specified if the file type is 'text'.  If a character, then it must be in quotes.
#' @param availmemGb a numeric value. It specifies the amount of available memory (in Gigabytes). 
#'         This should be set to be as large as possible for best performance.   
#' @param  quiet      a logical value. If set to \code{TRUE}, additional runtime output is printed. 
#' @details
#' 
#' \code{ReadMarker} can handle three different types of marker data; namely,
#' genotype data in a plain text file, PLINK ped files, and vcf files.  
#'  
#' \subsection{\strong{Reading in a plain text file containing the marker genotypes}}{
#' To load a text file that contains snp genotypes, run \code{ReadMarker} with \code{filename} set to the name of the file, 
#' and \code{AA}, \code{AB}, \code{BB} set to the corresponding genotype values. 
#' The genotype values in the text file can be numeric, character, or  a mix of both.
#'
#' We make the following assumptions
#' \itemize{
#' \item{The text file does not contain row or column headings}
#' \item{The file is allowed to contain missing genotypes that have been coded according to \code{missing}}
#' \item{Individuals are diploid}
#' \item{The rows of the text file are the individuals and the columns are the marker loci}
#' \item{The file is space separated}
#' \item{The mapping of the observed genotypes in the marker file to \code{AA}, \code{AB}, and \code{BB}, remains the same for all loci}
#' \item{Individuals are outbred when \code{AA}, \code{AB}, and \code{BB} are specified and 
#'  inbred when only \code{AA}, and \code{BB} are specified}
#' \item{For a text file, the same alphanumeric value is used for all missing marker genotypes. For a PLINK ped file, the missing allele is allowed to 
#' be '0' or '-'.} 
#'}
#'
#' For example, suppose we have a space separated text file with marker genotype data collected from five snp loci on three individuals
#' where the snp genotype AA has been coded 0, the snp genotype AB has been coded 1, the snp genotype BB has been coded 2, and missing genotypes are coded 
#' as 99
#' \tabular{ccccc}{
#'  0 \tab  1  \tab  2 \tab  0\tab   2 \cr
#'  1 \tab  1  \tab  0 \tab  2 \tab  0 \cr
#'  2 \tab  2  \tab  1 \tab  1 \tab  99
#'}
#' The file is called geno.txt and is located in the directory /my/dir/. 
#'
#'  
#' To load these data, we would use the command
#'
#' \preformatted{geno_obj <- ReadMarker(filename='/my/dir/geno.txt', AA=0, AB=1, BB=2, type='text', missing=99)}
#'
#' where the results from running the function are placed in \code{geno_obj}.
#'
#' As another example, suppose we have a space separated text file with marker genotype data collected from five snp loci on three individuals 
#' where the snp genotype AA has been coded a/a, the snp genotype AB has been coded a/b, and the snp genotype BB has been coded b/b
#' \tabular{c}{
#'  a/a a/b b/b a/a b/b \cr
#'  a/b a/b a/a b/b a/a \cr
#'  b/b b/b a/b a/b NA 
#'}
#' The file is called geno.txt and is located in the same directory from which R is being run (i.e. the working directory). 
#'
#' To load these data, we would use the command 
#'
#' \preformatted{geno_obj <- ReadMarker(filename='geno.txt', AA='a/a', AB='a/b', BB='b/b', 
#'                                        type='text', missing = 'NA')}
#'
#' where the results from running the function are placed in \code{geno_obj}.
#'}
#'
#' \subsection{\strong{Reading in a PLINK ped file}}{
#' PLINK is a well known toolkit for the analysis of genome-wide association data. See  
#' \url{https://www.cog-genomics.org/plink2}
#' for details. 
#'
#' Full details of PLINK ped files can be found \url{https://www.cog-genomics.org/plink/1.9/formats#ped}. Briefly, 
#' the PED file is a space delimited file (tabs are not allowed): the first six columns are mandatory:
#'
#'
#' \tabular{l}{
#'     Family ID  \cr
#'     Individual ID \cr
#'     Paternal ID \cr
#'     Maternal ID \cr
#'     Sex (1=male; 2=female; other=unknown) \cr
#'     Phenotype
#'}
#'
#'
#' Here, these columns can be any values since \code{ReadMarker} ignores these columns.  
#'
#' Genotypes (column 7 onwards) can be any character 
#' (e.g. 1,2,3,4 or A,C,G,T or anything else) except 0 which is, by default, the missing genotype character. All markers should be biallelic. 
#' All snps must have two alleles specified.  Missing alleles (i.e 0 or -) are allowed.
#' No column headings should be given. 
#' 
#' As an example, suppose we have data on three individuals  genotyped for four snp loci 
#' \tabular{cccccccccccccc}{
#'     FAM001 \tab 101  \tab  0    \tab 0   \tab   1  \tab 0 \tab  A  \tab  G  \tab C \tab C \tab C \tab G \tab A \tab A \cr
#'     FAM001 \tab 201  \tab  0    \tab 0   \tab   2  \tab 0 \tab  A  \tab  A  \tab C \tab T \tab G \tab G \tab T \tab A \cr 
#'     FAM001 \tab 300  \tab  101  \tab 201 \tab   2  \tab 0 \tab  G  \tab  A  \tab T \tab T \tab C \tab G \tab A \tab T 
#'}
#'
#' Then to load these data, we would use the command 
#'
#' \preformatted{geno_obj <- ReadMarker(filename='PLINK.ped', type='PLINK')}
#'
#' where \code{geno_obj} is used by \code{\link{AM}}, and the file PLINK.ped is located in the working directory (i.e. the directory from which R 
#' is being run). 
#'}
#'
#' \subsection{\strong{Reading in a vcf file}}{
#'  VCF is a tab separated text file containing  meta-information lines, a header line, and data 
#' lines. The data lines  contain information about a position in the genome.
#'
#' It is assumed that genotype information has been recorded on samples for each position. 
#'
#' Loci with more than two alleles will be removed automatically. 
#' 
#' Eagle will only accept a single (uncompressed) vcf file. If chromosomal information has been recorded in separate 
#' vcf files, these files need to be merged into a single vcf file.  This can be done by using the 
#'  BCFtools utility set with command line "bcftools concat".
#'}
#'
#'
#' @return  
#' To allow Eagle to handle data larger than the memory capacity of a machine, \code{ReadMarker} doesn't load 
#' the marker data into memory. Instead, it 
#' writes a reformatted version of the marker data, and its transpose, to the harddrive. These two files
#' are only temporary, being removed at the end of the R session. 
#' The object returned by
#' \code{ReadMarker} is a list object with the elements \code{tmpM} , \code{tmpMt}, and \code{dim_of_M}  
#' which is the full file name (name and path)  of the reformatted file for the marker  data,  the full file name of the reformatted file 
#' for the transpose of the marker  data,  and a 2 element vector with the first element the number of individuals and the second 
#' element the number of marker loci. 
#'
#' 
#'
#' @examples
#'   #--------------------------------
#'   #  Example 1
#'   #-------------------------------
#'   #
#'   # Read in the genotype data contained in the text file geno.txt
#'   #
#'   # The function system.file() gives the full file name (name + full path).
#'   complete.name <- system.file('extdata', 'geno.txt', package='Eagle')
#'   # 
#'   # The full path and name of the file is
#'   print(complete.name)
#'   
#'   # Here, 0 values are being treated as genotype AA,
#'   # 1 values are being treated as genotype AB, 
#'   # and 2 values are being treated as genotype BB. 
#'   # 4 gigabytes of memory has been specified. 
#'   # The file is space separated with the rows the individuals
#'   # and the columns the snp loci.
#'   geno_obj <- ReadMarker(filename=complete.name, type='text', AA=0, AB=1, BB=2, availmemGb=4) 
#'    
#'   # view list contents of geno_obj
#'   print(geno_obj)
#'
#'   #--------------------------------
#'   #  Example 2
#'   #-------------------------------
#'   #
#'   # Read in the allelic data contained in the PLINK ped file geno.ped
#'   #
#'   # The function system.file() gives the full file name (name + full path).
#'   complete.name <- system.file('extdata', 'geno.ped', package='Eagle')
#'
#'   # 
#'   # The full path and name of the file is
#'   print(complete.name)
#'   
#'   # Here,  the first 6 columns are being ignored and the allelic 
#'   # information in columns 7 -  10002 is being converted into a reformatted file. 
#'   # 4 gigabytes of memory has been specified. 
#'   # The file is space separated with the rows the individuals
#'   # and the columns the snp loci.
#'   geno_obj <- ReadMarker(filename=complete.name, type='PLINK', availmemGb=4) 
#'    
#'   # view list contents of geno_obj
#'   print(geno_obj)
#'
#'
#'
#'   #--------------------------------
#'   #  Example 3
#'   #-------------------------------
#'   #
#'   #
#'   # Read in the genotype data contained in the vcf file geno.vcf
#'   #
#'   # The function system.file() gives the full file name (name + full path).
#'   complete.name <- system.file('extdata', 'geno.vcf', package='Eagle')
#'   # 
#'   # The full path and name of the file is
#'   print(complete.name)
#'   
#'   # The file contains 5 marker loci recorded on 3 individuals
#'   # Two of the loci contain multiple alleles and are removed. 
#'   # A summary of the file is printed once the file has been read.
#'   geno_obj <- ReadMarker(filename=complete.name, type="vcf", availmemGb=4) 
#'    
#'   # view list contents of geno_obj
#'   print(geno_obj)
#'
ReadMarker <- function( filename=NULL, type='text', missing=NULL,
                           AA=NULL, AB=NULL, BB=NULL, 
                           availmemGb=16, 
                           quiet=TRUE ){


 if (nargs() == 0){
    ## checking that function has arguments
    message(" Please supply arguments to function \n")
    return(FALSE)
 }  else {
      ## read in either a text file or a PLINK file. The parameter type must be specified. Default it text file. 
   if(is.null(type)){
      message(" type must be set to \"text\" ,\"PLINK\", or \"vcf\". \n")
      message(" ReadMarker has terminated with errors.")
      return(FALSE)
   }
   if(!(type=="text" || type=="PLINK" || type=="vcf") ){
       message(" type must be set to \"text\" ,\"PLINK\", or \"vcf\". \n")
      message(" ReadMarker has terminated with errors")
      return(FALSE)
   }


    ## ------   PLINK ped file -------------------
    if (type=="PLINK"){
       
       ## checking if a PLINK file has been specified. 
       if (is.null(filename)){
            message(" The name of the PLINK ped file is missing.")
            message(" ReadMarker has terminated with errors.")
            return(FALSE)
       }
       if (!file.exists(fullpath(filename) )){
            message(" The PLINK ped file ", filename, " could not be found. ")
            message(" ReadMarker has terminated with errors ")
            return(FALSE)
       }

       ## Rcpp function to get dimensions of PLINK ped  file
       dim_of_M <- getRowColumn(fname=fullpath(filename))
       dim_of_M[2] <- (dim_of_M[2] - 6)/2  ## adjusting for PLINK structure
       ## Rcpp function to create binary packed Mt and M file (if needed)
       it_worked <- create.bin(file_genotype=fullpath(filename), type=type, availmemGb=availmemGb, dim_of_M=dim_of_M,  quiet=quiet  )
       if(!it_worked)
           return(FALSE) 

    if(.Platform$OS.type == "unix") {
       tmpM <- paste(tempdir(), "/", "M.bin", sep="")
       tmpMt <- paste(tempdir(), "/", "Mt.bin", sep="")
     } else {
       tmpM <- paste(tempdir()  , "\\", "M.bin", sep="")
       tmpMt <- paste(tempdir() , "\\", "Mt.bin", sep="")
     }




   }  
   if (type == "text"){

      ## ------------  text file -----------------------
      ## Assuming a text file that may be comma separated with numeric genotypes that need to be mapped onto AA, AB, and BB. 
      ## check of parameters
      error.code <- check.inputs(file_genotype=filename, availmemGb=availmemGb)
      if(error.code){
          message(" ReadMarker has terminated with errors.")
          return(FALSE)
       }
 ## Has AA, AB, BB been assigned character values
  if(is.null(AA) ||  is.null(BB))
  { 
     message("Error: The function parameters AA and BB must be assigned a numeric or character value since a text file is being assumed. \n")
     message(" Type help(ReadMarker) for help on how to read in the marker data. \n")
     message(" ReadMarker has terminated with errors")
     return(FALSE)
  }

  ## if there are no hets. 
  if(is.null(AB))
     AB <- "NA"  ## no hets 




  genofile <- fullpath( filename) 



  ## Rcpp function to get dimensions of ASCII genotype file
  message(" Getting number of individuals and snp from file ... ")
  dim_of_M <- getRowColumn(fname=genofile)

  ## Rcpp function to create ascii  M and Mt file from 
  message(" Beginning creation of reformatted file ... ")
  message(" Reading marker data ... \n")
  it_worked <- create.bin(file_genotype=genofile, type=type, AA=as.character(AA), AB=as.character(AB), BB=as.character(BB), 
              availmemGb=availmemGb, dim_of_M=dim_of_M, quiet=quiet, missing=missing  )

    if(!it_worked)   ## error has occurred. 
       return(FALSE)



     if(.Platform$OS.type == "unix") {
       tmpM <- paste(tempdir() , "/", "M.bin", sep="")
       tmpMt <- paste(tempdir() , "/", "Mt.bin", sep="")
     } else {
       tmpM <- paste(tempdir() , "\\", "M.bin", sep="")
       tmpMt <- paste(tempdir() , "\\", "Mt.bin", sep="")
     }


 
  }  ## end if else nargs()==1  (PLINK case)


  if (type == "vcf"){
     geno <- ReadVCF( filename=filename, availmemGb=availmemGb, quiet=TRUE )

  }




  if(type != "vcf"){
  geno <- list("tmpM"=tmpM, "tmpMt"=tmpMt,
               "dim_of_M" = dim_of_M,       
               "dim_of_Mt" = c(dim_of_M[2], dim_of_M[1]),
               "availmemGb" = availmemGb )

  if(.Platform$OS.type == "unix") {
       RDatafile <- paste(tempdir() , "/", "M.RData", sep="")
  } else {
       RDatafile <- paste( tempdir() , "\\", "M.RData", sep="")
  }




  save(geno, file=RDatafile)

}

  ## create M.Rdata file in current directory
  return(geno)

  } ## end if else nargs()==0

}  ## end function call ReadMarker



