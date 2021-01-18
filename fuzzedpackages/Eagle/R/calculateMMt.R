

calculateMMt <- function(geno=NULL, availmemGb, ncpu, selected_loci=NA, dim_of_M=NULL, quiet = TRUE )
{  
 ## internal function to AM
 ## R interface to Rcpp code to calculate M %*% t(M)
 ## Args
 ##      geno        absolute path + file name of binary packed M file
 ##      availmemGb    amount of memory in Gbytes available for creation of MMt
 ##      ncpu    number of cores for matrix operations
 ##      selectedloci an integer vector that gives the column number (0- L-1 ) of the loci that
 ##                   have been selected to act as fixed QTL effects in the model. 
 ##      dim_of_M    numeric vector with the row, column numbers of M. 
  #------------------------------------------
  # ascii file about to be overwritten
  #------------------------------------------
   
   
  if(!file.exists(geno)){
    message(" Error: The binary packed file ", geno, " cannot be found.\n")
    message(" calculateMMt has terminated with errors.")
    return(NULL)
   }
  if(!any(is.na(selected_loci))) selected_loci <- selected_loci-1
  MMt <- calculateMMt_rcpp( f_name=geno, selected_loci = selected_loci,
                               max_memory_in_Gbytes=availmemGb, num_cores=ncpu,
                               dims= dim_of_M, quiet = quiet, message=message) 
  return(MMt)
   
}  ## end function



   
   

