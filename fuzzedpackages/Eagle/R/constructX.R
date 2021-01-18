constructX <- function(Zmat=NULL, fnameMt=NULL, currentX=NULL, loci_indx=NULL,
                        dim_of_Mt=NULL,
                        map=NULL)
  {
    ## internal function for AM
    ## R function to construct the design matrix X
    ## Args
    ##   currentX    current model matrix
    ##   loci        the marker loci to be included as fixed QTL effects (additive model)

   if(is.na(loci_indx))
   {
     return(currentX)
   } else {
#       genodat <- extract_geno(fnameM=fnameM, colnum=loci_indx,
#                           availmemGb=availmemGb, dim_of_M=dim_of_M)

       genodat <- extract_geno_Mt(fnameMt=fnameMt, colnum=loci_indx,
                           dim_of_Mt=dim_of_Mt)
      if(is.null(Zmat)){
         newX <- cbind(currentX, genodat)
      } else {
         newX <- cbind(currentX, Zmat %*% genodat)
      }
      colnames(newX) <- c(colnames(currentX), as.character(map[[1]][loci_indx])) ## adding col names to new X  
      return(newX)
   }
  }



