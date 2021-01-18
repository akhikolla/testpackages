
extract_geno_Mt <- function(fnameMt=NULL, colnum=NULL, dim_of_Mt=NULL)
  { ## internal function for AM
    ## Rcpp function to extra a row of genotypes from bin file Mt 

    selected_locus <- colnum - 1  ## to be consistent with C++'s indexing starting from 0

    genodata <- extract_geno_Mt_rcpp(f_name=fnameMt,
                              selected_locus=selected_locus, dims=dim_of_Mt)


    return(genodata)

  }



