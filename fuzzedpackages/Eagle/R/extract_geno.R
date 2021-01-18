
extract_geno <- function(fnameM=NULL, colnum=NULL, availmemGb=8,
                          dim_of_M=NULL,
                          selected_locus=NA)
  { ## internal function for AM
    ## Rcpp function to extra a column of genotypes from ascii file M

    selected_locus <- colnum - 1  ## to be consistent with C++'s indexing starting from 0

    genodata <- extract_geno_rcpp(f_name_ascii=fnameM,
                               max_memory_in_Gbytes = availmemGb,
                              selected_locus=selected_locus, dims=dim_of_M)

    return(genodata)

  }



