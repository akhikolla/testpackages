.calcMMt <- function(geno,  ncpu, selected_loci, quiet)
  {
    ## internal function: used only in AM  and SummaryAM
    ## calculates M %*% t(M) via C++ for out of memory calculation
    MMt <- calculateMMt(geno=geno[["tmpM"]], availmemGb=geno[["availmemGb"]],
                           ncpu=ncpu,
                           dim_of_M = geno[["dim_of_M"]],
                           selected_loci=selected_loci, quiet = quiet )
    gc()


    ## Trick for dealing with singular MMt due to collinearity
    MMt <- MMt/max(MMt) + diag(0.95, nrow(MMt))
    return(MMt)
  }



