.form_results <- function(traitname, trait, selected_loci,   fformula, indxNA_pheno, 
                           ncpu,  availmemGb, quiet,  extBIC, lambda, geno, pheno, map, Zmat, outlierstat )
{
  ## internal function - used by AM for forming the results object
  if (length(selected_loci) > 1){
   sigres <- list(traitname=traitname, 
                    trait=trait,
                    fformula = fformula,
                    indxNA_pheno = indxNA_pheno,
                    Mrk=map[[1]][selected_loci],
                    Chr=map[[2]][selected_loci],
                    Pos=map[[3]][selected_loci],
                    Indx=selected_loci,
                    ncpu=ncpu,
                    ngpu=computer$ngpu,
                    availmemGb=availmemGb,
                    quiet=quiet,
                    extBIC=extBIC,
                    lambda=lambda,
                    geno=geno,
                    pheno=pheno,
                    map=map,
                    Zmat=Zmat,
                    outlierstat=outlierstat) 
  } else {
   sigres <- list(traitname=traitname,
                    trait=trait,
                    fformula = fformula,
                    indxNA_pheno = indxNA_pheno,
                    Mrk=NA,
                    Chr=NA,
                    Pos=NA,
                    Indx=selected_loci,
                    ncpu=ncpu,
                    ngpu=computer$ngpu,
                    availmemGb=availmemGb,
                    quiet=quiet,
                    extBIC=extBIC,
                    lambda=lambda,
                    geno=geno,
                    pheno=pheno,
                    map=map,
                    Zmat=Zmat, 
                    outlierstat=outlierstat) 
}

return(sigres)
}



