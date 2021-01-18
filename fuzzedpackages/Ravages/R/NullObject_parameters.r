NullObject.parameters <- function(pheno, RVAT, pheno.type, ref.level, data, formula){
  if(missing(data)) data <- NULL
  if(missing(formula)) formula <- NULL
  if(!(pheno.type %in% c("categorial", "continuous"))) stop("'pheno.type' should be 'categorial' or 'continuous'")
  if(!(RVAT %in% c("burden", "SKAT"))) stop ("'RVAT' should be 'burden' or 'SKAT'")
  if(pheno.type == "categorial"){
    if (!is.factor(pheno))  stop("'pheno' is not a factor")
    pheno <- droplevels(pheno)
    if(RVAT == "burden"){
      if(is.null(ref.level)) stop("'ref.level' should be specified for the multinomial regression")
      NullObject.params <- burden.NullObject(group = pheno, ref.level = ref.level, data = data, formula = formula)
    }
    if(RVAT == "SKAT") NullObject.params <- SKAT.NullObject(group = pheno, data = data, formula = formula)
  }
  if(pheno.type == "continuous"){
    if(RVAT == "burden") NullObject.params <- burden.NullObject.continuous(pheno = pheno, data = data, formula = formula)
    if(RVAT == "SKAT") NullObject.params <- SKAT.NullObject.continuous(pheno = pheno, data = data, formula = formula)
  }
  tt <- list(pheno.type = pheno.type)
  return(c(NullObject.params, tt))
}
