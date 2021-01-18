GRR.matrix <- function(genes.maf=Kryukov, n.case.groups=2, GRR=c("SKAT", "constant", "variable"), GRR.value, GRR.function, GRR.multiplicative.factor=2, select.gene){
  ##Select MAF from the file given by the user  
  if(nlevels(genes.maf$gene)>1){
    if(missing(select.gene)){
      warning("More than one gene in the file, only the first one is used")
      select.gene <- levels(genes.maf$gene)[[1]]
    }
    pop.maf <- subset(genes.maf, genes.maf$gene %in% select.gene)$maf
  }else{
    pop.maf <- genes.maf$maf
  }
  n.variants <- length(pop.maf)

  ##Check on multiplicative values
  if(is.null(GRR.multiplicative.factor) & n.case.groups>1) stop("Needs 'GRR.multiplicative.factor'")
  if(!is.null(GRR.multiplicative.factor)){ 
    if(n.case.groups==1){
       warning("Only one group of cases, 'GRR.multiplicative.factor' is ignored")
    }else{
      if(length(GRR.multiplicative.factor)!=(n.case.groups-1)) stop("Wrong number of multiplicative factors")
    }
  }

  GRR <- match.arg(GRR)
                      
  ##Same GRR for all variants
  if(GRR=="constant"){
    if(missing(GRR.value)) stop("Needs 'GRR.value' for variants")
    if(length(GRR.value)>1){
      warning("Only one GRR value needed, only the first value will be used")
      GRR.value <- GRR.value[1]
    }
    #If only one group
    if(n.case.groups==1){
      GRR.matrix <- matrix(rep(GRR.value, n.variants),nrow=n.case.groups)
    }else{
      #Need n.groups-1 multiplicative factor: multiplication between group of case 1 and other groups of cases
      GRR.matrix <- matrix(rbind(rep(GRR.value, n.variants), t(sapply(GRR.multiplicative.factor, function(z) z*rep(GRR.value, n.variants)))), nrow=n.case.groups)
    }
  }
                                                                                                  
  if(GRR=="SKAT"){
    if(is.null(pop.maf)) stop("Needs MAF in the population to compute GRR")
    if(n.case.groups==1){
      GRR.matrix <- matrix(exp(0.402*abs(log10(pop.maf))), nrow=n.case.groups)
    }else{
      GRR.matrix <- matrix(rbind(exp(0.402*abs(log10(pop.maf))), t(sapply(GRR.multiplicative.factor, function(z) z*exp(0.402*abs(log10(pop.maf)))))), nrow=n.case.groups)
    }
  }
                                                                              
  if(GRR=="variable"){
    if(is.null(pop.maf)) stop("Needs MAF in the population to compute GRR")
    if(missing(GRR.function)) stop("Needs a function to compute GRR")
    if(n.case.groups==1){
      GRR.matrix <- matrix(do.call(GRR.function, list(pop.maf)),nrow=n.case.groups, byrow=FALSE)
    }else{
      GRR.matrix <- matrix(rbind(do.call(GRR.function, list(pop.maf)), t(sapply(GRR.multiplicative.factor, function(z) z*do.call(GRR.function, list(pop.maf))))), nrow=n.case.groups, byrow=FALSE)
    }
  }
  return(GRR.matrix)
}
