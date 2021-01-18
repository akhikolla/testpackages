##random.bed.matrix with GRR
random.bed.matrix <- function(genes.maf = Kryukov, size, prev, replicates, 
                              GRR.matrix.del, GRR.matrix.pro, p.causal = 0.5, p.protect = 0, 
                              same.variant=FALSE, 
                              genetic.model=c("general", "multiplicative", "dominant", "recessive"), select.gene) {
  
  if (nlevels(genes.maf$gene) > 1){ 
    if(missing(select.gene)){
      warning("More than one gene in the file, only the first one is used")
      select.gene <- levels(genes.maf$gene)[[1]]
    }
    pop.maf <- subset(genes.maf, genes.maf$gene %in% select.gene)$maf
  }else{
    pop.maf <- genes.maf$maf
  }
  
  ##Check GRR
  if (!is.list(GRR.matrix.del)) {
    if (is.matrix(GRR.matrix.del)) {
      GRR.matrix.del <- list(GRR.matrix.del)
    }else{
      stop("GRR.matrix.del should be a list or a matrix")
    }
  GRR.het <- GRR.matrix.del[[1]]
  }

  genetic.model <- match.arg(genetic.model)

  if (length(GRR.matrix.del) == 1) {
    if (genetic.model == "general") {
      stop("Needs two GRR matrices in the general model")
    }else{
      GRR.homo.alt <- NULL
    }
  }else{
    if (genetic.model == "general") {
      GRR.homo.alt <- GRR.matrix.del[[2]]
    }else{
      warning("Only one GRR matrix needed for this model, only the first one is used")
      GRR.homo.alt <- NULL
    }
  }
  
  ##Same for protective
  if (!missing(GRR.matrix.pro)) {
    if (!is.list(GRR.matrix.pro)) {
      if (is.matrix(GRR.matrix.pro)) {
        GRR.matrix.pro <- list(GRR.matrix.pro)
      }else{
        stop("GRR.matrix.pro should be a list or a matrix")
      }
    }
    GRR.het.pro <- GRR.matrix.pro[[1]]
    if (length(GRR.matrix.pro) == 1) {
      if (genetic.model == "general") {
        stop("Needs two GRR matrices in the general model")
      }else{
        GRR.homo.alt.pro <- NULL
      }
    }else{
      if (genetic.model == "general") {
        GRR.homo.alt.pro <- GRR.matrix.pro[[2]]
      }else{
        warning("Only one GRR matrix needed for this model, only the first one is used")
        GRR.homo.alt.pro <- NULL
      }
    }
  }else{
    GRR.het.pro <- 1/GRR.het
    GRR.homo.alt.pro <- 1/GRR.homo.alt
  }
  
  ##Check on GRR values
  if (any(GRR.het < 1) | any(GRR.homo.alt < 1)) stop("Matrix of deleterious GRR has GRR values lower than 1")
  if (!missing(GRR.matrix.pro)) {
    if (any(GRR.het.pro > 1) | any(GRR.homo.alt.pro > 1)) stop("Matrix of protective GRR has GRR values greater than 1")
  }

  ##Order arguments
  GRR.pars <- list(OR.del = GRR.het, OR.pro = GRR.het.pro, p.causal = p.causal, prob.pro = p.protect)
  if (!is.null(GRR.homo.alt)) {
    GRR.homo.alt.pars <- list(OR.del = GRR.homo.alt, OR.pro = GRR.homo.alt.pro)
  }else{
    GRR.homo.alt.pars <- NULL
  }

  ##Check dimensions
  if(!is.null(GRR.homo.alt.pars)){
    if(nrow(GRR.het) != nrow(GRR.homo.alt.pars$OR.del) |  ncol(GRR.het) != ncol(GRR.homo.alt.pars$OR.del)) stop("GRR.het and GRR.homo.alt have different dimensions")
  }
  
  ##Choose the OR function
  if(same.variant){
    variant.function <- OR.matrix.same.variant
  }else{
    variant.function <- OR.matrix
  }
  
  GRR.pars$n.variants <- length(pop.maf)
  nb_snps <- GRR.pars$n.variants * replicates
  nb_inds <- sum(size)
  x <- new.bed.matrix(nb_inds, nb_snps);
  for(b in 1:replicates) {
    GRR.het <- do.call( variant.function, GRR.pars)
    if(!is.null(GRR.homo.alt.pars)){
      #Give GRR>1 for the same variants as GRR.het
      for(i in 1:length(prev)){
        GRR.homo.alt[i, which(GRR.het[i,]==1)] <- 1
        GRR.homo.alt[i, which(GRR.het[i,]<1)] <- GRR.homo.alt.pars$OR.pro[i, which(GRR.het[i,]<1)]
      }
    }else{
      GRR.homo.alt=NULL
    }    
    MAFS <- genotypic.freq(genes.maf=genes.maf, GRR.het=GRR.het, GRR.homo.alt=GRR.homo.alt, prev=prev, select.gene=select.gene, genetic.model=genetic.model)
  #Check if problems with model
    if(any(MAFS$freq.homo.ref[1,]>1 | MAFS$freq.het[1,]<0 | MAFS$freq.homo.alt[1,]<0)) stop("Impossible genetic model, please change your parametrization")
    .Call("oz_random_filling_bed_matrix_noHW", PACKAGE = "Ravages", x@bed, MAFS$freq.homo.ref, MAFS$freq.het, size, (b-1)*GRR.pars$n.variants)
  }
  x@ped$pheno <- factor(rep.int( 1:length(size) - 1, size))
  x@snps$genomic.region <- factor( rep( sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = GRR.pars$n.variants) )
  x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  x <- set.stats(x, verbose = FALSE)
  x
}

