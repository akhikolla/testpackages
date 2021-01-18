#We only accept matrices of good dimensions
genotypic.freq <- function(genes.maf = Kryukov, GRR.het, GRR.homo.alt, prev, genetic.model=c("general", "multiplicative", "dominant", "recessive"), select.gene) {
  
  genetic.model <- match.arg(genetic.model)
  #Test if a good genetic model is given
  genetic.model <- match.arg(genetic.model)

  #Selection of maf
  if (nlevels(genes.maf$gene) > 1) {
    if(missing(select.gene)){
      warning("More than one gene in the file, only the first one is used")
      select.gene <- levels(genes.maf$gene)[[1]]
    }
    pop.maf <- subset(genes.maf, genes.maf$gene %in% select.gene)$maf
  }else{
    pop.maf <- genes.maf$maf
  }

  #test dimensions of GRR.het
  if(nrow(GRR.het) != length(prev) | ncol(GRR.het) != length(pop.maf)) 
    stop("GRR.het dimensions mismatch")
    
  #Test on GRR.homo.alt only for general model , if not general model: only first GRR.het matrix used
  if(genetic.model=="general"){
    if(missing(GRR.homo.alt) | is.null(GRR.homo.alt)){
      stop("general model needs two GRR values per group")
    }else{
      if(nrow(GRR.het)!=nrow(GRR.homo.alt) | ncol(GRR.het)!=ncol(GRR.homo.alt)){
        stop("GRR.het and GRR.homo.alt have different dimensions")
      }
    }
  }else{
    if(!missing(GRR.homo.alt)){
      if(!is.null(GRR.homo.alt) ){
      warning("Only one GRR matrix needed for this model, only the first one is used")
      }
    }
  }

  
      
  #Creates matrix for each model
  if(genetic.model=="multiplicative"){
    GRR.homo.alt <- GRR.het**2
  }
  
  if(genetic.model=="dominant"){
    GRR.homo.alt <- GRR.het
  }
  
  #GRR.het=1 for heterozygous if genetic.model=recessive
  if(genetic.model=="recessive"){
    GRR.homo.alt <- GRR.het
    GRR.het <- matrix(rep(1, ncol(GRR.het)*nrow(GRR.het)), nrow=nrow(GRR.het)) 
  }    

  if(any(prev < 0) | sum(prev) > 1)
    stop("Unappropriate prev values")

  #Frequencies calculation
  homo.ref <- het <- homo.alt <- matrix(rep(0, ncol(GRR.het)*(nrow(GRR.het)+1)), nrow=nrow(GRR.het)+1)
  p.c <- matrix(rep(0,ncol(GRR.het)*nrow(GRR.het)), nrow=nrow(GRR.het))
  p.t <- numeric(ncol(GRR.het))
  for (i in 1:ncol(GRR.het)){
    freq.case <- p.case(pop.maf[i], GRR.het[,i], GRR.homo.alt[,i])
    freq.controls <- p.tem.GRR(pop.maf[i], GRR.het[,i], GRR.homo.alt[,i], prev=prev)

    homo.ref[,i] <- c(freq.controls[3], freq.case[,3])
    het[,i] <- c(freq.controls[2], freq.case[,2])
    homo.alt[,i] <- c(freq.controls[1], freq.case[,1])
  }
  
  if(nrow(homo.ref) == 2) 
    rownames(homo.ref) <- rownames(het) <- rownames(homo.alt) <- c("controls","cases")
  else 
    rownames(homo.ref) <- rownames(het) <- rownames(homo.alt) <- c("controls", sprintf("cases_%d", 1:(nrow(homo.ref)-1)))
  return(list(freq.homo.ref = homo.ref, freq.het = het, freq.homo.alt = homo.alt))
} 
  
  
p.case <- function(p, GRR, GRR.2){
  r <- matrix( 0.0, ncol = 3, nrow = length(GRR))
  r[,1] <- (GRR.2*p**2)/(GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2)     # freq homo alt
  r[,2] <- (GRR*2*p*(1-p))/(GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2)  # freq het
  r[,3] <- 1-r[,1]-r[,2]                                            # freq homo ref
  return(r)
}

p.tem.GRR <- function(p, GRR, GRR.2, prev){
  f <- prev / (GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2) #Frequence des cas chez les homo de ref
  r <- numeric(3)
  r[1] <- (p**2 * (1-sum(f*GRR.2))) / (1-sum(prev))           # freq homo alt
  r[2] <- (2*p*(1-p) * (1-sum(f*GRR))) / (1-sum(prev))        # freq het
  r[3] <- 1-r[1]-r[2]                                             # freq homo ref
  return(r)
}



