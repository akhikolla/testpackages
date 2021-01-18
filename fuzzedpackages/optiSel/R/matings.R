

"matings" <- function(phen, Kin, alpha=1, ub.n=NA, max=FALSE, solver="default", ...){
  if(!("Indiv" %in% colnames(phen))){stop("Column 'Indiv' is missing in phen.")}
  if(!("Sex"   %in% colnames(phen))){stop("Column 'Sex' is missing in phen.")}
  if(!("n"     %in% colnames(phen))){stop("Column 'n' is missing in phen.")}
  
  phen <- checkphen(phen, columns=c("Indiv", "Sex"), quiet=TRUE, na.Sex=FALSE) 
  if(!is.numeric(phen$n) | any(is.na(phen$n)) | any(phen$n<0) | any(floor(phen$n)!= phen$n)){
    stop("column 'n' must be numeric with non-negative integer values.")
  }
  if(sum(phen$n[phen$Sex=="male"])!=sum(phen$n[phen$Sex=="female"])){
    stop("Total numbers of matings for males amnd females must be equal.")
  }

  if("herd" %in% colnames(phen)){
    phen$herd <- as.character(phen$herd)
    phen$herd[phen$Sex=="male"] <- NA
    if(any(is.na(phen$herd[phen$Sex=="female"]))){
      stop("Column 'herd' contains NA for some females.")
    }
  }
  
  phen <- phen[phen$n>0,]
  Sire <- phen$Indiv[phen$Sex=="male"]
  Dam  <- phen$Indiv[phen$Sex=="female"]

  if(is.null(rownames(Kin))||is.null(colnames(Kin))){
    stop("Row names and column names of Matrix 'Kin' must be the individual IDs.")
  }

  if(!all(Sire %in% rownames(Kin))){
    stop("All male candidates must appear in the row names of matrix Kin.")
  }
  
  if(!all(Dam %in% colnames(Kin))){
    stop("All female candidates must appear in the column names of matrix Kin.")
  }
  
  if(alpha<1){
    if("herd" %in% colnames(phen)){
      phen$herd <- as.character(phen$herd)
      ub.herd   <- alpha * tapply(phen$n, phen$herd,sum)
    }else{
      stop("Column 'herd' is needed if alpha<1.")
    }
  }

  Kin <- Kin[Sire, Dam]
  Zeros <- 0*Kin
  
  rhsM <- phen$n[phen$Indiv %in% Sire]
  ConM <- NULL
  for(k in 1:length(Sire)){
    Con <- Zeros
    Con[k,] <- 1
    ConM <- rbind(ConM, c(Con))
  }
  
  rhsF <- phen$n[phen$Indiv %in% Dam]
  ConF <- NULL
  for(k in 1:length(Dam)){
    Con <- Zeros
    Con[,k] <- 1
    ConF <- rbind(ConF, c(Con))
  }

  if(alpha<1){
    herds <- setdiff(unique(phen$herd), NA)
    rhsH  <- rep(ub.herd[as.character(herds)], each=length(Sire))
    ConH  <- NULL
    for(l in herds){
      for(k in 1:length(Sire)){
        Con <- Zeros
        Con[k,] <- phen[Dam,"herd"] %in% l  
        ConH <- rbind(ConH, c(Con))
      }
    }
  }
  
  nVar <- nrow(Kin)*ncol(Kin)
  
  G <- NULL
  h <- NULL
  
  if(identical(solver, "default")){
    G <- rbind(G, -diag(nVar))
    h <- c(h, rep(0, nVar))
  }
  
  if(identical(solver, "default") & !is.na(ub.n)){
    G <- rbind(G, diag(nVar))
    h <- c(h, rep(ub.n, nVar))
  }
  
  if(alpha<1){
    G <- rbind(G,  ConH)
    h <- c(h, rhsH)
  }
  
  A <- rbind(ConF, ConM)
  b <- c(rhsF, rhsM)
  
  if(identical(solver, "default")){
    opt <- list(...)
    if("maxit"        %in% names(opt)){opt$maxit        <- as.integer(opt$maxit)}
    if("verbose"      %in% names(opt)){opt$verbose      <- as.integer(opt$verbose)}
    if("mi_max_iters" %in% names(opt)){opt$mi_max_iters <- as.integer(opt$mi_max_iters)}
    opt <- do.call(ecos.control, opt)
    
    dims <- list(l=length(h), q=NULL, e=0L)
    A    <- Matrix(A, sparse=TRUE)
    G    <- Matrix(G, sparse=TRUE)
    sig  <- ifelse(max,-1,1)
    
    res    <- ECOS_csolve(c=sig*c(Kin), G=G, h=h, dims=dims, A=A, b=b, int_vars=as.integer(1:nVar), control=opt)
    Mating <- round(res$x, 0)
    info   <- res$infostring
    objval <- sum(c(Kin)*res$x)/sum(res$x)
  }else{ #use Rsymphony_solve_LP
    if(is.na(ub.n)){
      bounds <- NULL
    }else{
      bounds <- list(upper=list(ind=1:(length(Sire)*length(Dam)), val=rep(ub.n, length(Sire)*length(Dam))))
    }
    
    Dir    <- c(rep("==", length(b)), rep("<=", length(h)))
    res    <- solver(obj=c(Kin), mat=rbind(A, G), dir=Dir, rhs=c(b, h), types="I", bounds=bounds, max=max, ...)
    Mating <- res$solution
    if(res$status==0L){
      info  <- "Optimum solution found"
    }else{
      info <- "No solution found"
    }
    objval <- sum(c(Kin)*res$solution)/sum(res$solution)
  }
  
  
  Mating <- matrix(Mating, nrow=nrow(Zeros), ncol=ncol(Zeros))
  Mating <- as.data.table(Mating)
  colnames(Mating) <- Dam
  Mating$Sire <- Sire
  Matings <- melt(Mating, id.vars="Sire", variable.name="Dam", value.name="n")
  Matings$Dam <- as.character(Matings$Dam)
  Matings <- Matings[Matings$n>0,]
  setDF(Matings)
  if(alpha<1){
    herds <- phen$herd
    names(herds) <- phen$Indiv
    Matings$herd <- herds[Matings$Dam]
    }
  attributes(Matings)$objval <- objval
  attributes(Matings)$info <- info
  cat(paste(info,"\n"))
  Matings
}