
"lc4cand" <- function(phen, bc, condf, Traits, considerSexes){
  phen[is.na(phen$Sex),"Sex"] <- "undefined"
  Breeds <- names(considerSexes)
  
  A    <- do.call(rbind, lapply(Breeds, function(x){1*(phen$Breed==x)}))
  rownames(A) <- paste0("bc.", Breeds)
  use  <- rep(TRUE, nrow(A))
  dir  <- rep("==", nrow(A))
  val  <- bc[Breeds]
  name <- rownames(A)
  
  if(any(considerSexes)){
    sexes <- 2*(phen$Sex=="female")-1
    A2    <- do.call(rbind, lapply(Breeds[considerSexes], function(x){(phen$Breed==x)*sexes}))
    rownames(A2) <- paste0("scd.", Breeds[considerSexes])
    A     <- rbind(A, A2)
    use   <- c(use, rep(TRUE, nrow(A2)))
    dir   <- c(dir, rep("==", nrow(A2)))
    val   <- c(val, rep(0,    nrow(A2)))
    name  <- c(name, rownames(A2))
  }

  if(length(Traits)>0){
    alpha <- ifelse(Traits %in% condf$var, condf[Traits, "val"], -Inf)
    alpha <- matrix(alpha, nrow=nrow(phen), ncol=length(Traits), byrow=TRUE)
    A2    <- t(phen[, Traits] - alpha)
    A2[is.na(A2)] <- 0
    A     <- rbind(A, A2)
    use   <- c(use,  Traits %in% condf$var)
    dir   <- c(dir,  ifelse(Traits %in% condf$var, condf[Traits, "dir"], ">="))
    val   <- c(val,  rep(0,length(Traits)))
    name  <- c(name, Traits)
  }
  
  lc <- lincon(A=A,dir=dir, val=val, id=phen$Indiv, use=use, name=name)
  
}