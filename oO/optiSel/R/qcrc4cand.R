
"qcrc4cand" <- function(kinship, condf, bc){
  qrc <- kinship
  for(i in seq_along(qrc)){
    qrc[[i]]$dir    <- "<="
    qrc[[i]]$use    <- qrc[[i]]$name %in% condf$var[!condf$isLin]
    qrc[[i]]$val    <- ifelse(qrc[[i]]$use, condf$val[condf$var==qrc[[i]]$name], NA)
    class(qrc[[i]]) <- str_replace(class(qrc[[i]]), "Fun", "Con")
    if(class(qrc[[i]])=="quadCon" && (qrc[[i]]$breed!="across breeds")){
      qrc[[i]]$val <- (qrc[[i]]$val-qrc[[i]]$d)*(bc[qrc[[i]]$breed])^2
      qrc[[i]]$d   <- 0
      qrc[[i]]$a   <- qrc[[i]]$a * bc[qrc[[i]]$breed]
    }
    if(class(qrc[[i]])=="ratioCon" && (qrc[[i]]$breed!="across breeds")){
      qrc[[i]]$d1 <- qrc[[i]]$d1*(bc[qrc[[i]]$breed])^2
      qrc[[i]]$d2 <- qrc[[i]]$d2*(bc[qrc[[i]]$breed])^2
      qrc[[i]]$a1 <- qrc[[i]]$a1* bc[qrc[[i]]$breed]
      qrc[[i]]$a2 <- qrc[[i]]$a2* bc[qrc[[i]]$breed]
    }
    
  }
  return(qrc)
}