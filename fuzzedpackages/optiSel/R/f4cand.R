

"f4cand" <- function(cand, method, Traits, bc){
  optiParam <- str_sub(method, 5, -1)
  if(optiParam %in% Traits){
    a <- cand$phen[[optiParam]]
    a[is.na(a)] <- 0
    f <- linfun(a=a, id=cand$phen$Indiv, name=optiParam)
  }else{
    f <- cand$kinship[[optiParam]]
    if(f$breed=="across breeds"){
      if("ratioFun" %in% class(f)){
        f$d1 <- 0
        f$d2 <- 0
        f$a1 <- 0*f$a1
        f$a2 <- 0*f$a2
      }
      if("quadFun" %in% class(f)){
        f$d <- 0
        f$a <- 0*f$a
      }      
    }else{
      if("ratioFun" %in% class(f)){
        f$d1 <- f$d1*(bc[f$breed])^2
        f$d2 <- f$d2*(bc[f$breed])^2
        f$a1 <- f$a1* bc[f$breed]
        f$a2 <- f$a2* bc[f$breed]
      }
      if("quadFun" %in% class(f)){
        f$d <- f$d*(bc[f$breed])^2
        f$a <- f$a*(bc[f$breed])
      }
    }
  }
  return(f)
}