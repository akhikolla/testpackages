.overcat <- function(msg, prevmsglength) {
  if (prevmsglength>0) {cat("\r")}    
  cat(msg)
  return(nchar(msg))
}

.ULI <- function(...) {
  redondGeo <- cbind(...) ## always a matrix
  if (ncol(redondGeo)==0L) return(rep(1L,nrow(redondGeo))) ## .generateInitPhi with constant predictor led here
  if (nrow(redondGeo)==1L) return(1L) ##  trivial case where the forllowingcode fails
  # redondFac <- apply(redondGeo,1,paste,collapse=" ") ## always characters whatever the number of columns 
  # redondFac <- as.integer(as.factor(redondFac)) ## as.factor effectively distinguishes unique character strings 
  # uniqueFac <- unique(redondFac) ## seems to preserve order ## unique(<integer>) has unambiguous behaviour
  # sapply(lapply(redondFac, `==`, uniqueFac ), which)
  redondFac <- apply(redondGeo,2L,factor,labels="") # not cute use of labels... 
  redondFac <- apply(redondFac,1L,paste,collapse="_") ## paste factors
  redondFac <- as.character(factor(redondFac))
  uniqueFac <- unique(redondFac) ## seems to preserve order ## unique(<integer>) has unambiguous behaviour
  uniqueIdx <- seq(length(uniqueFac))
  names(uniqueIdx) <- uniqueFac
  return(uniqueIdx[redondFac])
}

projpath <- function() {
  fn <- get("getActiveProject",envir = asNamespace("rstudioapi"))
  fn()
}