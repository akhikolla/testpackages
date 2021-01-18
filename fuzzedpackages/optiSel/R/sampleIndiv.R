
"sampleIndiv" <- function(Pedig, from="Born", each=100){
  IDs  <- split(Pedig$Indiv, Pedig[[from]])
  keep <- unlist(mapply(sample, x=IDs, size=pmin(lapply(IDs,length), each), SIMPLIFY=FALSE))
  names(keep)<-NULL
  keep
}
