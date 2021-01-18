getRSLocation <- function(rs, species=NULL){
#  rm(list=ls())
#rs="rs14225498"
#species="Gallus_gallus"
#url="jul2016.archive.ensembl.org"
  
  if(is.null(species)) stop("Currently you need to provide a species")
    

  ensemblReturn <- tryCatch(readLines(paste("https://www.ensembl.org/",species,"/Variation/Explore?v=",rs,sep="")),
                            error = function(e) conditionMessage(e))
  
  tmp <- ensemblReturn[grepl(rs,ensemblReturn)]
  tmp <- tmp[2]
  
  out <- strsplit(strsplit(tmp,"r=")[[1]][2],";v=")[[1]][1]

  out <- strsplit(out,":")[[1]]
  out <- list(Chromosome=out[1], BP=out[2])
  out
}