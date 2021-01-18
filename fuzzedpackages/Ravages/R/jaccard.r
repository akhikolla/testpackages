
Jaccard <- function(x, maf.threshold = 0.01) {
  if(maf.threshold > 0.5) 
    stop("MAF threshold can't be this high")
  
  w <- (x@snps$maf <= maf.threshold)
  inverse <- (x@p >= 1-maf.threshold)
 
  Ja <- .Call('oz_jaccard', PACKAGE = "Ravages", x@bed, w, inverse )
  J <- Ja$A/Ja$B
  # les individus qui ne portent aucun génotype rare 
  # ont un coeff d'identité nul avec tout le monde (y compris eux-même)
  J[Ja$B == 0] <- 0
  rownames(J) <- colnames(J) <- x@ped$id
  J
}
