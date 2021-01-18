 .calc_extBIC_MLE <- function(trait=NULL, currentX=NULL, MMt=NULL,  geno=NULL, Zmat=NULL, numberSNPselected=0, quiet=TRUE, lambda=NULL, eig.L=NULL)
 {
   ## internal function: used by AM 
   ## smallest extBIC and BIC is best
   ## internal function: use in AM only


   res_p <- emma.MLE(y=trait, X= currentX , K=MMt, Z=Zmat, llim=-100,ulim=100, 
                             eig.L=eig.L)
   BIC <- -2 * res_p$ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component


  # calculate lambda
  if(is.null(lambda)){
     lambda <- 1
     # found this to be anti-conservative when sample size is small
     #lambda <- log(geno$dim_of_M[2])/log(length(trait))
     #lambda <- 1-(1/(2*lambda))
   } ## outer if

  if(length(lambda)==1){
   extBIC <- BIC + 2 * lambda   *lchoose(geno$dim_of_M[2], numberSNPselected)  
   if(!quiet){
      cat(" Lambda  = ", lambda, "\n")
   }
   } else {
     extBIC <- rep(NA, length(lambda))
     for(ii in 1:length(lambda)){
        extBIC[ii] <- BIC + 2 * lambda[ii]  *lchoose(geno$dim_of_M[2], numberSNPselected)
     }

   }


    return(extBIC)
 }


