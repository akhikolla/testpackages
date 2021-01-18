  .calcVC <- function(trait, Zmat, currentX, MMt, eigenMMt, eig.L=NULL )
  {
    ## internal function: used by AM and SummaryAM
    ## perform likelihood ratio test for variance component Var_g
    # res_full <- emma.REMLE(y=trait, X= currentX , Z=Zmat, K=MMt, ngpu=ngpu)


 
    res_full <- emma.MLE(y=trait, X= currentX , Z=Zmat, K=MMt, eig.L=eig.L   )
    return(list("vg"=res_full$vg, "ve"=res_full$ve, "ML"=res_full$ML   ))

  }


