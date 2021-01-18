.build_design_matrix <- function(pheno=NULL,   fformula=NULL, quiet=TRUE, indxNA_pheno=NULL)
{

   ## assign model matrix X
   if(is.null(fformula))
   {  ## trait + intercept being fitted only
      Xmat <- matrix(data=1, nrow=nrow(pheno), ncol=1)
      colnames(Xmat) <- "(Intercept)"  ## to be consistent with model.matrix
   } else {

      # dealing with missing covariates
      fvars <- all.vars(fformula)  # list of fixed effects
      for (ii in fvars){ 
        # dealing with missingness in the covaries - setting to mean
        indx <- which(is.na(pheno[,ii]))
        if (length(indx) > 0){
           if(!is.factor(pheno[, ii])){            
              # this is a covariate. Replace missing value with mean
              m <- mean(pheno[, ii], na.rm=TRUE)
              pheno[indx, ii] <- m
           }
        } ## end if length
      } ## end for


      current.na.action <- options('na.action')[[1]]
      options(na.action='na.pass')   # this leaves NA's in the model matrix, preserving its correct dim
      Xmat <- model.matrix(fformula, data=pheno)
      options(na.action=current.na.action)

     # at this stage, if there are any NA's, it is due to a factor having NA values. 
     # replace all NAs with 0 values now
     indx <- which(is.na(Xmat), arr.ind=TRUE)
     if ( dim(indx)[1] > 0){
       Xmat[indx] <- 0
     }  ## end if dim
   }

   # dealing with missing trait values 
   # add extra factor manually
   if (!is.null(indxNA_pheno)){
     D <- diag(length(indxNA_pheno)) # creates a matrix, even if it is of length 1
     Zero <- matrix(data=0, nrow=nrow(Xmat)-length(indxNA_pheno) , ncol=ncol(D))  
     extrXmat <- rbind(D, Zero)
     # need to put rows in correct order to match missingness patern in trait
     indx <- 1:nrow(Xmat)
     indx <- indx[-indxNA_pheno]  # remove rows with missing data
     indx <- c(indxNA_pheno, indx)  # added back in but indx now out of order
     indx <- order(indx)
     extrXmat <- as.matrix(extrXmat[indx,])  # puts in correct order
     colnames(extrXmat) <- paste0("mv",1:length(indxNA_pheno)) # called mv? for missing value
     nms <- c(colnames(Xmat), colnames(extrXmat))
     Xmat <- cbind(Xmat, extrXmat )
     colnames(Xmat) <- nms
   }



 if (!quiet ){
   message("Dimension of design matrix, before addition of marker fixed effects is ", nrow(Xmat), " rows and ", ncol(Xmat), " columns.\n")
 }
if(!is.matrix(Xmat))
   Xmat <- matrix(data=Xmat, ncol=1)

## check that matrix doesn't all contain the same value
indx <- NULL
if (ncol(Xmat) > 1 ){
  for(ii in 2:ncol(Xmat)){
    # first column has intercept
    u <- length(unique(Xmat[, ii]))
    if(u == 1){
      indx <- c(indx, ii)
    }

  }
  if(length(indx) > 0)
     Xmat <- Xmat[, -indx]
}  


  return(Xmat)
}


