#' @title Summary of multiple locus association mapping results
#' @description    A summary function that provides additional information on the significant 
#'     marker-trait associations found by \code{\link{AM}}
#' @param  AMobj  the (list) object obtained from running \code{\link{AM}}. 
#' @details
#'
#' \code{SummaryAM} produces three  tables, an overall summary table, a table of the 
#'  SNP names and positions, and a table  of results with 
#' the  p-value for each 
#' fixed effect in the final model.  
#' @examples
#'  \dontrun{
#'   # Since the following code takes longer than 5 seconds to run, it has been tagged as dontrun. 
#'   # However, the code can be run by the user. 
#'   #
#'
#'   #---------------
#'   # read the map 
#'   #---------------
#'   #
#'   # File is a plain space separated text file with the first row 
#'   # the column headings
#'   complete.name <- system.file('extdata', 'map.txt', 
#'                                    package='Eagle')
#'   map_obj <- ReadMap(filename=complete.name) 
#'
#'  # to look at the first few rows of the map file
#'  head(map_obj)
#'
#'   #------------------
#'   # read marker data
#'   #------------------
#'   # Reading in a PLINK ped file 
#'   # and setting the available memory on the machine for the reading of the data to 8 gigabytes
#'   complete.name <- system.file('extdata', 'geno.ped', 
#'                                      package='Eagle')
#'   geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
#'  
#'   #----------------------
#'   # read phenotype data
#'   #-----------------------
#'
#'   # Read in a plain text file with data on a single trait and two fixed effects
#'   # The first row of the text file contains the column names y, cov1, and cov2. 
#'   complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
#'   
#'   pheno_obj <- ReadPheno(filename=complete.name)
#'            
#'   #-------------------------------------------------------
#'   # Perform multiple-locus genome-wide association mapping 
#'   #-------------------------------------------------------                   
#'   res <- AM(trait = 'y',
#'                            fformula=c("cov1 + cov2"),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj)
#'
#'   #-----------------------------------------
#'   # Produce additional summary information 
#'   #------------------------------------------
#'
#'   SummaryAM(AMobj=res)
#'  }
#'
#' 
#' 
#' @seealso \code{\link{AM}}
#'
SummaryAM <- function(AMobj=NULL)
{



 if(is.null(AMobj)){
    message(" SummaryAM function requires AMobj object to be specified. This object is obtained by running AM().")
    return(NULL)
    }
 if(!is.list(AMobj)){
    message(" SummaryAM function requires AMobj object to be a list object.")
    return(NULL)
   }




## Table 1
# create list with necessary components 
lst <- list( "ncpu" = AMobj[["ncpu"]],
             "memory" = AMobj[["availmemGb"]],
             "numsnp"=nrow(AMobj[["map"]]),
            "numsamples"=nrow(AMobj[["pheno"]]),
            "traitname"=AMobj[["traitname"]], 
            "fformula" = AMobj[["fformula"]],
            "nummissingpheno" = AMobj[["indxNA_pheno"]],
            "numSigSNP"=AMobj[["Mrk"]],
            "lambda" = AMobj[["lambda"]] 
         )



if (is.null(lst[["fformula"]])){
 lst[["fformula"]] <- "intercept only"
} else {
  lst[["fformula"]] <- as.character(AMobj[["fformula"]])[2]
} 

if (is.null(lst[["nummissingpheno"]])){
lst[["nummissingpheno"]] <- 0
} else {
lst[["nummissingpheno"]] <- length(lst[["nummissingpheno"]])
}
if (length(lst[["numSigSNP"]]) == 1) {
  lst[["numSigSNP"]] <- 0
} else {
  lst[["numSigSNP"]] <- length(lst[["numSigSNP"]]) - 1
}








message("\n\n Table 1: Summary Information \n   ")
message(  sprintf("%50s", "--------------------------------------------------------" ))
message( sprintf("%-40s  %-10s", "Number cpu: ", lst[["ncpu"]] ))
message( sprintf("%-40s  %-10s", "Max memory (Gb): ", lst[["memory"]] ))
message( sprintf("%-40s  %-10s", "Number of samples: ", lst[["numsamples"]] ))
message( sprintf("%-40s  %-10s", "Number of snp: ", lst[["numsnp"]] ))
message( sprintf("%-40s  %-10s", "Trait name: ", lst[["traitname"]] ))
message( sprintf("%-40s  %-30s", "Fixed model: ", lst[["fformula"]] ))
message( sprintf("%-40s  %-10s", "Number samples missing obs:", lst[["nummissingpheno"]] ))
message( sprintf("%-40s  %-10s", "Number significant snp-trait assocs:", lst[["numSigSNP"]] ))
message( sprintf("%-40s  %-4s", "Lambda value for extBIC: ", round(lst[["lambda"]],2)  ))
message(  sprintf("%50s", "--------------------------------------------------------" ))
message("\n\n")

# create data frame of summary information for use in shiny app
infodf <- data.frame("description"= c("Number cpu", "Max memory (Gb)", "Number of samples", "Number of snp", 
                                    "Trait name", "Fixed model", "Number samples missing obs", 
                                    "Number significant snp-trait assocs", "Lambda value for extBIC"  )   ,
                   "value" = c(lst[["ncpu"]], lst[["memory"]], lst[["numsamples"]], 
                     lst[["numsnp"]], lst[["traitname"]], lst[["fformula"]], lst[["nummissingpheno"]], 
                    lst[["numSigSNP"]], round(lst[["lambda"]],2)   )
         )




## Table Findings 
message("\n\n Table 2: Findings \n   ")

if ( length(AMobj$Mrk)==1){
   message(" No findings ... ")
   df_findings = NA
   message("\n\n")
} else {
  message(sprintf("%22s   %11s   %10s   %10s",  "SNP",   "Chr", "Position" , "Col index"))
message(  sprintf("%60s", "------------------------------------------------------------------" ))
  for(ii in 2:length(AMobj$Mrk ) )
  {
      message(sprintf("%22s    %10s   %10.2f   %10.0f",
            AMobj$Mrk[ii],  AMobj$Chr[ii], AMobj$Pos[ii], AMobj$Indx[ii]))
  }  ## end for ii
message(  sprintf("%60s", "------------------------------------------------------------------" ))
 message("\n\n\n")
  df_findings <- data.frame(SNP=AMobj$Mrk[-1], Chr=AMobj$Chr[-1], 
                            Pos=AMobj$Pos[-1], ColIndx=AMobj$Indx[-1])

}

  ## create table for use by Shiny

## Table III
 ## check to make sure that null model is not being supplied
 if (length(AMobj$Mrk)==1){
   message(" No significant marker-trait associations have been found by AM. \n")
   message(" No p-values to report \n")
   return(invisible(lst) )
 }



  ## build environmental effects design matrix
  baseX <- .build_design_matrix(pheno=AMobj$pheno,  
                                    fformula=AMobj$fformula,
                                   quiet=AMobj$quiet, indxNA_pheno=AMobj$indxNA_pheno)
  ## add genetic marker effects 
  fullX <- baseX


  for(loc in AMobj$Indx){
           
           fullX <- constructX(Zmat=AMobj$Zmat, fnameMt=AMobj$geno[["tmpMt"]], 
                              currentX=fullX, loci_indx=loc,
                               dim_of_Mt=AMobj$geno[["dim_of_Mt"]],
                                map=AMobj$map)
  }  ## end for loc
  cnames <- colnames(fullX)  ## names of all the variables in the model


  ## calculate MMt
  MMt <- .calcMMt(AMobj$geno,  AMobj$ncpu, AMobj$Indx, AMobj$quiet)

  ## calculate variance components of LMM

   Args <- list(y= AMobj$trait  , X= fullX , Z=AMobj$Zmat, K=MMt)

  #eR <- emma.REMLE(y=AMobj$trait, X= fullX , K=MMt, llim=-100,ulim=100)
  eR <-  do.call(emma.MLE, Args)






 ## calculating p values of fixed marker effect via Wald statistic
 mrks <- AMobj$Mrk[-1]  ## its -1 to remove the NA for the null model 
 pval <- vector("numeric", length(colnames(fullX)) )


 H <- calculateH(MMt=MMt, varE=eR$ve, varG=eR$vg, Zmat=AMobj$Zmat)

# H <-  eR$vg * MMt + eR$ve * diag(1, nrow(MMt))
 Hinv <- try(chol2inv(chol(H)))
 if (class(Hinv) ==  "try-error"){
    message(" Inverse for H matrix has failed. ")
    message(" Error in SummaryAM function has occurred. ")
    return(FALSE)
 }





 beta <- try(chol2inv(chol( t(fullX) %*% Hinv %*% fullX)) %*% t(fullX) %*% Hinv %*% matrix(data=AMobj$trait ,ncol=1)   )
 if (class(beta) ==  "try-error"){
    message(" Matrix inverse for beta vector has failed. ")
    message(" Error in SummaryAM function has occurred. ")
    return(FALSE)
 }



df <- rep(1, length(cnames))
pval <- rep(NA, length(cnames))
W <- rep(NA, length(cnames))

for(ii in cnames){
   indx <- which(cnames==ii)
   L <- matrix(data=rep(0, length(indx)*ncol(fullX)), byrow=TRUE, nrow=length(indx) )  # r x p matrix
   LL <- diag(length(indx))
   L[,indx] <- LL



  W[which(ii==cnames)]  <- t(L %*% beta) %*%
            chol2inv(chol( L %*% chol2inv(chol(t(fullX) %*% Hinv %*% fullX)) %*% t(L) )) %*%
            (L %*% beta)
 pval[which(ii==cnames)] <- 1 - pchisq(W[which(ii==cnames)], length(indx))  ## its not -1 here because fullX already has 1
                                                      ## less factor level
 }




# create newvarnames that does not have mv? covariates (if any)
newvarnames <- cnames
if (length(AMobj$indxNA_pheno)> 0){
  nms <- paste0("mv", 1:length(AMobj$indxNA_pheno))
  indx <- match(nms, cnames)
   newvarnames <- cnames[-indx] 

}

message("\n\n Table 3: Size and Significance of Effects in Final Model \n   ")

  message(sprintf("%22s %11s %6s   %10s   %13s", "", "Effect Size",   "Df", "Wald statstic" , "Pr(Chisq)"))
message(  sprintf("%70s", "----------------------------------------------------------------------------" ))
  for(ii in newvarnames )
  {
      indx <- which(cnames == ii)
      message(sprintf("%20s    %10.2f %6i   %13.2f       %.3E",
         ii,   beta[indx], df[indx], W[indx], pval[indx ]))
  }  ## end for ii
message(  sprintf("%70s", "----------------------------------------------------------------------------" ))
 message("\n\n\n")



df_size <- data.frame("Effects"=cnames, "Size"=as.character(round(beta,2)),  "Df"=as.character(df),   
                      "Wald statistic"=as.character(round(W,2)),        
                      "Pr(Chisq)"=pval, check.names=FALSE)


  res <- list()
  res[["pvalue"]] <- pval
  res[["size"]] <- df_size
  res[["Waldstat"]] <- W
  res[["df"]] <- df
  res[["summarylist"]] <- infodf
  res[["findings"]] <- df_findings

  
  return(invisible(res))
}
