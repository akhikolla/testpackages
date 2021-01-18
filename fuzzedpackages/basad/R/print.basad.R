### print.basad.R
### The default print method for the basad
###
###
### Author: Qingyan Xiang


print.basad <- function(x, ...){

#### Check that object is compatible
  if (!inherits(x, "basad"))
     stop("This function only works for objects of class 'basad'")

#### extract summary data
  verboseList = x$verbose
  
#### --------------------------------------------------------------
###	Terminal Output
### ---------------------------------------------------------------

### basad output
    
cat("----------------------------------", "\n")
cat("Sample size                      :", verboseList[[1]], "\n" )
cat("Dimension                        :", verboseList[[2]], "\n" )
cat("Burn-in length                   :", verboseList[[3]], "\n" )
cat("Iteration length                 :", verboseList[[4]], "\n" )
cat("Block updating split sizes       :", verboseList[[5]], "\n" )
cat("Alternative fast sampling        :", verboseList[[6]], "\n" )
cat("Model selection criteria         :", verboseList[[7]], "\n" )
cat("\n\n")
cat("-----> Selected variables:\n")
print(round( x$select.var, 4 ))
cat("----------------------------------", "\n")

}
