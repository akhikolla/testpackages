##' Joint modeling of longitudinal continuous data and competing risks
##' @title Joint Modelling for Continuous outcomes 
##' @param p  The dimension of fixed effects (include intercept) in yfile.
##' @param yfile Y matrix for longitudinal measurements in long format. For example, for a subject with n measurements, there should be n rows for this subject. The # of rows in y matrix is the total number of measurements for all subjects in the study. The columns in Y should start with the longitudinal outcome (column 1), the covariates for the random effects, and then the covariates for the fixed effects.
##' @param cfile C matrix for competing risks failure time data. Each subject has one data entry, so the number of rows equals to the number of subjects. The survival / censoring time is included in the first column, and the failure type coded as 0 (censored events), 1 (risk 1), or 2 (risk 2) is given in the second column. Two competing risks are assumed. The covariates are included in the third column and on.
##' @param mfile M vector to indicate the number of longitudinal measurements per subject. The number of rows equals to the number of subjects.
##' @param point Quadrature points used in the EM procedure.Default is 20.
##' @param maxiterations Maximum values of iterations. Default is 100000.
##' @param do.trace Print detailed information of each iteration. Default is false, i.e., not to print the iteration details.
##' @param type_file Types of inputs. Default is true, i.e.  data files with headers. If set to "F", inputs are changed to data matrixes or data.frames (with headers)
##' @return Object of class \code{JMcmprsk} with elements
##'   \tabular{ll}{
##'       \code{vcmatrix}    \tab  The variance-covariance matrix for all the parameters. The parameters are in the order: \eqn{\beta}, \eqn{\sigma^2}, \eqn{\gamma}, \eqn{\nu}, and \eqn{\Sigma}. The elements in \eqn{\Sigma} are output in the order along the main diagonal line, then the second main diagonal line, and so on. \cr
##'       \code{betas} \tab The point  estimates of \eqn{\beta}. \cr
##'       \code{se_betas} \tab The standard error estimate of \eqn{\beta}. \cr
##'       \code{gamma_matrix} \tab  The point  estimate of \eqn{\gamma}. \cr
##'       \code{se_gamma_matrix}   \tab  The standard error estimate of \eqn{\gamma}. \cr
##'       \code{v_estimate} \tab The point  estimate of \eqn{\nu}. \cr
##'       \code{se_v_estimate}    \tab The standard error estimate of \eqn{\nu}. \cr
##'       \code{sigma2_val}     \tab  The point estimate of \eqn{\sigma^2}.\cr
##'       \code{se_sigma2_val}     \tab  The standard error estimate of \eqn{\sigma^2}.\cr
##'       \code{sigma_matrix}     \tab The point estimate of \eqn{\Sigma} (only the upper triangle portion of the matrix is output).\cr
##'       \code{se_sigma}     \tab The standard error estimate of \eqn{\Sigma}.The standard errors are given in this order: main diagonal, the second main diagonal, and so on. \cr
##'       \code{loglike}     \tab Log Likelihood.\cr
##'   }
##'   
##' @examples
##' # A toy example on a dataset called from file paths
##' require(JMcmprsk)
##' set.seed(123)
##' yfile=system.file("extdata", "jmcsimy.txt", package = "JMcmprsk")
##' cfile=system.file("extdata", "jmcsimc.txt", package = "JMcmprsk")
##' mfile=system.file("extdata", "jmcsimm.txt", package = "JMcmprsk")
##' jmc_0fit = jmc_0(p=4, yfile, cfile, mfile, point=6, do.trace = FALSE)

##' \dontrun{
##' # A toy example on data frames/matrices
##' require(JMcmprsk)
##' set.seed(123)
##' data(lung)
##' lungY <- lung[, c(2:11)]
##' lungC <- unique(lung[, c(1, 12, 13, 6:10)])
##' lungC <- lungC[, -1]
##' lungM <- data.frame(table(lung$ID))
##' lungM <- as.data.frame(lungM[, 2])
##' res1=jmc_0(p=8, lungY, lungC, lungM, point=20, do.trace = FALSE, type_file = FALSE)
##' res1
##' }


##' @references
##' \itemize{
##' \item Elashoff, Robert M., Gang Li, and Ning Li. "A joint model for longitudinal measurements and survival data in the presence of multiple failure types." Biometrics 64.3 (2008): 762-771.
##' }
##' @seealso \code{\link{jmo}}
##' @export
jmc_0<- function (p,yfile,cfile,mfile,point=20,maxiterations=100000,do.trace=FALSE,type_file=TRUE)
{
  if (do.trace) { 
    trace=1;
  }else{
    trace=0;
  }
  
  
  #Gaussian-Hermite quadrature nodes and weights
  #The dimension of xs/ws is half of the point value since they are symmetric
  
  gq_vals <- statmod::gauss.quad(n = point, kind = "hermite")
  
  xs <- gq_vals$nodes[(point / 2 + 1) : point]
  
  ws <- gq_vals$weights[(point / 2 + 1) : point]
  
  if (type_file){
  # store header names for future useage
  ydata=read.table(yfile,header = TRUE)
  ynames=colnames(ydata)
  yfile=tempfile(pattern = "", fileext = ".txt")
  writenh(ydata,yfile)
  
  cdata=read.table(cfile,header = TRUE)
  cnames=colnames(cdata)
  cfile=tempfile(pattern = "", fileext = ".txt")
  writenh(cdata,cfile)
  
  
  ydim=dim(ydata)
  # number of observations in study is equals to the #of rows in Y matrix
  n1=ydim[1];
  # dim of fixed effects plus dim of random effects should be 
  # total the column of y -the survival time column 1
  p1a=ydim[2]-1-p;
  
  if((p<1)|(p1a<1)){
    stop("Possibe wrong dimension of fixed effects in Y!")
  }
  
  if (p1a > 3) {
    stop("Maximum of 3 random effects are allowed. Please reconsider the random effect covariates you need!")
  }
  
  cdim=dim(cdata);
  
  ##Check the completeness of data
  if (sum(complete.cases(ydata)) != ydim[1]) {
    stop("Missing values detected! Please make sure your longitudinal data is complete!")
  }
  
  if (sum(complete.cases(cdata)) != cdim[1]) {
    stop("Missing values detected! Please make sure your survival data is complete!")
  }
 
  # number of subjects in study is equals to the #of rows in C matrix
  k=cdim[1];
  #
  p2=cdim[2]-2;

  maxl=max(read.table(mfile));
  myresult=jmc_main(k,n1, p,p2, maxl,p1a, maxiterations, point,xs,ws, yfile,cfile,mfile,trace)
  
  PropComp <- round(table(cdata[, 2])/k * 100, 2)
  
  }else{
    ynames=colnames(yfile)
    yfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(yfile,yfilenew)
    
    cnames=colnames(cfile)
    cfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(cfile,cfilenew)
    
    mfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(mfile,mfilenew)
    
    ydim=dim(yfile)
    # number of observations in study is equals to the #of rows in Y matrix
    n1=ydim[1];
    # dim of fixed effects plus dim of random effects should be 
    # total the column of y -the survival time column 1
    p1a=ydim[2]-1-p;
    
    if((p<1)|(p1a<1)){
      stop("Possibe wrong dimension of fixed effects in Y!")
    }
    
    if (p1a > 3) {
      stop("Maximum of 3 random effects are allowed. Please reconsider the random effect covariates you need!")
    }
    
    cdim=dim(cfile);
    
    ##Check the completeness of data
    if (sum(complete.cases(yfile)) != ydim[1]) {
      stop("Missing values detected! Please make sure your longitudinal data is complete!")
    }
    
    if (sum(complete.cases(cfile)) != cdim[1]) {
      stop("Missing values detected! Please make sure your survival data is complete!")
    }
    
    # number of subjects in study is equals to the #of rows in C matrix
    k=cdim[1];
    #
    p2=cdim[2]-2;
    
    maxl=max(mfile);
  
    
  myresult=jmc_main(k,n1, p,p2, maxl,p1a, maxiterations, point,xs,ws, yfilenew,cfilenew,mfilenew,trace)  
    
    PropComp <- round(table(cfile[, 2])/k * 100, 2)
    
  }
  
  

 


  myresult$type="jmc";

  #names
  names(myresult$betas)=ynames[(ydim[2]-p+1):ydim[2]]

  colnames(myresult$gamma_matrix)=cnames[3:(p2+3-1)]

  myresult$k=k
  
  ##create longitudinal submodel formula
  #ynames <- colnames(ydata)
  LongOut <- ynames[1]
  LongX <- paste0(ynames[(3+p1a):length(ynames)], collapse = "+")
  FunCall_long <- as.formula(paste(LongOut, LongX, sep = "~"))
  
  ##create survival submodel formula
  #cnames <- colnames(cdata)
  SurvOut <- paste0("Surv(", cnames[1], ",", cnames[2], ")")
  SurvX <- paste0(cnames[-(1:2)], collapse = "+")
  FunCall_survival <- as.formula(paste(SurvOut, SurvX, sep = "~"))
  
  
  DataPath <- NULL
  SummaryInfo <- list(k, n1, PropComp, FunCall_long, FunCall_survival, DataPath)
  names(SummaryInfo) <- c("NumSub", "Numobs", "PropComp", 
                          "LongitudinalSubmodel", "SurvivalSubmodel", "DataPath")
  
  myresult$SummaryInfo <- SummaryInfo
  myresult$point <- point
  class(myresult) <- "JMcmprsk"
  
  mycall=match.call()
  myresult$call=mycall
  
  return (myresult)
}




