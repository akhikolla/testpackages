##' Joint modeling of longitudinal ordinal data and competing risks
##' @title Joint Modelling for Ordinal outcomes 
##' @param p  The dimension of proportional odds covariates (not including intercept) in yfile.
##' @param s  The dimension of non-proportional odds covariates in yfile.
##' @param yfile Y matrix for longitudinal measurements in long format. For example, for a subject with n measurements, there are n rows for this subject. The # of rows in y matrix is the total number of measurements for all subjects. The columns in Y are ordered this way: the longitudinal outcome (column 1), then the covariates for random effects, and lastly, the covariates for fixed effects (no intercept).
##' @param cfile C matrix for competing risks failure time data. Each subject has one data entry, so the number of rows equals to the number of subjects. The survival / censoring time is included in the first column, and the failure type coded as 0 (censored events), 1 (risk 1), or 2 (risk 2) is given in the second column. Two competing risks are assumed. The covariates are included in the third column and on.
##' @param mfile M vector to indicate the number of longitudinal measurements per subject. The number of rows equals to the number of subjects.
##' @param point Quadrature points used in the EM procedure. Default is 20.
##' @param maxiterations Maximum values of iterations. Default is 100000.
##' @param do.trace Print detailed information of each iteration. Default is false, not to print the iteration details.
##' @param type_file Types of inputs. Default is true, i.e.  data files with headers. If set to "F", inputs are changed to data matrixes or data.frames (with headers)
##' @return Object of class \code{JMcmprsk} with elements
##'   \tabular{ll}{
##'       \code{vcmatrix}    \tab  The variance-covariance matrix for all the parameters. The parameters are in the order: \eqn{\beta}, \eqn{\alpha}, \eqn{\theta}, \eqn{\gamma}, \eqn{\nu}, and \eqn{\Sigma}. The elements in \eqn{\Sigma} are output in the order along the main diagonal line, then the second main diagonal line, and so on. \cr
##'       \code{betas} \tab The point  estimates of \eqn{\beta}. \cr
##'       \code{se_betas} \tab The standard error estimate of \eqn{\beta}. \cr
##'       \code{alphamatrix} \tab The point  estimates of \eqn{\alpha}. \cr
##'       \code{se_alphas} \tab The standard error estimate of \eqn{\alpha}. \cr
##'       \code{theta} \tab The point  estimates of \eqn{\theta}. \cr 
##'       \code{se_theta} \tab The standard error estimate of \eqn{\theta}. \cr
##'       \code{gamma_matrix} \tab  The point  estimate of \eqn{\gamma}. \cr
##'       \code{se_gamma_matrix}   \tab  The standard error estimate of \eqn{\gamma}. \cr
##'       \code{v_estimate} \tab The point  estimate of \eqn{\nu}. \cr
##'       \code{se_v_estimate}    \tab The standard error estimate of \eqn{\nu}. \cr
##'       \code{sigma_matrix}     \tab The point estimate of \eqn{\Sigma} (only the upper triangle portion of the matrix is output).\cr
##'       \code{se_sigma}     \tab The standard error estimate of \eqn{\Sigma}.The standard errors are given in this order: main diagonal, the second main diagonal, and so on. \cr
##'       \code{loglike}     \tab Log Likelihood.\cr
##'   }
##' @examples
##'require(JMcmprsk)
##'set.seed(123)
##'
##'# A toy example on a dataset called from file paths
##'yfn=system.file("extdata", "jmosimy.txt", package = "JMcmprsk")
##'cfn=system.file("extdata", "jmosimc.txt", package = "JMcmprsk")
##'mfn=system.file("extdata", "jmosimm.txt", package = "JMcmprsk")
##'fit <- jmo_0(p=3,s=1, yfn,cfn,mfn,point=6,do.trace = FALSE)
##'fit
##'
##'\dontrun{
##'# A toy example on a dataset called from data frame
##'data(ninds)
##'yread <- ninds[, c(2:14)]
##'mread <- as.data.frame(table(ninds$ID))
##'mread <- as.data.frame(mread[, 2])
##'cread <- ninds[, c(1, 15, 16, 6, 10:14)]
##'cread <- unique(cread)
##'cread <- cread[, -1]
##'jmofit=jmo_0(p=9,s=2, yread,cread,mread,point=6,do.trace = FALSE, type_file = FALSE)
##'jmofit
##'}

##' @references
##' \itemize{
##' \item Ning Li,Robert M. Elashoff,Gang Li and Jeffrey Saver. "Joint modeling of longitudinal ordinal data and competing risks survival times and analysis of the NINDS rt-PA stroke trial." Statistics in medicine 29.5 (2010): 546-557.
##' }
##' @seealso \code{\link{jmc_0}}
##' @export
jmo_0 <- function (p,s,yfile,cfile,mfile,point=20,maxiterations=100000,do.trace=FALSE,type_file=TRUE)
{
 
  # more error control here.
  
  
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
  
  # store header names for future useage
  if (type_file){
  ydata=read.table(yfile,header = TRUE)
  ynames=colnames(ydata)
  yfile=tempfile(pattern = "", fileext = ".txt")
  writenh(ydata,yfile)
 
  cdata=read.table(cfile,header = TRUE)
  cnames=colnames(cdata)
  cfile=tempfile(pattern = "", fileext = ".txt")
  writenh(cdata,cfile)

  ydim=dim(ydata)
  
  # number of observations in study is equals to the #of rows in Y matrix and delete header here
  n1=ydim[1];
  
  ydata[,1]=factor(ydata[,1])
  #generate data column names for further useage
  y_names=colnames(ydata)[1]
  fixed_col_names=paste(names(ydata[,(ncol(ydata)-p+1):ncol(ydata)]), collapse='+')
  initvalues=MASS::polr(as.formula(paste(y_names,"~",fixed_col_names)),data=ydata)
  betas=initvalues$coefficients
  thetas=initvalues$zeta
  
  # the levels of the first column of yfile
  K_num=length(unique(ydata[,1]))
  # dim of fixed effects plus dim of random effects should be 
  # total the column of y -the survival time column 1
  
  # dim of random effects
  p1a=ydim[2]-1-p-s;
  
  #if((p1<1)|(p1a<1)){
   # stop("Possibe wrong dimension of fixed effects in Y!")
  #}
  
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
  # The dimension of fixed effects in C
  p2=cdim[2]-2;
  
  PropComp <- round(table(cdata[, 2])/k * 100, 2)

  #The max number observations for a subject
  j_max=max(read.table(mfile));
  myresult=jmo_main(k,n1, p,p2, p1a,s,K_num, j_max, point,xs,ws,betas,thetas, maxiterations, yfile,cfile,mfile,trace)
 
  }else{
  #in this case yfile=ydata  
  ynames=colnames(yfile)
  yfilenew=tempfile(pattern = "", fileext = ".txt")
  writenh(yfile,yfilenew)
  cnames=colnames(cfile)
  cfilenew=tempfile(pattern = "", fileext = ".txt")
  writenh(cfile,cfilenew)
  
  mfilenew=tempfile(pattern = "", fileext = ".txt")
  writenh(mfile,mfilenew)
  
  
  
  ydim=dim(yfile)
  
  # number of observations in study is equals to the #of rows in Y matrix and delete header here
  n1=ydim[1];
  
  yfile[,1]=factor(yfile[,1])
  #generate data column names for further useage
  y_names=colnames(yfile)[1]
  fixed_col_names=paste(names(yfile[,(ncol(yfile)-p+1):ncol(yfile)]), collapse='+')
  initvalues=MASS::polr(as.formula(paste(y_names,"~",fixed_col_names)),data=yfile)
  betas=initvalues$coefficients
  thetas=initvalues$zeta
  
  # the levels of the first row of yfile
  K_num=length(unique(yfile[,1]))
  # dim of fixed effects plus dim of random effects should be 
  # total the column of y -the survival time column 1
  
  # dim of random effects
  p1a=ydim[2]-1-p-s;
  
  #if((p1<1)|(p1a<1)){
  # stop("Possibe wrong dimension of fixed effects in Y!")
  #}
  
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
  # The dimension of fixed effects in C
  p2=cdim[2]-2;
  
  PropComp <- round(table(cfile[, 2])/k * 100, 2)
  #The max number observations for a subject
  j_max=max(mfile);
  myresult=jmo_main(k,n1, p,p2, p1a,s,K_num, j_max, point,xs,ws,betas,thetas, maxiterations, yfilenew,cfilenew,mfilenew,trace)
  
}
  


  
  
#the program is estimating -beta, -alpha, and -bi in equation (1) of the stats #in med paper.when reporting the results, the sign of beta, alpha (which is #beta2 in the program), and rho_bu (or sigma_bu, which is the off-diagonal  #elements of the sig matrix) should be flipped (i.e., negative values should be #positive, and vice versa)
  #beta and alpha
  myresult$betas=-myresult$betas
  myresult$alphamatrix=-myresult$alphamatrix
  
  #names
  names(myresult$betas)=ynames[(ydim[2]-p+1):ydim[2]]
  colnames(myresult$alphamatrix)=ynames[3:(s+3-1)]
  colnames(myresult$gamma_matrix)=cnames[3:(p2+3-1)]
  
  
  #off-diagnoal elements
  
  for (i in 1:(dim(myresult$sigma_matrix)[1] - 1)) 
    for (j in (i +1):(dim(myresult$sigma_matrix)[2])) 
      {
      myresult$sigma_matrix[i,j]=-myresult$sigma_matrix[i,j]
    }
  
  myresult$k=k
  myresult$type="jmo";
  
  
  ##create longitudinal submodel formula
  #ynames <- colnames(ydata)
  LongOut <- ynames[1]
  LongP <- paste0(ynames[(ydim[2]-p+1):ydim[2]], collapse = "+")
  LongNP <- paste0(ynames[3:(s+3-1)], collapse = " + ")
  FunCall_long <- as.formula(paste(LongOut, LongP, sep = "~"))
  
  
  ##create survival submodel formula
  #cnames <- colnames(cdata)
  SurvOut <- paste0("Surv(", cnames[1], ",", cnames[2], ")")
  SurvX <- paste0(cnames[-(1:2)], collapse = "+")
  FunCall_survival <- as.formula(paste(SurvOut, SurvX, sep = "~"))
  
  DataPath <- NULL
  SummaryInfo <- list(k, n1, PropComp, FunCall_long, FunCall_survival, DataPath, LongNP)
  names(SummaryInfo) <- c("NumSub", "Numobs", "PropComp", 
                          "LongitudinalSubmodel", "SurvivalSubmodel", "DataPath", "LongNP")
  
  myresult$SummaryInfo <- SummaryInfo
  myresult$point <- point
  
  mycall=match.call()
  myresult$call=mycall
  
  class(myresult) <- "JMcmprsk"
  return (myresult)
}




