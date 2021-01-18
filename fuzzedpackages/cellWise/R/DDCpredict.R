
DDCpredict = function(Xnew,InitialDDC,DDCpars=NULL){
  #
  # This function performs the DDC prediction on Xnew.
  #
  # The inputs are:
  #
  # Xnew       : the new data (test data), which must be a 
  #              matrix or a data frame. 
  #              It must always be provided.
  # InitialDDC : the output of the DDCcore function on the 
  #              initial (training) dataset. Must be provided.
  # DDCpars    : the input options to be used for the prediction.
  #              By default the options of InitialDDC are used:
  
  if (is.null(DDCpars)) { 
    DDCpars = InitialDDC$DDCpars
  } else {
    if (!is.list(DDCpars)) {
      stop("DDCpars must be a list")
    }
    InitialDDC$DDCpars[names(DDCpars)] <- DDCpars 
    DDCpars <- InitialDDC$DDCpars
  }

  
  if (!"precScale" %in% names(DDCpars)) {
    DDCpars$precScale <- 1e-12
  }
  if (!"tolProb" %in% names(DDCpars)) {
    DDCpars$tolProb <- 0.99
  }
  if (!"corrlim" %in% names(DDCpars)) {
    DDCpars$corrlim <- 0.5
  }
  if (!"combinRule" %in% names(DDCpars)) {
    DDCpars$combinRule <- "wmean"
  }
  
  
  if (!"tolProbCell" %in% names(DDCpars)) {
    DDCpars$tolProbCell <- DDCpars$tolProb
  }
  if (!"tolProbRow" %in% names(DDCpars)) {
    DDCpars$tolProbRow <- DDCpars$tolProb
  }
  if (!"includeSelf" %in% names(DDCpars)) {
    DDCpars$includeSelf <- 1
  }
  
  #
  # The outputs are:
  #
  # DDCpars     : the options used in the call.
  # locX        : the locations of the columns, from InitialDDC.
  # scaleX      : the scales of the columns, from InitialDDC.
  # Z           : Xnew standardized by locX and scaleX
  # nbngbrs           : predictions use a combination of nbngbrs columns
  # ngbrs       : for each column, the list of its neighbors,
  #               from InitialDDC.
  # robcors     : for each column, the correlations with its 
  #               neighbors, from InitialDDC.
  # robslopes   : slopes to predict each column by its neighbors,
  #               from InitialDDC.
  # deshrinkage : for each column, its deshrinkage factor used
  #               in InitialDDC. 
  # Xest        : predicted values for every cell of Xnew.
  # scalestres  : scale estimate of the residuals (Xnew - Xest),
  #               from InitialDDC.
  # stdResid    : columnwise standardized residuals of Xnew.
  # indcells    : positions of cellwise outliers in Xnew.
  # Ti          : outlyingness of rows in Xnew.
  # medTi       : median of the Ti in InitialDDC.
  # madTi       : mad of the Ti in InitialDDC.  
  # indrows     : row numbers of the outlying rows in Xnew.
  # indNAs      : positions of the NAs in Xnew.
  # indall      : positions of NA's and outlying cells in Xnew.
  # Ximp        : Xnew where all cells in indall are imputed by their
  #               prediction.
  #  
  
  ## FUNCTION FOR PREDICTING VALUES IN A COLUMN:
  
  predictCol2 = function(colj,U,ngbrs,corrweight,robslopes,combinRule){
    # Predicts the values in column colj using the set 'ngbrs' of
    # columns of U, by applying the combination rule 'combinRule' whose
    # inputs are the weights in 'corrweight' and the slopes in 'robslopes'.
    #
    # Assumes that the first entry of colj is the number of the column.
    # Assumes the remainder of colj is a vector with same length as the 
    # remainder of the columns of U, and that all of these columns are
    # already centered. 
    j    = colj[1]
    colj = colj[-1] # take out first entry
    U    = matrix(U[-1,],ncol = dim(U)[2])   
    # take out first row 
    # matrix(): in case of vector
    contributors = (corrweight[j,] > 0)
    if(length(contributors) < 1){
      estcol = 0
    } else {
      ngb1    = ngbrs[j,contributors]
      slopes1 = robslopes[j,contributors]
      corrwt1 = corrweight[j,contributors]
      ZestAllh = t(t(matrix(U[,ngb1],nrow=dim(U)[1])) * slopes1) 
      # is n by k matrix
      # Predicts column j from each column h, using slope(Z_j ~ Z_h).
      # This array has the estimates from k variables.    
      if (combinRule =="wmean"){
        estcol = apply(ZestAllh,1,weighted.mean,w=corrwt1,na.rm=TRUE)
      } 
      if (combinRule =="wmedian"){
        estcol = apply(ZestAllh,1,weightedMedian,w=corrwt1,na.rm=TRUE)  
      }
      if (combinRule =="mean"){
        estcol = apply(ZestAllh,1,mean,na.rm=TRUE) 
      }
      if (combinRule =="median"){
        estcol = apply(ZestAllh,1,median,na.rm=TRUE)  
      } 
      estcol
    }    
  }
  
  ## Here the actual computations start.
  #
  # Retrieve parameters from the list:
  #
  # precScale   = DDCpars$precScale
  if (DDCpars$tolProb < 0.5) stop("tolProb must be >= 0.5")
  if (DDCpars$tolProb >= 1)  stop("tolProb must be < 1.0")
  qCell       = sqrt(qchisq(DDCpars$tolProbCell,1))
  # = cutoff used for flagging cells.
  qRow        = sqrt(qchisq(DDCpars$tolProbRow,1)) 
  # = cutoff used for flagging rows.
  includeSelf = DDCpars$includeSelf
  combinRule  = DDCpars$combinRule
  corrlim     = DDCpars$corrlim
  uniDetect  = limitFilt  # DDCpars$uniDetect
  
  # if(combinRule == "wmedian") { library(matrixStats) }    
  # library matrixStats contains the function weightedMedian
  
  # Turn data into a matrix
  if (is.data.frame(Xnew) | is.matrix(Xnew) | is.vector(Xnew)) { 
    Xnew = data.matrix(Xnew) } else {
      stop("Data matrix must be of class matrix or data.frame") }
  if( (nrow(Xnew) > ncol(Xnew)) & ncol(Xnew)==1 )  Xnew=t(Xnew)
  n = nrow(Xnew)
  d = ncol(Xnew)
  
  if(d < 2) stop(" DDCpredict needs at least 2 columns")
  
  wnq(paste(" The computation started at: ",date(),sep=""))
  
  #####################################
  ##    STEP 1: STANDARDIZE DATA     ##
  #####################################
  
  # Robust standardization
  locX   = InitialDDC$locX
  Z      = sweep(Xnew,2,locX)
  scaleX = InitialDDC$scaleX
  Z      = sweep(Z,2,scaleX,"/")
  
  #####################################
  ##    STEP 2: UNIVARIATE ANALYSIS  ##
  #####################################  
  
  indNAs      = which(is.na(Z))
  # Univariate analysis by column: replace outliers by NAs
  U           = uniDetect(Z,qCut=qCell)
  rownames(U) = rownames(Z)
  colnames(U) = colnames(Z)    
  UniIndex    = setdiff(which(is.na(U)),indNAs) 
  # does not include original missings
  
  #####################################################
  ##    STEP 3: CALCULATE CORRELATIONS AND SLOPES    ##
  #####################################################      
  
  # k = min(d-1,1000) # -1 since we do not have to compute the
  # # correlation of a column with itself
  nbngbrs = InitialDDC$nbngbrs
  
  # For each column j of U, take the k columns h != j of U that
  # it has the highest absolute correlation robCorr with:
  
  ngbrs     = InitialDDC$ngbrs      # identities of these columns
  robcors   = InitialDDC$robcors    # the correlations
  robslopes = InitialDDC$robslopes  # the slopes
  if(includeSelf) 
    corrweight = abs(robcors)       # should have no NAs
  if (corrlim > 0) { corrweight[corrweight < corrlim] = 0 }
  
  if(!includeSelf) colStandalone = which(rowSums(corrweight)==0)
  if(includeSelf) colStandalone = which(rowSums(corrweight[,-1])==0)
  colConnected  = which(!((1:d) %in% colStandalone))
  
  indexStandalone = col2cell(colStandalone,n=n)
  indexStandalone = indexStandalone[indexStandalone %in% UniIndex]
  # = list of flagged cells in standalone variables.
  
  numiter = 1 # increase this if you want to iterate
  for (iter in 1:numiter){
    
    ####################################
    ##    STEP 4 : ESTIMATE CELLS     ##
    ####################################     
    
    Zest = U # These values will remain for standalone columns.  
    
    # Estimation for connected variables:
    U = rbind(1:d,U)
    Zest[,colConnected] = apply(U[,colConnected],2,predictCol2,U=U,
                                ngbrs=ngbrs,corrweight=corrweight,
                                robslopes=robslopes,
                                combinRule=combinRule) 
    U = U[-1,]
    
    ####################################
    ##    STEP 5 : DESHRINKAGE        ##
    ####################################
    
    # Deshrinkage: rescale Zest[,j] 
    Zest[,colConnected] = t(InitialDDC$deshrinkage * 
                              t(Zest[,colConnected]))
    
    # Finally, all NAs are replaced by zeroes:
    Zest[is.na(Zest)] = 0
    
    ####################################
    ##    STEP 6 : FLAGGING CELLS     ##
    ####################################      
    
    # Compute cell residuals:
    Zres = Z-Zest # original minus estimated
    Zres[,colStandalone] = Z[,colStandalone] 
    
    Zres[,colConnected] = scale(matrix(Zres[,colConnected],
                                       ncol=length(colConnected)),
                                center=FALSE,InitialDDC$scalestres)    
    # where the generic R-function scale() scaled these residuals.
    # We don't have to scale the standalone columns, as they
    # were already standardized in the beginning.
    
    # Next, flag outlying cells by their large residuals:
    indcells = which(abs(Zres) > qCell) 
    # does not flag the NAs as cells
    U[indcells] = NA
    
  } # ends the iteration
  rm(U) # we no longer use it.
  
  indcells = setdiff(indcells,col2cell(colStandalone,n=n))
  # the indices of outlying cells in connected variables only
  indcells = unique(sort(c(indcells,indexStandalone))) 
  # are the indices of both types of outlying cells
  
  ####################################
  ##    STEP 7 : FLAGGING ROWS      ##
  ####################################      
  
  Ti      = integer(0)
  indrows = integer(0)
  indall  = indcells
  
  compT = function(rowi){ mean(pchisq(rowi^2,1),na.rm=T) - 0.5 }
  
  
  Ti = as.vector(apply(Zres,1,FUN=compT))    
  # calculate the test value (outlyingness) of each row:     
  Ti = scale(Ti, InitialDDC$medTi, InitialDDC$madTi)
  rownames(Ti) = rownames(Z)
  indrows = which(is.na(uniDetect(Ti,qCut=qRow)) & (Ti>0))
  indall  = unique(c(indcells,row2cell(indrows,n,d))) 
  
  
  #############################################
  ##    STEP 8: UNSTANDARDIZE AND IMPUTE     ##
  #############################################
  
  # compute Xest (storing in the existing matrix Zest to save space)
  Zest = sweep(sweep(Zest,2,scaleX,"*"),2,locX,"+")
  
  # compute Ximp (storing in the existing matrix Xnew to save space)
  Xnew[indcells]  = Zest[indcells] # imputes the outlying cells of Xnew
  Xnew[indNAs] = Zest[indNAs] # imputes the missing values of Xnew
  
  attr(Zres,"scaled:scale") = NULL
  attr(Ti,"scaled:center")  = NULL
  attr(Ti,"scaled:scale")   = NULL
  
  wnq(paste(" The computation ended at: ",date(),sep=""))
  wnq(" ")
  
  
  DDCpars <- DDCpars[!(names(DDCpars) %in% c("tolProbCell", "tolProbRow", 
                                             "includeSelf"))]
  return(list(DDCpars=DDCpars,
              locX=locX,
              scaleX=scaleX,
              Z=Z,  
              nbngbrs=nbngbrs,
              ngbrs=ngbrs,
              robcors=robcors,
              robslopes=robslopes,
              deshrinkage=InitialDDC$deshrinkage,
              Xest=Zest,
              scalestres=InitialDDC$scalestres,
              stdResid=Zres,
              indcells=indcells,
              Ti=Ti,
              medTi=InitialDDC$medTi,
              madTi=InitialDDC$madTi,
              indrows=indrows,
              indNAs=indNAs,
              indall=indall,              
              Ximp=Xnew))  
} # ends DDCpredict


## AUXILIARY FUNCTIONS:


limitFilt <- function(v,qCut) {
  # Detects outliers and sets them to NA.
  # Assumes that the data have already been standardized.
  vout <- v
  vout[(abs(v) > qCut)] <- NA
  return(vout)
}

wnq <- function(string,qwrite=0){ # auxiliary function
  # writes a line without quotes
  if(qwrite==1) write(noquote(string),file="",ncolumns=100)
}


pnq <- function(string,qwrite=0){ # auxiliary function
  # prints a line without quotes
  if(qwrite==1) print(noquote(string))
} 


col2cell = function(colNrs,n) {
  # Transforms column indices to cellwise indices.
  # Here colNrs is a vector with column numbers between 1 and d.
  cindex = t(matrix(rep((colNrs-1)*n,n),ncol=n,byrow=FALSE))
  cindex = cindex + seq(1,n,1) # contains the cell numbers
  return(as.vector(cindex))
}

row2cell = function(rowNrs,n,d) {
  # Transforms row indices to cellwise indices.
  # Here rowNrs is a vector with row numbers between 1 and n.
  as.vector(t(matrix(rep(rowNrs,d),ncol=d,byrow=FALSE))+
              seq(0,n*(d-1),n))
}