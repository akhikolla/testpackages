#' @include coclusterStrategy.R
#' @include optionclasses.R
#' 
NULL

#' Co-Clustering function.
#' 
#' This function performs Co-Clustering (simultaneous clustering of rows and columns )
#' for Contingency data-sets using latent block models.It can also be used to
#' perform semi-supervised co-clustering.  
#' 
#' @param data Input data as matrix (or list containing data matrix, numeric vector for row effects and numeric 
#'        vector column effects in case of contingency data with known row and column effects.)
#' @param semisupervised Boolean value specifying whether to perform semi-supervised
#' co-clustering or not. Make sure to provide row and/or column labels if
#' specified value is true. The default value is false.
#' @param rowlabels Integer Vector specifying the class of rows. The class number starts from zero. Provide -1 for unknown row class. 
#' @param collabels Integer Vector specifying the class of columns. The class number starts from zero. Provide -1 for unknown column class.
#' @param model This is the name of model. The following models exists for Poisson data:
#' \tabular{rlll}{
#'     pik_rhol_unknown(default) \tab contingency \tab unequal \tab N.A \cr
#'     pi_rho_unknown \tab contingency \tab equal \tab N.A \cr
#'     pik_rhol_known \tab contingency \tab unequal \tab N.A \cr
#'     pi_rho_known \tab contingency \tab equal \tab N.A \cr
#' }
#' 
#' @param nbcocluster Integer vector specifying the number of row and column clusters respectively.
#' @param strategy Object of class \code{\linkS4class{strategy}}.
#' @param nbCore number of thread to use (OpenMP must be available), 0 for all cores. Default value is 1.
#' 
#' @return Return an object of \code{\linkS4class{BinaryOptions}} or \code{\linkS4class{ContingencyOptions}}
#' or \code{\linkS4class{ContinuousOptions}} depending on whether the data-type is Binary, Contingency or Continuous
#' respectively.
#' 
# @export
#' 
# @exportPattern "^[[:alpha:]]+"
# @useDynLib RCocluster
#' 
#' @examples
#' 
#' ## Simple example with simulated contingency data
#' ## load data
#' data(contingencydataunknown)
#' ## usage of coclusterContingency function in its most simplest form
#' strategy = coclusterStrategy( nbinititerations = 5, nbxem = 2, nbiterations_int = 2
#'                             , nbiterationsxem = 10, nbiterationsXEM = 100, epsilonXEM=1e-5)
#' out<-coclusterContingency( contingencydataunknown, nbcocluster=c(2,3), strategy = strategy)
#' ## Summarize the output results
#' summary(out)
#' ## Plot the original and Co-clustered data 
#' plot(out)
#' 
coclusterContingency<-function( data, semisupervised = FALSE
                              , rowlabels = integer(0), collabels = integer(0)
                              , model = NULL, nbcocluster, strategy = coclusterStrategy()
													    , nbCore = 1) 
{
	#Check for data
	if(missing(data)){ stop("Data is missing.")}
  if(!is.list(data))
  {
    if(!is.matrix(data)) { stop("Data should be matrix.")}
    dimData = dim(data)
  }
  else
  {
    if(!is.matrix(data[[1]])) { stop("Data should be matrix.")}
    if(!is.numeric(data[[2]])||!is.numeric(data[[3]]))
    { stop("Row/Column effects should be numeric vectors.")}
    if(length(data[[2]])!=dim(data[[1]])[1]||length(data[[3]])!=dim(data[[1]])[2])
    { stop("Dimension mismatch in Row/column effects  and Data.")}
    dimData = dim(data[[1]])
  }
  
    #check for row and column labels
  if (semisupervised)
  {
    if(missing(rowlabels)&&missing(collabels))
      stop("Missing row and column labels. At-least one should be provided to perform semi-supervised Co-clustering.")
    if(!missing(rowlabels)&&!is.numeric(rowlabels))
      stop("Row labels should be a numeric vector.")
    if(!missing(collabels)&&!is.numeric(collabels))
      stop("Column labels should be a numeric vector.")
    
    if(missing(rowlabels))      rowlabels = rep(-1,dimData[1])
    else if(missing(collabels)) collabels = rep(-1,dimData[2])
    
    if(dimData[1]!=length(rowlabels))
      stop("rowlabels length does not match number of rows in data (also ensure to put -1 in unknown labels)")
    
    if(dimData[2]!=length(collabels))
      stop("collabels length does not match number of columns in data (also  ensure to put -1 in unknown labels)")
  }
	#check for number of coclusters
  if(missing(nbcocluster))     { stop("Mention number of CoClusters.")}
  if(dimData[1]<nbcocluster[1]){ stop("Number of Row clusters exceeds numbers of rows.")}
  if(dimData[2]<nbcocluster[2]){ stop("Number of Column clusters exceeds numbers of columns.")}
  if(nbcocluster[1]<1 || nbcocluster[2]<1) { stop("Number of cluster must be at least 1.")}

  #check for Algorithm name (and make it compatible with version 1)
	if(strategy@algo=="XEMStrategy")
  {
    warning("The algorithm 'XEMStrategy' is renamed as BEM!")
    strategy@algo == "BEM"
  }
  else if(strategy@algo == "XCEMStrategy")
  {
    warning("The algorithm 'XCEMStrategy' is renamed as BCEM!")
    strategy@algo = "BCEM"
  }
  else if(strategy@algo!="BEM" && strategy@algo!="BCEM" && strategy@algo!="BSEM" )
  {
    stop("Incorrect Algorithm, Valide algorithms are: BEM, BCEM, BSEM") 
  }
	#check for stopping criteria
	if(strategy@stopcriteria!="Parameter" && strategy@stopcriteria!="Likelihood")
		stop("Incorrect stopping criteria, Valid stopping criterians are: Parameter, Likelihood")
  
	#check for datatype and models and create input object to be passed in .Call function.
	if(is.null(model) && !is.list(data)){ model = "pik_rhol_unknown"}
	else if(is.null(model) && is.list(data)) {	model = "pik_rhol_known"}
	else
  {
    if(model!="pik_rhol_unknown" && model!="pik_rhol_known" && 
			model!="pi_rho_unknown" && model!="pi_rho_known")
    {
		  stop("Incorrect Model, Valid Contingency models are:pik_rhol_unknown, pik_rhol_known, pi_rho_unknown, pi_rho_known")
    }
  }
	if((model=="pi_rho_known"||model=="pik_rhol_known")&& (length(data)!=3))
	{stop("Missing Row/Column effects.")}
	
	if(length(strategy@initmethod)==0){ strategy@initmethod = "emInitStep"}
  ## {
  ##   if((model=="pi_rho_known"||model=="pik_rhol_known"))
  ##   { strategy@initmethod = "emInitStep"}
  ##   else
  ##   { strategy@initmethod = "emInitStep"}
  ## }
  ## else
  ## {
  ##    if(strategy@initmethod!="randomInit"&&(model=="pi_rho_known"||model=="pik_rhol_known"))
  ##    { stop("Incorrect initialization method, valid method(s) are: randomInit")}
  ##    else if((strategy@initmethod!="cemInitStep" && strategy@initmethod!="emInitStep")&&(model=="pi_rho_unknown"||model=="pik_rhol_unknown"))
  ##      stop("Incorrect initialization method, valid method(s) are: cemInitStep, emInitStep")
  ## }
	#  check nbCore
	if(!is.numeric(nbCore) && length(nbCore) != 1) stop("nbCore must be an integer")
	# check options
	if(!is.list(data))
        inpobj<-new( "ContingencyOptions", data = data
                   , semisupervised = semisupervised
                   , rowlabels = rowlabels, collabels = collabels
                   , datatype = "contingency", model = model, nbcocluster = nbcocluster, strategy = strategy)
  else
     inpobj<-new( "ContingencyOptions",data = data[[1]]
                , datamui=data[[2]],datanuj=data[[3]]
                , semisupervised = semisupervised
                , rowlabels = rowlabels, collabels = collabels
                , datatype = "contingency", model = model, nbcocluster = nbcocluster, strategy = strategy)

  .Call("CoClustmain",inpobj, nbCore,PACKAGE = "blockcluster")
  cat(inpobj@message,"\n")
  return(inpobj)
}

