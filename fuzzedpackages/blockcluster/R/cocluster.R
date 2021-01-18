#' @include coclusterStrategy.R
#' @include optionclasses.R
#' 
NULL

#' Co-Clustering function.
#' 
#' This function performs Co-Clustering (simultaneous clustering of rows and columns ) for Binary, Contingency
#' and Continuous data-sets using latent block models.It can also be used to perform semi-supervised co-clustering.  
#' 
#' @param data Input data as matrix (or list containing data matrix, numeric vector for row effects and numeric 
#'        vector column effects in case of contingency data with known row and column effects.)
#' @param datatype This is the type of data which can be "binary" , "contingency", "continuous" or "categorical".
#' @param semisupervised Boolean value specifying whether to perform semi-supervised co-clustering or not. Make sure to provide row and/or
#' column labels if specified value is true. The default value is false.
#' @param rowlabels Integer Vector specifying the class of rows. The class number starts from zero. Provide -1 for unknown row class. 
#' @param collabels Integer Vector specifying the class of columns. The class number starts from zero. Provide -1 for unknown column class.
#' @param model This is the name of model. The following models exists for various types of data:
#' \tabular{rlll}{
#'     Model  \tab Data-type \tab Proportions \tab Dispersion/Variance \cr
#'     pik_rhol_epsilonkl(Default) \tab binary \tab unequal \tab unequal \cr
#'     pik_rhol_epsilon \tab binary \tab unequal \tab equal \cr
#'     pi_rho_epsilonkl \tab binary \tab equal \tab unequal \cr
#'     pi_rho_epsilon \tab binary \tab equal \tab equal \cr
#'     pik_rhol_sigma2kl(Default) \tab continuous \tab unequal \tab unequal \cr
#'     pik_rhol_sigma \tab continuous \tab unequal \tab equal \cr
#'     pi_rho_sigma2kl \tab continuous \tab equal \tab unequal \cr
#'     pi_rho_sigma2 \tab continuous \tab equal \tab equal \cr
#'     pik_rhol_unknown(default) \tab contingency \tab unequal \tab N.A \cr
#'     pi_rho_unknown \tab contingency \tab equal \tab N.A \cr
#'     pik_rhol_known \tab contingency \tab unequal \tab N.A \cr
#'     pi_rho_known \tab contingency \tab equal \tab N.A \cr
#'     pik_rhol_multi \tab categorical \tab unequal \tab unequal \cr
#'     pi_rho_multi \tab categorical \tab equal \tab unequal \cr
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
# @useDynLib blockcluster
#' 
#' @examples
#' 
#' # Simple example with simulated binary data
#' #load data
#' data(binarydata)
#' #usage of cocluster function in its most simplest form
#' out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3))
#' #Summarize the output results
#' summary(out)
#' #Plot the original and Co-clustered data 
#' plot(out)
#' 
cocluster<-function( data, datatype
                   , semisupervised = FALSE
                   , rowlabels = integer(0), collabels = integer(0)
                   , model = NULL, nbcocluster, strategy = coclusterStrategy()
							     , nbCore =1
							     ) 
{
	#check for datatype and models and create input object to be passed in .Call function.
	if (missing(datatype))
	{ stop("In cocluster, mention datatype.")} 
	else
	{
    if(datatype == "binary")
	  {
      inpobj<-coclusterBinary( data, semisupervised, rowlabels, collabels
												     , model, nbcocluster, strategy, nbCore)
	  }
    else if(datatype == "continuous")
	  {
      inpobj<-coclusterContinuous( data, semisupervised, rowlabels, collabels
								                 , model, nbcocluster, strategy, nbCore) 
	  }	
    else if(datatype == "contingency")
	  {
		  inpobj<-coclusterContingency( data, semisupervised, rowlabels, collabels
									                , model, nbcocluster, strategy, nbCore)
	  }
    else if(datatype == "categorical")
	  {
		  inpobj<-coclusterCategorical( data, semisupervised, rowlabels, collabels
                                  , model, nbcocluster, strategy, nbCore)
	  }
  }
  return(inpobj)
}


