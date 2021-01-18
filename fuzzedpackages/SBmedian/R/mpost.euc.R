#' Median Posterior for Subset Posterior Samples in Euclidean Space
#' 
#' \code{mpost.euc} is a general framework to \emph{merge} multiple 
#' empirical measures \eqn{Q_1,Q_2,\ldots,Q_M \subset R^p} from independent subset of data by finding a median 
#' \deqn{\hat{Q} = \textrm{argmin}_Q \sum_{m=1}^M d(Q,Q_m)}
#' where \eqn{Q} is a weighted combination and \eqn{d(P_1,P_2)} is distance in RKHS between two empirical measures \eqn{P_1} and \eqn{P_2}. 
#' As in the references, we use RBF kernel with bandwidth parameter \eqn{\sigma}.
#' 
#' @param splist a list of length \eqn{M} containing vectors or matrices of univariate or multivariate subset posterior samples respectively. 
#' @param sigma bandwidth parameter for RBF kernel.
#' @param maxiter maximum number of iterations for Weiszfeld algorithm.
#' @param abstol stopping criterion for Weiszfeld algorithm.
#' @param show.progress a logical; \code{TRUE} to show iteration mark, \code{FALSE} otherwise.
#' 
#' @return a named list containing:
#' \describe{
#' \item{med.atoms}{a vector or matrix of all atoms aggregated.}
#' \item{med.weights}{a weight vector that sums to 1 corresponding to \code{med.atoms}.}
#' \item{weiszfeld.weights}{a weight for \eqn{M} subset posteriors.}
#' \item{weiszfeld.history}{updated parameter values. Each row is for iteration, while columns are weights corresponding to \code{weiszfeld.weights}.}
#' }
#' 
#' @examples 
#' ## Median Posteior from 2-D Gaussian Samples
#' #  Step 1. let's build a list of atoms whose numbers differ
#' set.seed(8128)                   # for reproducible results
#' mydata = list()
#' mydata[[1]] = cbind(rnorm(96, mean= 1), rnorm(96, mean= 1))
#' mydata[[2]] = cbind(rnorm(78, mean=-1), rnorm(78, mean= 0))
#' mydata[[3]] = cbind(rnorm(65, mean=-1), rnorm(65, mean= 1))
#' mydata[[4]] = cbind(rnorm(77, mean= 2), rnorm(77, mean=-1))
#' 
#' #  Step 2. Let's run the algorithm
#' myrun = mpost.euc(mydata, show.progress=TRUE)
#' 
#' #  Step 3. Visualize
#' #  3-1. show subset posterior samples
#' opar <- par(mfrow=c(2,3), no.readonly=TRUE)
#' for (i in 1:4){
#'   plot(mydata[[i]], cex=0.5, col=(i+1), pch=19, xlab="", ylab="", 
#'        main=paste("subset",i), xlim=c(-4,4), ylim=c(-3,3))
#' }
#' 
#' #  3-2. 250 median posterior samples via importance sampling
#' id250 = base::sample(1:nrow(myrun$med.atoms), 250, prob=myrun$med.weights, replace=TRUE)
#' sp250 = myrun$med.atoms[id250,]
#' plot(sp250, cex=0.5, pch=19, xlab="", ylab="", 
#'      xlim=c(-4,4), ylim=c(-3,3), main="median samples")
#' 
#' #  3-3. convergence over iterations
#' matplot(myrun$weiszfeld.history, xlab="iteration", ylab="value",
#'         type="b", main="convergence of weights")
#' par(opar)
#'         
#' @references 
#' \insertRef{minsker_scalable_2014}{SBmedian}
#' 
#' \insertRef{minsker_robust_2017}{SBmedian}
#' 
#' @export
mpost.euc <- function(splist, sigma = 0.1, maxiter = 121, abstol = 1e-6, show.progress = FALSE){
  ##-----------------------------------------------------------------------------------------------
  ## Check the input
  
  if ((!is.list(splist))||(length(splist)<2)){
    stop(" * mpost.euc : 'splist' should be a LIST of length larger than 1.")
  }
  if (inherits(splist[[1]], "vector")){
    check_vector(splist)
    for (i in 1:length(splist)){
      splist[[i]] = matrix(splist[[i]], ncol = 1) # transform to matrices
    }
    vflag = TRUE
  } else if (inherits(splist[[1]], "matrix")){
    check_matrix(splist)
    vflag = FALSE
  } else {
    stop(" * mpost.euc : elements in 'splist' should ALL be either VECTORS or MATRICES.")
  }
  
  ##-----------------------------------------------------------------------------------------------
  ## Run the main code
  mysigma    = as.double(sigma)
  mymaxiter  = round(maxiter)
  myabstol   = as.double(abstol)
  medMeasure = engine_main(splist, mysigma, mymaxiter, myabstol, show.progress, "mpost.euc")
  
  ##-----------------------------------------------------------------------------------------------
  ## Manipulating the returned output
  output = list()
  # 1. med.atoms : median atoms which are collection of all atoms
  if (vflag){
    output$med.atoms = as.vector(medMeasure$medianAtoms) # case : numbers in a vector
  } else {
    output$med.atoms = medMeasure$medianAtoms            # case : vectors row-stacked matrix
  }
  # 2. med.weights
  natoms = as.vector(medMeasure$natoms)
  med.weights = c()
  for (i in 1:length(natoms)){
    med.weights = c(med.weights, rep(1/natoms[i], natoms[i])*medMeasure$weiszfeldWts[i])
  }
  output$med.weights = med.weights
  # 3. weiszfeld.weights
  output$weiszfeld.weights = as.vector(medMeasure$weiszfeldWts)
  # 4. weiszfeld.history
  output$weiszfeld.history = t(medMeasure$historyWeiszfeldWts)
  rownames(output$weiszfeld.history) = paste("iteration",1:ncol(medMeasure$historyWeiszfeldWts))
  # and Return
  return(output)
}