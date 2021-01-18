#' Median Posterior for Subset Posterior Samples in SPD manifold
#' 
#' SPD manifold is a collection of matrices that are symmetric and positive-definite and 
#' it is well known that using Euclidean geometry for data on the manifold is rather inaccurate. 
#' Here, we propose a function for dealing with SPD matrices specifically where valid examples include 
#' full-rank covariance and precision matrices. Note that \eqn{N_M = \sum_{m=1}^M n_m}.
#' 
#' @param splist a list of length \eqn{M} containing \eqn{(p\times p)} matrix or 3d array of size \eqn{(p\times p\times n_m)} whose slices are SPD matrices from subset posterior samples respectively. 
#' @param sigma bandwidth parameter for RBF kernel.
#' @param maxiter maximum number of iterations for Weiszfeld algorithm.
#' @param abstol stopping criterion for Weiszfeld algorithm.
#' @param show.progress a logical; \code{TRUE} to show iteration mark, \code{FALSE} otherwise.
#' 
#' @return a named list containing:
#' \describe{
#' \item{med.atoms}{a \eqn{(p\times p\times N_M)} 3d array whose slices are atoms aggregated.}
#' \item{med.weights}{a weight vector that sums to 1 corresponding to \code{med.atoms}.}
#' \item{weiszfeld.weights}{a weight for \eqn{M} subset posteriors.}
#' \item{weiszfeld.history}{updated parameter values. Each row is for iteration, while columns are weights corresponding to \code{weiszfeld.weights}.}
#' }
#' 
#' @examples 
#' ## Median Posteior from 5-dimension Wishart distribution
#' ## Visualization will be performed for distribution of larget eigenvalue
#' ## where RED is for estimated density and BLUE is density from all samples.
#' 
#' #  Step 1. let's build a list of atoms whose numbers differ
#' set.seed(8128)                   # for reproducible results
#' mydata = list()
#' mydata[[1]] = stats::rWishart(96, df=10, Sigma=diag(5))
#' mydata[[2]] = stats::rWishart(78, df=10, Sigma=diag(5))
#' mydata[[3]] = stats::rWishart(65, df=10, Sigma=diag(5))
#' mydata[[4]] = stats::rWishart(77, df=10, Sigma=diag(5))
#' 
#' #  Step 2. Let's run the algorithm
#' myrun = mpost.spd(mydata, show.progress=TRUE)
#' 
#' #  Step 3. Compute largest eigenvalues for the samples
#' eig4 = list()
#' for (i in 1:4){
#'   spdmats = mydata[[i]]        # SPD atoms
#'   spdsize = dim(spdmats)[3]    # number of atoms
#'   eigvals = rep(0,spdsize)     # compute largest eigenvalues
#'   for (j in 1:spdsize){
#'     eigvals[j] = max(base::eigen(spdmats[,,j])$values)
#'   }
#'   eig4[[i]] = eigvals
#' }
#' eigA   = unlist(eig4)
#' eiglim = c(min(eigA), max(eigA))
#' 
#' #  Step 4. Visualize
#' #  4-1. show distribution of subset posterior samples' eigenvalues
#' opar <- par(mfrow=c(2,3), no.readonly=TRUE)
#' for (i in 1:4){
#'   hist(eig4[[i]], main=paste("subset", i), xlab="largest eigenvalues", 
#'        prob=TRUE, xlim=eiglim, ylim=c(0,0.1))
#'   lines(stats::density(eig4[[i]]), lwd=1, col="red")
#'   lines(stats::density(eigA),      lwd=1, col="blue")
#' }
#' 
#' #  4-2. 250 median posterior samples via importance sampling
#' id250 = base::sample(1:length(eigA), 250, prob=myrun$med.weights, replace=TRUE)
#' sp250 = eigA[id250]
#' hist(sp250, main="median samples", xlab="largest eigenvalues", 
#'      prob=TRUE, xlim=eiglim, ylim=c(0,0.1))
#' lines(stats::density(sp250), lwd=1, col="red")
#' lines(stats::density(eigA),  lwd=1, col="blue")
#' 
#' #  4-3. convergence over iterations
#' matplot(myrun$weiszfeld.history, xlab="iteration", ylab="value",
#'         type="b", main="convergence of weights")
#' par(opar)
#' 
#' @export
mpost.spd <- function(splist, sigma = 0.1, maxiter = 121, abstol = 1e-6, show.progress = FALSE){
  ##-----------------------------------------------------------------------------------------------
  ## Check the input
  if ((!is.list(splist))||(length(splist)<2)){
    stop(" * mpost.spd : 'splist' should be a LIST of length larger than 1.")
  }
  datlist = list()
  for (i in 1:length(splist)){
    tgtspd = splist[[i]]
    if ((length(dim(tgtspd))==2)&&(inherits(tgtspd, "matrix"))){ # single matrix
      datlist[[i]] = matrix(return_spd_single(tgtspd), nrow = 1)
    } else if ((length(dim(tgtspd))==3)&&(inherits(tgtspd, "array"))){
      datlist[[i]] = return_spd_multiple(tgtspd)
    } else {
      stop(paste0("* mpost.spd : ",i,"-th element in the given list should have either SPD matrix or 3d array of SPD matrices."))
    }
  }
  if (length(unique(unlist(lapply(datlist, ncol))))!=1){
    stop("* mpost.spd : every element in the given list should have same dimension.")
  }
  
  ##-----------------------------------------------------------------------------------------------
  ## Run the main code
  mysigma    = as.double(sigma)
  mymaxiter  = round(maxiter)
  myabstol   = as.double(abstol)
  medMeasure = engine_main(datlist, mysigma, mymaxiter, myabstol, show.progress, "mpost.spd")
  
  ##-----------------------------------------------------------------------------------------------
  ## Manipulating the returned output
  output = list()
  # 1. med.atoms : median atoms which are collection of all atoms
  p = sqrt(ncol(medMeasure$medianAtoms))
  N = nrow(medMeasure$medianAtoms)
  med.atoms = array(0,c(p,p,N))
  for (n in 1:N){
    tgt = matrix(medMeasure$medianAtoms[n,], nrow=p)
    tgt = (tgt + t(tgt))/2
    med.atoms[,,n] = expm::expm(tgt)
  }
  output$med.atoms = med.atoms
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


# auxiliary functions for SPD cases ---------------------------------------
#' @keywords internal
#' @noRd
return_spd_single <- function(mat){
  cond1 = (nrow(mat)==ncol(mat))
  cond2 = isSymmetric(mat)
  cond3 = (min(base::eigen(mat, symmetric = TRUE, only.values=TRUE)$values) > 0)
  if (cond1&&cond2&&cond3){
    return(as.vector(expm::logm(mat)))
  }
}
#' @keywords internal
#' @noRd
return_spd_multiple <- function(array3){
  N = dim(array3)[3]
  output = c()
  for (n in 1:N){
    output = rbind(output, return_spd_single(array3[,,n]))
  }
  return(output)
}