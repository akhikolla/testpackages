#' Multi Resolution Scanning for one-way ANDOVA using the multi-scale Beta-Binomial model
#'
#' This function executes the Multi Resolution Scanning algorithm to detect differences 
#' across the distributions of multiple groups having multiple replicates.
#'
#' @param X Matrix of the data. Each row represents an observation.
#' @param G Numeric vector of the group label of each observation. Labels are integers starting from 1. 
#' @param H Numeric vector of the replicate label of each observation. Labels are integers starting from 1. 
#' @param n_groups Number of groups.
#' @param n_subgroups Vector indicating the number of replicates for each grop.
#' @param Omega Matrix defining the vertices of the sample space. 
#' The \code{"default"} option defines a hyperrectangle containing all the data points.
#' Otherwise the user can define a matrix  where each row represents a dimension,  and the two columns contain the associated lower and upper limit.
#' @param K Depth of the tree. Default is \code{K = 6}, while the maximum is \code{K = 14}.
#' @param init_state Initial state of the hidden Markov process.
#' The three states are \emph{null}, \emph{altenrative} and \emph{prune}, respectively.
#' @param beta Spatial clustering parameter of the transition probability matrix. Default is \code{beta = 1.0}.
#' @param gamma Parameter of the transition probability matrix. Default is \code{gamma = 0.07}.
#' @param delta Parameter of the transition probability matrix. Default is \code{delta = 0.4}.
#' @param eta Parameter of the transition probability matrix. Default is \code{eta = 0.0}.
#' @param alpha Pseudo-counts of the Beta random probability assignments.
#' @param nu_vec The support of the discrete uniform prior on nu.
#' @param return_global_null Boolean indicating whether to return the marginal posterior probability of the global null.
#' @param return_tree Boolean indicating whether to return the posterior representative tree. 
#' @return An \code{mrs} object. 
#' @references Ma L. and Soriano J. (2016). 
#' Analysis of distributional variation through multi-scale Beta-Binomial modeling. 
#'  \emph{arXiv}. 
#'  \url{http://arxiv.org/abs/1604.01443}
#' @export
#' @examples
#' set.seed(1) 
#' n = 1000
#' M = 5
#' class_1 = sample(M, n, prob= 1:5, replace=TRUE  )
#' class_2 = sample(M, n, prob = 5:1, replace=TRUE )
#' 
#' Y_1 = rnorm(n, mean=class_1, sd = .2)
#' Y_2 = rnorm(n, mean=class_2, sd = .2)
#' 
#' X = matrix( c(Y_1, Y_2), ncol = 1)
#' G = c(rep(1,n),rep(2,n))
#' H = sample(3,2*n, replace = TRUE  )
#' 
#' ans = andova(X, G, H)
#' ans$PostGlobNull
#' plot1D(ans)
andova <- function(  X, 
                         G, 
                         H,
                         n_groups = length(unique(G)), 
                         n_subgroups = NULL,
                         Omega = "default", 
                         K = 6, 
                         init_state = c(0.8,0.2,0), 
                         beta = 1.0, 
                         gamma = 0.07, 
                         delta = 0.4,
                         eta = 0, 
                         alpha = 0.5,
                         nu_vec = 10^(seq(-1,4)),
                         return_global_null = TRUE,
                         return_tree = TRUE )
{
  X = as.matrix(X)
  if(Omega[1] == "default")
  {
    Omega = t(apply(X,2,range))
    Omega[,2] = Omega[,2]*1.0001
  }
  else
  {
    EmpOmega = t(apply(X,2,range))
    if( sum((EmpOmega[,1] >= Omega[,1]) + (EmpOmega[,2] < Omega[,2])) != 2*ncol(X)  )
    {
      print("ERROR: 'X' out of bounds, check 'Omega'")
      return(0);
    }
  }
  
  if( (K > 14) || (K <= 1) )
  {
    print("ERROR: 0 < K < 15")
    return(0);
  }
  
  if( min(G) < 1 )
  {
    print("ERROR: min(G) should be positive")
    return(0);
  }
  
  
  if(is.null(n_subgroups))
  {
    n_subgroups = rep(NA, n_groups)
    for(g in 1:n_groups)
    {
      n_subgroups[g] = max(  H[which(G==g)] )
    }
  }
  
  if( sum( init_state < 0 ) > 0 | sum( init_state ) != 1 )
  {
    print("ERROR: init_state should be non-negative and sum to 1")
    return(0);
  }
  
  if( alpha<=0 )
  {
    print("ERROR: alpha > 0")
    return(0);
  }
  
  if( beta < 1 | beta > 2 )
  {
    print("ERROR: 1 <= beta <= 2")
    return(0);    
  }
  
  if( gamma < 0 | gamma > 1 )
  {
    print("ERROR: 0 <= gamma <= 1")
    return(0);
  }
  
  if( delta < 0 | delta > 1 )
  {
    print("ERROR: 0 <= delta <= 1")
    return(0);
  }
  
  ans = fitMRSNESTEDcpp( X, 
                         G, 
                         H,
                         n_groups, 
                         n_subgroups,
                         init_state, 
                         Omega, 
                         nu_vec,
                         K,
                         alpha, 
                         beta, 
                         gamma,
                         delta,
                         eta, 
                         return_global_null, 
                         return_tree )
  
  if (return_tree) {
    ans$RepresentativeTree$EffectSizes = matrix( unlist(ans$RepresentativeTree$EffectSizes), 
                                                 nrow = length(ans$RepresentativeTree$Levels), byrow = TRUE)
    ans$RepresentativeTree$Regions = matrix( unlist(ans$RepresentativeTree$Regions), 
                                             nrow = length(ans$RepresentativeTree$Levels), byrow = TRUE)
  }

  
  colnames(ans$Data$X) = colnames(X)
  
  class(ans) = "mrs"
  return(ans)
}
