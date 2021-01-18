#' generateSimulationData
#'
#' @param sigma noise level - default is 1
#' @param sparsity percentage of relevent features - default is 0.02
#' @param seed Random seed - default is 1234
#'
#' @return Returns simulation data - a 21000x1000 sparse matrix with 4 clusters (2 clusters are n = 10,000 and 2 clusteraer n = 500)
#' @export
#'
#' @examples
#' X = generateSimulationData()
generateSimulationData = function( sigma = 1, sparsity = 0.02, seed = 1234 )
{
  # size of various clusters
  n1 = 500
  n2 = 500
  n3 = 300
  n4 = 300
  
  n = c( n1, n2, n3, n4 ) # Unequal Cluster size settings
  p  = 1000  # first 20 being informative and remaining ones are non-informative
  C = 4 # true numebr of clusters
  N = sum( n )
  
  # generate sparsity and matrix
  set.seed( seed )
  X = matrix( rnorm( N * p, 0, sigma ), nrow = N, ncol = p, byrow = TRUE )
  dsparse = 0.5
  p0 = dsparse * sparsity * p
  p1 = ( 1 - dsparse ) * sparsity * p
  
  #Cluster 1: N( 2.5, sigma )( 1-p0 ), N( 1.5, sigma )( p0-( p0 + p1 ) )
  X[1:n[1], 1:p0] = rnorm( n[1] * p0, 2.5, sigma )
  X[1:n[1], ( 1 + p0 ):( p0 + p1 )] = rnorm( n[1] * p1, 1.5, sigma )
  
  #Cluster 2: N( 0, sigma )( 1-10 ), N( 1.5, sigma ) ( 11-20 )
  X[( n[1] + 1 ):( n[1] + n[2] ),1:p0] = rnorm( n[2] * p0, 0, sigma )
  X[( n[1] + 1 ):( n[1] + n[2] ),( 1 + p0 ):( p0 + p1 )] = rnorm( n[2] * p1, 1.5, sigma )
  
  #Cluster 3: N( 0, sigma )( 1-10 ), N( -1.5,sigma )( 11-20 )
  X[( n[1] + n[2] + 1 ):( n[1] + n[2] + n[3] ),1:p0] = rnorm( n[3] * p0,  0, sigma )
  X[( n[1] + n[2] + 1 ):( n[1] + n[2] + n[3] ),( 1 + p0 ):( p0 + p1 )] = rnorm( n[3] * p1, -1.5, sigma )
  #
  #Cluster 4: N( -2.5,sigma )( 1-10 ), N( -1.5, sigma )( 11-20 )
  X[( n[1] + n[2] + n[3] + 1 ):( n[1] + n[2] + n[3] + n[4] ),1:p0] = rnorm( n[4] * p0, -2.5, sigma )
  X[( n[1] + n[2] + n[3] + 1 ):( n[1] + n[2] + n[3] + n[4] ),( 1 + p0 ):( p0 + p1 )] = rnorm( n[4] * p1, -1.5, sigma )
  
  Y = c( rep( 1, n1 ), rep( 2, n2 ), rep( 3, n3 ), rep( 4, n4 ) )
  
  return( list( X = X, Y = Y ) )
}
