# calculate contingency table
f_nij = function( v1, v2, l1, l2 )
{
  m = matrix( 0, l1, l2 )
  # creates contingency table n( i,j )=t( i,j )
  for ( i in 1:length( v1 ) )
  {
    m[v1[i], v2[i]] = m[v1[i], v2[i]] + 1
  }
  return( m )
}

# calculate estimate mutual information
f_emi = function( s1, s2, l1, l2, n )
{
  s_emi = 0
  for( i in 1:l1 )
  {
    for ( j in 1:l2 ){
      min_nij = max( 1, s1[i] + s2[j] - n )
      max_nij = min( s1[i], s2[j] )
      n.ij = seq( min_nij, max_nij )   #sequence of consecutive numbers
      t1= ( n.ij / n ) * log( ( n.ij * n ) / ( s1[i] * s2[j] ) )
      t2 = exp( lfactorial( s1[i] ) + lfactorial( s2[j] ) + lfactorial( n - s1[i] ) +
                  lfactorial( n - s2[j] ) - lfactorial( n ) - lfactorial( n.ij ) - lfactorial( s1[i] - n.ij ) -
                  lfactorial( s2[j] - n.ij ) - lfactorial( n - s1[i] - s2[j] + n.ij ) )
      emi = sum( t1 * t2 )
      s_emi = s_emi + emi
    }
  }
  return( s_emi )
}


#' Calculate ami and nmi
#'
#' @param v1 first cluster vector
#' @param v2 second cluster vector
#'
#' @return Returns ami and nmi values (a value bewteen 0 and 1) that measures how close the classification between the two clusters is. A value closer to 1 means the clusters are more similar across v1 and v2.
#' @export
#'
#' @examples
#' cluster1 <- sample(1:5, size = 10, replace = TRUE)
#' cluster2 <- sample(1:2, size = 10, replace = TRUE)
#' f_rez( cluster1, cluster2 )
f_rez = function( v1, v2 )
{
  # get number of occurance of each argument
  s1 = tabulate( v1 );
  s2 = tabulate( v2 );
  
  # get length variables
  l1 = length( s1 )
  l2 = length( s2 )
  N = length( v1 )
  
  # get contingency table n( i,j )=t( i,j ). this would be equivalent with table( v1,v2 )
  tij = f_nij( v1, v2, l1, l2 )   #
  
  #function for Mutual Information from package infotheo
  mi = infotheo::mutinformation( v1, v2 )
  
  # sum over logs of tabulations
  h1 = -sum( s1 * log( s1 / N ) ) / N
  h2 = -sum( s2 * log( s2 / N ) ) / N
  
  # calculate nmi/emi
  nmi = mi / max( h1, h2 )        # NMI Normalized MI
  nmi = mi / sqrt( h1 * h2 )
  emi = f_emi( s1, s2, l1, l2, N ) # EMI Expected MI
  ami = ( mi - emi ) / ( sqrt( h1 * h2 ) - emi )
  
  # erturn nmi and ami
  return( list( nmi = nmi,ami = ami ) )
}
