cleanFindCC = function( temp )
{
  for ( i in 1:length( temp ) )
  {
    temp[[i]] = temp[[i]][, 1]
  }

  return( temp )
}


### STEP 1
# initial clustering
initialClustering = function( X, num_cut, seed = 1 )
{
  # get dimensions
  p = ncol( X )
  N = nrow( X )

  # sample data
  set.seed( seed )
  fold_ids = sample( rep( seq_len( num_cut ), length.out = N ) )

  # if not doing parallel: d is the number of clusters, CC is an empty list
  CC = list() # if not doing parallel

  # get clusters in each cut
  for (i in 1:num_cut) 
  {
    # get relevent indices
    temp_index = which( fold_ids == i ) # get the index
    Z1 = X[temp_index, ]

    # get scaled cross product
    GG1 = tcrossprod( Z1 ) / p

    # remove the diagonal and add as a new column at the end
    gg_wodiag = GG1 - diag( diag( GG1 ) )
    cut = length( temp_index ) - 1
    GG_new = cbind( gg_wodiag + diag( colSums( gg_wodiag ) / ( cut-1 ) ), diag( GG1 ) )

    # calling actual clustering
    MclustGG1 = Mclust( GG_new, modelNames = "VVI", verbose = F )

    # if not doing parallel
    CC_i = getCCmatrix_c( MclustGG1$classification, temp_index, MclustGG1$G )
    CC_i = cleanFindCC( CC_i )
    CC = c( CC, CC_i )
  }

  return( CC = CC )
}
# 
# # initial clustering
# initialClusteringparallel = function( X, num_cut, seed = 1 )
# {
#   # get dimensions
#   p = ncol( X )
#   N = nrow( X )
#   
#   # sample data
#   set.seed( seed )
#   fold_ids = sample( rep( seq_len( num_cut ), length.out = N ) )
#   
#   # register parallel session with total number of cores - 2
#   # doParallel::registerDoParallel( cores = parallel::detectCores() - 2 )
#   
#   chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#   
#   if (nzchar(chk) && chk == "TRUE") {
#     # use 2 cores in CRAN/Travis/AppVeyor
#     num_workers <- 2L
#   } else {
#     # use all cores in devtools::test()
#     num_workers <- parallel::detectCores()
#   }
#   
#   cl = parallel::makeCluster(num_workers)
#   doParallel::registerDoParallel(cl)
#   
#   CC_loop <- foreach( i = 1:num_cut, .combine='c', .multicombine=TRUE, .packages = c( "mclust" ) ) %dopar% 
#   {
#     # get relevent indices
#     temp_index = which( fold_ids == i ) # get the index
#     Z1 = X[temp_index, ]
#     
#     # get scaled cross product
#     GG1 = tcrossprod( Z1 ) / p
#     
#     # remove the diagonal and add as a new column at the end
#     gg_wodiag = GG1 - diag( diag( GG1 ) )
#     cut = length( temp_index ) - 1
#     GG_new = cbind( gg_wodiag + diag( colSums( gg_wodiag ) / ( cut-1 ) ), diag( GG1 ) )
#     
#     # calling actual clustering
#     MclustGG1 = Mclust( GG_new, modelNames = "VVI", verbose = F )
#     return_information = list( MclustGG1$classification, temp_index, MclustGG1$G )
#   }
#   
#   parallel::stopCluster(cl)
#   
#   CC = list()
#   for ( i in seq( 1, length( CC_loop ), 3 ) )
#   {
#     CC_i = getCCmatrix_c( CC_loop[[ i ]], CC_loop[[ i + 1 ]], CC_loop[[ i + 2 ]] )
#     CC_i = cleanFindCC( CC_i )
# 
#     CC = c( CC, CC_i )
#   }
#   
#   return( CC = CC )
# }

### step 2
# get the matrix of means
getMeansMatrix = function( X, CC )
{
  d = length( CC )
  p = ncol( X )

  Cmeans = matrix( 0, nrow = d, ncol = p )
  # get matrix of means
  Cmeans = getMatrixMeans_c( CC, X, d )

  return( Cmeans )
}


### step 3
secondClustering = function( Cmeans )
{
  # get dimensions
  n = nrow( Cmeans )
  p = ncol( Cmeans )

  # reset for notational match with original code
  ZZ = Cmeans

  # getting ZZ'
  GG = tcrossprod( ZZ ) / p

  # set up new matrix for RJ algorithm
  gg = GG
  gg_wodiag = gg - diag( diag( gg ) )
  gg_wdiag = cbind( gg_wodiag, diag( gg ) )
  GG_new = cbind( gg_wodiag + diag( colSums( gg_wodiag ) / ( n-1 ) ), diag( gg ) )

  # do second clustering
  Clust.J = Mclust( GG_new, modelNames = "VVI", verbose = F )

  return( Clust.J )
}

### step 4 - assign each data point a group
assignGroups = function( N, G, classification, CC )
{
  Group = assignGroups_c( N, G, classification, CC )
  return( Group[,1] )
}

### step 5 - final clustering
finalClustering = function( G, Group, X )
{
  N = nrow( X )
  p = ncol( X )

  Lmark = getFinalMeans_c( G, Group, X )

  LL = tcrossprod( X, Lmark ) / p

  newRJ = Mclust( LL, modelNames = "VVI", G, verbose = F )

  return( newRJ )
}


## Actual function called by user if RJ scale
RJclust_backend = function( X, num_cut, seed = seed )
{
  # 1- initia lclustering of cuts of matrix
  # CC = initialClusteringparallel( X, num_cut, seed = seed )
  CC = initialClustering( X, num_cut, seed = seed )
  
  # 2- get first dxp matrix of means
  Cmeans = getMeansMatrix( X, CC )

  # 3- get second clustering based on means
  clustering = secondClustering( Cmeans )

  # 4 - assign each data point a group
  G = as.numeric( clustering[6] )
  # classification = clustering[14] # from original code, but should be 15 - classification index
  classification = clustering[15] # this is the classification
  classification = as.vector( unlist( classification[1] ) )
  
  Group = assignGroups( nrow( X ), G, classification, CC )

  # 5 - final clustering
  RJ  = finalClustering( G, Group, X )

  return( RJ )
}


# RJ function with no scaling
RJclust_noscale = function( Z, seed = 1 )
{
  set.seed( seed )
  p_temp = ncol( Z )
  n_temp = nrow( Z )
  
  # RJ steps
  gg = tcrossprod( Z )/p_temp
  gg_wodiag = gg - diag( diag( gg ) )
  # gg_wdiag = cbind( gg_wodiag, diag( gg ) )
  GG_new = cbind( gg_wodiag + diag( colSums( gg_wodiag ) / ( n_temp - 1 ) ), diag( gg ) )
  
  # par( mfrow = c( 1,1 ) )
  # image( rotate( GG_new ), col = gray.colors( 33 ) )
  
  step_one_clust = Mclust( GG_new, verbose = FALSE, G = 10 )
  return( step_one_clust )
}


#' RJclust
#'
#' @param X Data input 
#' @param num_cut Number of cuts for RJ algorithm (suggested sqrt(p))
#' @param seed Seed (defalt = 1)
#'
#' @return Returns RJ algorithm result
#' @export
#'
#' @examples
#' X = generateSimulationData()
#' X = X$X
#' clust = RJclust(X, 100)
RJclust = function( X, num_cut = NULL, seed = 1 )
{
  # check that data is a matrix
  if ( !is.matrix( X ) )
  {
    stop( "Data must be in a matrix form" )
  }

  # check that num_cut is not too big
  if ( num_cut >= nrow( X ) )
  {
    stop( "num_cut must be < n" )
  }
  
  if ( num_cut >= nrow( X ) / 4 )
  {
    warning( "RJclust will preform beter with a num_cut that divides the data into larger chunks" )
  }
  
  if ( is.null( num_cut ) & nrow( X ) > 5000 )
  {
    warning( "RJclust will preform better with the scaled version, try passing in a num_cut value")
  }
  
  # if there is a num_cut indicated, run RJ_scale, otherwise run with no scale
  if ( !is.null( num_cut ) )
  {
    to_return = RJclust_backend( X, num_cut, seed ) 
  } else {
    to_return = RJclust_noscale( X, seed )
  }
  # run RJ algorithm
  return( to_return )
}


