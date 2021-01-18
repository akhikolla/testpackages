# Simulate a mixed data set
# @param models a vector of models to simulate (note that proportion free or
# proportion equal models will be treated in the same mannear)
# @param params parameters to used, for each model, the list provides a vector
# with the number of variables and the proportion of missing values
# @param nbSample size of the sample to simulate (not used if z is given)
# @param nbCluster number of class to simulate (not used if prop is given)
# @param z vector giving the classes of each sample
# @param prop vector with the probabilities of each classes (not used if z is given)
#
# @return two lists
simulMixedData <- function(models, params, nbSample, nbCluster=2, z = NULL, prop = NULL)
{
  # check
  if (!is.vector(models) < 2) { stop("models has to be a vector")}
  if (!is.list(params) < 2) { stop("params has to be a list")}
  if (length(models) != length(params)) { stop("params and models must have the same size")}
  if (nbCluster < 2) { stop("The number of clusters must be greater or equal to 2")}
  
  # check if we have to compute z or it is already given
  if (is.null((z)))
  {
    if (nbSample < nbCluster) { stop("The number of sample must be greater to the number of class")}
    if(is.null(prop)) { prop=rep(1/nbCluster, length.out=nbCluster); }
    z <- sample.int(1:nbCluster, size=nbSample, replace = TRUE, prob= prop);
  }
  else
  { nbSample <- length(z);}
  
  res <- vector("list", length = length(models)+1);
  names(res) <- c("z", models);
  res$z <- z;
  for (i in 1:length(models))
  {
    # get current model and parameters
    if (is.list(models)) { model <- models[[i]];}
    else                 { model <- models[i];}
    param <- params[[i]];
    if (clusterValidCategoricalNames(model))
    {
      # simulate Categorical
      all = c( "categorical_pk_pjk", "categorical_pk_pk", "categorical_p_pjk", "categorical_p_pk")
      
    }
    else if (clusterValidDiagGaussianNames(models[i]))
    {
      # simulate Gaussian
      all = c( "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_pk_sk", "gaussian_pk_s"
          , "gaussian_p_sjk", "gaussian_p_sj", "gaussian_p_sk", "gaussian_p_s")
    }
    else if (clusterValidGammaNames(models[i]))
    {
      # simulate Gamma
      all = c( "gamma_p_ajk_bjk",  "gamma_p_ajk_bk",  "gamma_p_ajk_bj",  "gamma_p_ajk_b"
          , "gamma_p_ak_bjk",  "gamma_p_ak_bk",  "gamma_p_ak_bj",  "gamma_p_ak_b"
          , "gamma_p_aj_bjk", "gamma_p_aj_bk"
          , "gamma_p_a_bjk", "gamma_p_a_bk"
          , "gamma_pk_ajk_bjk", "gamma_pk_ajk_bk", "gamma_pk_ajk_bj", "gamma_pk_ajk_b"
          , "gamma_pk_ak_bjk", "gamma_pk_ak_bk", "gamma_pk_ak_bj", "gamma_pk_ak_b"
          , "gamma_pk_aj_bjk", "gamma_pk_aj_bk"
          , "gamma_pk_a_bjk", "gamma_pk_a_bk"
      )
    }
    else if (clusterValidPoissonNames(models[i]))
    {
      # simulate Poisson
    all = c( "poisson_pk_ljk", "poisson_pk_lk", "poisson_pk_ljlk", "poisson_p_ljk", "poisson_p_lk", "poisson_p_ljlk")
  }
    else
    {
      stop("Invalid model name")
    }
  }
  res
}
