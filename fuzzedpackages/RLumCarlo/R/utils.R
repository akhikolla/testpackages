## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Internal, non-exported package functions
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Combine Model Output Array
#'
#' @param ... Arrays to be combined
#'
#' @author Johannes Friedrich, University of Bayreuth (Germany)
#'
#' @noRd
comb_array <- function(...) abind::abind(..., along = 3)

#' Create RLumCarlo Model Output List
#'
#' @param signal [numeric] (**required**): signal vector
#'
#' @param time [numeric] (**required**): time vector
#'
#' @param model [character] (*with default*): the name of the model, the functions tries
#' to set this automatically.
#'
#' @author Sebastian Kreutzer, Geography & Earth Sciences, Aberystwyth University (United Kingdom)
#'
#' @md
#' @noRd
.return_ModelOutput <- function(signal, time, model = as.character(sys.call(which = -1))[1]){
  list <- list(signal = signal, time = time)
  class(list) <- "RLumCarlo_Model_Output"
  attr(list, "model") <- if(!is.null(model)) model else NA_character_

  ## return
  return(list)
}

#' @title Register Multi-Core back end (helper function)
#'
#' @param method [character] (*with default*): Sequential `'seq'` or parallel `'par'`processing. In
#' the parallel mode the function tries to run the simulation on multiple CPU cores (if available) with
#' a positive effect on the computation time.
#'
#' @param cores [integer] (*with default*): allows to specify the number of used cores
#'
#'@md
#'@noRd
.registerClusters <- function(method, cores = parallel::detectCores(), verbose = FALSE){
  ## check the method parameter
  if(!method %in% c("par", "seq"))
    stop(paste0("[",as.character(sys.call(which = -1))[1],"()] Allowed keywords for 'method' are either 'par' or 'seq'!"),
             call. = FALSE)

  ##get number of cores
  if(is.na(cores) || is.null(cores) || !is.numeric(cores) || cores == 1)
    method <- "seq"

  if(method != "par"){
    cl <- parallel::makeCluster(1)
    doParallel::registerDoParallel(cl)
    ##ensures that we do not have any particular problems
    foreach::registerDoSEQ()

  } else {
    ##we never use all cores, this is not nice
    cl <- parallel::makeCluster(cores-1)
    doParallel::registerDoParallel(cl)

  }

  ##provide a feedback
  if(verbose) print(cl)

  return(cl)
}


#'@title Distribute electrons over cluster (helper function)
#'
#'@description Once the cluster system is created, the number of
#'total electrons needs to be distributed over the cluster, according
#'the number of cluster groups. This function does the job for it.
#'
#'@param clusters [list] (**required**): output from [create_Clusters]
#'
#'@param N_system [numeric] (**required**): total number of electrons in the
#'system created by [create_Clusters]
#'
#'@author Sebastian Kreutzer, Geography & Earth Sciences, Aberystwyth University (United Kingdom)
#'
#'@md
#'@noRd
.distribute_electrons <- function(clusters, N_system){
  ## get number of elements per group using a lookup table
  cl_groups <- table(clusters$cl_groups)[clusters$cl_groups]

  ## get number of cluster groups
  n_group <- max(clusters$cl_groups)

  ## generate data.fame with electrons per cluster
  df <- data.frame(
    GROUP = clusters$cl_groups,
    N_TOTAL = N_system[1],
    e_in_GROUP = N_system[1] / n_group,
    cl_in_GROUP = as.numeric(cl_groups),
    e_in_cluster = round((N_system / n_group) / as.numeric(cl_groups), 0)
  )

  ## distribute missing or superfluous electrons
  ## background: the electron number can be only an integer value,
  ## this is quick an dirty and sometimes the number of electrons
  ## differ by one
  if((e_diff <- N_system[1] - sum(df$e_in_cluster)) != 0) {
    e_count <- e_diff / abs(e_diff)
    id <- sample(1:nrow(df), size = abs(e_diff), replace = FALSE)
    df[["e_in_cluster"]][id] <- df[["e_in_cluster"]][id] + e_count
    df[["e_in_cluster"]][df[["e_in_cluster"]][id] < 0] <- 0
  }

  return(df)
}

