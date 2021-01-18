
#' @title Add new constraint to constraint table
#' @description Add a new constraint to the mapping between a given macro dataset (class "macroACS")
#' and a matching micro dataset (class "micro_synthetic). May be called repeatedly to create a 
#' set of constraints.
#' @param attr_name The name of the attribute, or variable, that you wish to constrain.
#' @param attr_totals A named integer vector of counts per level of the new constraining attribute.
#' @param micro_data The micro dataset, of class \code{"micro_synthetic"}, for which you wish to
#' add a constraint.
#' @param constraint_list A \code{list} of prior constraints on the same dataset which you wish to
#' add to. Defaults to \code{NULL} (ie. the default is that this is the first constraint.)
#' @return A list of constraints.
#' @examples \dontrun{
#' ## assumes that you have a micro_synthetic dataset named test_micro and attribute counts
#' ## named a,e,g respectively 
#' c_list <- add_constraint(attr_name= "age", attr_totals= a, micro_data= test_micro)
#' c_list <- add_constraint(attr_name= "edu_attain", attr_totals= e, micro_data= test_micro,
#'                         constraint_list= c_list)
#' c_list <- add_constraint(attr_name= "gender", attr_totals= g, micro_data= test_micro,
#'                          constraint_list= c_list)
#' }
#' @export
add_constraint <- function(attr_name= "variable", attr_totals, micro_data,
                           constraint_list= NULL) {
  # 00. error checking
  #------------------------------------
  if (missing(attr_name) | missing(attr_totals)) 
    stop("attr_name and attr_totals must be supplied.")
  if (!is.micro_synthetic(micro_data)) stop("Input appropriate micro_data -- use class 'micro_synthetic'.")
  if (!is.character(attr_name)) stop("attr_name must be a string.")
  if (!exists(attr_name, as.environment(micro_data)))
    stop(paste("variable", attr_name, "is not in micro_data.", sep= " "))
  if (!is.numeric(attr_totals) | !all(attr_totals %% 1 == 0) | is.null(names(attr_totals))) 
    stop("attr_totals must a named numeric vector of integers.")
  if (!all(names(attr_totals) == levels(get(attr_name, as.environment(micro_data)))))
    stop("names of attr_totals must match levels and order of the associated variable in micro_data.")
  
  # 01. add constraint
  #------------------------------------
  if (is.null(constraint_list)) {
    constraint_list <- list(attr_totals)
    names(constraint_list)[1] <- attr_name
  } else {
    # check that attr_totals match prior constraints
    cur_cnt <- sum(constraint_list[[1]])
    new_cnt <- sum(attr_totals)
    if (cur_cnt != new_cnt) {
      stop ("current constraint total count does NOT match prior constraints.")
    } else {
      # add new constraint
      constraint_list[[length(constraint_list) + 1]] <- attr_totals
      names(constraint_list)[length(constraint_list)] <- attr_name
    }
  }
  # 02. return
  return(constraint_list)
}

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#' @title Calculate the total absolute error (TAE) between sample data and constraints.
#' @description Calculates the total absolute error (TAE) between sample micro data and constraining
#' totals from the matching macro data. Allows for updating of prior TAE instead of re-calculating
#' to improve speed in iterating. The updating feature is particularly helpful for optimizing
#' micro data fitting via simulated annealing (see \code{\link{optimize_microdata}}).
#' @param sample_data A \code{data.frame} with attributes matching \code{constraint_list}.
#' @param constraint_list A \code{list} of constraints. See \code{\link{add_constraint}}.
#' @param prior_sample_totals An optional \code{list} containing attribute counts of a prior sample 
#' corresponding to the constraint list. Defaults to \code{NULL}.
#' @param dropped_obs_totals An optional \code{list} containing attribute counts from the dropped 
#' observations in a prior sample. Defaults to \code{NULL}.
#' @param new_obs An optional \code{data.frame} containing new observations with attributes matching 
#' those in \code{sample_data}, \code{constraint_list}, and \code{prior_sample_totals}. Defaults 
#' to \code{NULL}.
#' @examples \dontrun{
#' ## assumes that you have a micro_synthetic dataset named test_micro and attribute count
#' ## named g respectively 
#' c_list <- add_constraint(attr_name= "gender", attr_totals= g, micro_data= test_micro,
#'             constraint_list= c_list)
#' calculate_TAE(test_micro, c_list)
#' }
#' @export
calculate_TAE <- function(sample_data, constraint_list, 
                          prior_sample_totals= NULL, dropped_obs_totals= NULL, new_obs= NULL) {
  ## 01. error checking
  #------------------------------------
  if (is.null(prior_sample_totals)) { # can leave sample data as promise (un-eval'd) if iterating
    if (!is.data.frame(sample_data)) stop("sample_data must be a data.frame")
    cnt <- sum(constraint_list[[1]])
    if (nrow(sample_data) != cnt) 
      stop("sample_data observation count does not match macro observation count.")  
  }
  
  if (!is.null(prior_sample_totals)) {
    if (is.null(dropped_obs_totals) | is.null(new_obs))
      stop("if prior_sample_totals is supplied, dropped_obs_totals and new_obs must also be supplied.")
    if (!is.list(prior_sample_totals) | !all(unlist(lapply(prior_sample_totals, is.numeric)))) 
      stop("prior_sample_totals must be a list of numerics.")
  }
  if (!is.null(dropped_obs_totals)) {
    if (is.null(prior_sample_totals) | is.null(new_obs))
      stop("if dropped_obs_totals is supplied, prior_sample_totals and new_obs must also be supplied.")
    if (!is.list(dropped_obs_totals) | !all(unlist(lapply(dropped_obs_totals, is.numeric)))) 
      stop("dropped_obs_totals must be a list of numerics.")
  }
  if (!is.null(new_obs)) {
    if (is.null(prior_sample_totals) | is.null(dropped_obs_totals))
      stop("if new_obs is supplied, prior_sample_totals and dropped_obs_totals must also be supplied.")
    if (!is.data.frame(new_obs)) 
      stop("new_obs must be a data.frame")
  } 
  
  ## 02. calculations
  #------------------------------------
  con_names <- names(constraint_list)
  if (is.null(prior_sample_totals)) { 
    # get sample totals, calculate attribute TAE. return total TAE and sample counts
    samp_totals <- sample_totals(attr_names= con_names, df= sample_data)
    tae <- tae_mapply(samp_totals, constraint_list)
    return(list(tae= sum(tae), samp_totals= samp_totals))
  } else {
    # update: attribute totals, attribute tae, total TAE. return
    new_obs_totals <- sample_totals(attr_names= con_names, df= new_obs)
    samp_totals <- mapply(
      function(prior_cnts, dropped_cnts, new_cnts) {
        return(prior_cnts - dropped_cnts + new_cnts)
      }, prior_cnts= prior_sample_totals, dropped_cnts= dropped_obs_totals,
       new_cnts= new_obs_totals)
    tae <- tae_mapply(samp_totals, constraint_list)
    return(list(tae= sum(tae), samp_totals= samp_totals))
  }
}

# Save typing duplication.
# helper function for calculating attribute counts. 
sample_totals <- function(attr_names, df) {
  lapply(attr_names, function(nm, df) {
    table(get(nm, as.environment(df)))
  }, df= df)
}

# Save typing duplication.
# helper function for calculating attribute TAE
# combines sample totals w/ constraints; then calcs TAE
tae_mapply <- function(samples, constraints) {
  mapply(function(samp_totals, constraint_totals) {
    d <- cbind(samp_totals, constraint_totals)
    return(sum(abs(d[,1] - d[,2])))
  }, samp_totals= samples, constraint_totals= constraints) 
}
 
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#' @title Optimize the selection of a micro data population.
#' @description Optimize the candidate micro dataset such that the lowest loss against the 
#' macro dataset constraints is obtained. Loss is defined here as total absolute error (TAE)
#' and constraints are defined by the \code{constraint_list}. Optimization is done by
#' simulated annealing--see details.
#' 
#' @section Details:
#' Spatial microsimulation involves the study of individual-level phenomena within a specified set of
#' geographies in which these individuals act. It involves the creation of synthetic data to model,
#' via simulation, these phenomena. As a first step to simulation, an appropriate micro-level 
#' (ie. individual) dataset must be generated. This function creates such appropriate micro-level
#' datasets given a set of candidate observations and macro-level constraints.
#' 
#' Optimization is done via simulated annealing, where we wish to minimize the total absolute error
#' (TAE) between the micro-data and the macro-constraints. The annealing procedure is controlled by 
#' the parameters \code{tolerance}, \code{resample_size}, \code{p_accept}, and 
#' \code{max_iter}. Specifically, \code{tolerance} indicates the maximum allowable TAE between the
#' output micro-data and the macro-constraints within a given \code{max_iter} allowable iterations 
#' to converge. \code{resample_size} and \code{p_accept} control movement about the candidate space. 
#' Specfically, \code{resample_size} controls the jump size between neighboring 
#' candidates and \code{p_accept} controls the hill-climbing rate for exiting local minima.
#' 
#' Please see the references for a more detailed discussion of the simulated annealing procedure.
#' 
#' @param micro_data A \code{data.frame} of micro data observations. 
#' @param prob_name It is assumed that observations are weighted and do not have an equal probability
#' of occurance. This string specifies the variable within \code{micro_data} that contains the probability
#' of selection.
#' @param constraint_list A \code{list} of constraining macro data attributes. See \code{\link{add_constraint}}
#' @param tolerance An integer giving the maximum acceptable loss (TAE), enabling early stopping.
#' Defaults to a misclassification rate of 1 individual per 1,000 per constraint. 
#' @param resample_size An integer controlling the rate of movement about the candidate space. 
#' Specifically, it specifies the number of observations to change between iterations. Defaults to 
#' \code{0.5\%} the number of observations.
#' @param p_accept The acceptance probability for the Metropolis acceptance criteria.
#' @param max_iter The maximum number of allowable iterations. Defaults to \code{10000L}
#' @param seed A seed for reproducibility. See \code{\link[base]{set.seed}}
#' @param verbose Logical. Do you wish to see verbose output? Defaults to \code{TRUE}
#'
#' @references Ingber, Lester. "Very fast simulated re-annealing." Mathematical and computer 
#' modelling 12.8 (1989): 967-973.
#' @references Metropolis, Nicholas, et al. "Equation of state calculations by fast computing 
#' machines." The journal of chemical physics 21.6 (1953): 1087-1092.
#' @references Szu, Harold, and Ralph Hartley. "Fast simulated annealing." Physics letters A 122.3 
#' (1987): 157-162.
#' @examples \dontrun{
#' ## assumes you have micro_synthetic object named test_micro and constraint_list named c_list
#' opt_data <- optimize_microdata(test_micro, "p", c_list, max_iter= 10, resample_size= 500, 
#'               p_accept= 0.01, verbose= FALSE)
#' }
#' @export
optimize_microdata <- function(micro_data, prob_name= "p", constraint_list, 
                               tolerance= round(sum(constraint_list[[1]]) / 2000 * length(constraint_list), 0),
                               resample_size= min(sum(constraint_list[[1]]), max(500, round(sum(constraint_list[[1]]) * .005, 0))), 
                               p_accept= 0.40, max_iter= 10000L, 
                               seed= sample.int(10000L, size=1, replace=FALSE),
                               verbose= TRUE) {
  ## 01. error checking
  #------------------------------------
  if (!is.micro_synthetic(micro_data)) stop("micro_data must be of class micro_synthetic.")
  if (!exists(prob_name, as.environment(micro_data))) {
    message("Probability vector not found in micro_data. \n Uniform probabilities used for sampling (ie- SRS).")
    micro_data[prob_name] <- 1/nrow(micro_data)
  }
  if (is.null(names(constraint_list)) | !is.list(constraint_list) | 
      !all(unlist(lapply(constraint_list, is.numeric))) |
      !all(sapply(names(constraint_list), function(n, df) {
        exists(n, as.environment(df))}, df= micro_data)))
    stop("constraint_list must be a named list of numeric vectors with corresponding attributes in micro_data.")
  if ((tolerance %% 1 != 0) | tolerance < 0) stop("tolerance must be specified as an integer >= 0.")
  if ((resample_size %% 1 != 0) | resample_size < 1) stop("resample_size must be specified as a positive integer.")
  if (!is.numeric(p_accept) | p_accept <= 0 | p_accept >= 1) stop("p_accept must be numeric in (0,1).")
  if ((max_iter %% 1 != 0) | max_iter < 1) stop("max_iter must be an integer.")

  ## create output structure for proposal and accepted TAE
  tae_path <- matrix(NA, nrow= max_iter, ncol= 2, dimnames= list(1:max_iter, c("proposal TAE", "current TAE")))
    
  ## 02. Take initial sample / iteration
  #------------------------------------
  mc <- match.call()
  set.seed(seed)
  sz <- sum(constraint_list[[1]])
  iter <- 1L
  
  cur_samp <- sample_micro(micro_data, sz, prob_name)
  tae_0 <- calculate_TAE(sample_data= cur_samp, constraint_list,
                         prior_sample_totals= NULL, dropped_obs_totals= NULL, new_obs= NULL)
  # add TAE to tracker
  tae_path[iter, ] <- rep(tae_0[[1]], 2)
  
  if (verbose) {
    cat("Iteration", iter, ": TAE =", sprintf("%.0f", tae_0[[1]]), "... \n")
  }
  
  # check if we got lucky
  if (tae_0$tae < tolerance) {
    tae_path <- tae_path[1:iter,]
    # got lucky, return:
    return(list(best_fit= cur_samp, tae= tae_0[[1]], call= mc, p_accept= p_accept, 
                iter= iter, max_iter= max_iter, tae_path= tae_path, seed= seed))
  } else {
    iter <- iter + 1L
  ## 03. Anneal to convergence
  #------------------------------------
    con_names <- names(constraint_list)
    resample_size <- min(resample_size, round(sz * .05,0)) # check for small geogs, never more than 5%
    
    # set cooling schedule; eg. T = {T_1, ..., T_{max_iter} }
    cool_rt <- p_accept * exp(-1/20 * seq(length.out= max_iter) / length(constraint_list))
    
    repeat {
      ##  (A) drop obs, grab new ones
      if (iter < 100) { # initially, take larger jumps
        drop_ind <- sample.int(nrow(cur_samp), size= sz * 0.1, replace=FALSE)
        new_obs  <- sample_micro(micro_data, sz * 0.1, prob_name)
      } else if (iter %% 500 == 0) { # every 500 iterations, take a big jump in the sample space
        drop_ind <- sample.int(nrow(cur_samp), size= round(sz * 0.2, 0), replace=FALSE)
        new_obs  <- sample_micro(micro_data, round(sz * 0.2, 0), prob_name)  
      } else {
        drop_ind <- sample.int(nrow(cur_samp), size= resample_size, replace=FALSE)
        new_obs  <- sample_micro(micro_data, resample_size, prob_name)  
      }
      
      ## (B) calculate new TAE
      drop_totals <- sample_totals(con_names, cur_samp[drop_ind,])
      tae_1 <- calculate_TAE(cur_samp, constraint_list,
                             prior_sample_totals= tae_0[[2]],
                             dropped_obs_totals= drop_totals,
                             new_obs= new_obs)
      
      if (verbose) {
        cat("Iteration", iter, ": Current TAE =", sprintf("%.0f", tae_0[[1]]), "\n Sample TAE =",
            sprintf("%.0f", tae_1[[1]]),"... \n")
      }
      
      if (tae_1[[1]] < tae_0[[1]]) {
        cur_samp <- rbind(cur_samp[-drop_ind, ], new_obs) # P(Accept | \delta E <0) == 1
        tae_0 <- tae_1
      } else { # P(Accept | \delta E > 0, T_k) \propto d_tae * U(0,1)
        tae_rel <- tae_1[[1]] / tae_0[[1]]
        if(stats::runif(1, 0, 1 * tae_rel) < cool_rt[iter]) {
          cur_samp <- rbind(cur_samp[-drop_ind, ], new_obs)
          tae_0 <- tae_1
        } # else -- stays the same
      }
      
      # add TAE to tracker
      tae_path[iter,] <- c(tae_1[[1]], tae_0[[1]]) 
      
      ## (C) check to exit
      if (tae_0[[1]] < tolerance | iter >= max_iter) {
        tae_path <- tae_path[1:iter,]
        return(list(best_fit= cur_samp, tae= tae_0[[1]], call= mc, p_accept= p_accept, 
                    iter= iter, max_iter= max_iter, tae_path= tae_path, seed= seed))
      } else { # (d) update for next iteration
        iter <- iter + 1
      }
    }
  }
}

# helper function to save typing
sample_micro <- function(df, size, prob_name) {
  if (!data.table::is.data.table(df)) {data.table::setDT(df)}
  data.table::data.table(df[sample.int(nrow(df), size= size, replace= TRUE, prob= df[[prob_name]]),
                -which(names(df) == prob_name), with= FALSE])
}