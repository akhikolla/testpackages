
#' @title Add new constraint to a set of geographies
#' @description Add a new constraint to the mapping between a a set of macro datasets and a matching 
#' set of micro dataset (supplied as class 'macro_micro'). May be called repeatedly to create a 
#' set of constraints across the sub-geographies.
#' @param attr_name The name of the attribute, or variable, that you wish to constrain.
#' @param attr_total_list A list of named integer vectors containing counts per level of the new 
#' constraining attribute for each geography.
#' @param macro_micro The geographical dataset of macro and micro data. Should be of class 
#' \code{"macro_micro"}.
#' @param constraint_list_list A \code{list} of lists containing prior constraints on the same dataset 
#' for which you wish to add to. Defaults to \code{NULL} (ie. the default is that this is the first 
#' constraint.)
#' @return A list of constraint lists.
#' @seealso \code{\link{add_constraint}}
#' @export
#' 
#' @examples \dontrun{
#' # assumes that micro_synthetic already exists in your environment
#' 
#' # 1. build constraints for gender and age
#' g <- all_geog_constraint_gender(micro_synthetic, method= "macro.table")
#' 
#' a <- all_geog_constraint_age(micro_synthetic, method= "macro.table")
#' 
#' # 2. bind constraints to geographies and macro-data
#' cll <- all_geogs_add_constraint(attr_name= "age", attr_total_list= a, 
#'           macro_micro= micro_synthetic)
#' cll <- all_geogs_add_constraint(attr_name= "gender", attr_total_list= g, 
#'           macro_micro= micro_synthetic, constraint_list_list= cll)
#' 
#' }
all_geogs_add_constraint <- function(attr_name= "variable", attr_total_list, macro_micro,
                                     constraint_list_list= NULL) {
  # 00. error checking
  #------------------------------------
  if (missing(attr_name)) 
    stop("attr_name must be supplied.")
  if (!is.list(macro_micro) | !is.synthACS(macro_micro))
    stop("macro_micro must be supplied as a class 'synthACS' object.")
  if (!is.character(attr_name)) stop("attr_name must be a string.")
  if (!all(unlist(lapply(macro_micro, function(l, nm) {exists(nm, as.environment(l[[2]]))}, nm= attr_name))))
    stop(paste("variable", attr_name, "is not contained in all elements of macro_micro.", sep= " "))
  if (is.list(attr_total_list) &
      (!all(unlist(lapply(attr_total_list, is.numeric))) | 
      !all(unlist(lapply(attr_total_list, function(l) all(l %% 1 == 0)))) |
      any(unlist(lapply(attr_total_list, function(l) is.null(names(l)))))))
    stop("attr_totals_list must a list of named numeric integer vectors.")
  
  if ( # similar check to \code{\link{add_constraint}} but across lists
    !all(mapply(function(mm, attr_name, attr_totals) {
    all(names(attr_totals) == levels(get(attr_name, as.environment(mm[[2]]))))
  }, mm= macro_micro, attr_name= rep(attr_name, length(attr_total_list)), 
     attr_totals= attr_total_list, SIMPLIFY= TRUE))
     ) stop("names of attr_totals must match levels and order of the associated variable in micro_data.")
  
  # 01. add constraint -- wrapper
  #------------------------------------
  micro_datas <- lapply(macro_micro, "[[", 2)
  
  if (!is.null(constraint_list_list)) {
    constraint_wrap <- mapply(add_constraint, 
                              attr_name= rep(attr_name, length(micro_datas)),
                              attr_totals= attr_total_list,
                              micro_data= micro_datas, 
                              constraint_list= constraint_list_list, SIMPLIFY= FALSE,
                              USE.NAMES= FALSE)  
  } else {
    constraint_wrap <- mapply(add_constraint, 
                              attr_name= rep(attr_name, length(micro_datas)),
                              attr_totals= attr_total_list,
                              micro_data= micro_datas, 
                              SIMPLIFY= FALSE, USE.NAMES= FALSE)
  }
  
  # 02. return
  return(constraint_wrap)
}


#' @title Optimize the selection of a micro data population for a set of geographies.
#' @description Optimize the candidate micro datasets such that the lowest loss against the 
#' macro dataset constraints are obtained. Loss is defined here as total absolute error (TAE)
#' and constraints are defined by the \code{constraint_list_list}. Optimization is done by
#' simulated annealing and geographies are run in parallel.
#' 
#' @param macro_micro The geographical dataset of macro and micro data. Should be of class 
#' \code{"macro_micro"}.
#' @param prob_name It is assumed that observations are weighted and do not have an equal probability
#' of occurance. This string specifies the variable within each dataset that contains the probability
#' of selection.
#' @param constraint_list_list A list of constraint lists. See \code{\link{add_constraint}}, 
#' \code{\link{all_geogs_add_constraint}}
#' @param p_accept The acceptance probability for the Metropolis acceptance criteria.
#' @param max_iter The maximum number of allowable iterations. Defaults to \code{10000L}
#' @param seed A seed for reproducibility. See \code{\link[base]{set.seed}}
#' @param leave_cores An \code{integer} for the number of cores you wish to leave open for other
#' processing.
#' @param verbose Logical. Do you wish to see verbose output? Defaults to \code{TRUE}
#' @seealso \code{\link{optimize_microdata}}
#' @export
#' 
#' @examples \dontrun{
#'  # assumes that micro_synthetic and cll already exist in your environment
#'  # see: examples for derive_synth_datasets() and all_geogs_add_constraint()
#'  optimized_la <- all_geog_optimize_microdata(micro_synthetic, prob_name= "p", 
#'      constraint_list_list= cll, p_accept= 0.01, max_iter= 1000L)
#' }
all_geog_optimize_microdata <- function(macro_micro, prob_name= "p", constraint_list_list, 
                                        p_accept= 0.40, max_iter= 10000L,
                                        seed= sample.int(10000L, size=1, replace=FALSE),
                                        leave_cores= 1L,
                                        verbose= TRUE) {
  
  # 01. error checking
  #------------------------------------
  if (!is.list(macro_micro) | !is.synthACS(macro_micro))
    stop("macro_micro must be supplied as a class 'synthACS' object.")
  if (!is.numeric(p_accept) | p_accept <= 0 | p_accept >= 1) stop("p_accept must be numeric in (0,1).")
  if ((max_iter %% 1 != 0) | max_iter < 1) stop("max_iter must be an integer.")
  if ((leave_cores %% 1 != 0) | leave_cores < 0) stop("leave_cores must be a non-negative integer.")
  
  # 02. wrap optimize micro in parallel
  #------------------------------------
  mc <- match.call()
  micro_datas <- lapply(macro_micro, "[[", 2)
  
  if (verbose) message("Beginning parallel optimization...")
  
  nnodes <- min(parallel::detectCores() - leave_cores, length(micro_datas))
  if (.Platform$OS.type != "unix") {cl <- parallel::makeCluster(nnodes, type= "PSOCK")}
  else {cl <- parallel::makeCluster(nnodes, type= "FORK")}
  
  parallel::clusterExport(cl, "data.table", envir= as.environment("package:data.table"))
  
  geography_anneal <- parallel::clusterMap(cl, RECYCLE= TRUE, SIMPLIFY= FALSE, .scheduling= "dynamic",
                                 fun= optimize_microdata, 
                                 micro_data= micro_datas, prob_name= prob_name,
                                 constraint_list= constraint_list_list,
                                 p_accept= p_accept, max_iter= max_iter,
                                 seed= seed, verbose= FALSE)
  
  parallel::stopCluster(cl)
  if (verbose) message("... Optimization complete")
  
  # 03. return
  #------------------------------------
  best_fits <- lapply(geography_anneal, function(l) return(l[["best_fit"]]))
  taes <- lapply(geography_anneal, function(l) return(l[["tae"]]))
  iters <- lapply(geography_anneal, function(l) return(l[["iter"]]))
  tae_paths <- lapply(geography_anneal, function(l) return(l[["tae_path"]]))
  
  smsm <- list(best_fit= best_fits, tae= taes, call= mc, p_accept= p_accept, 
               iter= iters, max_iter= max_iter, tae_paths= tae_paths, seed= seed, 
               D= length(constraint_list_list[[1]]))
  # add class
  class(smsm) <- "smsm_set"
  return(smsm)
}







#' @title Add a new attribute to a set (ie list) of synthetic_micro datasets
#' @description Add a new attribute to a set (ie list) of synthetic_micro datasets using conditional 
#' relationships between the new attribute and existing attributes (eg. wage rate conditioned on age 
#' and education level). The same attribute is added to *each* synthetic_micro dataset, where each
#' dataset is supplied a distinct relationship for attribute creation.
#' @param df_list A \code{list} of R objects each of class "synthetic_micro". 
#' @param prob_name A string specifying the column name of each \code{data.frame} in \code{df_list} 
#' containing the probabilities for each synthetic observation.
#' @param attr_name A string specifying the desired name of the new attribute to be added to the data.
#' @param conditional_vars An character vector specifying the existing variables, if any, on which 
#' the new attribute (variable) is to be conditioned on for each dataset. Variables must be specified 
#' in order. Defaults to \code{NULL} ie- an unconditional new attribute.
#' @param st_list A \code{list} of equal length to \code{df_list}. Each element of \code{st_list} is 
#' a \code{data.frame} symbol table with N + 2 columns. The last two columns must be:
#' 1. A vector containing the new attribute counts or percentages; 2. is a vector of the new attribute 
#' levels. The first N columns must match the conditioning scheme imposed by the variables in 
#' \code{conditional_vars}. See \code{\link{synthetic_new_attribute}} and examples.
#' @param leave_cores An \code{integer} for the number of cores you wish to leave open for other
#' processing.
#' @return A list of new synthetic_micro datasets each with class "synthetic_micro".
#' @seealso \code{\link{synthetic_new_attribute}}
#' @examples \dontrun{
#'  set.seed(567L)
#'  df <- data.frame(gender= factor(sample(c("male", "female"), size= 100, replace= TRUE)),
#'                  age= factor(sample(1:5, size= 100, replace= TRUE)),
#'                  pov= factor(sample(c("lt_pov", "gt_eq_pov"),
#'                                     size= 100, replace= TRUE, prob= c(.15,.85))),
#'                  p= runif(100))
#' df$p <- df$p / sum(df$p)
#' class(df) <- c("data.frame", "micro_synthetic")
#' 
#' # and example test elements
#' cond_v <- c("gender", "pov")
#' levels <- c("employed", "unemp", "not_in_LF")
#' sym_tbl <- data.frame(gender= rep(rep(c("male", "female"), each= 3), 2),
#'                       pov= rep(c("lt_pov", "gt_eq_pov"), each= 6),
#'                       cnts= c(52, 8, 268, 72, 12, 228, 1338, 93, 297, 921, 105, 554),
#'                       lvls= rep(levels, 4))
#' 
#' 
#' 
#' df_list <- replicate(10, df, simplify= FALSE)
#' st_list <- replicate(10, sym_tbl, simplify= FALSE)
#' 
#' # run
#' library(parallel)
#' syn <- all_geog_synthetic_new_attribute(df_list, prob_name= "p", attr_name= "variable",
#'                                         conditional_vars= cond_v,st_list= st_list)
#' }
#' @export
all_geog_synthetic_new_attribute <- function(df_list, prob_name= "p",
                                             attr_name= "variable",
                                             conditional_vars= NULL,
                                             st_list= NULL,
                                             leave_cores= 1L) {
 
  # 01. error checking
  #------------------------------------
  if (!is.null(st_list)) {
    if (length(df_list) != length(st_list))
      stop("when conditioning, st_list and df_list must have equal lengths")
  }
  if ((leave_cores %% 1 != 0) | leave_cores < 0) stop("leave_cores must be a non-negative integer.")
 
  # 02. wrap synthetic_new_attribute in parallel
  #------------------------------------
  len <- length(df_list)
  
  if (!is.null(conditional_vars)) { # need to replicate() these for clusterMap 
    conditional_vars <- replicate(len, conditional_vars, simplify= FALSE)
  }
  
  nnodes <- min(parallel::detectCores() - leave_cores, len)
  if (.Platform$OS.type != "unix") {cl <- parallel::makeCluster(nnodes, type= "PSOCK")}
  else {cl <- parallel::makeCluster(nnodes, type= "FORK")}
  
  # to allow simplified testing (using non synthACS class objects)
  if (is.synthACS(df_list)) {df_list2 <- lapply(df_list, function(l) return(l[[2]]))} 
  else {df_list2 <- df_list}
  
  parallel::clusterExport(cl, "setnames", envir= as.environment("package:data.table")) 
  
  synthetic_data <- parallel::clusterMap(cl, RECYCLE= TRUE, SIMPLIFY= FALSE, .scheduling= "dynamic",
                                fun= synthetic_new_attribute,
                                df= df_list2,
                                prob_name= prob_name, attr_name= attr_name,
                                conditional_vars= conditional_vars,
                                sym_tbl= st_list)
  
  parallel::stopCluster(cl)
  # 03. return
  #------------------------------------
  df_list <- mapply(function(x,y) {return(list(macro_constraints= x, synthetic_micro= y))},
                    x= lapply(df_list, function(l) return(l[[1]])),
                    y= synthetic_data, 
                    SIMPLIFY = FALSE)
  
  class(df_list) <- c("synthACS", "list")
  return(df_list) 
}