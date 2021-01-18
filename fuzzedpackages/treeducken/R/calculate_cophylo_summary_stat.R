#' @describeIn cophy_summary_stat Calculates the summary statistics for one index of the list of cophylogenetic objects
#'
#' @return A vector consisting of (in order) cospeciations, host speciations, host extinctions, symbiont speciations, symbiont extinctions, parafit statistic, and parafit p-value
#' @export
cophy_summary_stat_by_indx <- function(cophy_obj, cophy_obj_indx){
  if(cophy_obj_indx < 1)
    stop("'cophy_obj_indx' must be greater than 0")
  if(!is.numeric(cophy_obj_indx))
    stop("'cophy_obj_indx' must be a number.")
  if(!is.null(cophy_obj[[cophy_obj_indx]]$event_history)){
    cat("Calculating for replicate ", cophy_obj_indx, "\n")
    events <- cophy_obj[[cophy_obj_indx]]$event_history
    event_types <- levels(as.factor(events$Event_Type))
    num_each_event <- tabulate(as.factor(events$Event_Type))
    names(num_each_event) <- event_types
    cospeciations <- num_each_event["C"]
    host_speciations <- num_each_event["HG"]
    host_extinctions <- num_each_event["HL"]
    symbiont_speciations <- num_each_event["SG"]
    symbiont_extinctions <- num_each_event["SL"]
  }
  else{
    cospeciations <- 0
    host_speciations <- 0
    host_extinctions <- 0
    symbiont_speciations <- 0
    symbiont_extinctions <- 0
  }
  host_expansions <- NULL # not sure why I don't have a specific tag for this...
  # if(length(cophy_obj[[cophy_obj_indx]]$host_tree$tip.label) < 3 ||
  #    length(cophy_obj[[cophy_obj_indx]]$symb_tree$tip.label) < 3)
  # {
  #   parafits <- NA
  #   parafit_test <- NA
  # }
  # else
  # {
    parafits <- treeducken::parafit_stat(cophy_obj[[cophy_obj_indx]]$host_tree,
                                         cophy_obj[[cophy_obj_indx]]$symb_tree,
                                         cophy_obj[[cophy_obj_indx]]$association_mat)
    parafit_test <- treeducken::parafit_test(cophy_obj[[cophy_obj_indx]]$host_tree,
                                             cophy_obj[[cophy_obj_indx]]$symb_tree,
                                             cophy_obj[[cophy_obj_indx]]$association_mat,
                                             parafits)
    # cophy_eigen <- treeducken::almost_parafit_stat(cophy_obj[[cophy_obj_indx]]$host_tree,
    #                                   cophy_obj[[cophy_obj_indx]]$symb_tree,
    #                                   cophy_obj[[cophy_obj_indx]]$association_mat)
 # }
  c(cospeciations,
    host_speciations,
    host_extinctions,
    symbiont_speciations,
    symbiont_extinctions,
    #host_expansions,
    parafits,
    parafit_test)
   # cophy_eigen)
}
#' Calculates summary statistics for cophylogenetic objects
#'
#' @description For cophylogenetic objects produced in treeducken via `sim_cophylo_bdp`, calculates the numbers of different events of interest. In addition, calculates and tests the ParaFit test.
#'
#' @param cophy_obj The cophylogenetic object produced via `sim_cophylo_bdp`
#' @param cophy_obj_indx The index with `cophy_obj` for `cophylo_summary_stat_by_indx`
#'
#' @return A dataframe containing statistics relevant to cophylogenetic analysis
#' @examples
#' host_mu <- 0.5 # death rate
#' host_lambda <- 2.0 # birth rate
#' numb_replicates <- 1
#' time <- 1.0
#' symb_mu <- 0.2
#' symb_lambda <- 0.4
#' host_shift_rate <- 0.0
#' cosp_rate <- 2.0
#'
#' cophy_pair <- sim_cophylo_bdp(hbr = host_lambda,
#'                            hdr = host_mu,
#'                            cosp_rate = cosp_rate,
#'                            host_exp_rate = host_shift_rate,
#'                            sdr = symb_mu,
#'                            sbr = symb_lambda,
#'                            numbsim = numb_replicates,
#'                            time_to_sim = time)
#' summary_stats <- cophy_summary_stat(cophy_pair)
#' @export
cophy_summary_stat <- function(cophy_obj) {
    if(class(cophy_obj) != "multiCophy") {
      if(class(cophy_obj) == "cophy") {
        mult_cophy_obj <- list(cophy_obj)
        class(mult_cophy_obj) <- "multiCophylo"
        stat_df <- data.frame(matrix(0, nrow = 1, ncol = 7))
        stat_df[1,] <- treeducken::cophy_summary_stat_by_indx(mult_cophy_obj, 1)
      }
      else
        stop("'cophy_obj' must be an object of class 'multiCophylo")
    }
    else{
      num_cophy_obj <- length(cophy_obj)
      stat_df <- data.frame(matrix(0, nrow = num_cophy_obj, ncol = 7))
      for(i in 1:num_cophy_obj){
        stat_df[i,] <- treeducken::cophy_summary_stat_by_indx(cophy_obj, i)
      }
    }
    colnames(stat_df) <- c("Cospeciations",
                           "Host_Speciations",
                           "Host_Extinctions",
                           "Symbiont_Speciations",
                           "Symbiont_Extinctions",
                        #   "Host Expansions",
                            "Parafit_Stat",
                            "Parafit_P-value")
 ##                           "Cophylo_Eigen")
    if(!is.null(cophy_obj$event_history)){
      stat_df[which(is.na(stat_df[, 1:5]), arr.ind = TRUE)] <- 0
    }
    stat_df
}
#' Calculate the ParafitGlobal statistic on 2 trees and their association matrix
#'
#' @description Calculate the ParafitGlobal Statistic to be used in the hypothesis test described in Legendre et al. (2002).
#' The null hypothesis of this test being that the evolution of the two trees together with their associations at the present
#' have been independent.
#' @param host_tr The host tree of class "phy"
#' @param symb_tr The symbiont tree of class "phy
#' @param assoc_mat Association matrix between the extant tips of `host_tr` and `symb_tr`
#' @details
#' `parafit_stat` drops any non-extant tips from the tree. Then the phylogenetic distance matrix is obtained for both host and symbiont tree.
#' Next the principal coordinates are found for the host and symbiont distance matrices before these PCoA vectors are used in the following
#' matrix multiplication following Legendre et al. (2002): D = H t(A) A. The trace is then found of this to get our ParaFitGlobal Statistic.
#'
#' The test function `parafit_test` performs a row-wise permutation of the association matrix as described in Legendre et al. 2002. This is
#' performed a number of times set by the user (default is 999) and a p-value is output.
#'
#' The value from this is input into the test function. Note that this gives only the raw statistic unlike `ape::parafit`. That is the
#' only reason it is implemented here in treeducken (similar to `treeducken::cherries`).
#' @examples
#' tr_pair <- sim_cophylo_bdp(hbr=0.1,
#'                           hdr=0.05,
#'                           sdr=0.1,
#'                           host_exp_rate=0.4,
#'                           sbr = 0.05,
#'                           cosp_rate = 1.0,
#'                           numbsim = 1,
#'                           time_to_sim = 1)
#' # maybe we are interested in only cophylogenetic object 1
#' ht <- tr_pair[[1]]$host_tree
#' st <- tr_pair[[1]]$symb_tree
#' A <- tr_pair[[1]]$association_mat
#' pfs <- parafit_stat(host_tr = ht, symb_tr = st, assoc_mat = A)
#'
#' parafit_test(ht, st, A, pfs, reps = 19)
#' @seealso parafit_test
#' @references
#' Legendre, P., Y. Desdevises and E. Bazin. 2002. A statistical test for host-parasite coevolution. Systematic Biology, 51(2), 217–234.
#' @export
parafit_stat <- function(host_tr, symb_tr, assoc_mat){
    if(class(host_tr) != "phylo"){
      stop("'host_tr' must be an object of class 'phylo'")
    }

    if(class(symb_tr) != "phylo"){
      stop("'symb_tr' must be an object of class 'phylo'")
    }
    if(!("matrix" %in% class(assoc_mat)))
      stop("'assoc_mat' must be an object of class 'matrix'")
    host_tree <- treeducken::drop_extinct(host_tr, tol= 0.001)
    symb_tree <- treeducken::drop_extinct(symb_tr, tol = 0.001)

    if(length(host_tree$tip.label) < 3)
    {
      warning("'host_tr' must be a tree with more than 2 extant tips to calculate the parafit stat returning NA.\n")
      return(NA)
    }
    if(length(symb_tree$tip.label) < 3)
    {
      warning("'symb_tr' must be a tree with more than 2 extant tips to calculate the parafit stat returning NA.")
      return(NA)
    }
    if(length(host_tree$tip.label) != ncol(assoc_mat))
      stop("'assoc_mat' must have the same number of columns as extant tips in 'host_tr'. It does not.")
    if(length(symb_tree$tip.label) != nrow(assoc_mat))
      stop("'assoc_mat' must have the same number of rows as extant tips in 'symb_tr'. It does not.")
    H <- ape::cophenetic.phylo(host_tree)
    S <- ape::cophenetic.phylo(symb_tree)
    H_pcoas <- ape::pcoa(H)
    S_pcoas <- ape::pcoa(S)
# remember the weirdness with pca$values corresponding to cospeciations
    D <- t(H_pcoas$vectors) %*% t(assoc_mat) %*% S_pcoas$vectors
    sum(diag(D)^2)
}
#' @describeIn parafit_stat Perform ParaFit Hypothesis Test
#' @param D the statistic calculated using `parafit_stat`
#' @param reps Number of permutations to perform on the association matrix for the hypothesis test
#' @return A p-value for the hypothesis test described above
#' @export
parafit_test <- function(host_tr, symb_tr, assoc_mat, D, reps = 99){
    if(!is.numeric(D))
    {
      if(is.na(D))
        return(NA)
      else
        stop("'D' must be a number. You seem to have input nan")
    }
    if(!is.numeric(reps))
      stop("'reps' must be a number. Whatever you put in is not a number.")
    if(class(host_tr) != "phylo"){
      stop("'host_tr' must be a binary phylogenetic tree")
    }
    if(class(symb_tr) != "phylo"){
      stop("'symb_tr' must be a binary phylogenetic tree")
    }
    if(!("matrix" %in% class(assoc_mat)))
      stop("'assoc_mat' must be a matrix")
    null_dist <- vector(length = reps)
    for(i in 1:reps){
        shuffled_A <- t(apply(assoc_mat, 1, sample))
        null_dist[i] <- parafit_stat(host_tr, symb_tr, shuffled_A)
    }
    null_dist <- append(null_dist, D)
    length(null_dist[null_dist >= D]) / (reps + 1)
}
