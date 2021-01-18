#' Drops extinct tips from tree
#' @author LJ Harmon, and JW Brown
#' This is a direct port of the geiger function, I import it here for convenience.
#' This code is copied under  GPL 3 license. 
#' 
#' @param phy a 'phylo' class object
#' @param tol tolerance in decimal values for branch lengths
#' @return A 'phylo' class object with extinct tips removed
#' @references 
#' Pennell M, Eastman J, Slater G, Brown J, Uyeda J, Fitzjohn R, Alfaro M, Harmon L (2014). “geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees.” Bioinformatics, 30, 2216-2218
#' @examples
#' mu <- 0.5 # death rate
#' lambda <- 2.0 # birth rate
#' numb_replicates <- 10
#' numb_extant_tips <- 4
#' # simulate trees under the GSA so first simulates a tree with
#' # numb_extant_tips * 100 tips counting each time we have a tree with 10 tips
#' # then randomly picks one of those trees
#'
#' tree_list <- sim_sptree_bdp(sbr = lambda,
#'                 sdr = mu,
#'                 numbsim = numb_replicates,
#'                 n_tips = numb_extant_tips)
#' pruned <- drop_extinct(tree_list[[1]])
#' 
#' @export
drop_extinct <- function(phy, tol = NULL) {

    if (!"phylo" %in% class(phy)) {
        stop("'phy' is not of class 'phylo'.");
    }
    if (is.null(phy$edge.length)) {
        stop("'phy' does not have branch lengths.");
    }
    if (is.null(tol)) {
        tol <- min(phy$edge.length) / 1000;
    }
    aa <- treeducken::is_extinct(phy = phy, tol = tol);
    if (length(aa) > 0) {
        #cat("Dropping", length(aa), "taxa:", aa, "\n", sep=" ");
        phy <- ape::drop.tip(phy, aa);
    }
    return(phy);
}
#' Identify extinct tips from tree
#' 
#' 
#' This is a direct port of the geiger function, I import it here for convenience.
#' This code is copied under  GPL 3 license.
#' 
#' @param phy a 'phylo' class object
#' @param tol tolerance in decimal values for branch lengths
#' @return A list of the tips that are extinct
#' @references 
#' Pennell M, Eastman J, Slater G, Brown J, Uyeda J, Fitzjohn R, Alfaro M, Harmon L (2014). “geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees.” Bioinformatics, 30, 2216-2218
#' @examples
#' mu <- 0.5 # death rate
#' lambda <- 2.0 # birth rate
#' numb_replicates <- 10
#' numb_extant_tips <- 4
#' # simulate trees under the GSA so first simulates a tree with
#' # numb_extant_tips * 100 tips counting each time we have a tree with 10 tips
#' # then randomly picks one of those trees
#'
#' tree_list <- sim_sptree_bdp(sbr = lambda,
#'                 sdr = mu,
#'                 numbsim = numb_replicates,
#'                 n_tips = numb_extant_tips)
#' is_extinct(tree_list[[1]])
#' @export
is_extinct <- function (phy, tol = NULL) {
    if (!"phylo" %in% class(phy)) {
        stop("'phy' is not of class 'phylo'.");
    }
    if (is.null(phy$edge.length)) {
        stop("'phy' does not have branch lengths.");
    }
    if (is.null(tol)) {
        tol <- min(phy$edge.length) / 1000;
    }
   # phy <- reorder(phy);
    xx <- numeric(ape::Ntip(phy) + phy$Nnode);
    for (i in seq_len(length(phy$edge[,1]))) {
        xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i];
    }
    aa <- max(xx[1:ape::Ntip(phy)]) - xx[1:ape::Ntip(phy)] > tol;
    if (any(aa)) {
        return(phy$tip.label[which(aa)]);
    } else {
        return(NULL);
    }
}
