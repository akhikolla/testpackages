#' Calculate summary statistics for gene trees
#'
#' @description  Calculates summary statistics including Colless' statistic, gamma statistic of the locus tree input as an index as part of a list, gamma statistic of gene tree, sackin statistic, cherry statistic, and time to most recent common ancestor
#'
#' @param container_tree_gene_tree_obj Locus tree object obtain from `sim_locustree_genetree_mlc`
#' @param container_tree_indx Index of locus tree object of interest
#'
#' @return Dataframe with summary statistics for each gene tree
#' @examples
#' # first simulate a species tree
#' mu <- 0.5
#' lambda <- 1.0
#' nt <- 6
#' tr <- sim_sptree_bdp(sbr = lambda, sdr = mu, numbsim = 1, n_tips = nt)
#' # for a locus tree with 100 genes sampled per locus tree
#' gentrees <- sim_multispecies_coal(tr[[1]],
#'                                   ne = 10000,
#'                                   num_sampled_individuals = 1,
#'                                   num_genes = 100)
#'
#' gt_df <- genetree_summary_stat(gentrees, container_tree_indx = 1)
#' @export
genetree_summary_stat <- function(container_tree_gene_tree_obj, container_tree_indx){

    genetrees <- container_tree_gene_tree_obj[[container_tree_indx]]$gene.trees
    locus_tree <- container_tree_gene_tree_obj[[container_tree_indx]]$container.tree
    num_genetrees <- length(genetrees)
    colless <- vector(length = num_genetrees)
    gamma_locus <- rep(ape::gammaStat(locus_tree), times = length(genetrees))
    gamma <- vector(length = num_genetrees)
    sackin <- vector(length = num_genetrees)
    tmrca <- vector(length = num_genetrees)
    cherries <- vector(length = num_genetrees)
    for (i in 1:num_genetrees) {
        colless[i] <- apTreeshape::colless(apTreeshape::as.treeshape.phylo(genetrees[[i]]))
        gamma[i] <- ape::gammaStat(genetrees[[i]])
        cherries[i] <- treeducken::count_cherries(genetrees[[i]])
        sackin[i] <- apTreeshape::sackin(apTreeshape::as.treeshape.phylo(genetrees[[i]]))
        tmrca[i] <- max(ape::node.depth.edgelength(genetrees[[i]]))
    }
    data.frame(colless, sackin, tmrca, gamma_locus, gamma, cherries)
}
#' Calculate cherry statistic for gene-trees
#' @author Emmanuel Paradis
#' @description Calculate cherry statistic according to the definition given in  McKenzie and Steel 2000 (see below for reference)
#' @param tree an object of class "phylo"
#' @return The value fo cherries on a tree
#' @details This calculates the value for the cherry test statistic on a rooted tree. Note that this does not perform the actual
#' hypothesis test against Yule or uniform tree models.
#'
#' @examples
#' # first simulate a species tree
#' mu <- 0.5
#' lambda <- 1.0
#' nt <- 6
#' tr <- sim_sptree_bdp(sbr = lambda, sdr = mu, numbsim = 1, n_tips = nt)
#' treeducken::count_cherries(tr[[1]])
#' # to do the hypothesis test you can use the ape version of this function
#' ape::cherry(tr[[1]])
#' @references
#' McKenzie, A. and Steel, M. (2000) Distributions of cherries for two models of trees. Mathematical Biosciences, 164, 81–92.
# this is copied from ape::cherry directly, I have copied it here because
# I wanted a return value of the statistic rather than the test itself
# copyright Emmanuel Paradis
# used under GNU public license
#' @export
count_cherries <- function(tree){
    n <- length(tree$tip.label)
    sum(tabulate(tree$edge[, 1][tree$edge[, 2] <= n]) == 2)
}