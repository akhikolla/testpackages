#' Simulates multi-locus coalescent on a given locus tree
#'
#' @description separates a locus tree into loci broken up by
#' duplications  and simulates the coalescent on each loci.
#'
#' @param locus_tree a locus tree from `sim_locustree_bdp` of class `phy`
#' @param effective_pop_size the effective population size
#' @param generation_time unit time per generation (default 1 year per generation)
#' @param mutation_rate number of mutations per unit time
#' @param num_reps number of coalescent simulations per locus
#' @return A list of list of gene trees of length `num_reps` simulated along each locus.
#' The first member of the list is the parent tree, all others are child trees
#'
#' @details
#' This simulation follows the algorithm given in Rasmussen and Kellis 2012.
#' The locus tree is scaled into coalescent units prior to being used.
#' The generation_time parameter default assumes 1 generation
#' per year if the units of the tree are in millions of years.
#' The mutation_rate parameter is by default set to 1e-6 mutations per year 
#' (this is totally arbitrary).
#' Also note that the return type is a list of many trees so for
#' sufficiently complicated locus trees with `num_reps` set to a larger
#' value may slow things considerably so use with caution. 
#'
#' @examples
#' # first simulate a species tree
#' mu <- 0.5
#' lambda <- 1.0
#' nt <- 6
#' tr <- sim_sptree_bdp(sbr = lambda, sdr = mu, numbsim = 1, n_tips = nt)
#' # for a locus tree with 100 genes sampled per locus tree
#' gene_br <- 0.1
#' gene_dr <- 0.02
#' transfer_rate <- 0.2
#' locus_tree <- sim_locustree_bdp(species_tree = tr[[1]],
#'                   gbr = gene_br,
#'                   gdr = gene_dr,
#'                   lgtr = transfer_rate,
#'                   num_loci = 1)
#' effect_popsize <- 1e6
#' gene_tree_obj <- sim_multilocus_coal(locus_tree[[1]],
#'                                       effect_popsize,
#'                                       num_reps = 20)
sim_multilocus_coal <- function(locus_tree,
                                effective_pop_size,
                                generation_time = 1,
                                mutation_rate = 1e-6,
                                num_reps) {
    if(effective_pop_size <= 0) {
        stop("'effective_pop_size' must be a strictly positive number")
    }
    if(generation_time <= 0) {
        stop("'generation_time' must be a strictly positive number")
    }
    if(num_reps < 1) {
        stop("'effective_pop_size' must be  at leat 1")
    }
    if(class(locus_tree) != "phylo") {
        stop("'locus_tree' must be an object of class 'phylo")
    }
    # first rescale the locus tree to be not in time but in coalescent units
    # tree is in time units so to get to generations we want

    ipp <- 1

    # split tree into subtrees at each duplication
    locus_trees_by_dup <- seperate_into_loci(locus_tree)
    if(class(locus_trees_by_dup) == "phylo") {
        message("This is a locus tree with only one loci")
        big_tree <- treeducken::sim_multispecies_coal(locus_tree,
                                           ne = effective_pop_size,
                                           mutation_rate = mutation_rate,
                                           generation_time = generation_time,
                                           num_sampled_individuals = ipp,
                                           rescale = TRUE,
                                           num_genes = num_reps)
        gt_list <- list("parent_tree" = big_tree[[1]]$gene.trees, "child_trees" = NULL)
        class(gt_list) <- "mlc_genetrees"
        return(gt_list)
    }
    mlc_df <- list(length = length(locus_trees_by_dup))
    # loop through the locus trees simulate using MSC then prune the loci
    # that have that locus as a subtree to one tip (to simulate that
    # that coalescent
    # event has occurred and won't affect the sims)xz
    little_trees <- list(length = length(num_reps))
    for(i in seq_len(length(locus_trees_by_dup) - 1)) {
        # add artificial root for treeducken

        if(is.null(locus_trees_by_dup[[i]]$root.edge)) {
             locus_trees_by_dup[[i]]$root.edge <- 0.0
        }
        mlc_df[[i]]  <- treeducken::sim_multispecies_coal(
                                                locus_trees_by_dup[[i]],
                                                ne = effective_pop_size,
                                                mutation_rate = mutation_rate,
                                                generation_time = generation_time,
                                                num_sampled_individuals = ipp,
                                                rescale = TRUE,
                                                num_genes = num_reps)
        little_trees[[i]] <- mlc_df[[i]][[1]]$gene.trees
        # collapse those loci on all
        locus_trees_by_dup <- collapse_locus_subtree(
                                            locus_trees_by_dup,
                                                    locus_trees_by_dup)
    }
    locus_trees_by_dup[[length(locus_trees_by_dup)]]$root.edge <- locus_tree$root.edge
    big_tree <- treeducken::sim_multispecies_coal(
                                      locus_trees_by_dup[[length(locus_trees_by_dup)]],
                                      ne = effective_pop_size,
                                      mutation_rate = mutation_rate,
                                      generation_time = generation_time,
                                      num_sampled_individuals = ipp,
                                      num_genes = num_reps)
    gt_list <- list("parent_tree" = big_tree[[1]]$gene.trees, "child_trees" = little_trees)
    class(gt_list) <- "mlc_genetrees"
    gt_list
}
#' Separate a locus tree into loci
#'
#' @details This separates loci based on node labels "D[A-Z]". This is intended
#' to be used internally, but should work with other trees where duplications
#' are marked similarly.
#' 
#' @param locus_tree tree of type `phy`
#' @return list of subtrees (with `locus_tree at the end`)
#' @examples
#' # first simulate a species tree
#' mu <- 0.5
#' lambda <- 1.0
#' nt <- 6
#' tr <- sim_sptree_bdp(sbr = lambda, sdr = mu, numbsim = 1, n_tips = nt)
#' # for a locus tree with 100 genes sampled per locus tree
#' gene_br <- 0.1
#' gene_dr <- 0.02
#' transfer_rate <- 0.0
#' locus_tree <- sim_locustree_bdp(species_tree = tr[[1]],
#'                   gbr = gene_br,
#'                   gdr = gene_dr,
#'                   lgtr = transfer_rate,
#'                   num_loci = 1)
#' locus_tree_subtrees <- seperate_into_loci(locus_tree[[1]])
#' @export
seperate_into_loci <- function(locus_tree) {
    if(class(locus_tree) != "phylo") {
        stop("'locus_tree' must be an object of class 'phylo")
    }
    if(!(any(grep("D[A-Z]", locus_tree$node.label)))) {
        return(locus_tree)
    }

    # split tree into subtrees at each duplication
    loc_trees_subtree_all <- ape::subtrees(locus_tree)
    ul_trees <- unlist(loc_trees_subtree_all, recursive = FALSE)
    node_labels <- ul_trees[which(names(ul_trees) == "node.label")]
    indices_of_dup <- which(node_labels$node.label != "")
    rev(loc_trees_subtree_all[indices_of_dup])

}
#' Collapse a clade into a single tip
#'
#' @details Takes a clade as input and collapses that clade to one tip
#' in all trees in `list_of_subtrees`.
#' @param list_of_subtrees a list of type `multiPhylo`
#' @param locus_to_collapse a subtree found within a subst of `list_of_subtrees`
#' @return multiPhy (list of trees) with the subtree in question collapse
#' @examples
#' lambda <- 1.0
#' mu <- 0.2
#' nt <- 10
#' trees <- sim_sptree_bdp(sbr = lambda, sdr = mu, numbsim = 1, n_tips = nt)
#' subtrees_of_trees <- ape::subtrees(trees[[1]])
#' st_of_interest <- subtrees_of_trees[[1]]
#' collapse_st_of_interest <- collapse_locus_subtree(trees, st_of_interest)
#' @export
collapse_locus_subtree <- function(list_of_subtrees,
                                   locus_to_collapse) {
    tip_labels_subtrees <- get_tip_labels_tree_list(list_of_subtrees)
    tips_to_remove <- locus_to_collapse$tip.label[-1]
        # figure out which loci have those tips

    tree_indices_to_drop <- lapply(tip_labels_subtrees,
                                   function(x) x %in% tips_to_remove)


    tree_indices_to_keep <- which(
                                unlist(
                                    lapply(tree_indices_to_drop, all))
                                        , useNames = FALSE)
    names(tree_indices_to_keep) <- NULL
    tree_indices_to_drop <- which(
                                unlist(
                                    lapply(tree_indices_to_drop, any
                            , useNames = FALSE)))
    names(tree_indices_to_drop) <- NULL
    drop_from_drop <- intersect(tree_indices_to_drop, tree_indices_to_keep)
    if(length(drop_from_drop) != 0)
        tree_indices_to_drop <- tree_indices_to_drop[-drop_from_drop]
    trees_to_prune <- list_of_subtrees[tree_indices_to_drop]
    pruned_trees <- lapply(trees_to_prune,
                                function(x) ape::drop.tip(x,
                                                         tip = tips_to_remove,
                                                         rooted = TRUE))
    list_of_subtrees[tree_indices_to_drop] <- pruned_trees
    list_of_subtrees
}

#' Get all the tip labels of a `multiPhylo` object
#' @details Retrieves the member "tip.label" from each tree in multi_tree
#' @param multi_tree an object of class `multiPhylo`
#' @return a list of the same length as `multi_tree` with only the tip labels
#' @examples
#' mu <- 0.5
#' lambda <- 1.0
#' nt <- 6
#' tr <- sim_sptree_bdp(sbr = lambda, sdr = mu, numbsim = 5, n_tips = nt)
#' tips_of_tr <- get_tip_labels_tree_list(tr)
get_tip_labels_tree_list <- function(multi_tree) {

    ul_multi_tree <- unlist(multi_tree, recursive = FALSE)
    ul_multi_tree[which(names(ul_multi_tree) == "tip.label")]
}

#' Retrieve all gene trees of the parent tree from a list generated from sim_multilocus_coal
#' @param gene_tree_list A list of length 2: "parent_tree" and "child_trees" both of which are of class "multiPhylo"
#' @return A `multiPhylo` object of only the gene trees generated on the parent subtree
#' @examples
#' #' # first simulate a species tree
#' mu <- 0.5
#' lambda <- 1.0
#' nt <- 6
#' tr <- sim_sptree_bdp(sbr = lambda, sdr = mu, numbsim = 1, n_tips = nt)
#' # for a locus tree with 100 genes sampled per locus tree
#' gene_br <- 0.1
#' gene_dr <- 0.02
#' transfer_rate <- 0.2
#' locus_tree <- sim_locustree_bdp(species_tree = tr[[1]],
#'                   gbr = gene_br,
#'                   gdr = gene_dr,
#'                   lgtr = transfer_rate,
#'                   num_loci = 1)
#' effect_popsize <- 1e6
#' gene_tree_obj <- sim_multilocus_coal(locus_tree[[1]],
#'                                       effect_popsize,
#'                                       num_reps = 2)
#' parent_trees <- retrieve_parent_genetrees(gene_tree_obj)
retrieve_parent_genetrees <- function(gene_tree_list) {
    if(class(gene_tree_list) != "mlc_genetrees")
        stop("list is not in the correct format, did you not use `sim_multilocus_coal`?")
    return(gene_tree_list$parent_tree)
}
#' Retrieves the gene trees of the child subtrees 
#' @param gene_tree_list A list of length 2: "parent_tree" and "child_trees" both of which are of class "multiPhylo"
#' @describeIn retrieve_parent_genetrees Returns a list of objects of class `multiPhylo`
retrieve_child_genetrees <- function(gene_tree_list) {
    if(class(gene_tree_list) != "mlc_genetrees")
        stop("list is not in the correct format, did you not use `sim_multilocus_coal`?")
    return(gene_tree_list$child_trees)
}