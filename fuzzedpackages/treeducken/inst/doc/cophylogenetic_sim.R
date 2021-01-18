## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------

library(treeducken)
set.seed(42)

exp_tips_host <- calculate_expected_leaves_sptree(1.0, 0.5, 2.0)
# on average we will get about 5 tips
# maybe we want something more like 8, so let's decrease the extinction rate
exp_tips_host <- calculate_expected_leaves_sptree(1.0, 0.3, 2.0)
# that looks about right, let's assume there are no host speciations that are
# not cospeciations
h_lambda <- 0.0
h_mu <- 0.3
c_lambda <- 1.0
# now we need to worry about the symbiont tree since only cospeciations are
# occurring
# we can use the same math as for the locus tree simulator
exp_tips_symb <- calculate_expected_leaves_locustree(t = 2.0,
                                                     dup_rate = 0.2,
                                                     loss_rate = 0.1, 8)
# a little less than 10 symbiont tips on average
s_lambda <- 0.2
s_mu <- 0.1
# let's assume that when symbionts speciate that they always inherit their
# ancestor's host repertoire
s_her <- 0.0

## -----------------------------------------------------------------------------
# simulate 10 paired trees with our set parameters.
host_symb_sets <- sim_cophylo_bdp(hbr = h_lambda,
                hdr = h_mu,
                sbr = s_lambda,
                cosp_rate = c_lambda,
                sdr = s_mu,
                host_exp_rate = s_her,
                time_to_sim = 2.0,
                numbsim = 10)
print(host_symb_sets, details = TRUE)

## -----------------------------------------------------------------------------
# # plot one or two of them
tree_set_of_interest <- host_symb_sets[[1]]
print(tree_set_of_interest)

## -----------------------------------------------------------------------------
plot(tree_set_of_interest, col = "red", lty = "dotted")

## -----------------------------------------------------------------------------
host_tree <- tree_set_of_interest$host_tree
symb_tree <- tree_set_of_interest$symb_tree
a <- tree_set_of_interest$association_mat
d <- parafit_stat(host_tr = host_tree, symb_tr = symb_tree, assoc_mat = a)
parafit_test(host_tr = host_tree,
             symb_tr = symb_tree,
             assoc_mat = a,
             D = d,
             reps = 99)

## -----------------------------------------------------------------------------
# of course that bit is maybe a little too verbose
# we may be interested in a lot of trees
# we can instead use cophylo_summary_stat
host_symb_summary_df <- cophy_summary_stat(host_symb_sets)
host_symb_summary_df

## -----------------------------------------------------------------------------
host_tree_locus_trees <- sim_locustree_bdp(host_tree,
                                           gbr = 0.4,
                                           gdr = 0.2,
                                           lgtr = 0.0,
                                           num_loci = 10)
host_tree_locus_trees

## -----------------------------------------------------------------------------
plot(host_tree_locus_trees)

## -----------------------------------------------------------------------------
symb_tree_locus_trees <- sim_locustree_bdp(symb_tree,
                                           gbr = 0.2,
                                           gdr = 0.1,
                                           lgtr = 0.1,
                                           num_loci = 10)
str(symb_tree_locus_trees[[1]])

## -----------------------------------------------------------------------------
host_locus_tree <- host_tree_locus_trees[[3]]
# randomly choose an effective popsize 
popsize <- 1e6
host_loci_gene_trees <- sim_multilocus_coal(host_locus_tree,
                                            effective_pop_size = popsize,
                                            num_reps = 100)

## -----------------------------------------------------------------------------
gopher_lice_map <- read.table(system.file("extdata",
                                          "gopher_lice_mapping.txt",
                                          package = "treeducken"),
                              stringsAsFactors = FALSE, header = TRUE)

gopher_lice_assoc_matrix <- convert_assoc_table_to_matrix(gopher_lice_map)
gopher_tree <- ape::read.nexus(system.file("extdata",
                                           "gophers_bd.tre",
                                           package = "treeducken"))
lice_tree <- ape::read.nexus(system.file("extdata",
                                         "lice_bd.tre",
                                         package = "treeducken"))
gopher_lice_cophylo <- convert_to_cophy(hostTree = gopher_tree,
                                          symbTree = lice_tree,
                                          assocMat = gopher_lice_assoc_matrix)
print(gopher_lice_cophylo)
cophy_summary_stat(gopher_lice_cophylo)

plot(gopher_lice_cophylo,
     fsize = 0.5,
     show_tip_label = FALSE,
     gap = 1,
     col = "purple",
     lty = "dashed")

