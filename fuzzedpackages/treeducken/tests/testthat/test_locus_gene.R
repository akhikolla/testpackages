# test species tree output is a list of trees with correct length
sim_test_spt_loct <- function(gbr, gdr, lgtr, numLoci, species_tree_len = 5.0){
    tr <- sim_sptree_bdp_time(0.1, 0.05, 1, species_tree_len)
    sim_locustree_bdp(tr[[1]],
                      gbr = gbr,
                      gdr = gdr,
                      lgtr = lgtr,
                      num_loci = numLoci)

}

sim_test_spt_loct <- function(gbr, gdr, lgtr, numLoci, species_tree_len = 5.0){
    tr <- sim_sptree_bdp_time(0.1, 0.05, 1, species_tree_len)
    sim_locustree_bdp(tr[[1]],
                      gbr = gbr,
                      gdr = gdr,
                      lgtr = lgtr,
                      num_loci = numLoci)

}

sim_test_spt_loct_equality <- function(gene_birth = 0.0, gene_death = 0.0, transfers = 0.0, numLoci){
    tr <- sim_sptree_bdp_time(0.1, 0.05, 1, 5.0)
    loctr <- sim_locustree_bdp(tr[[1]],
                               gbr = gene_birth,
                               gdr = gene_death,
                               lgtr = transfers,
                               num_loci = numLoci)
    all(sapply(loctr, ape::all.equal.phylo, current = tr[[1]], use.tip.label = FALSE))
}

test_that("sim_locustree_bdp produces the right number of trees", {
    expect_equal(length(sim_test_spt_loct(0.1, 0.05, 0.05, 10)), 10)
    expect_equal(length(sim_test_spt_loct(0.1, 0.05, 0.05, 5)), 5)
    expect_equal(length(sim_test_spt_loct(0.1, 0.05, 0.05, 50)), 50)
})



# test that input species tree is the same as locus tree if every param is set
# 0.0
test_that("sim_locustree_bdp returns tree concordant with species tree when locus tree parameters are 0.0",{
    expect_true(sim_test_spt_loct_equality(numLoci = 20))
})


get_length_tree <- function(tr){
    max(ape::node.depth.edgelength(tr)) + tr$root.edge
}

get_all_tree_lengths <- function(multiTree){
    min(sapply(multiTree, get_length_tree))
}
