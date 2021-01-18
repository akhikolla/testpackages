# library(testthat)
# library(AnaCoDa)
# 
# context("Gene")
# 
# test_that("general gene functions", {
#   expect_equal(testGene(), 0)
# })
# 
# g <- new(Gene, "ATGCTCATTCTCACTGCTGCCTCGTAG", "2", "New Test Gene")
# test_that("get AA Count", {
#   expect_equal(g$getAACount("M"), 1)
#   expect_equal(g$getAACount("L"), 2)
#   expect_equal(g$getAACount("I"), 1)
#   expect_equal(g$getAACount("T"), 1)
#   expect_equal(g$getAACount("A"), 2)
#   expect_equal(g$getAACount("S"), 1)
#   expect_equal(g$getAACount("X"), 1)
#   expect_equal(g$getAACount("G"), 0)
#   
#   #Checking invalid cases
#   expect_equal(g$getAACount("g"), 0)
#   expect_equal(g$getAACount("AA"), 0)
# })
# 
# test_that("get Codon Counts", {
#   expect_equal(g$getCodonCount("ATG"), 1)
#   expect_equal(g$getCodonCount("CTC"), 2)
#   expect_equal(g$getCodonCount("ATT"), 1)
#   expect_equal(g$getCodonCount("ACT"), 1)
#   expect_equal(g$getCodonCount("GCT"), 1)
#   expect_equal(g$getCodonCount("GCC"), 1)
#   expect_equal(g$getCodonCount("TCG"), 1)
#   expect_equal(g$getCodonCount("TAG"), 1)
#   expect_equal(g$getCodonCount("AAA"), 0)
#   
#   #Checking invalid cases
#   expect_equal(g$getCodonCount("atg"), 0)
#   expect_equal(g$getCodonCount("ATGG"), 0)
# })
# 
# test_that("get RFP Value", {
#   expect_equal(g$getSumRFPCountForCodon("TGC", 1), 0)
#   expect_equal(g$getSumRFPCountForCodon("CAC", 1), 0)
#   expect_equal(g$getSumRFPCountForCodon("GTG", 1), 0)
#   expect_equal(g$getSumRFPCountForCodon("TCC", 1), 0)
#   
#   #Checking invalid cases
#   expect_equal(g$getSumRFPCountForCodon("atg", 1), 0)
#   expect_equal(g$getSumRFPCountForCodon("ATGG", 1), 0)
# })
# 
# test_that("get Codon Positions", {
#   expect_equal(g$getCodonPositions("ATG"), c(0))
#   expect_equal(g$getCodonPositions("CTC"), c(1,3))
#   expect_equal(g$getCodonPositions("ATT"), c(2))
#   expect_equal(g$getCodonPositions("ACT"), c(4))
#   expect_equal(g$getCodonPositions("GCT"), c(5))
#   expect_equal(g$getCodonPositions("GCC"), c(6))
#   expect_equal(g$getCodonPositions("TCG"), c(7))
#   expect_equal(g$getCodonPositions("TAG"), c(8))
#   expect_equal(g$getCodonPositions("GTG"), numeric(0))
#   
#   #Checking invalid cases
#   expect_equal(g$getCodonPositions("atg"), numeric(0))
#   expect_equal(g$getCodonPositions("ATGG"), numeric(0))
# })
