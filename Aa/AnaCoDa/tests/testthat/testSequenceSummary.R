library(testthat)
library(AnaCoDa)

context("Sequence Summary")

test_that("general sequence summary functions", {
  expect_equal(testSequenceSummary(), 0)
})

#Used to test AAToCodonRange, but that is now returned by
#reference and cannot deal with that in R. Since this will
#change with codonTable being completed, the tests have been
#ommitted for the time being.

#this function will be moved to CodonTable when it is finished.
test_that("AA To Codon", {
  expect_equal(AAToCodon("A", FALSE), c("GCA", "GCC", "GCG", "GCT"))
  expect_equal(AAToCodon("C", FALSE), c("TGC", "TGT"))
  expect_equal(AAToCodon("D", FALSE), c("GAC", "GAT"))
  expect_equal(AAToCodon("E", FALSE), c("GAA", "GAG"))
  expect_equal(AAToCodon("F", FALSE), c("TTC", "TTT"))
  expect_equal(AAToCodon("G", FALSE), c("GGA", "GGC", "GGG", "GGT"))
  expect_equal(AAToCodon("H", FALSE), c("CAC", "CAT"))
  expect_equal(AAToCodon("I", FALSE), c("ATA", "ATC", "ATT"))
  expect_equal(AAToCodon("K", FALSE), c("AAA", "AAG"))
  expect_equal(AAToCodon("L", FALSE), c("CTA", "CTC", "CTG", "CTT", "TTA", "TTG"))
  expect_equal(AAToCodon("M", FALSE), c("ATG"))
  expect_equal(AAToCodon("N", FALSE), c("AAC", "AAT"))
  expect_equal(AAToCodon("P", FALSE), c("CCA", "CCC", "CCG", "CCT"))
  expect_equal(AAToCodon("Q", FALSE), c("CAA", "CAG"))
  expect_equal(AAToCodon("R", FALSE), c("AGA", "AGG", "CGA", "CGC", "CGG", "CGT"))
  expect_equal(AAToCodon("S", FALSE), c("TCA", "TCC", "TCG", "TCT"))
  expect_equal(AAToCodon("T", FALSE), c("ACA", "ACC", "ACG", "ACT"))
  expect_equal(AAToCodon("V", FALSE), c("GTA", "GTC", "GTG", "GTT"))
  expect_equal(AAToCodon("W", FALSE), c("TGG"))
  expect_equal(AAToCodon("Y", FALSE), c("TAC", "TAT"))
  expect_equal(AAToCodon("Z", FALSE), c("AGC", "AGT"))
  expect_equal(AAToCodon("X", FALSE), c("TAA", "TAG", "TGA"))
  expect_equal(AAToCodon("A", TRUE), c("GCA", "GCC", "GCG"))
  expect_equal(AAToCodon("C", TRUE), c("TGC"))
  expect_equal(AAToCodon("D", TRUE), c("GAC"))
  expect_equal(AAToCodon("E", TRUE), c("GAA"))
  expect_equal(AAToCodon("F", TRUE), c("TTC"))
  expect_equal(AAToCodon("G", TRUE), c("GGA", "GGC", "GGG"))
  expect_equal(AAToCodon("H", TRUE), c("CAC"))
  expect_equal(AAToCodon("I", TRUE), c("ATA", "ATC"))
  expect_equal(AAToCodon("K", TRUE), c("AAA"))
  expect_equal(AAToCodon("L", TRUE), c("CTA", "CTC", "CTG", "CTT", "TTA"))
  expect_equal(AAToCodon("M", TRUE), character(0))
  expect_equal(AAToCodon("N", TRUE), c("AAC"))
  expect_equal(AAToCodon("P", TRUE), c("CCA", "CCC", "CCG"))
  expect_equal(AAToCodon("Q", TRUE), c("CAA"))
  expect_equal(AAToCodon("R", TRUE), c("AGA", "AGG", "CGA", "CGC", "CGG"))
  expect_equal(AAToCodon("S", TRUE), c("TCA", "TCC", "TCG"))
  expect_equal(AAToCodon("T", TRUE), c("ACA", "ACC", "ACG"))
  expect_equal(AAToCodon("V", TRUE), c("GTA", "GTC", "GTG"))
  expect_equal(AAToCodon("W", TRUE), character(0))
  expect_equal(AAToCodon("Y", TRUE), c("TAC"))
  expect_equal(AAToCodon("Z", TRUE), c("AGC"))
  expect_equal(AAToCodon("X", TRUE), character(0))
})

test_that("Amino Acid Vector", {
  expect_equal(aminoAcids(), c(
  "A",
  "C",
  "D",
  "E",
  "F",
  "G",
  "H",
  "I",
  "K",
  "L",
  "M",
  "N",
  "P",
  "Q",
  "R",
  "S",
  "T",
  "V",
  "W",
  "Y",
  "Z",
  "X"))
})

test_that("Codon Vector", {
  expect_equal(codons(), c(
    "GCA", 
    "GCC",
    "GCG",
    "GCT",
    "TGC",
    "TGT",
    "GAC",
    "GAT",
    "GAA",
    "GAG",
    "TTC",
    "TTT",
    "GGA",
    "GGC",
    "GGG",
    "GGT",
    "CAC",
    "CAT",
    "ATA",
    "ATC",
    "ATT",
    "AAA",
    "AAG",
    "CTA",
    "CTC",
    "CTG",
    "CTT",
    "TTA",
    "TTG",
    "ATG",
    "AAC",
    "AAT",
    "CCA",
    "CCC",
    "CCG",
    "CCT",
    "CAA",
    "CAG",
    "AGA",
    "AGG",
    "CGA",
    "CGC",
    "CGG",
    "CGT",
    "TCA",
    "TCC",
    "TCG",
    "TCT",
    "ACA",
    "ACC",
    "ACG",
    "ACT",
    "GTA",
    "GTC",
    "GTG",
    "GTT",
    "TGG",
    "TAC",
    "TAT",
    "AGC",
    "AGT",
    "TAA",
    "TAG",
    "TGA"
  ))
})
