cat("tests_linkage.R:\n")

locinfile <- genepopExample('sample.txt')
outfile <- test_LD(locinfile,"sample.txt.DIS")
testthat::expect_equal(readLines(outfile)[68],
             "last pop        ADH thre ADH-5      0.334362     0.000890   166751")
