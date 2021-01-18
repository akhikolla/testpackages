cat("tests_private_all.R:\n")

locinfile <- genepopExample('sample.txt')
outfile <- Nm_private(locinfile,"sample.txt.PRI")
testthat::expect_equal(readLines(outfile)[15],
             "Number of migrants after correction for size= 0.255145")
outfile <- basic_info(locinfile,"sample.txt.INF")
testthat::expect_equal(readLines(outfile)[798], "last pop   0.750 0.250 4     ")
