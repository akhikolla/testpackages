cat("test_Fis.R:\n")

locinfile <- genepopExample('sample.txt')
outfile <- genedivFis(locinfile,outputFile = "sample.txt.DIV")
testthat::expect_equal(readLines(outfile)[130], "                      0.6736     -0.0722")
outfile <- genedivFis(locinfile, sizes=TRUE
                      , outputFile = "sample.txt.MSD")
testthat::expect_equal(readLines(outfile)[141], "                      4.0556      0.4658")
