cat("test_Fst.R:\n")


locinfile <- genepopExample('sample.txt')
outfile <- Fst(locinfile,outputFile="sample.txt.FST")
testthat::expect_equal(readLines(outfile)[141],
             "           All:  0.0225      0.1921      0.2103")
outfile <- Fst(locinfile,pairs=TRUE,outputFile="sample.txt.ST2")
testthat::expect_equal(readLines(outfile)[68], "4      0.1514  0.0469  0.2688 ")
outfile <- Fst(locinfile,sizes=TRUE,outputFile="sample.txt.RHO")
testthat::expect_equal(readLines(outfile)[153],
             "           All:  0.2777      0.3095      0.5013")
outfile <- Fst(locinfile,sizes=TRUE,pairs=TRUE,outputFile="sample.txt.ST2")
testthat::expect_equal(readLines(outfile)[80], "4      0.6250  0.2408  0.2812 ")
