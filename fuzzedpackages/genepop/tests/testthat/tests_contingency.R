cat("tests_contingency.R:\n")

locinfile <- genepopExample('structest.txt')
outfile <- struc(locinfile) 
testthat::expect_equal(readLines(outfile)[34],
                       "P value=0.598782   S.E=0.006606 (243439 switches)")
