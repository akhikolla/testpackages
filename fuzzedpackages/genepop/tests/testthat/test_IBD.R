cat("test_IBD.R:\n")

locinfile <- genepopExample('w2.txt')
outfile <- ibd(locinfile,"w2.txt.ISO", geographicScale = "Log", statistic="e")
if (sessionInfo()[["R.version"]][["arch"]]=="i386") { # i386 vs x64
  testthat::expect_equal(readLines(outfile)[229], "0.01514 [ -0.0143952 , 0.0374583 ]") 
} else testthat::expect_equal(readLines(outfile)[229], "0.01514 [ -0.0143947 , 0.0374584 ]")


locinfile <- genepopExample('PEL1600withCoord.txt')
outfile <- ibd(locinfile,"PEL1600withCoord.ISO", statistic = "SingleGeneDiv",
               geographicScale = "1D")
nums <- as.numeric(unlist(strsplit(readLines(outfile)[59], "[^0-9e.-]+")))
testthat::expect_equal(nums,c(2.92606e-06, 6.28251e-07, 6.25044e-06))

