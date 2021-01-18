context("Test 2: more testing filtering")

require(fastcmh)



test_that("only one interval after filtering", {
        
#        filename <- "../../inst/extdata/unfilteredtest2.csv"
        filename <- file.path(system.file(package="fastcmh"), "extdata", "unfilteredtest2.csv")
        filename <- fixSeparator(filename)
        df <- read.csv(filename)
        tau <- df$tau
        l <- df$l
        pvalue <- df$pvalue

        df <- data.frame(tau=tau, l=l, pvalue=pvalue)
        filtered <- cpp_test_filtering(df)

        expect_equal(doAnyOverlap(filtered), FALSE)

#        cat("\n")
#        print(filtered)
        
        })
