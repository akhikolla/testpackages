context("Test 3: fastcmh on sample data")

require(fastcmh)



test_that("fastcmh find significant intervals", {
        
#        folder <- "../../inst/extdata/"
        folder <- file.path(system.file(package="fastcmh"), "extdata")
        folder <- fixSeparator(folder)

        mylist <- runfastcmh(folder)
        sig <- mylist$sig
        soln <- data.frame(start=c(100), end=c(103), pvalue=c(4.83e-13))


        expect_equal(sig$start, soln$start)
        expect_equal(sig$end, soln$end)
        expect_equal(sig$pvalue, soln$pvalue)

#        cat("\n")
#        print(filtered)
        
        })
