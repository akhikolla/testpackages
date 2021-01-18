context("Test 1: testing filtering")

require(fastcmh)

# Example 1:
#
# raw:
# [100, 102], [101, 103], [98, 101]
#
# filtered: 
# [101, 103]
# (has smallest p-value)

test_that("only one interval after filtering", {
        tau <- c(100, 101, 98)
        l <- c(2, 3, 4)
        pvalue <- c(1e-5, 1e-6, 1e-2)
        filteredSoln <- data.frame(start=c(101), end=c(103), pvalue=c(1e-6) )

        df <- data.frame(tau=tau, l=l, pvalue=pvalue)
        filtered <- cpp_test_filtering(df)

        expect_equal(filtered$start[1], filteredSoln$start[1])
        expect_equal(filtered$end[1], filteredSoln$end[1])
        expect_equal(filtered$pvalue[1], filteredSoln$pvalue[1])

        
        })



# Example 2:
#
# raw:
# [10, 12], [21, 23]
#
# filtered: 
# [10, 12], [21, 23]
# (non-overlapping)

test_that("two non-overlapping intervals", {
        tau <- c(10, 21)
        l <- c(2, 3)
        pvalue <- c(1e-5, 1e-6)
        filteredSoln <- data.frame(start=c(10, 21), end=c(11, 23), pvalue=c(1e-5, 1e-6) )

        df <- data.frame(tau=tau, l=l, pvalue=pvalue)
        filtered <- cpp_test_filtering(df)

        expect_equal(filtered$start, filteredSoln$start)
        expect_equal(filtered$end, filteredSoln$end)
        expect_equal(filtered$pvalue, filteredSoln$pvalue)

        
        })



# Example 3:
# (result from dataset)
#
# raw: (all + 500,000)
# [780, 781], [801, 802], [794, 796], [800, 802], [799, 802]
#
# filtered: 
# [780, 781], [794, 796], [800, 802]
# (first two non-overlapping, third interval has the smallest p-value
#  out of the three overlapping)

test_that("some overlapping, others non-overlapping, multiple intervals after filtering", {
        tau <- c(505780, 505801, 505794, 505800, 505799)
        l <-   c(2, 2, 3, 3, 4)
        pvalue <- c(2.31e-10, 4.32e-9, 4.19e-10, 1.66e-10, 3.11e-10)

        filteredSoln <- data.frame(start=c(505780, 505794, 505800), 
                                   end=  c(505781, 505796, 505802), 
                                   pvalue=c(2.31e-10, 4.19e-10, 1.66e-10))

        df <- data.frame(tau=tau, l=l, pvalue=pvalue)
        filtered <- cpp_test_filtering(df)

        expect_equal(filtered$start, filteredSoln$start)
        expect_equal(filtered$end, filteredSoln$end)
        expect_equal(filtered$pvalue, filteredSoln$pvalue)

        })



# Example 4:
# (checking edge case 1)
#
# [10, 11], [13, 15], [9, 11]
#
# (first and third become [10, 11] - 12 is not in any interval
#
# filtered: 
# [10, 11], [13, 15]
# (third has largest p-value)

test_that("edge case: one empty space between overlapping intervals", {
        tau <- c(10, 13, 9)
        l <- c(2, 3, 3)
        pvalue <- c(1e-5, 1e-6, 1e-4)
        filteredSoln <- data.frame(start=c(10, 13), end=c(11, 15), pvalue=c(1e-5, 1e-6) )

        df <- data.frame(tau=tau, l=l, pvalue=pvalue)
        filtered <- cpp_test_filtering(df)

        expect_equal(filtered$start, filteredSoln$start)
        expect_equal(filtered$end, filteredSoln$end)
        expect_equal(filtered$pvalue, filteredSoln$pvalue)
        
        })



# Example 5:
# (checking edge case 2)
#
# [10, 11], [13, 15], [9, 12]
#
# (first and third become [10, 12]  cluster, now adjacent to [13, 15],
#  BUT they do not overlap, just adjacent, so considered separate clusters)
#
# filtered: 
# [10, 11], [13, 15]
# (third again has largest p-value)
test_that("edge case: no empty space between overlapping intervals, but adjacent", {
        tau <- c(10, 13, 9)
        l <- c(2, 3, 4)
        pvalue <- c(1e-5, 1e-6, 1e-4)
        filteredSoln <- data.frame(start=c(10, 13), end=c(11, 15), pvalue=c(1e-5, 1e-6) )

        df <- data.frame(tau=tau, l=l, pvalue=pvalue)
        filtered <- cpp_test_filtering(df)

        expect_equal(filtered$start, filteredSoln$start)
        expect_equal(filtered$end, filteredSoln$end)
        expect_equal(filtered$pvalue, filteredSoln$pvalue)

        })



# Example 6:
# (checking edge case 3)
#
# [10, 11], [13, 15], [9, 13]
#
# (first and third become [10, 13]  cluster, now OVERLAPS with [13, 15],
#  EVEN though [10, 11] and [9, 13] do not overlap, the interval [9, 13]
#  'bridges' them, and so now there is only one cluster, and result is [13, 15],
#  which has the smallest p-value)
#
# filtered: 
# [13, 15]
test_that("edge case: a not-very-significant interval causes two nearby intervals to overlap", {
        tau <- c(10, 13, 9)
        l <- c(2, 3, 5)
        pvalue <- c(1e-5, 1e-6, 1e-2)
        filteredSoln <- data.frame(start=c(13), end=c(15), pvalue=c(1e-6) )

        df <- data.frame(tau=tau, l=l, pvalue=pvalue)
        filtered <- cpp_test_filtering(df)

        expect_equal(filtered$start, filteredSoln$start)
        expect_equal(filtered$end, filteredSoln$end)
        expect_equal(filtered$pvalue, filteredSoln$pvalue)

#        cat("\n")
#        print(filtered)
        
        })
