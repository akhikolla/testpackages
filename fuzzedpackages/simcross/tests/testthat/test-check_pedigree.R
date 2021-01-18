
context("check_pedigree")

test_that("simple tests of check_pedigree", {

    tab <- sim_ril_pedigree(7)
    expect_true(check_pedigree(tab))

    # not in correct order
    tab2 <- tab[c(9, 1:8, 10:nrow(tab)),]
    expect_error(check_pedigree(tab2))

    # male mom
    tab3 <- tab
    tab3[1,4] <- 1
    expect_error(check_pedigree(tab3))

    # female dad
    tab4 <- tab
    tab4[2,4] <- 0
    expect_error(check_pedigree(tab4))

})

