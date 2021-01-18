
context("sim_ril_pedigree")

test_that("result of sim_ril_ped passes check_pedigree", {

    # parents must have length (2, 4, 8)
    expect_error(sim_ril_pedigree(parents=1:6))
    expect_error(sim_ril_pedigree(parents=NULL))
    expect_error(sim_ril_pedigree(parents=1))

    # results should pass check_pedigree
    for(i in c(2,4,8)) {
        expect_true(check_pedigree(sim_ril_pedigree(parents=1:i, selfing=FALSE)))
        expect_true(check_pedigree(sim_ril_pedigree(parents=1:i, selfing=TRUE),
                                   ignore_sex=TRUE))
    }

})
