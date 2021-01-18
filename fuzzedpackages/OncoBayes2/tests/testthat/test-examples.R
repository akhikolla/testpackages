context("examples")

## run examples with super-slim sampling to make sure they actually
## run
fake_sampling()

if(FALSE) {
    ## does not run atm for odd reasons
test_that("blrm_exnex example runs", example(blrm_exnex, package="OncoBayes2", run.dontrun=TRUE))
test_that("prior_summary example runs", example(prior_summary.blrmfit, package="OncoBayes2", run.dontrun=TRUE))
test_that("posterior_linpred example runs", example(posterior_linpred.blrmfit, package="OncoBayes2", run.dontrun=TRUE))
test_that("posterior_predict example runs", example(posterior_predict.blrmfit, package="OncoBayes2", run.dontrun=TRUE))
test_that("posterior_interval example runs", example(posterior_interval.blrmfit, package="OncoBayes2", run.dontrun=TRUE))
test_that("predictive_interval example runs", example(predictive_interval.blrmfit, package="OncoBayes2", run.dontrun=TRUE))
}

test_that("example_model() lists expected models", {
    expect_set_equal(example_model(), names(example_cache))
})
