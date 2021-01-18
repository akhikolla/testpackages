library(checkmate)
library(assertthat)
library(Formula)
library(abind)
library(rstan)

set_sampling_default  <- function(iter, warmup, chains, cores=1) {
    options(OncoBayes2.MC.iter=iter, OncoBayes2.MC.warmup=warmup, OncoBayes2.MC.chains=chains, mc.cores=cores)
}

very_fast_sampling <- function() {
    set_sampling_default(300, 150, 1, 1)
}

fake_sampling <- function() {
    set_sampling_default(4, 2, 1, 1)
}

default_sampling <- function() {
    set_sampling_default(NULL, NULL, NULL, NULL)
}

run_example <- function(example) {
    env <- new.env()
    suppressWarnings(example_model(example, env, silent=TRUE))
    invisible(env)
}

fake_sampling()

## set up slim sampling in case we are on CRAN
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    very_fast_sampling()
}

