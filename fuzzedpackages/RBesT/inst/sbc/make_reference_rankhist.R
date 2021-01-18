library(pkgbuild)
## ensure that the current dev version of RBesT is loaded
pkgbuild::compile_dll("../..")

library(assertthat)
library(batchtools)
library(dplyr)
library(knitr)
source("sbc_tools.R")
set.seed(453453)

#' according to the docs this speeds up the reduce step
options(batchtools.progress = FALSE)

registry_tmp <- Sys.getenv("TMP_NFS", tempdir())

#' Evaluate dense and sparse data-scenario

#' Dense: 10 trials with 40 entries each
dense_data  <- list(group=rep(1:10, each=40))
#' Sparse: 2 trials with 40 entries each
sparse_data <- list(group=rep(1:2,  each=40))

reg <- makeExperimentRegistry(
    file.dir = tempfile("sbc_", registry_tmp),
    ## use the default configured batchtools configuration batchtools.conf
    ## - found in the environemnt variable R_BATCHTOOLS_SEARCH_PATH
    ## - current working directory
    ## - $(HOME)/.batchtools.conf
    seed = 47845854,
    ## our worker functions and package loading
    source="sbc_tools.R"
)

## resources of each job: Less than 55min, 2000MB RAM and 1 cores
job_resources <- list(walltime=55, memory=2000, ncpus=1, max.concurrent.jobs=500)

if(FALSE) {
    ## for debugging here
    removeProblems("dense")
    removeProblems("sparse")
}

addProblem("dense",
           data=dense_data,
           fun=simulate_fake,
           seed=2345,
           cache=FALSE
           )

addProblem("sparse",
           data=sparse_data,
           fun=simulate_fake,
           seed=2346,
           cache=FALSE
           )

addAlgorithm("RBesT", fit_rbest)

## family, mean_mu, sd_mu, sd_tau, samp_sd
scenarios <- data.frame(
    family=c("binomial", "gaussian", "poisson"),
    mean_mu=c(-1, 0, 0),
    sd_mu=c(1),
    sd_tau=c(rep(0.5, 3), rep(1, 3)),
    samp_sd=c(1),
    stringsAsFactors=FALSE)

pdes <- list(sparse=scenarios, dense=scenarios)
ades <- list(RBesT=data.frame())

#' Add the defined problems and analysis methods to the registry and
#' set the number of replications:
S  <- 1E4L
addExperiments(pdes, ades, repls=S)

summarizeExperiments()

if(FALSE) {
    ## used for debugging
    job1 <- testJob(1)
    job1

    job2 <- testJob(6)
    job3 <- testJob(11)

    job <- makeJob(1)

    debug(fit_rbest)

    res  <- fit_rbest(dense_data, job, job$instance )
    res

    data <- dense_data
    instance  <- job$instance

    job1

}


#'
#' Chunk the jobs into 500 chunks to run
#'
ids <- getJobTable()
ids <- ids[, chunk:=chunk(job.id, 500)]

#' Once things run fine let's submit this work to the cluster.
auto_submit(ids, reg, job_resources)

#' Ensure that no error occured
assert_that(nrow(findErrors()) == 0)

#' Collect results.
calibration_data <- ijoin(
    ## grab job parameters
    unwrap(getJobPars()),
    unwrap(reduceResultsDataTable())
)

calibration_data[,algorithm:=NULL]

## collect sampler diagnostics
sampler_diagnostics <- calibration_data %>%
    group_by(family, problem, sd_tau) %>%
    summarize(N=n(),
              total_divergent=sum(n_divergent),
              min_ess=min(min_Neff),
              max_Rhat=max(max_Rhat),
              total_large_Rhat=sum(max_Rhat > 1.2),
              min_lp_ess_bulk=min(lp_ess_bulk),
              min_lp_ess_tail=min(lp_ess_tail))


cat("\nSampler diagnostics:\n\n")
kable(sampler_diagnostics, digits=3)
cat("\n")

if(sum(sampler_diagnostics$total_divergent) != 0) {
    warning("There were some divergent transitions!")
}
if(any(sampler_diagnostics$max_Rhat > 1.2) ) {
    warning("There were some parameters with large Rhat!")
}

#' Bin raw data as used in the analysis.
scale64  <- scale_ranks(1024, 2^4)
B <- 1024L / 2^4
calibration_data_binned <- calibration_data[, scale64(.SD), by=c("problem", "family", "sd_tau")]

#' Save as data.frame to avoid data.table dependency.
calibration_data <- as.data.frame(calibration_data)
calibration_data_binned <- as.data.frame(calibration_data_binned)

#' Further identification and verification data of run
git_hash <- system2("git", c("rev-parse", "HEAD"), stdout=TRUE)
created <- Sys.time()
created_str <- format(created, "%F %T %Z", tz="UTC")

calibration <- list(raw=calibration_data,
                    data=calibration_data_binned,
                    sampler_diagnostics = sampler_diagnostics,
                    S=S,
                    B=B,
                    git_hash=git_hash,
                    created=created)

saveRDS(calibration, file="calibration.rds")

library(tools)
md5 <- md5sum("calibration.rds")
cat(paste0("Created:  ", created_str, "\ngit hash: ", git_hash, "\nMD5:      ", md5, "\n"), file="calibration.md5")

#'
#' Summarize execution time
#'
job_report <- unwrap(getJobTable())
units(job_report$time.running)  <- "mins"

chunk_cols  <- c("job.id", "chunk")
job_report  <- job_report[ids[, ..chunk_cols], on="job.id", nomatch=0]

runtime_by_problem_family  <- job_report %>%
    group_by(family, problem) %>%
    summarize(total=sum(time.running), mean=mean(time.running), max=max(time.running))

runtime_by_problem  <- job_report %>%
    group_by(problem) %>%
    summarize(total=sum(time.running), mean=mean(time.running), max=max(time.running))

runtime  <- job_report %>%
    group_by(family) %>%
    summarize(total=sum(time.running), mean=mean(time.running), max=max(time.running))

runtime_by_problem_chunk  <- job_report %>%
    group_by(problem, chunk) %>%
    summarize(chunk_total=sum(time.running)) %>%
    summarize(total=sum(chunk_total), mean=mean(chunk_total), max=max(chunk_total))

cat("Summary on job runtime on cluster:\n\n")

cat("\nRuntime by problem and chunk:\n")
kable(runtime_by_problem_chunk, digits=2)

cat("\nRuntime by problem and family:\n")
kable(runtime_by_problem_family, digits=2)

cat("\nRuntime by family:\n")
kable(runtime, digits=2)

cat("\nRuntime by problem:\n")
kable(runtime_by_problem, digits=2)

duration_by_problem <- job_report %>%
    group_by(problem) %>%
        summarize(total=difftime(max(done), min(submitted), units="mins"))

duration <- job_report %>%
    summarize(total=difftime(max(done), min(submitted), units="mins"))

cat("\nDuration by problem:\n")
kable(duration_by_problem, digits=2)

cat("\nDuration:\n")
kable(duration, digits=2)

#' Cleanup
removeRegistry(0)

#' Session info
sessionInfo()
