library(pkgbuild)
## ensure that the current dev version of OncoBayes2 is loaded
pkgbuild::compile_dll("../..")

library(knitr)
library(batchtools)
set.seed(453453)
## load utilities and current dev version of OncoBayes2
source("sbc_tools.R")

#' according to the docs this speeds up the reduce step
options(batchtools.progress = FALSE)

registry_tmp <- Sys.getenv("TMP_NFS", tempdir())

reg <- makeExperimentRegistry(
    file.dir = tempfile("sbc-", registry_tmp),
    ## use the default configured batchtools configuration batchtools.conf
    ## - found in the environemnt variable R_BATCHTOOLS_SEARCH_PATH
    ## - current working directory
    ## - $(HOME)/.batchtools.conf
    ## conf.file = NA,
    seed = 47845854,
    ## our worker functions and package loading
    source="sbc_tools.R")

## resources of each job: Less than 55min, 2000MB RAM and 2 cores
job_resources <- list(walltime=55, memory=2000, ncpus=2, max.concurrent.jobs=250)

if(FALSE) {
  ## for debugging here
  removeProblems("combo2_EX")
  removeProblems("combo2_NEX")
  removeProblems("combo2_EXNEX")
  removeProblems("base")
}

#' Evaluate dense and sparse data-scenario
source("sbc_example_models.R")

## family, mean_mu, sd_mu, sd_tau, samp_sd
scenarios <- data.frame(model = names(example_models),
                        stringsAsFactors=FALSE)

base_data  <- list(models = example_models)

addProblem("warmup_base",
           data = base_data,
           fun = simulate_fake,
           seed = 2345,
           ## caching speeds up the reduce step
           cache = TRUE
)

addAlgorithm("OncoBayes2", fit_exnex)


pdes_warmup <- list(warmup_base = scenarios)
ades <- list(OncoBayes2 = data.frame())

#' Add the defined problems and analysis methods to the registry and
#' set the number of replications:
S <- 1E4

## leads to 50 warmup infos per condition based on 100 chains
S_warmup <- 50

S_final <- S - S_warmup
addExperiments(pdes_warmup, ades, repls = S_warmup)

summarizeExperiments()

#'
#' Run first batch where we learn the warmup info
#'

#'
#' Number of jobs per chunk
#'
chunk_size <- 60

ids_warmup <- unwrap(getJobPars())

#'
#' run smaller chunk sizes at the beginning to quickly have the
#' desired warmup info
#'
ids_warmup[,chunk:=chunk(job.id, chunk.size=chunk_size/2)]

## submitJobs(warmup_ids)
## print(getStatus())

auto_submit(ids_warmup, reg, job_resources)

#'
#' Collect warmup info
#'
warmup_info  <- reduceResultsList(
    fun=function(run) {
        run[c("stepsize", "inv_metric", "draw")]
    }
)

warmup_info_by_model <- lapply(split(ids_warmup$job.id,  ids_warmup$model), function(mj) list(warmup_info=warmup_info[mj]))

#'
#' Schedule remaining runs with available warmup info
#'
base_data_final <- base_data
base_data_final$models  <- modifyList(base_data_final$models, warmup_info_by_model)

addProblem("base",
           data = base_data_final,
           fun = simulate_fake,
           seed = 2345 + nrow(ids_warmup),
           ## caching speeds up the reduce step
           cache = TRUE
)

pdes <- list(base = scenarios)

addExperiments(pdes, ades, repls = S_final)

summarizeExperiments()


#'
#' Chunk the jobs into packs of size 20 chunks to run
#'
ids <- getJobTable()
ids_main <- ids[problem=="base"]
ids_main[, chunk:=chunk(job.id, chunk.size=chunk_size)]

num_jobs <- nrow(ids)

#' Once things run fine let's submit this work to the cluster.
##submitJobs(ids, job_resources)
##print(getStatus())

#' This function deals with unstable nodes in the cluster
auto_submit(ids_main, reg, job_resources)

#' Ensure that no error occured
assert_that(nrow(findErrors()) == 0)

#' Collect results.
calibration_data <- ijoin(
  ## grab job parameters
  unwrap(getJobPars()),
  unwrap(reduceResultsDataTable(fun=function(x) c(rank=as.list(x$rank),
                                                  list(n_divergent=x$n_divergent,
                                                       min_Neff=x$min_Neff,
                                                       max_Rhat=x$max_Rhat,
                                                       lp_ess_bulk = x$lp_ess_bulk,
                                                       lp_ess_tail = x$lp_ess_tail
                                                       )) ))
)

## check that indeed all jobs have finished
assert_that(nrow(calibration_data) == num_jobs)

## collect sampler diagnostics
sampler_diagnostics <- calibration_data %>%
    group_by(model, problem) %>%
    summarize(N=n(),
              total_divergent=sum(n_divergent),
              min_ess=min(min_Neff),
              max_Rhat=max(max_Rhat),
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

# there is only one algorithm. remove that column.
calibration_data <- calibration_data %>% select(-algorithm, -problem)

#' Bin raw data as used in the analysis.
B <- 1024L / 2^5

rank_params <- names(calibration_data)[grepl(names(calibration_data), pattern = "rank")]

# calibration_data_binned <- calibration_data[, scale64(.SD), by=c("problem", "model", params)]

calibration_data_binned <- calibration_data %>%
  mutate_at(.vars = rank_params, .funs = function(x) ceiling((x + 1) / (1024 / B) - 1))

names(calibration_data_binned) <- gsub(
  names(calibration_data_binned),
  pattern = "rank[.]",
  replacement = ""
)

params <- gsub(rank_params, pattern = "rank[.]", replacement = "")

calibration_binned <- calibration_data_binned %>%
  dplyr::select(-job.id, - n_divergent, - min_Neff) %>%
  tidyr::gather(key = "param", value = "bin", - model) %>%
  group_by(model, param, bin) %>%
  tally() %>%
  right_join(
    expand.grid(
      model = unique(calibration_data_binned$model),
      param = params,
      bin = 0:(B - 1),
      stringsAsFactors = FALSE
    ),
    c("model", "param", "bin")
  ) %>%
  replace_na(list(n = 0)) %>%
  arrange(model, param, bin) %>%
  spread(key = param, value = n)



#' Save as data.frame to avoid data.table dependency.
calibration_data <- as.data.frame(calibration_data)
calibration_binned <- as.data.frame(calibration_binned)

#' Further identification and verification data of run
git_hash <- system2("git", c("rev-parse", "HEAD"), stdout=TRUE)
created <- Sys.time()
created_str <- format(created, "%F %T %Z", tz="UTC")

calibration <- list(raw = calibration_data,
                    data = calibration_binned,
                    sampler_diagnostics = sampler_diagnostics,
                    S = S,
                    B = B,
                    git_hash = git_hash,
                    created = created)

saveRDS(calibration, file = "calibration.rds")

library(tools)
md5 <- md5sum("calibration.rds")
cat(paste0("Created:  ", created_str, "\ngit hash: ", git_hash, "\nMD5:      ", md5, "\n"),
    file="calibration.md5")


#'
#' Summarize execution time
#'
job_report <- unwrap(getJobTable())
units(job_report$time.running)  <- "mins"

chunk_cols  <- c("job.id", "chunk")
job_report  <- rbind(
    job_report[ids_warmup[, ..chunk_cols], on="job.id", nomatch=0],
    job_report[ids_main[, ..chunk_cols], on="job.id", nomatch=0]
)

runtime_by_problem_model  <- job_report %>%
    group_by(model, problem) %>%
    summarize(total=sum(time.running), mean=mean(time.running), max=max(time.running))

runtime_by_problem  <- job_report %>%
    group_by(problem) %>%
    summarize(total=sum(time.running), mean=mean(time.running), max=max(time.running))

runtime  <- job_report %>%
    group_by(model) %>%
    summarize(total=sum(time.running), mean=mean(time.running), max=max(time.running))

runtime_by_problem_chunk  <- job_report %>%
    group_by(problem, chunk) %>%
    summarize(chunk_total=sum(time.running)) %>%
    summarize(total=sum(chunk_total), mean=mean(chunk_total), max=max(chunk_total))

cat("Summary on job runtime on cluster:\n\n")

cat("\nRuntime by problem and chunk:\n")
kable(runtime_by_problem_chunk, digits=2)

cat("\nRuntime by problem and model:\n")
kable(runtime_by_problem_model, digits=2)

cat("\nRuntime by model:\n")
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
