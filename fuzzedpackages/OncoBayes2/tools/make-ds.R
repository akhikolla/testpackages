# Dataset preparation
# -------------------
#

make_ds <- function() {

    examples  <- new.env()
    flist <- list.files("data-raw", pattern = "*[.]R")
    for(f in flist) source(file.path("data-raw", f), examples)

    ex <- names(examples)
    for(d in ex[!ex %in% c("df", "group_defs", "Combo3Data")]){
      assign(d, examples[[d]], envir = .GlobalEnv)
    }

    use_data(hist_combo3, # testdata2, testdata3,
             hist_combo2,
             codata_combo2,
             hist_SA,
             drug_info_combo2,
             dose_info_combo2, # Add further data sets here (separated by comma)
             
             overwrite=TRUE)
  # use_data will make a file called testdata.rda available in data/

}


# Pepare some internal example data sets

make_internal_ds <- function() {

    example_directory <- "man-roxygen"
    all_examples  <- list.files(example_directory,
                                pattern="^example-.*.R$", full.names=TRUE)
    example_cache <<- list()
    for(ex in all_examples) {
        id <- gsub("\\.R$", "", gsub("^example-", "", basename(ex)))
        cat("Saving example", id, "\n")
        example_cache[[id]] <<- new.env()
        ex_str <- readLines(ex)
        ex_str <- gsub("^#\'", "", ex_str)
        ## filter out @examples, dontrun and trailing }
        ex_str <- grep("@examples", ex_str, value=TRUE, invert=TRUE)
        ex_str <- grep("dontrun", ex_str, value=TRUE, invert=TRUE)
        ex_str <- grep(" \\}$", ex_str, value=TRUE, invert=TRUE)
        example_cache[[id]] <<- ex_str
    }

    calibration  <- readRDS("inst/sbc/calibration.rds")

    calibration_meta <- calibration[c("S", "B", "git_hash", "created")]
    calibration_data <- calibration$data

    ## reshape calibration data
    calibration_data <- calibration$data %>%
        tidyr::gather(key = "param", value = "count", - model, -bin) %>%
        arrange(param, bin) %>%
        group_by(model, param) %>%
        mutate(ecdf = cumsum(count) / calibration$S,
               ecdf_ref = (bin + 1) / calibration$B,
               allzero = all(ecdf == 0)) %>%
        filter(!allzero) %>%
        select(-allzero, -ecdf, -ecdf_ref) %>%
        as.data.frame()

    calibration_md5  <- strsplit(readLines("inst/sbc/calibration.md5"), ": +")
    vals  <- sapply(calibration_md5, function(x) { x[[2]] } )
    keys  <- sapply(calibration_md5, function(x) { x[[1]] } )
    names(vals)  <- keys
    calibration_meta["MD5"]  <- vals["MD5"]

    pkg_create_date  <- Sys.time()

    use_data(
        calibration_data,
        calibration_meta,
        pkg_create_date,
        example_cache,
        internal=TRUE, ## Suppress external use: Data will be accessible to
                       ## the project code, but not to the users
        overwrite=TRUE)

}


library(devtools)
library(dplyr)
library(tibble)
library(tidyr)

make_ds()
make_internal_ds()
