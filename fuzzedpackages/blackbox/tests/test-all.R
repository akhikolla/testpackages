if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") { ## not on CRAN
  if(require("testthat", quietly = TRUE)) {
    pkg   <- "blackbox"
    require(pkg, character.only=TRUE, quietly=TRUE)
    ## test_package(pkg) ## for an installed package
    if (interactive()) {
      if (FALSE) { ## tests not included in package (using unpublished data, etc.)
        oldpath <- getwd() 
        ptpath <- paste0(projpath(),"/package/tests_private/")
        subpaths <- dir(ptpath,full.names = TRUE)
        virgin_opts <- blackbox.options()
        priv_timings <- list()
        # the migraine scripts have rm(list=ls()) hence we copy some variable in the options 'keep_over_blackbox_tests' and 'blackbox_path'
        options(keep_over_blackbox_tests=mget(c("oldpath","subpaths","virgin_opts","priv_timings"),.GlobalEnv))
        for (st in getOption("keep_over_blackbox_tests")$subpaths) {
          options("blackbox_path"=st)
          setwd(st)
          blackbox.options(getOption("keep_over_blackbox_tests")$virgin_opts)
          priv_timing <- system.time(source("migraine_1.R"))
          keep_over_blackbox_tests <- getOption("keep_over_blackbox_tests")
          keep_over_blackbox_tests$priv_timings[[getOption("blackbox_path")]] <- priv_timing
          options("keep_over_blackbox_tests"=keep_over_blackbox_tests)
        }
        setwd(getOption("keep_over_blackbox_tests")$oldpath)
        #priv_testfiles <- dir(ptpath,pattern="*.R",full.names = TRUE)
        #priv_timings <- t(sapply(priv_testfiles, function(fich){system.time(source(fich))}))
        print(getOption("keep_over_blackbox_tests")$priv_timings)
      }
    } else {    
      report <- test_check(pkg) ## for R CMD check ## report is NULL...
      print(warnings()) # TODO? catch most of these by expect_warning(..)
    }  
  } else {
    cat( "package 'testthat' not available, cannot run unit tests\n" )
  }
}
