# original from https://github.com/ropensci/rfishbase/blob/master/tests/testthat/helper_rfishbase.R
# Tests that contain this function will not be run if any of these conditions fail:
needs_api <- function() {
  skip_on_cran()  
}

# checks if fishbase is not working due to website down or internet connection
not_working <- function(url = "https://www.fishbase.us"){
  test <- try(suppressWarnings(readLines(url, n = 1)), silent = TRUE)
  inherits(test, "try-error")
}

# skip test if fishbase site cannot be reached
check_api <- function() {
  if (not_working()) {
    skip("API not available")
  }
}
