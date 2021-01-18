
# @description This function confirms that 
# @param d A 'macroACS' class object to check
# @section Details:
# The 'macroACS' class should be a list of length 6 as noted below:
  # endyear: an integer scalar
  # span: an integer scalar in {1,3,5}
  # estimates: a list of data.frame's containing only numeric elements
  # standard_error: a list of data.frame's containing only numeric elements
  # geography: a data.frame describing the geography from \code{\link[acs]{geo.make}}
  # geo_title: a formall class 'geo' object from \code{\link[acs]{geo.make}}
confirm_macroACS_class <- function(d) {
  # 1. confirm correct class structure
  testthat::expect_true(is.list(d))
  testthat::expect_true(length(d) == 6)
  testthat::expect_equal(names(d), c("endyear", "span", "estimates", "standard_error", "geography",
                                "geo_title"))
  testthat::expect_true(is.numeric(d$endyear))
  testthat::expect_true(is.numeric(d$span))
  testthat::expect_true(d$span %in% c(1L, 3L, 5L))
  testthat::expect_true(all(unlist(lapply(d$estimates, is.data.frame))))
  testthat::expect_true(all(unlist(lapply(d$standard_error, is.data.frame))))
  testthat::expect_true(is.list(d$geo_title))
  testthat::expect_true(is.data.frame(d$geography))
  
  # 2. confirm conforming data elements
  testthat::expect_true(all(unlist(
    mapply(function(est, se) {all.equal(dim(est), dim(se))}, 
           est= d$estimates, se= d$standard_error, SIMPLIFY = FALSE)
    )))
  testthat::expect_true(all(unlist(
    mapply(function(est, se) {all.equal(names(est), names(se))}, 
           est= d$estimates, se= d$standard_error, SIMPLIFY = FALSE)
  )))
  testthat::expect_true(all(unlist(lapply(d$estimates, function(l) {
    all(unlist(lapply(l, is.numeric)))
  }))))
  testthat::expect_true(all(unlist(lapply(d$standard_error, function(l) {
    all(unlist(lapply(l, is.numeric)))
  }))))
  
  # 3. Confirm sorted
  alpha <- order(d$geography$NAME)
  alpha_e <- order(rownames(d$estimates[[1]]))
  alpha_se <- order(rownames(d$standard_error[[1]]))
  nr <- length(d$geography$NAME)
  testthat::expect_equal(alpha, 1:nr)
  testthat::expect_equal(alpha, alpha_e)
  testthat::expect_equal(alpha, alpha_se)
}