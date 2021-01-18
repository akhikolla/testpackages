# Script to perform the Operation Qualification
# ---------------------------------------------

library(testthat)
library(OncoBayes2)

## enforce that all tests are run
Sys.setenv(NOT_CRAN="true")

cat("TEST RUN DATE:", date(), "\n")

cat("TESTING PACKAGE:\n")
print(packageDescription("OncoBayes2"))

cat("RUNNING PACKAGE TESTS:\n")
# Run each section separately to get subsequent numbering per section
# of the TAP reporter; execution order must be aligned with steps described 
# in the vignette
for(test in c("blrm_exnex", "blrm_trial", "examples", "posterior", "sbc")) {
    test_package("OncoBayes2", filter=test, reporter="tap")
}

# Finally run all tests once more, but with the stop reporter. This
# ensures that the last line of this script is only displayed if and
# only if all tests run successful
test_package("OncoBayes2", reporter="stop")

cat("\n\nR SESSION INFO:\n")

print(sessionInfo())

cat("\nTEST FINISH DATE:", date(), "\n")
cat("\n\nALL TESTS SUCCESSFUL\n")

