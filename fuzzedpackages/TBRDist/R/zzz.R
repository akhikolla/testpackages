#' @useDynLib TBRDist
.onUnload <- function (libpath) {
  library.dynam.unload("TBRDist", libpath)
}

# Suppress "NOTE: Nothing imported from Rdpack":
#' @importFrom Rdpack reprompt
NULL

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you checked the output of TBRDist.Rmd?",
    "Have you updated the version number in NEWS.md & DESCRIPTION?"
  )
}

# Additional tests:
#
# spell_check()
# pkgdown::build_reference_index()
# check_win_devel(); rhub::check_for_cran()
# codemetar::write_codemeta()
# # revdepcheck::revdep_check()
