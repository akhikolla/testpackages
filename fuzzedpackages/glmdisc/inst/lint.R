# Lint
if (!require(styler, quietly = TRUE)) {
     install.packages("styler")
}
docOrg <- lapply(dir("R", pattern="*.R", recursive = TRUE, full.names = TRUE), readLines)
styler::style_pkg()
docNew <- lapply(dir("R", pattern="*.R", recursive = TRUE, full.names = TRUE), readLines)
# Test
if (!identical(docOrg, docNew)) {
     stop("Lint was not performed, try running styler::style_pkg().")
}