# Package-wide global variables
.dmbcEnv <- new.env()

#' @importFrom utils globalVariables
#' @importFrom tools file_path_as_absolute
.onLoad <- function(lib, pkg){
  # needed to avoid annoying notes in R CMD CHECK
  # (see https://github.com/tidyverse/magrittr/issues/29)
  if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(".", "S", "cl", "italic", "lbl", "p_i", "p_j", "p", "G", "DCIC"))
  }
  .dmbcEnv$path.to.me <- tools::file_path_as_absolute(lib)
  .dmbcEnv$nlog.double.eps <- -log(.Machine[["double.eps"]])
  .dmbcEnv$allowedfamilies <- c("binomial", "multinomial", "gaussian", "poisson")
  .dmbcEnv$current_p <- 2
  .dmbcEnv$current_G <- 3
  .dmbcEnv$current_family <- "binomial"
}

.onAttach <- function(lib, pkg) {
  packageStartupMessage(sprintf("Package %s (%s) loaded.\nTo cite, type citation(\"%s\")",
    pkg, utils::packageDescription(pkg)$Version, pkg))
}

.onUnload <- function(lib) {
	library.dynam.unload("dmbc", lib)
}
