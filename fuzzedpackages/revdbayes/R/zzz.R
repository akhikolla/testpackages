.onUnload <- function (libpath) {
  library.dynam.unload("revdbayes", libpath)
}

.onAttach <- function(...) {
  ggplot2::theme_set(bayesplot::theme_default())
}
