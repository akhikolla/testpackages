#' @importFrom gdtools raster_view
NULL

svglite_path <- function(...) {
  file.path("R", "svglite", ...)
}

for (file in list.files(svglite_path(), pattern = "*.R")) {
  source(svglite_path(file), local = TRUE)
}
