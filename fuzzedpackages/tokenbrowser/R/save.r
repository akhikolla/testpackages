#' Wrap html body in the template and save
#'
#' @param data     The html body data
#' @param template The html header/footer template
#' @param filename The name of the file to save the html. Default is a temp file
#'
#' @return The (local) url to the html file
save_html <- function(data, template, filename=NULL) {
  if (is.null(filename)) {
    filename = tempfile("tokenbrowser_", fileext = ".html")
    #message("Writing html to ", filename)
  }
  sink(filename)
  cat(template$header)
  cat(data)
  cat(template$footer)
  sink()
  filename
}
