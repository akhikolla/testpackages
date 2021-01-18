#' Build a Book without Underscores
#' 
#' @description 
#' Since the use of underscores ('_') is not permitted when streaming 
#' \strong{bookdown} documents via 
#' \href{https://pages.github.com/}{GitHub Pages}, this wrapper function serves 
#' to remove any unwanted underscores from subfolders and link \code{.html} 
#' documents created by \code{\link[bookdown]{render_book}}. 
#'  
#' @param output_dir Output directory as \code{character}. 
#' @param ... Arguments passed to \code{\link[bookdown]{render_book}}.
#' 
#' @seealso \code{\link[bookdown]{render_book}}.
#' 
#' @note 
#' While all remaining arguments passed to \code{\link[bookdown]{render_book}} 
#' remain untouched, and hence, their specification is freely up to the user, 
#' the default value of 'output_dir' is explicitly set to \code{"book"} here. If 
#' this were not the case (i.e. if the default value were used), the output 
#' document would be created in \code{"_book"} which is not desirable for 
#' obvious reasons.
#' 
#' @author Florian Detsch
#' 
# @examples 
# \dontrun{
# buildBook(input = "index.Rmd"
#           , output_format = "bookdown::gitbook"
#           , output_dir = "book")
# }
#' 
#' @export buildBook
#' @name buildBook
buildBook = function(output_dir = "book", ...) {
  
  if (is.null(output_dir))
    stop("Output directory must be other than 'NULL'.\n")
    
  ## if a previous build failed, '_main.Rmd' needs to be removed manually
  if (file.exists("_main.Rmd")) 
    jnk = file.remove("_main.Rmd")
  
  bookdown::render_book(output_dir = output_dir, ...)
  
  ## remove leading underscore ("_*") from image links 
  html = list.files(output_dir, pattern = ".html$", full.names = TRUE)
  
  for (i in html) {
    lns = readLines(i)
    file.remove(i)
    
    lns = gsub("_main_files/", "main_files/", lns)
    writeLines(lns, i)
  }
  
  ## remove leading underscore ("_*") from image subfolder
  mf = file.path(output_dir, "main_files")
  
  if (dir.exists(mf)) {
    jnk = suppressWarnings(
      try(unlink(mf, recursive = TRUE), silent = TRUE)
    )
    
    if (jnk != 0 | inherits(jnk, "try-error"))
      stop("Something went wrong when trying to delete ", mf, ".\n")
  }
  
  jnk = file.rename(file.path(output_dir, "_main_files")
                    , file.path(output_dir, "main_files"))

  cat("Process finished successfully.\n")
}

