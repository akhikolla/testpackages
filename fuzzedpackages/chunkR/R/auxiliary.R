
#' chunkR devel
#' @description The function opens the chunkR-devel web site:
#' https://github.com/leandroroser/chunkR
#' @export

chunkR_devel <- function(){
  cat("Opening link: https://leandroroser.github.io/chunkR\n")
  browseURL("https://github.com/leandroroser/chunkR")
}

is_meta <- function(X) {
  meta <- c("\\.", "\\\\", "\\|", "\\[", "\\]", "\\{", "\\}", 
            "\\(", "\\)", "\\^", "\\*", "\\?", "\\+", "\\$")
  any(meta %in% paste("\\", X, sep = ""))
}