#' Set working directory dependent on current OS
#' 
#' @description
#' Similar to \code{\link{setwd}}, this function sets the working directory to a
#' user-defined path. Rather than supplying a single 'dir' argument, however, 
#' both an OS-sensitive path to the desired hard disk partition and, optionally, 
#' an extension of this file path are required.  
#' 
#' @param lin,win Absolute file paths to the Linux and Windows partition as 
#' \code{character}.
#' @param ext Optional file path extension as \code{character} that will be 
#' added to 'lin' or 'win' after automatic OS determination.
#' 
#' @author 
#' Florian Detsch
#' 
#' @seealso
#' \code{\link{setwd}}, \code{\link{switch}}
#' 
#' @examples
#' \dontrun{
#' # desired partition
#' setwdOS()
#' 
#' # including file path extension
#' setwdOS(ext = "kilimanjaro/nubiscope")
#' }
#' 
#' @export setwdOS
#' @name setwdOS
setwdOS <- function(lin = "/media/permanent/", 
                    win = "C:/", 
                    ext = NULL) {
  
  ## determine os
  ch_dir_os <- switch(Sys.info()[["sysname"]], 
                      "Linux" = lin, 
                      "Windows" = win)
  
  ## setwd
  if (is.null(ext))
    setwd(ch_dir_os)
  else 
    setwd(paste(ch_dir_os, ext, sep = "/"))
  
  return(invisible())
}