skip_on_macOS <- function()
  skip_if(Sys.info()["sysname"] == "Darwin")
