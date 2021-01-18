## preliminary code

# .printeSaddleVersion <- function()
# { library(help=synlik)$info[[1]] -> version
#   version <- version[pmatch("Version",version)]
#   um <- strsplit(version," ")[[1]]
#   version <- um[nchar(um)>0][2]
#   cat(paste("This is \"esaddle\" ",version,"\n"),sep="")
# }
# 
# 
# .onAttach <- function(...) { 
#   .printeSaddleVersion()
# }

# .onUnload <- function(libpath) library.dynam.unload("esaddle", libpath)
# 
# .onLoad <- function(lib,pkg) {
#    library.dynam("esaddle", pkg, lib)
# }
