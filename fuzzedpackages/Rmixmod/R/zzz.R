.onAttach <-function(lib,pkg){
  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage("Rmixmod v. ", as.character(ver), " / URI: www.mixmod.org")
}
