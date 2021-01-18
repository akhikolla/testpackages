
.onLoad <- function(libname, pkgname) {
# 	loadRcppModules()
	loadModule(module="flan_module",what=TRUE,loadNow=TRUE)
}
