ManifoldOptim.system.file <- function(...){
	tools::file_path_as_absolute( base::system.file( ..., package = "ManifoldOptim" ) )
}

ManifoldOptim.LdFlags <- function() {
	ManifoldOptim.system.file("libs/ManifoldOptim.so")
}

LdFlags <- function() {
	cat(ManifoldOptim.LdFlags())
}

