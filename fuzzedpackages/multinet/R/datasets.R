loadModule("multinet",TRUE)

# Functions loading some well-known small datasets
# Larger datasets can be downloaded from the book webpage: multilayer.it.uu.se

ml_aucs <- function() {
	read_ml(system.file("extdata", "aucs.mpx", package="multinet"),"AUCS")
}

ml_bankwiring <- function() {
    read_ml(system.file("extdata", "bankwiring.mpx", package="multinet"),"Bank Wiring")
}

ml_florentine <- function() {
	read_ml(system.file("extdata", "florentine.mpx", package="multinet"),"Florentine families")
}

ml_monastery <- function() {
    read_ml(system.file("extdata", "monastery.mpx", package="multinet"),"Monastery")
}

ml_tailorshop <- function() {
    read_ml(system.file("extdata", "tailorshop.mpx", package="multinet"),"Tailorshop")
}

ml_toy <- function() {
    read_ml(system.file("extdata", "book.mpx", package="multinet"),"toy")
}
