
fertilizers <- function() {
	f <- system.file("extdata/fertilizers.csv", package="Rquefts")
	utils::read.csv(f)
}



nutrientRates <- function(supply, treatment) {
#	result <- matrix(nrow=ncol(supply)-1, ncol=ncol(treatment)-1)
	NPK <- supply[, c("N", "P", "K")] * treatment / 100
	apply(NPK, 2, sum)
}

