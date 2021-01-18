#' @title Simulate character displacement data
#' @description Simulates data under a Brownian motion or character displacement model
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param n.steps Number of time steps the for the simulation (default = 1000 time steps).
#' @param sigma The value of Brownian variance in the simulation
#' @param a The strength of competition between inter-specific lineages
#' @param ntraits Number of traits to be simulated
#' @param sympatry an optional matrix giving the time that each pair of species starts to interact
#' @param allopatry an optional matrix giving the times when species stop interacting
#' @param trait.lim an optional parameter that puts limits on the available trait-space, preventing trait values with magnitude greater than the value of lim
#' @return A list containing the the simulated data (tval) showing the sigma, a, mean gap and gap standard deviation. Additionally, if used, the user input sympatry (symp) and/or allopatry (allo) matrices
#' @useDynLib motmot
#' @importFrom Rcpp sourceCpp
#' @seealso \code{\link{chr.disp.param}}, \code{\link{chr.disp.lrt}}
#' @references Clarke M, Thomas GH, Freckleton RP. 2017. Trait evolution in adaptive radiations: modelling and measuring interspecific competition on phylogenies. The American Naturalist. 189, 121-137.
#' @author Magnus Clarke and Mark Puttick
#' @examples
#' ## import finch data form Clarke et al. (2017)
#' data(finches)
#' emp.tree <- finch.tree
#' emp.data <- finch.data
#' ## simulate small amount of data 
#' ## (example only - many more datasets are required for accuracy)
#' sim.data <- chr.disp.sim(emp.tree, n.steps=100,
#' sigma=1, a=2, ntraits=1, sympatry=NA, allopatry=NA, trait.lim=NA)
#' @export

chr.disp.sim <- function(phy, n.steps=1000, sigma=1, a=0, ntraits=1, sympatry=NA, allopatry=NA, trait.lim=NA) {
    
    phy.names <- phy$tip.label
    phy.chr.disp <- checktree(phy)
    num_tips <- length(phy.chr.disp$data_order)
    splitting_nodes <- phy.chr.disp$splitting_nodes - 1
    times <- phy.chr.disp$times
    tval <- rep(0, num_tips * ntraits)
    dt <- nodeTimes(phy)[1,1] / n.steps

    symp <- vectortime(sympatry, 0, num_tips)
    allo <- vectortime(allopatry, 99e9, num_tips )
    if(is.na(trait.lim)) trait.lim <- 9e99

    result <- .C("pathsim", ntip=as.integer(num_tips), dt=as.double(dt), rate = as.double(sigma^2), a=as.double(a), r_intervals=as.double(times), splitters=as.integer(splitting_nodes), tval = as.double(tval), ntraits=as.integer(ntraits), symp=as.double(symp), allo=as.double(allo), lim=as.numeric(trait.lim))

    result$tval <- as.matrix(reorder_data(phy.chr.disp, result$tval, ntraits))
    rownames(result$tval) <- phy.names
    if(all(result$symp == 0)) result$symp <- NULL
    if(all(result$allo == 0)) result$allo <- NULL
    result
}
