getBlocksLine <- function(mask, nblock){
    nvertex <- sum(mask == 1)
    lapply(1:nblock, function(x) seq(x, nvertex, by=nblock))
}

getBlocksPlane <- function(mask, nblock){
    ind <- which(mask == 1, arr.ind=T)
    oe1 <- ind[,1] %% 2
    oe2 <- ind[,2] %% 2
    if(nblock==2){
        mark <- abs(oe1 - oe2)
        lapply(1:2, function(x) which(mark == (x-1)))
    }
    else
        lapply(1:4, function(x) which((oe1 + oe2*2) == (x-1)))
}

getBlocksCube <- function(mask, nblock){
    ind <- which(mask == 1, arr.ind=T)
    x.mark <- ind[,1] %% 2
    y.mark <- ind[,2] %% 2
    z.mark <- ind[,3] %% 2

    if(nblock==2){
        BW <- ifelse(rowSums(cbind(x.mark, y.mark, z.mark)) %% 2 == 1, 1, 0)
        lapply(1:2, function(x) which(BW == (x-1)))
    }

    else
        lapply(1:8, function(x) which((x.mark + y.mark*2 + z.mark*4) == (x-1)))
}

#' Get Blocks of a Graph
#' 
#'   Obtain blocks of vertices of a 1D, 2D, or 3D graph, in order to use
#'   the conditional independence to speed up the simulation (chequerboard idea). 
#'
#' @param mask a vector, matrix, or 3D array specifying vertices of a graph. Vertices of value 1 are within the graph and 0 are not.
#' @param nblock a scalar specifying the number of blocks. For a 2D graph \code{nblock} could be either 2 or 4, and for a 3D graph \code{nblock} could be either 2 or 8.
#'
#' @return A list with the number of components equal to \code{nblock}. Each component consists of vertices within the same block.
#' @export
#' @details The vertices within each block are mutually independent given the vertices in other blocks. Some blocks could be empty.
#' @references
#'     Wilkinson, D. J. (2005)
#'     "Parallel Bayesian Computation"
#'     \cite{Handbook of Parallel Computing and Statistics}, pp. 481-512
#'     \emph{Marcel Dekker/CRC Press}
#' @keywords spatial
#' @usage getBlocks(mask, nblock)
#' @examples 
#'   #Example 1: split a line into 2 blocks
#'   getBlocks(mask=c(1,1,1,1,0,0,1,1,0), nblock=2)
#'   
#'   #Example 2: split a 4*4 2D graph into 4 blocks in order
#'   #           to use the chequerboard idea for a neighbourhood structure
#'   #           corresponding to the second-order Markov random field.
#'   getBlocks(mask=matrix(1, nrow=4, ncol=4), nblock=4)
#'   
#'   #Example 3: split a 3*3*3 3D graph into 8 blocks
#'   #           in order to use the chequerboard idea for a neighbourhood
#'   #           structure based on the 18 neighbors definition, where the
#'   #           neighbors of a vertex comprise its available
#'   #           adjacencies sharing the same edges or faces.
#'   mask <- array(1, dim=rep(3,3))
#'   getBlocks(mask, nblock=8)
getBlocks <- function(mask, nblock){
    if(!(is.vector(mask) || is.matrix(mask) || length(dim(mask))==3))
        stop("The graph has to be in 1D, 2D, or 3D.")
    if(nblock <= 0 || floor(nblock) != nblock)
        stop("The number of blocks has to be and integer bigger than 0")
    if(length(mask==1) < nblock)
        stop("The number of blocks has to be less than
             the number of vertices inside the mask.")
    if(is.matrix(mask) && !nblock %in% c(2, 4))
        stop("For a plane, the number of blocks has to be equal to 2 or 4.")
    if(length(dim(mask))==3 && !nblock %in% c(2, 8))
        stop("For a cube, the number of blocks has to be equal to 2 or 8.")

    if(is.vector(mask))
        getBlocksLine(mask, nblock)

    else if(is.matrix(mask))
        getBlocksPlane(mask, nblock)

    else if(length(dim(mask))==3)
        getBlocksCube(mask, nblock)
}
