pre <- function(mask){
    n <- sum(mask==1)
    if(is.vector(mask)){
        mask <- c(0, mask, 0)
        mask.new <- cumsum(mask)
        mask.new[mask==0] <- n + 1
    }
    else if(is.matrix(mask)){
        nr <- nrow(mask)
        nc <- ncol(mask)
        mask.new0 <- matrix(0, nr + 2, nc + 2)
        mask.new0[2:(nr+1), 2:(nc+1)] <- mask
        mask.new <- cumsum(mask.new0)
        mask.new[mask.new0==0] <- n + 1
    }
    else{
        nr <- dim(mask)[1]
        nc <- dim(mask)[2]
        nz <- dim(mask)[3]
        mask.new0 <- array(0, dim=c(nr + 2, nc + 2, nz + 2))
        mask.new0[2:(nr+1), 2:(nc+1), 2:(nz + 1)] <- mask
        mask.new <- cumsum(mask.new0)
        mask.new[mask.new0==0] <- n + 1
    }
    list(mask=mask.new, n=n)
}

get2NeighborsLine <- function(mask, n){
    focus <- which(mask <= n)
    cbind(mask[focus-1], mask[focus+1])
}

get2NeighborsPlane <- function(mask, n, nr, neiStruc){
    neighbors <- vector("list", 4)
    focus <- which(mask <= n)
    if(neiStruc[1]==T)
        neighbors[[1]] <- cbind(mask[focus-1], mask[focus+1])
    if(neiStruc[2]==T)
        neighbors[[2]] <- cbind(mask[focus-nr], mask[focus+nr])
    if(neiStruc[3]==T)
        neighbors[[3]] <- cbind(mask[focus-nr-1], mask[focus+nr+1])
    if(neiStruc[4]==T)
        neighbors[[4]] <- cbind(mask[focus-nr+1], mask[focus+nr-1])
    neighbors
}

get2NeighborsCube <- function(mask, n, nr, nc, neiStruc){
    neighbors <- vector("list", 12)
    focus <- which(mask <= n)
    n12 <- nr*nc
    if(neiStruc[1]==T)
        neighbors[[1]] <- cbind(mask[focus-1], mask[focus+1])
    if(neiStruc[2]==T)
        neighbors[[2]] <- cbind(mask[focus-nr], mask[focus+nr])
    if(neiStruc[3]==T)
        neighbors[[3]] <- cbind(mask[focus-nr-1], mask[focus+nr+1])
    if(neiStruc[4]==T)
        neighbors[[4]] <- cbind(mask[focus-nr+1], mask[focus+nr-1])

    if(neiStruc[5]==T)
        neighbors[[5]] <- cbind(mask[focus-1], mask[focus+1])
    if(neiStruc[6]==T)
        neighbors[[6]] <- cbind(mask[focus-n12], mask[focus+n12])
    if(neiStruc[7]==T)
        neighbors[[7]] <- cbind(mask[focus-n12-1], mask[focus+n12+1])
    if(neiStruc[8]==T)
        neighbors[[8]] <- cbind(mask[focus-n12+1], mask[focus+n12-1])

    if(neiStruc[9]==T)
        neighbors[[9]] <- cbind(mask[focus-nr], mask[focus+nr])
    if(neiStruc[10]==T)
        neighbors[[10]] <- cbind(mask[focus-n12], mask[focus+n12])
    if(neiStruc[11]==T)
        neighbors[[11]] <- cbind(mask[focus-n12-nr], mask[focus+n12+nr])
    if(neiStruc[12]==T)
        neighbors[[12]] <- cbind(mask[focus-n12+nr], mask[focus+n12-nr])


    neighbors
}

get2NeighborsCube3DDiad <- function(mask){
    expand <- matrix(c(-1, -1, -1,
                       -1, +1, +1,
                       -1, +1, -1,
                       -1, -1, +1,
                       +1, -1, -1,
                       +1, +1, +1,
                       +1, +1, -1,
                       +1, -1, +1), ncol=8)
    browser()
    n <- sum(mask)
    maskn <- array(cumsum(mask), dim=dim(mask))
    maskn[mask==0] <- n + 1
    focus <- which(mask==1, arr.ind=TRUE)
    focus <- t(focus)
    neighbors <- do.call(cbind, lapply(1:8, function(i)
                         maskn[t(focus + expand[,i])])) 
    
}

getNeighborsB <- function(neiStruc, neighbors2){
    nvertex <- nrow(neighbors2)
    neighbors <- matrix(nvertex+1, nrow=nvertex+1, ncol=neiStruc)
    neighbors[1:nvertex,1:2] <- neighbors2
    if(neiStruc > 2)
        for(i in 1:(neiStruc/2-1)){
            neighbors[,i*2+1] <- neighbors[neighbors[,(i-1)*2+1],1]
            neighbors[,i*2+2] <- neighbors[neighbors[,(i-1)*2+2],2]
        }
    neighbors[-(nvertex+1),]
}

getNeighborsLine <- function(mask, neiStruc){
    r <- pre(mask)
    neighbors2 <- get2NeighborsLine(r$mask, r$n)
    getNeighborsB(neiStruc, neighbors2)

}

getNeighborsPlane <- function(mask, neiStruc){
    r <- pre(mask)
    nr <- nrow(mask) + 2
    neighbors2 <- get2NeighborsPlane(r$mask, r$n, nr, neiStruc > 0)
    if(all(neiStruc <= 2))
        neighbors <- do.call(cbind, neighbors2)
    else{
        neighbors <- NULL
        for(i in 1:4){
            if(!is.null(neighbors2[[i]])){
                new <- getNeighborsB(neiStruc[i], neighbors2[[i]])
                neighbors <- cbind(neighbors, new)
            }
        }
    }
    neighbors
}

getNeighborsCube <- function(mask, neiStruc){
    r <- pre(mask)
    nr <- dim(mask)[1] + 2
    nc <- dim(mask)[2] + 2
    neiStruc <- as.vector(t(neiStruc))
    neighbors2 <- get2NeighborsCube(r$mask, r$n, nr, nc, neiStruc > 0)
    if(all(neiStruc <= 2))
        neighbors <- do.call(cbind, neighbors2)
    else{
        neighbors <- NULL
        for(i in 1:12){
            if(!is.null(neighbors2[[i]])){
                new <- getNeighborsB(neiStruc[i], neighbors2[[i]])
                neighbors <- cbind(neighbors, new)
            }
        }
    }
    neighbors
}


#' Get Neighbours of All Vertices of a Graph
#' 
#' Obtain neighbours of vertices of a 1D, 2D, or 3D graph.
#'
#' @param mask a vector, matrix, or 3D array specifying vertices within a graph. Vertices of value 1 are within the graph and 0 are not.
#' @param neiStruc a scalar, vector of four components, or \eqn{3\times4} matrix
#' corresponding to 1D, 2D, or 3D graphs. It gives the definition of
#' neighbours of a graph.
#' All components of \code{neiStruc} should be positive (\eqn{\ge 0})
#' even numbers.
#' For 1D graphs, \code{neiStruc} gives
#' the number of neighbours of each vertex.
#' For 2D graphs, \code{neiStruc}[1] specifies
#' the number of neighbours on vertical direction, \code{neiStruc}[2]
#' horizontal direction, \code{neiStruc}[3] north-west (NW) to south-east (SE)
#' diagonal direction, and \code{neiStruc}[4] south-west (SW) to
#' north-east (NE) diagonal direction. For 3D
#' graphs, the first row of \code{neiStruc} specifies the number of neighbours on
#' vertical direction, horizontal direction and two diagonal directions from
#' the 1-2 perspective, the second row the 1-3 perspective, and the
#' third row the 2-3 perspective. The index to
#' perspectives is represented with the leftmost subscript of the
#' array being the smallest.
#'
#' @return   A matrix with each row giving the neighbours of a vertex.
#' The number of the rows is equal to the number of 
#' vertices within the graph and the number or columns
#' is the number of neighbours of each vertex.
#' 
#' For a 1D graph, if each vertex has two neighbours,
#' The first column are the neighbours on the left-hand side of
#' corresponding vertices and the second column the right-hand side.
#' For the vertices on boundaries, missing neighbours are represented by
#' the number of vertices within a graph plus 1.
#' When \code{neiStruc} is bigger than 2, The first two columns
#' are the same as when \code{neiStruc} is equal to 2; the third column
#' are the neighbours on the left-hand side of the vertices on the first column;
#' the forth column are the neighbours on the right-hand side of the vertices
#' on the second column, and so on and so forth. And again for the
#' vertices on boundaries, their missing neighbours are represented by
#' the number of vertices within a graph plus 1.
#' 
#' For a 2D graph, the index to vertices is column-wised. For each
#' vertex, the order of neighbours are as follows. First are those
#' on the vertical direction, second the horizontal
#' direction, third the NW to SE diagonal
#' direction, and forth the SW to NE diagonal
#' direction. For each direction, the neighbours of every vertex
#' are arranged in the same way as in a 1D graph.
#' 
#' For a 3D graph, the index to vertices is that
#' the leftmost subscript of the array moves the fastest.  For each
#' vertex, the neighbours from the 1-2 perspective
#' appear first and then the 1-3 perspective and finally
#' the 2-3 perspective. For each perspective, the neighbours are arranged
#' in the same way as in a 2D graph.
#'
#' @export
#' @details   There could be more than one way to define the same 3D neighbourhood
#' structure for a graph (see Example 3 for illustration). 
#' @keywords spatial
#' @references 
#'  Winkler, G. (2003)
#'  "Image Analysis, Random Fields and Markov Chain Monte Carlo Methods: A Mathematical Introduction" (2nd ed.)
#'  \emph{Springer-Verlag}
#'  
#'  Feng, D. (2008)
#'  "Bayesian Hidden Markov Normal Mixture Models with Application to MRI Tissue Classification"
#'  \emph{Ph. D. Dissertation, The University of Iowa} 
#' @usage getNeighbors(mask, neiStruc)
#' @examples
#'   #Example 1: get all neighbours of a 1D graph.
#'   mask <- c(0,0,rep(1,4),0,1,1,0,0,1,1,1)
#'   getNeighbors(mask, neiStruc=2)
#'   
#'   #Example 2: get all neighbours of a 2D graph based on neighbourhood structure
#'   #           corresponding to the second-order Markov random field.
#'   mask <- matrix(1, nrow=2, ncol=3)
#'   getNeighbors(mask, neiStruc=c(2,2,2,2))
#'   
#'   #Example 3: get all neighbours of a 3D graph based on 6 neighbours structure
#'   #           where the neighbours of a vertex comprise its available
#'   #           N,S,E,W, upper and lower adjacencies. To achieve it, there
#'   #           are several ways, including the two below.
#'   mask <- array(1, dim=rep(3,3))
#'   n61 <- matrix(c(2,2,0,0,
#'                   0,2,0,0,
#'                   0,0,0,0), nrow=3, byrow=TRUE)
#'   n62 <- matrix(c(2,0,0,0,
#'                   0,2,0,0,
#'                   2,0,0,0), nrow=3, byrow=TRUE)
#'   n1 <- getNeighbors(mask, neiStruc=n61)
#'   n2 <- getNeighbors(mask, neiStruc=n62)
#'   n1 <- apply(n1, 1, sort)
#'   n2 <- apply(n2, 1, sort)
#'   all(n1==n2)
#'   
#'   #Example 4: get all neighbours of a 3D graph based on 18 neighbours structure
#'   #           where the neighbours of a vertex comprise its available
#'   #           adjacencies sharing the same edges or faces.
#'   #           To achieve it, there are several ways, including the one below.
#'   
#'   n18 <- matrix(c(2,2,2,2,
#'                   0,2,2,2,
#'                   0,0,2,2), nrow=3, byrow=TRUE)  
#'   mask <- array(1, dim=rep(3,3))
#'   getNeighbors(mask, neiStruc=n18)
getNeighbors <- function(mask, neiStruc){
    if(!(is.vector(mask) || is.matrix(mask) || length(dim(mask))==3))
        stop("The graph has to be in 1D, 2D, or 3D.")

    if(length(mask) < 3)
        stop("There are at least 3 vertices.")
    if(is.vector(mask))
        if(!is.vector(neiStruc) || length(neiStruc) != 1)
            stop("For a line, you need to tell the number of neighbors
                 on only one direction.")
    if(is.matrix(mask))
        if(!is.vector(neiStruc) || length(neiStruc)!=4)
            stop("For a plane, you need to tell the number of neighbors
                 on four directions.")
    if(length(dim(mask))==3)
        if(!is.matrix(neiStruc) || nrow(neiStruc)!=3 || ncol(neiStruc)!=4)
            stop("For a cube, you need to tell the number of neighbors
                 on four directions from all three perspective.")

    if(! all(neiStruc >= 0))
        stop("Number of neighbors have to be >= 0.")
    if(! all(neiStruc %% 2 == 0))
        stop("Number of neighbors have to be all even.")

    if(is.vector(mask))
        neighbors <- getNeighborsLine(mask, neiStruc)

    if(is.matrix(mask))
        neighbors <- getNeighborsPlane(mask, neiStruc)

    if(length(dim(mask))==3)
        neighbors <- getNeighborsCube(mask, neiStruc)

    neighbors
}
