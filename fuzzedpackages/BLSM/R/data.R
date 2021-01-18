#' Example Adjacency Matrix
#' 
#' Adjacency matrix of a 10 nodes random network for testing purposes
#' 
#' @name example_adjacency_matrix
#' @docType data
#' @format A binary adjacency matrix representing links between nodes.
"example_adjacency_matrix"

#' Example Weights Matrix
#' 
#' "BLSM weights" matrix of a 10 nodes random network for testing purposes
#' 
#' @name example_weights_matrix
#' @docType data
#' @format A matrix containing positive weights for all pairs of nodes. 
#' 
#' Given a couple of nodes, a weight expresses the importance of the distance between the 
#' coordinates associated to the two nodes in the latent space in terms of the overall likelihood of the graph. 
#' For this reason, even missing links must have a coefficient, otherwise the relative positioning of disconnected nodes
#' would have no effect at all on the graph likelihood.
#'  
#' The exact probability equation is described in \link[BLSM]{BLSM}, as well as the notation used.
#' 
#' A few examples:
#' \itemize{
#' \item for unweighted networks, the "BLSM weights" matrix has all the values set to 1. 
#' \item if two nodes share a strong connection, then
#' the weight coefficient should be greater than 1 so that their positions in the latent space will be closer than they would be in an unweighted framework. 
#' \item if two nodes share a weak connection, a coefficient smaller than 1 will allow the latent coordinates to be pretty far from each other even though the nodes are connected. 
#' }
"example_weights_matrix"

#' Example BLSM object
#' 
#' BLSM object obtained by applying the Procrustean version of the latent space model to the unweighted network 
#' whose adjacency matrix is \link[BLSM]{example_adjacency_matrix}. Further details concerning the 
#' simulation are contained in the BLSM object itself.
#'  
#' @name example_blsm_obj
#' @docType data
#' @format BLSM object (\code{blsm_obj}), i.e. a list containing:
#' \itemize{
#' \item \code{Alpha }{\eqn{\alpha} values from the sampled iterations}
#' \item \code{Likelihood }{Log-likelihood values from the sampled iterations}
#' \item \code{Iterations }{Latent space coordinates from the sampled iterations. Latent positions are stored in a
#' 3D array whose dimensions are given by (1: number of nodes, 2: space dimensionality, 3: number of iterations).
#' In the non-Procrustean framework the latent distances are given instead of the positions: another 3D array is returned, whose dimensions
#' are given by (1: number of nodes, 2: number of nodes, 3: number of iterations). The command needed in order to get the average values over the iterations for
#' either the positions or the distances is \code{rowMeans(blsm_obj$Iterations, dims=2)} (see example below).}
#' \item \code{StartingPositions }{Latent space coordinates right after the initialization step. In the non-Procrustean framework starting distances are given instead.}
#' \item \code{Matrix }{Original matrices of the network (adjacency and BLSM weights)}
#' \item \code{Parameters }{List of parameters specified during the call to \link[BLSM]{estimate_latent_positions}}
#' }
"example_blsm_obj"

