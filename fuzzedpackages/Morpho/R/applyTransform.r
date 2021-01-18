#' apply affine transformation to data
#'
#' apply affine transformation to data
#' @param x matrix or mesh3d
#' @param trafo 4x4 transformation matrix or an object of class "tpsCoeff"
#' @param inverse logical: if TRUE, the inverse of the transformation is applied (for TPS coefficients have to be recomputed)
#' @param threads threads to be used for parallel execution in tps deformation.
#' @param ... additional arguments, currently not used.
#' @return the transformed object
#' @examples
#' data(boneData)
#' rot <- rotonto(boneLM[,,1],boneLM[,,2])
#' trafo <- getTrafo4x4(rot)
#' boneLM2trafo <- applyTransform(boneLM[,,2],trafo)
#' @rdname applyTransform
#' @seealso \code{\link{rotonto}, link{rotmesh.onto}, \link{tps3d}, \link{computeTransform}}
#' @export
applyTransform <- function(x,trafo,...)UseMethod("applyTransform")

#' @rdname applyTransform
#' @export
applyTransform.matrix <- function(x,trafo,inverse=FALSE,threads=1,...) {
    if (is.matrix(trafo)) {
        if (ncol(trafo) == 3 && ncol(x) ==3)
            trafo <- mat3x3tomat4x4(trafo)
        if (inverse)
            trafo <- solve(trafo)
        out <-homg2mat(trafo%*%mat2homg(x))
    } else if (inherits(trafo,"tpsCoeff")) {
        if (ncol(trafo$refmat) != ncol(x))
            stop("TPS must be computed from control points of the same dimensionality")
        if (inverse)
            trafo <- computeTransform(trafo$refmat,trafo$tarmat,type="tps")
        out <- .fx(trafo$refmat,x,trafo$coeff,threads=threads)
    }
    return(out)
}

#' @rdname applyTransform
#' @export
applyTransform.mesh3d <- function(x,trafo,inverse=FALSE,threads=1,...) {
    x$vb[1:3,] <- t(applyTransform(t(x$vb[1:3,]) ,trafo,inverse = inverse,threads=threads))
    ## case affine transformation
    reflect <- FALSE
    if (is.matrix(trafo)) {
        if (det(trafo) < 0) 
            reflect <- TRUE
    } else { ##case transform is tps
        if (det(computeTransform(trafo$refmat,trafo$tarmat,reflection = T)) < 0)
            reflect <- TRUE
    }
        if (reflect) {
            x <- invertFaces(x)
            message("faces' orientation has been inverted")
        }
    
    ## update normals
    if (!is.null(x$normals)) {
        if (!is.null(x$it) || !is.matrix(trafo)) {
            x <- vcgUpdateNormals(x,silent=TRUE)
        } else {
            ntrafo <- trafo
            ntrafo[1:3,4] <- 0
            x$normals[1:3,] <- t(applyTransform(t(x$normals[1:3,]),ntrafo))
        }
    }
    return(x)
 }

#' @rdname applyTransform
#' @export
applyTransform.default <- function(x,trafo,inverse=FALSE,threads=1,...) {
    x <- t(x)
    if (is.matrix(trafo)) {
        if (ncol(trafo) == 3 && ncol(x) == 3)
            trafo <- mat3x3tomat4x4(trafo)
        if (inverse)
            trafo <- solve(trafo)
        out <- homg2mat(trafo%*%mat2homg(x))
    } else if (inherits(trafo,"tpsCoeff")) {
        if (ncol(trafo$refmat) != ncol(x))
            stop("TPS must be computed from control points of the same dimensionality")
        if (inverse)
            trafo <- computeTransform(trafo$refmat,trafo$tarmat,type="tps")
        out <- .fx(trafo$refmat,x,trafo$coeff,threads=threads)
    }
    return(out)
}

mat3x3tomat4x4 <- function(x) {
    n <- ncol(x)
    x <- rbind(cbind(x,0),0);x[n+1,n+1] <-1
    return(x)
}
