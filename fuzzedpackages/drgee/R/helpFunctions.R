## helpFunctions.R

## Just a wrapper for centering a design matrix
## around the cluster mean defined by a cluster
## identification variable id
centerX <- function(x, id)
{
    return( .Call(center, x, id) )
}

## Just a wrapper of nleqslv in the R package nleqslv
## Used to solve a system of non-linear equations
## To use another equation solver,
## one can write another wrapper with the same input
## arguments and send this function as an input argument
findRoots <-
    function(beta.init, eq.func, d.eq.func = NULL, arg.list, ...){

        optim.object <- do.call(nleqslv, c(list(x = beta.init,
                                                fn = eq.func,
                                                jac = d.eq.func,
                                                arg.list = arg.list),
                                           list(...)))

        beta.hat <- optim.object$x
        
        if (optim.object$termcd > 2) {
            warning(paste("\nnleqslv: ", optim.object$message))
        }

        return( list(roots = beta.hat, optim.object = optim.object))

    }




