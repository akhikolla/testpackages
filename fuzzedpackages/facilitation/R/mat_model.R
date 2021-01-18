#' matrix population model 
#' 
#' Produces the Matrix Population Model matrix for a continuous time structured
#' population model, to be applied in a linear ODE. Unlike \code{\link{mat.model.base}}, only
#' works with a single population. If only the number of stages is provided, returns a ramdom
#' population matrix.
#' @param n The number of life stages. Default is 3.
#' @param Ds An n-array with death rates for each stage.
#' @param Gs An (n-1)-array with growth rates for each stage but the last.
#' @param Rs Either a single reproduction rate for the oldest stage, or an n-array of reproduction rates for each stage.
#' @examples
#' mat <- mat.model.base(5)
#' mat2 <- mat.model.base(3,c(1,2,3),c(10,10),100)
#' @export
#' @import stats
mat.model.base  <- function(n=3,Ds=runif(n,0,5),Gs=runif(n-1,0,5),Rs=runif(n,0,5)){
    if(n==1){
        Rs-Ds
    }
    else{
        if(length(Rs)==1){Rs <- c(rep(0,n-1),Rs)}
        Gs[n] <- 0
        M <- diag(-Ds-Gs) + diag(Gs)[c(n,1:(n-1)),]
        M[1,] <- M[1,] + Rs
        M
    }
}

#' matrix population model 
#' 
#' Produces the Matrix Population Model matrix for a continuous time structured
#' population model, to be applied in a linear ODE. If there is more than one population,
#' returns a list of matrices, or one block-diagonal matrix created by the combination.
#' @param data Either the result of a simulation, to extract the parameters from, or a
#' data.frame containing the parameters.
#' @param ns an array of numbers of stages. Use when \code{data} is a data.frame and the is
#' more than one population.
#' @param combine.matrices Logical. Combine the matrices into a single, multi-population matrix?
#' @examples
#' # example 1
#' mat.model(create.parameters(n=4))
#' 
#' # example 2 
#' data(malthusian)
#' mat.model(malthusian)
#' 
#' # example 3
#' data(twospecies)
#' mat.model(twospecies,combine.matrices=TRUE)
#' @export
#' @importFrom Matrix bdiag
mat.model <- function(data, ns, combine.matrices=FALSE){
    if(class(data)=="data.frame"){
        rates<-data
        if(missing(ns)){
            n<-1
        }
        else {
            n <- length(ns)
        }
    }
    else {
        rates<-data$param
        n<-data$num.pop
        ns<-data$num.stages
    }

    if(n==1){
        ns<-nrow(rates)
        mat.model.base(ns,rates$D,rates$G,rates$R)
    }
    else{
        nstarts<-c(1,sapply(1:(n-1),function(i)sum(ns[1:i])+1))
        Ms <- lapply(1:n,function(i){
                            r<-rates[1:ns[i]+nstarts[i]-1,]
                            mat.model.base(ns[i],r$D,r$G,r$R)
                        }
                    )
        if(combine.matrices){
            Matrix::bdiag(Ms)
        }
        else {
            Ms
        }
    }
}

#' solution.matrix 
#' 
#' The \code{solution.matrix} function returns the solution to a linear ODE of the form P' = MP,
#' which is merely P(t) = exp(Mt)p0 where p0 is the initial condition
#' @param p0 initial condition, as an array
#' @param M a square matrix with as many rows as P0
#' @param times an array containing the times in which to calculate the solution
#' @export
#' @importFrom Matrix expm
#' @examples
#' mat <- mat.model.base(5)
#' solution.matrix(c(1,0,0,0,0),mat)
solution.matrix <- function(p0, M, times = c(1:10)){
    expm <- function(M) as.matrix(Matrix::expm(M))

    S <- matrix(nrow=nrow(M),ncol=length(times))
    for(i in 1:length(times)){
        S[,i] <- expm(M*times[i]) %*% p0
    }
    colnames(S) <- times
    t(S)
}

#' Limiting Rate
#' 
#' This function returns the real dominant eigenvalue of a Matrix Population Model matrix. 
#' That is a real number that corresponds to the per-capita growth rate that a population
#' approaches as time passes, in a model with no interactions.
#'
#' A structured population can grow at exactly this rate if the
#' distribution between stages corresponds exactly to the distribution of the dominant
#' eigenvector. The models that can be simulated by this package are of a class that always
#' has a real dominant eigenvector. Note that these are continuous-time models, in which r >
#' 0 means the population will grow, and r < 0 means it will decrease. This function doesn't
#' throw errors, instead it returns 'NA'.
#' 
#' @param mat a square matrix
#' @export
#' @examples
#' mat <- mat.model.base(5)
#' limiting.rate(mat)
limiting.rate <- function(mat){tryCatch(max(Re(eigen(mat,symmetric=F)$values)),error=function(e) NA)}
