initparam <- function(obs, model, matT){
  lambda <- beta <- sigma <- list()
  prop <- runif(model@G)
  prop <-  prop / sum(prop)
  for (g in 1:model@G){
    if (model@spatial){
      lambda[[g]] <- matrix(0, model@K, 4)
      colnames(lambda[[g]]) <- c("cste", "row", "col", "time")      
    }else{
      lambda[[g]] <- matrix(0, model@K, 2)
      colnames(lambda[[g]]) <- c("cste", "time")
    }
    rownames(lambda[[g]]) <- paste("polynom", 1:model@K,sep=".")
    beta[[g]] <- matrix(0, model@K, model@Q+1)
    for (k in 1:model@K) beta[[g]][k,] <- coefficients(lm(as.numeric(obs@x[sample(1:obs@n,1),sample(1:obs@JJ,1),]) ~ matT+0))
  }
  sigma <- matrix(max(var(obs@x)), model@G, model@K)
  rownames(sigma) <- names(lambda) <- names(beta) <- names(prop) <- paste("compo", 1:model@G, sep=".")
  colnames(sigma) <- paste("polynom", 1:model@K, sep=".")
  return(new("STCparam", proportions=prop, lambda=lambda, beta=beta, sigma=sigma))
}

# Probacond <- function(obs, model, param, newtool, toollogistic){
#   output <- list(condintra=list(), condintramargin=list(), logcondinter=matrix(0, obs@n, model@G, byrow = TRUE))
#   for (g in 1:model@G){
#     output$condintra[[g]] <-  list()
#     output$condintramargin[[g]] <- matrix(0, obs@JJ*obs@TT, obs@n)
#     poidspolynom <- matrix(0, obs@TT*obs@JJ, model@K)
#     for (k in 1:model@K) poidspolynom[,k] <- toollogistic %*% param@lambda[[g]][k,]
#     poidspolynom<- exp(sweep(poidspolynom,1,apply(poidspolynom,1,max), "-"))
#     poidspolynom <- poidspolynom / rowSums(poidspolynom)
#     
#     for (k in 1:model@K){
#       output$condintra[[g]][[k]] <- matrix(NA, nrow(poidspolynom), ncol(obs@x))
#       for (loc in 1:nrow(poidspolynom)){
#         output$condintra[[g]][[k]][loc,] <- poidspolynom[loc,k] * dnorm(obs@x[loc,], mean = sum(newtool[loc,] * param@beta[[g]][k,]) , sd = sqrt(param@sigma[g,k]))
#       } 
#       output$condintramargin[[g]] <- output$condintramargin[[g]] + output$condintra[[g]][[k]]
#     }  
#     output$logcondinter[,g] <- colSums(log(output$condintramargin[[g]])) + log(param@proportions[g]) 
#   }
#   return(output)
# }


TuneOutput <- function(output){
  output@param@proportions <- as.numeric(output@param@proportions)
  for (g in 1:output@model@G){
    output@param@lambda[[g]] <- matrix(output@param@lambda[[g]], output@model@K, 4)
    colnames(output@param@lambda[[g]]) <- c("intercept", "coordinate1", "coordinate2", "time")
    output@param@beta[[g]] <- matrix(output@param@beta[[g]], output@model@K, output@model@Q+1)
    colnames(output@param@beta[[g]]) <- paste("deg.", rep(0:output@model@Q), sep="")
    rownames(output@param@lambda[[g]]) <- rownames(output@param@beta[[g]]) <- paste("polynom", 1:output@model@K, sep=".")
  }
  names(output@param@beta) <- names(output@param@lambda) <- paste("component", 1:output@model@G, sep=".")
  for (g in 1:output@model@G){
    for (k in 1:output@model@K) names(output@param@beta[[g]][k,]) <- paste("coeff", 0:(length(output@param@beta[[g]][k,])-1), sep=".")
  }  
#   output@partitions@hardseg <- list()
#   for (g in 1:output@model@G){
#     output@partitions@hardseg[[g]] <- matrix(apply(exp(toollogistic %*% t(output@param@lambda[[g]])), 1, which.max), output@data@JJ, output@data@TT)
#     colnames(output@partitions@hardseg[[g]]) <- paste("time", 1:output@data@TT, sep=".")
#     rownames(output@partitions@hardseg[[g]]) <- paste("site", 1:output@data@JJ, sep=".")
#   }
#   names(output@partitions@hardseg) <- paste("comp", 1:output@model@G, sep=".")

   output@partitions@hardind <- apply(output@partitions@fuzzyind, 1, which.max)
   names(output@partitions@hardind) <- rownames(output@partitions@fuzzyind) <- colnames(output@data@x)
   colnames(output@partitions@fuzzyind) <- paste("comp", 1:output@model@G, sep=".")
   
  output@criteria@AIC <- output@criteria@loglike - output@model@nbparam
  output@criteria@BIC <- output@criteria@loglike - 0.5 * output@model@nbparam * log(output@data@n)
  rownames(output@param@sigma) <- paste("component", 1:output@model@G, sep=".")
  colnames(output@param@sigma) <- paste("polynom", 1:output@model@K, sep=".")
  entrop <- 0
  for (g in 1:output@model@G) entrop <- entrop + sum(log(output@partitions@fuzzyind[which(output@partitions@hardind==g),g]))
  output@criteria@ICL <- output@criteria@BIC + entrop
  return(output)
}