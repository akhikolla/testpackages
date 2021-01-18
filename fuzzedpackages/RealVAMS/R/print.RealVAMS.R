print.RealVAMS <-
function (x, ...) 
{
    cat("Object of class 'RealVAMS'. Use generic functions 'coef','plot,'summary' on this object.", "\n") 
    cat("Object contains elements:\n")
    print(names(x))
    cat("See the documentation for the RealVAMS function for a description of these objects\n")
    #cat("Number of iterations: ",x$iter,"\n",sep="")
    #cat("Log-likelihood: ",x$loglik,"\n",sep="")
    cat("Coefficients:\n")
    print(x$parameters[,1])
  #  cat("Teacher Effects\n")
   # print(x$teach.effects)
    cat("\n")
    
    
    #rchol <- try(chol(R_inv))
    #yhat.s <- try(as.vector(rchol %*% (yhat)))
    #sresid <- try(as.vector(rchol %*% Y - yhat.s))
    
}
