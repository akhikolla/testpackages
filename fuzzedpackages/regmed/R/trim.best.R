
trim.best <-
function(obj, lambda=0.0, mediator.epsilon=0.0001){
  
    if(class(obj) != "regmed.grid") {
      stop("input not regmed.grid class")
    }

    
    fit.best <- getFit.regmed.grid(obj)
      
    keep.med<-(abs(fit.best$alpha*fit.best$beta) >= mediator.epsilon)

    ## if best fitting model has no mediators, return fit
    if(all(!keep.med)){
        return(fit.best)
    }

    
    ## otherwise reduce to selected mediators and refit 
    keep.med<-paste(paste0("\"",rownames(keep.med)[keep.med],"\""),collapse=",")
    flasso.txt <- paste0(", frac.lasso = ", as.character(fit.best$frac.lasso))
    med.txt <- paste0(",mediator=",as.character(fit.best$call)[names(fit.best$call)=="mediator"],"[,c(",keep.med,"),drop=FALSE]")

    out<-eval(parse(text=paste0("update(fit.best,lambda=",lambda, med.txt, flasso.txt, ")")))
      
    return(out)
}
