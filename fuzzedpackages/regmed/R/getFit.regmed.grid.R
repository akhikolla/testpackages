
getFit.regmed.grid <- function(obj,lambda=NULL){

    ## pull out best model based on minimum BIC, if ties use
    ## most parsimonuous model (largest lambda) ###

    if(is.null(lambda)) lambda<-max(obj$grid.data$lambda[obj$grid.data$bic==min(obj$grid.data$bic,na.rm=T)])
    
    zed.n<-names(obj$call)[-1]
    zed.v<-as.character(obj$call)[-1]
    zed.v<-zed.v[zed.n %in% names(formals(regmed.fit))]
    zed.n<-zed.n[zed.n %in% names(formals(regmed.fit))]
    
 
    zed.v[zed.n == "frac.lasso"] <- as.character(obj$frac.lasso)
    
    new.call<-eval(parse(text=paste0("call(\"regmed.fit\",",paste(paste0(zed.n,"=quote(",zed.v,")"),collapse=","),",lambda=quote(",lambda,"))")))
    
    out<-c(obj$fit.list[which(obj$grid.data$lambda==lambda)][[1]],list(MedCov=obj$MedCov,call=new.call))
    
    out$frac.lasso <- obj$frac.lasso
    class(out)<-"regmed"
    return(out)
}
