
EMVSplot <- function(result, 
                     plot_type = c("both", "reg", "gf"), 
                     omit.zeroes = FALSE){

	betas<-result$betas
	intersects<-result$intersects
	v0s<-result$v0
	
	if(result$log_v0 == F){
	
	if(result$independent == F){
	  log_post<-result$log_g_function
	  
	  if (plot_type=="both"){
	    par(mfrow=c(1,2))
	  }
	  
	  if (plot_type=="both"| plot_type=="reg"){
	    
	    select<-apply(result$prob_inclusion,2,function(x){as.numeric(x>0.5)})
	    
	    if (omit.zeroes){betas[select==0]=0}
	    
	    matplot(v0s,betas,xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,col="grey",lty=2,type="l")
	    matpoints(v0s,betas*select,xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,col=4,lty=2,pch=19)
	    matpoints(v0s,betas*(1-select),xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,col=2,lty=2,pch=19)
	    title("EMVS Regularization Plot")
	    
	    par(xpd=T)
	    labels=paste("X",1:ncol(betas),sep="")
	    labels[select[length(v0s),]==0]<-""
	    text(max(v0s)*(1.1),betas[length(v0s),],labels=labels,col=4)
	  }
	  
	  if (plot_type=="both"| plot_type=="gf"){
	    
	    plot(as.numeric(log_post)~v0s,pch=4,type="b",col=2,lwd=1,xlab=expression(v[0]),ylab=
	           expression(log(g(gamma))))
	    
	    title(expression(Log(g(gamma))))
	  }
	}

	if(result$independent == T) {
	  par(mfrow = c(1, 1))
	  select<-apply(result$prob_inclusion,2,function(x){as.numeric(x>0.5)})
	  
	  if (omit.zeroes){betas[select==0]=0}
	  
	  matplot(v0s,betas,xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,col="grey",lty=2,type="l")
	  matpoints(v0s,betas*select,xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,col=4,lty=2,pch=19)
	  matpoints(v0s,betas*(1-select),xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,col=2,lty=2,pch=19)
	  title("EMVS Regularization Plot")
	  
	  par(xpd=T)
	  labels=paste("X",1:ncol(betas),sep="")
	  labels[select[length(v0s),]==0]<-""
	  text(max(v0s)*(1.1),betas[length(v0s),],labels=labels,col=4)
	  
	}
	} else {
	  if(result$independent == F){
	    log_post<-result$log_g_function
	    
	    if (plot_type=="both"){
	      par(mfrow=c(1,2))
	    }
	    
	    if (plot_type=="both"| plot_type=="reg"){
	      
	      select<-apply(result$prob_inclusion,2,function(x){as.numeric(x>0.5)})
	      
	      if (omit.zeroes){betas[select==0]=0}
	      
	      matplot(log(v0s),betas,xlab=expression(log(v[0])),ylab=expression(hat(beta)),lwd=1,col="grey",lty=2,type="l")
	      matpoints(log(v0s),betas*select,xlab=expression(log(v[0])),ylab=expression(hat(beta)),lwd=1,col=4,lty=2,pch=19)
	      matpoints(log(v0s),betas*(1-select),xlab=expression(log(v[0])),ylab=expression(hat(beta)),lwd=1,col=2,lty=2,pch=19)
	      title("EMVS Regularization Plot")
	      
	      par(xpd=T)
	      labels=paste("X",1:ncol(betas),sep="")
	      labels[select[length(v0s),]==0]<-""
	      text(max(v0s)*(1.1),betas[length(v0s),],labels=labels,col=4)
	    }
	    
	    if (plot_type=="both"| plot_type=="gf"){
	      
	      plot(as.numeric(log_post)~v0s,pch=4,type="b",col=2,lwd=1,xlab=expression(v[0]),ylab=
	             expression(log(g(gamma))))
	      
	      title(expression(Log(g(gamma))))
	    }
	  }
	  
	  if(result$independent == T) {
	    par(mfrow = c(1, 1))
	    select<-apply(result$prob_inclusion,2,function(x){as.numeric(x>0.5)})
	    
	    if (omit.zeroes){betas[select==0]=0}
	    
	    matplot(log(v0s),betas,xlab=expression(log(v[0])),ylab=expression(hat(beta)),lwd=1,col="grey",lty=2,type="l")
	    matpoints(log(v0s),betas*select,xlab=expression(log(v[0])),ylab=expression(hat(beta)),lwd=1,col=4,lty=2,pch=19)
	    matpoints(log(v0s),betas*(1-select),xlab=expression(log(v[0])),ylab=expression(hat(beta)),lwd=1,col=2,lty=2,pch=19)
	    title("EMVS Regularization Plot")
	    
	    par(xpd=T)
	    labels=paste("X",1:ncol(betas),sep="")
	    labels[select[length(v0s),]==0]<-""
	    text(max(log(v0s))+0.75,betas[length(v0s),],labels=labels,col=4)
	    
	  }
	}
	
}


