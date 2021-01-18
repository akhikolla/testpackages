summary.splmm <- function(object,...)
 {
    if (object$stopped)
      {
      cat("Algorithm stopped due to","\n")
      cat("sparsity assumption","\n")
      cat("no summary is available.","\n")
      cat("try to increase lambda!","\n")
      } else {
        
    # 1. part
    table1 <- round(c(object$aic,object$bic,object$bbic,object$ebic,object$logLik,object$deviance,object$objective),1)
    names(table1) <- c("AIC","BIC","BBIC","EBIC","logLik","deviance","objective")
    cat("Model fitted by",object$method,"for","lambda1 =",round(object$lambda1,3),",lambda2 =",round(object$lambda2,3),":\n")
    print(table1)
    
    # 2. part
    cat("","\n")
    cat("Fixed effects:","\n")
    penalty.b <- rep("",length(object$coefficients))
    penalty.b[object$nonpen.b] <- "(n)"
    table3 <- data.frame(object$coefficients,penalty.b)
    
    xnames <- colnames(object$data$x)
    
    if (is.null(xnames))
    {
      rownames(table3) <- c("(Intercept)",paste("X", 2:length(object$coefficients),sep=""))
    } else
    {
      rownames(table3) <- c("(Intercept)",colnames(object$data$x)[-1]) 	
    }
    table3 <- table3[object$coefficients!=0,]
    colnames(table3) <- c("Estimate"," ")
    cat("Number of nonzero fixed effects=",sum(object$coefficients!=0),"\n")
    print(table3)

    # 3. part
    cat("","\n")
    cat("Random effects:","\n")
    penalty.L <- rep("",nrow(object$D))
    penalty.L[object$nonpen.L] <- "(n)"
    penalty.L = c(penalty.L,"")
    
    tableVar <- c(diag(object$D),object$sigma^2)
    tableSd <- c(sqrt(diag(object$D)),object$sigma)
    matrixVar <- data.frame(tableVar,tableSd, penalty.L)
    colnames(matrixVar) <- c("Variance","Std.Dev.","")
    

    znames <- colnames(object$data$z)

    if (is.null(znames))
      {
        rownames(matrixVar) <-c("(Intercept)",paste("Z", 2:nrow(object$D),sep=""),"")
      } else
      {
        rownames(matrixVar) <- c("(Intercept)",colnames(object$data$z)[-1],"sigma") 	
      }
    
    matrixVar <- matrixVar[diag(object$D)!=0,]
    cat("Number of nonzero random effects=",sum(diag(object$D)!=0),"\n")
    print(matrixVar)
    
    

    # 4th part
    cat("","\n")
    cat("Number of iterations:",object$counter,"\n")
  }
 }

print.splmm <- function(x,...) {
    # 0. part
    print(x$call)
  
  # 2. part
  cat("","\n")
  cat("Fixed effects:","\n")
  cat("Number of nonzero fixed effects=",sum(x$coefficients!=0),"\n")
  


  
    # 3. part
    cat("","\n")
    cat("Random effects:","\n")
    cat("Number of nonzero random effects=",sum(diag(x$D)!=0),"\n")
    
    

}

utils::globalVariables(c("lambda1", "Value", "Criteria","lambda2"))

plot.splmm <- function(x,...){
  
  # check if the object is splmm.tuning
  if(class(x)!="splmm.tuning") stop("input is not splmm.tuning object.")
  
  if(x$lam1.tuning&!x$lam2.tuning){
    BIC = cbind.data.frame(x$lam1.seq,rep("BIC",length(x$BIC.lam1)),x$BIC.lam1)
    names(BIC) = c("lambda1","Criteria","Value")
    AIC = cbind.data.frame(x$lam1.seq,rep("AIC",length(x$AIC.lam1)),x$AIC.lam1)
    names(AIC) = c("lambda1","Criteria","Value")
    BBIC = cbind.data.frame(x$lam1.seq,rep("BBIC",length(x$BBIC.lam1)),x$BBIC.lam1)
    names(BBIC) = c("lambda1","Criteria","Value")
    EBIC = cbind.data.frame(x$lam1.seq,rep("EIC",length(x$EBIC.lam1)),x$EBIC.lam1)
    names(EBIC) = c("lambda1","Criteria","Value")
    plot.df = rbind.data.frame(BIC,AIC,BBIC,EBIC)
    plot.df$Value = round(plot.df$Value,2)
    
    ggplot(plot.df, aes(x = lambda1, y = Value, group = Criteria)) +
      geom_line(aes(linetype = Criteria, color = Criteria))+
      geom_point(aes(color = Criteria))+
      ggtitle("Tuning for lambda 1")
    
  }else if(!x$lam1.tuning&x$lam2.tuning){
    
    BIC = cbind.data.frame(x$lam2.seq,rep("BIC",length(x$BIC.lam2)),x$BIC.lam2)
    names(BIC) = c("lambda2","Criteria","Value")
    AIC = cbind.data.frame(x$lam2.seq,rep("AIC",length(x$AIC.lam2)),x$AIC.lam2)
    names(AIC) = c("lambda2","Criteria","Value")
    BBIC = cbind.data.frame(x$lam2.seq,rep("BBIC",length(x$BBIC.lam2)),x$BBIC.lam2)
    names(BBIC) = c("lambda1","Criteria","Value")
    EBIC = cbind.data.frame(x$lam2.seq,rep("EIC",length(x$EBIC.lam2)),x$EBIC.lam2)
    names(EBIC) = c("lambda1","Criteria","Value")
    plot.df = rbind.data.frame(BIC,AIC,BBIC,EBIC)
    plot.df$Value = round(plot.df$Value,2)
    
    
    ggplot(plot.df, aes(x = lambda2, y = Value, group = Criteria)) +
      geom_line(aes(linetype = Criteria, color = Criteria))+
      geom_point(aes(color = Criteria))+
      ggtitle("Tuning for lambda 2")
    
  }else if(x$lam1.tuning&x$lam2.tuning){
    BIC1 = cbind.data.frame(x$lam1.seq,rep("BIC",length(x$BIC.lam1)),x$BIC.lam1)
    names(BIC1) = c("lambda1","Criteria","Value")
    AIC1 = cbind.data.frame(x$lam1.seq,rep("AIC",length(x$AIC.lam1)),x$AIC.lam1)
    names(AIC1) = c("lambda1","Criteria","Value")
    BBIC1 = cbind.data.frame(x$lam1.seq,rep("BBIC",length(x$BBIC.lam1)),x$BBIC.lam1)
    names(BBIC1) = c("lambda1","Criteria","Value")
    EBIC1 = cbind.data.frame(x$lam1.seq,rep("EBIC",length(x$EBIC.lam1)),x$EBIC.lam1)
    names(EBIC1) = c("lambda1","Criteria","Value")
    plot.df1 = rbind.data.frame(BIC1,AIC1,BBIC1,EBIC1)
    plot.df1$Value = round(plot.df1$Value,2)
    
    BIC2 = cbind.data.frame(x$lam2.seq,rep("BIC",length(x$BIC.lam2)),x$BIC.lam2)
    names(BIC2) = c("lambda2","Criteria","Value")
    AIC2 = cbind.data.frame(x$lam2.seq,rep("AIC",length(x$AIC.lam2)),x$AIC.lam2)
    names(AIC2) = c("lambda2","Criteria","Value")
    BBIC2 = cbind.data.frame(x$lam2.seq,rep("BBIC",length(x$BBIC.lam2)),x$BBIC.lam2)
    names(BBIC2) = c("lambda2","Criteria","Value")
    EBIC2 = cbind.data.frame(x$lam2.seq,rep("EBIC",length(x$EBIC.lam2)),x$EBIC.lam2)
    names(EBIC2) = c("lambda2","Criteria","Value")
    plot.df2 = rbind.data.frame(BIC2,AIC2,BBIC2,EBIC2)
    plot.df2$Value = round(plot.df2$Value,2)
    
    
    p1 <- ggplot(plot.df1, aes(x = lambda1, y = Value, group = Criteria)) +
      geom_line(aes(linetype = Criteria, color = Criteria))+
      geom_point(aes(color = Criteria))+
      ggtitle("Tuning for lambda 1")
    
    p2 <- ggplot(plot.df2, aes(x = lambda2, y = Value, group = Criteria)) +
      geom_line(aes(linetype = Criteria, color = Criteria))+
      geom_point(aes(color = Criteria))+
      ggtitle("Tuning for lambda 2")
    
    grid.arrange(p1,p2,ncol=2)
  }
    
    
  
  
}

plot3D.splmm <- function(x,criteria=c("BIC","AIC","BBIC","EBIC"),type=c("line","surface"),...){
  # check if the object is splmm.tuning
  if(class(x)!="splmm.tuning") stop("input is not splmm.tuning object.")
  if(!x$lam1.tuning|!x$lam2.tuning) stop("input needs to be splmm.tuning object for both lam1 and lam2.")
  
  criteria <- match.arg(criteria)
  type <- match.arg(type)
  
  if(type=="line"){
    if(criteria=="BIC"){
      BIC.df = cbind.data.frame(expand.grid(x$lam1.seq,x$lam2.seq),as.vector(x$fit.BIC))
      names(BIC.df) = c("lambda 1","lambda 2","Value")
      
      scatter3D(x=BIC.df$`lambda 1`, y=BIC.df$`lambda 2`, z=BIC.df$Value, phi = 0, bty = "g",  type = "h", 
                ticktype = "detailed", pch = 19, cex = 0.5, xlab="lambda 1",ylab="lambda 2",zlab="",main="BIC value")
    }else if(criteria=="AIC"){
      AIC.df = cbind.data.frame(expand.grid(x$lam1.seq,x$lam2.seq),as.vector(x$fit.AIC))
      names(AIC.df) = c("lambda 1","lambda 2","Value")
      
      scatter3D(x=AIC.df$`lambda 1`, y=AIC.df$`lambda 2`, z=AIC.df$Value, phi = 0, bty = "g",  type = "h", 
                ticktype = "detailed", pch = 19, cex = 0.5, xlab="lambda 1",ylab="lambda 2",zlab="",main="AIC value")
      
    }else if(criteria=="BBIC"){
      BBIC.df = cbind.data.frame(expand.grid(x$lam1.seq,x$lam2.seq),as.vector(x$fit.BBIC))
      names(BBIC.df) = c("lambda 1","lambda 2","Value")
      
      scatter3D(x=BBIC.df$`lambda 1`, y=BBIC.df$`lambda 2`, z=BBIC.df$Value, phi = 0, bty = "g",  type = "h", 
                ticktype = "detailed", pch = 19, cex = 0.5, xlab="lambda 1",ylab="lambda 2",zlab="",main="BBIC value")
    }else if(criteria=="EBIC"){
      EBIC.df = cbind.data.frame(expand.grid(x$lam1.seq,x$lam2.seq),as.vector(x$fit.EBIC))
      names(EBIC.df) = c("lambda 1","lambda 2","Value")
      
      scatter3D(x=EBIC.df$`lambda 1`, y=EBIC.df$`lambda 2`, z=EBIC.df$Value, phi = 0, bty = "g",  type = "h", 
                ticktype = "detailed", pch = 19, cex = 0.5, xlab="lambda 1",ylab="lambda 2",zlab="",main="EBIC value")
    }
    
  }else if(type=="surface"){
    if(criteria=="BIC"){
      persp3D(x = x$lam1.seq, y = x$lam2.seq, 
              z = x$fit.BIC, phi = 0, bty = "g", type = "h", 
              ticktype = "detailed", pch = 19, cex = 0.5, 
              xlab = "lambda 1", ylab = "lambda 2", 
              zlab = "", main = "BIC value")
    }else if(criteria=="AIC"){
      persp3D(x = x$lam1.seq, y = x$lam2.seq, 
              z = x$fit.AIC, phi = 0, bty = "g", type = "h", 
              ticktype = "detailed", pch = 19, cex = 0.5, 
              xlab = "lambda 1", ylab = "lambda 2", 
              zlab = "", main = "AIC value")
      
    }else if(criteria=="BBIC"){
      persp3D(x = x$lam1.seq, y = x$lam2.seq, 
              z = x$fit.BBIC, phi = 0, bty = "g", type = "h", 
              ticktype = "detailed", pch = 19, cex = 0.5, 
              xlab = "lambda 1", ylab = "lambda 2", 
              zlab = "", main = "BBIC value")
    }else if(criteria=="EBIC"){
      persp3D(x = x$lam1.seq, y = x$lam2.seq, 
              z = x$fit.EBIC, phi = 0, bty = "g", type = "h", 
              ticktype = "detailed", pch = 19, cex = 0.5, 
              xlab = "lambda 1", ylab = "lambda 2", 
              zlab = "", main = "EBIC value")
    }
    
  }
  
  
}
#plot.splmm <- function(x,...) {

#par(mfrow=c(3,2))

## Tukey-Anscombe plot

#plot(x$residuals~x$fitted.values,col=x$data$grp,main="Tukey-Anscombe Plot",xlab="fitted values",ylab="raw residuals")
#abline(h=0,col="grey")

## QQ-Plot of the residuals

#qqnorm(x$residuals,col=x$data$grp,main="QQ-Plot of the residuals")
#qqline(x$residuals)

## QQ-plot of the random effects
#ranef <- unlist(x$random,recursive=FALSE)
#ranef.sd <- rep(sqrt(diag(x$Psi)),length(levels(x$data$grp)))
#ranef.st <- ranef/ranef.sd
#qqnorm(ranef.st,main="QQ-Plot of the standardized random effects",col=rep(1:dim(x$data$z)[[2]],length(levels(x$data$grp))))
#qqline(ranef.st)

## boxplots of the random effects
#boxplot(ranef~as.factor(rep(1:dim(x$data$z)[[2]],length(levels(x$data$grp)))),border=1:dim(x$data$z)[[2]],
#        main="Random effects by effects",xlab="random effects index")
#abline(h=0,col="grey",lty=2)

## fitted vs. observed y
#plot(x$fitted,x$data$y,xlab="fitted values",ylab="observed values",main="response ~ fitted")
#abline(a=0,b=1)

## histogram of the fixed effects
#hist(x$coefficients,main=paste("|active set| = ",sum(x$coefficients!=0)),xlab="fixed effects")
#rug(x$coefficients,lwd=3)

#}
