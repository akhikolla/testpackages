print.bess.one=function(x, ...)
{
  if(x$family=="bess_gaussian") print(c(Df=sum(x$beta!=0),MSE=x$mse,AIC=x$AIC,BIC=x$BIC,EBIC=x$EBIC))else
    print(c(Df=sum(x$beta!=0),Dev=x$deviance,AIC=x$AIC,BIC=x$BIC,EBIC=x$EBIC))
}

