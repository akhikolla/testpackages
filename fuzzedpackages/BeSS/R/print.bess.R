print.bess=function(x, ...)
{
  if(x$family=="bess_gaussian") print(cbind(Df=x$s.list,MSE=x$mse,AIC=x$AIC,BIC=x$BIC,EBIC=x$EBIC))else
    print(cbind(Df=x$s.list,Dev=x$deviance,AIC=x$AIC,BIC=x$BIC,EBIC=x$EBIC))
}

