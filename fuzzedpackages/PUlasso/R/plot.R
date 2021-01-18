#'@import ggplot2
#'@export
plot.cvPUfit<-function(x,...){
  fit =x
  dt=data.frame(logLambda = log(fit$lambda), cvm = fit$cvm, cvlo = fit$cvm-fit$cvsd,cvup = fit$cvm+fit$cvsd )
  logL.min = log(fit$lambda.min); logL.1se = log(fit$lambda.1se)
  ggplot(dt,aes_string('logLambda','cvm'))+
    theme_classic(base_size = 20)+
    geom_vline(aes_string(xintercept='logL.min'),linetype="dashed",color="darkgrey")+
    geom_vline(aes_string(xintercept='logL.1se'),linetype="dashed",color="darkgrey")+
    geom_errorbar(aes_string(ymin='cvlo',ymax='cvup'),color="grey")+geom_point(size=2)+
    ylab("Deviance")+xlab("log(Lambda)")
}
