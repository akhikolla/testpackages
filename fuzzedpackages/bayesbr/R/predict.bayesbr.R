#'@title Prediction Method for \code{bayesbr} Objects
#'@name predict.bayesbr
#'@aliases predict.bayesbr
#'@description A function that informs various types of prediction through a beta regression by the Bayesian view.
#'@usage
#'\method{predict}{bayesbr}(object, newdata = NULL, type = c("response", "link",
#'"precision", "variance", "quantile"), na.action=c("exclude",
#'"replace"),at = 0.5,...)
#'@param object an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@param newdata A data frame in which to look for variables with which to predict. If omitted, the original observations are used.
#'@param type A character containing the types of predictions: if type is \emph{"response"} the function will calculate fitted values for theta; if type is \emph{"link"} the function will calculate the linear predictor for theta and zeta;if type is \emph{"precision"} the function will calculate fitted values for zeta parameter;if type is \emph{"variance"} the function will calculate fitted variances of response; if type is \emph{"quantile"} the function will calculate fitted quantiles of theta values;
#'@param na.action Characters provided or treatment used in NA values. If \code{na.action} is equal to exclude (default value), the row containing the NA will be excluded in all variables of the model. If \code{na.action} is equal to replace, the row containing the NA will be replaced by the average of the variable in all variables of the model.
#'@param at numeric vector indicating the quantiles that will be informed by the function (only if \code{type = "quantile"}). Its default is 0.5.
#'@param ... further arguments passed to or from other methods.
#'@seealso \code{\link{bayesbr}},\code{\link{residuals.bayesbr}},\code{\link{pmse}}
#'@examples
#'data("CarTask", package = "bayesbr")
#'
#'bbr = bayesbr(probability~ NFCCscale,data=CarTask,
#'              iter = 100, mean_betas = c(1,1.2))
#'
#'predict(bbr, type = "response")
#'predict(bbr, type = "link")
#'predict(bbr, type = "precision")
#'predict(bbr, type = "variance")
#'predict(bbr, type = "quantile", at = c(0.25, 0.5, 0.75))
#'
#'
#'df = data.frame(NFCCscale = rnorm(10,4,1.4))
#'
#'predict(bbr, newdata = df, type = "response")
#'predict(bbr, newdata = df, type = "link")
#'predict(bbr, newdata = df, type = "precision")
#'predict(bbr, newdata = df, type = "variance")
#'predict(bbr, newdata = df, type = "quantile", at = c(0.25, 0.5, 0.75))
#' @export
predict.bayesbr = function(object,newdata = NULL,type = c("response","link", "precision",
                                      "variance","quantile"),
                           na.action=c("exclude","replace"),
                           at = 0.5,...){
 type = match.arg(type)
 if(!is.null(newdata)){
   newdata = data.frame(newdata)
   na.action = match.arg(na.action)
   variaveis = colnames(newdata)
   n = nrow(newdata)
   x_names = object$info$names$names_x
   x_names = x_names[x_names != "(Intercept)"]
   w_names = object$info$names$names_w
   w_names = w_names[w_names != "(Intercept)"]
   p = object$info$p
   q = object$info$q
   formula = object$formula
   formula = Formula(as.formula(paste('~',as.character(formula)[3])))
   model_frame <- model.frame(formula, data = newdata)

   X = as.data.frame(model.matrix(formula, data = model_frame, rhs = 1))
   names_x = colnames(X)
   W = as.data.frame(model.matrix(formula, data = model_frame, rhs = 2))
   names_w = colnames(W)

   if(is.null(names_x)){
      X = NULL
   }
   else{X = as.matrix(X)}
   if(is.null(names_w)){
      W = NULL
   }
   else{W = as.matrix(W)}

   if(na.action == "exclude"){
     model = data.frame(cbind(X,W))
     na_values = which(is.na(model), arr.ind=TRUE)
     if(nrow(na_values)>0){
       model = drop_na(model)
       aux1 = p
       aux2 = 1+p
       aux3 = p+q
       if(p>0){
       X = model[,(1:aux1)]
       }
       if(q>0){
         W = model[,(aux2:aux3)]
       }
     }
   }
   if(na.action == "replace"){
       model = cbind(X,W)
       na_values = which(is.na(model), arr.ind=TRUE)
       if(nrow(na_values)>0){
         for(i in 1:nrow(na_values)){
           row = na_values[i,1]
           col = na_values[i,2]
           mean = mean(model[,col],na.rm=TRUE)
           model[row,col] = mean
         }
         aux1 = p
         aux2 = 1+p
         aux3 = p+q
         if(p>0){
           X = model[,(1:aux1)]
         }
         if(q>0){
           W = model[,(aux2:aux3)]
         }
       }
   }
   linear_predictor_theta = NULL
   linear_predictor_zeta = NULL
   if(p>0){
     mat_betas = c()
     for(i in 1:p){
       betas = object$info$samples$beta
       mat_betas = rbind(mat_betas,betas[[i]])
     }
     linear_predictor_theta = X %*% mat_betas
     linear_predictor_theta = round(linear_predictor_theta,5)
     theta = exp(linear_predictor_theta)/(1+exp(linear_predictor_theta))
   }
   else{
     theta = object$info$samples$theta$theta
     theta = matrix(theta,nrow = n, ncol=length(theta),byrow = T)
   }

   if(q>0){
     mat_gammas = c()
     for(i in 1:q){
       gammas = object$info$samples$gamma
       mat_gammas = rbind(mat_gammas,gammas[[i]])
     }
     linear_predictor_zeta = W %*% mat_gammas
     linear_predictor_zeta = round(linear_predictor_zeta,5)
     zeta = exp(linear_predictor_zeta)
   }
   else{
     zeta = object$info$samples$zeta$zeta
     zeta = matrix(theta,nrow=n,ncol=length(zeta),byrow = T)
   }

   mean_theta = c()
   mean_zeta = c()
   for(i in 1:n){
      mean_theta = c(mean_theta,mean(theta[i,]))
      mean_zeta = c(mean_zeta,mean(zeta[i,]))
   }
   names(mean_theta) = 1:n
   names(mean_zeta) = 1:n

   if(type=="response"){
      return(round(mean_theta,5))
  }
   else if(type == "link"){
     lista = list()
     if(is.null(linear_predictor_theta)){
       odds = mean_theta/(1-mean_theta)
       pred_linear = log(odds)
       pred_linear = round(pred_linear,5)
       names(pred_linear) = 1:n
       lista$mean = pred_linear
     }
     else{
       pred_linear = rowMeans(linear_predictor_theta)
       pred_linear = round(pred_linear,5)
       names(pred_linear) = 1:n
       lista$mean = pred_linear
     }

     if(is.null(linear_predictor_zeta)){
       pred_linear = log(mean_zeta)
       pred_linear = round(pred_linear,5)
       names(pred_linear) = 1:n
       lista$precision = pred_linear
     }
     else{
       pred_linear = rowMeans(linear_predictor_zeta)
       pred_linear = round(pred_linear,5)
       names(pred_linear) = 1:n
       lista$precision = pred_linear
     }
     return(lista)
   }
   else if(type == "precision"){
      return(round(mean_zeta,5))
   }
   else if(type == "quantile"){
     quantis = c()
     for (i in 1:n) {
       v_theta = theta[i,]

       quantis = rbind(quantis,round(quantile(v_theta,at),5))
       }
     rownames(quantis) = 1:n
     colnames(quantis) = paste0(at*100,"%")
     return(quantis)
   }
   else if(type=="variance"){
     variance = mean_theta*(1-mean_theta)/(1+mean_zeta)
     variance = round(variance,5)
     names(variance) = 1:n
     return(variance)
   }
 }
 else{
   n = object$info$n
   fitted = object$fitted.values
   theta = object$info$samples$theta
   zeta = object$info$samples$zeta


   if(length(fitted) == 1){
     fitted = rep(fitted,n)
   }
   names(fitted) = 1:n

   precision = c()
   if(length(zeta)==1){
     v_zeta = zeta$zeta
     meanv_theta = c(round(mean(v_zeta),5))
     if(length(meanv_theta) == 1){
       meanv_theta = rep(meanv_theta,n)
     }
     precision = meanv_theta
   }
   else{
     for (i in 1:n) {
       v_zeta = c()
       aux = paste0('zeta[',i,']')
       v_zeta = zeta[[aux]]
       precision = c(precision,round(mean(v_zeta),5))
     }
   }
   names(precision) = 1:n



   if(type=="response"){
     return(fitted)
   }
   else if(type == "link"){
      lista = list()
      odds = fitted/(1-fitted)
     pred_linear = log(odds)
     pred_linear = round(pred_linear,5)
     names(pred_linear) = 1:n
     lista$mean = pred_linear

     pred_linear = log(precision)
     pred_linear = round(pred_linear,5)
     names(pred_linear) = 1:n
     lista$precision = pred_linear
     return(lista)
   }
   else if(type == "precision"){

      return(precision)
   }
   else if(type == "quantile"){
     quantis = c()
     if(length(theta) == 1){
        quantis = rbind(quantis,round(quantile(theta$theta,at),5))
        rownames(quantis) = "theta"
      }
     else{
       aux_v = c()
       for (i in 1:n) {
         v_theta = c()
         aux = paste0('theta[',i,']')
         aux_v = c(aux_v,aux)
         v_theta = theta[[aux]]
         quantis = rbind(quantis,round(quantile(v_theta,at),5))
       }
       rownames(quantis) = aux_v
     }
     colnames(quantis) = paste0(at*100,"%")
     return(quantis)
   }
   else if(type=="variance"){
     variance = fitted*(1-fitted)/(1+precision)
     variance = round(variance,5)
     names(variance) = 1:n
     return(variance)
   }
 }

}
