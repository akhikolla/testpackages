plotl0 <- function(re){
 if (length(re)==3){  y = re$y
  beta = re$beta
  n = length(y)
  plot(1:n, y, xlab = "index", ylab = "value", main = "L0 fitted coefficient")
  points(beta, type = "l", col = "red")}
 if(length(re)==5)
 {
   x = re$df
   y = re$y
   df = re$df
   y = matrix(rep(y, length(df)), ncol = length(df))
   beta.all = re$beta.all
   res = (y-beta.all)^2
   mse = colMeans(res)
   plot(x, mse, type = "l", xlab  = "number of knots", ylab  = "MSE", col = "red")
 }
  
  
}
