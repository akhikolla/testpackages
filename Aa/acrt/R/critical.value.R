#               
#    Copyright (C) 2016  David Preinerstorfer
#    david.preinerstorfer@econ.au.dk
#
#    This file is a part of acrt.
#
#    acrt is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

critical.value <- function(alpha,                 
               ar.order.max,          
               bandwidth,
               ker,                   
               R,                 
               X,
               N0,
               N1,      
               N2,
               Mp,
               M1,
               M2,
               Eicker = FALSE,          
               opt.method.1 = "Nelder-Mead",
               opt.method.2 = "Nelder-Mead",  
               control.1 = list("reltol" = N1^(-.5), "maxit" = dim(X)[1]*20),
               control.2 = list("reltol" = N2^(-.5), "maxit" = dim(X)[1]*30),   
               cores = 1, 
               margin = rep(1, length = ar.order.max)
               ){

#compute dimension variables

n <- dim(X)[1]
if(ar.order.max >= n){
 ar.order.max <- n-1
}
k <- dim(X)[2]
q <- dim(R)[1]

#run input checks

check.alpha(alpha)
check.bandwidth(bandwidth)
check.ker(ker)
check.Eicker(Eicker)
check.X.R.order(X, R, ar.order.max)
check.N.M(N0, N1, N2, Mp, M1, M2)
check.cores(cores)
if( ar.order.max > 0 ){
check.margin(margin, ar.order.max)
}
check.method(opt.method.1, opt.method.2)

#compute auxiliary quantities

qrX <- qr(X)
Bmat <- Bfactor.matrix(qrX, R)
Coefpremult <- solve(t(qr.R(qrX)) %*% qr.R(qrX))%*%t(qr.R(qrX))%*%t(qr.Q(qrX))
Wmat <- wm(n, bandwidth, ker)

#if ar.order.max = 0 no optimization is needed:

if( ar.order.max == 0 ) {
y <- matrix(rnorm(N2 * n), ncol = N2, nrow = n)
vals <- 
 .Call('acrt_testvals', PACKAGE = 'acrt', y, Coefpremult, 0, 
       X, Wmat, Bmat, R, n, dim(y)[2], q, cores, Eicker)
return(list("critical.value" = quantile(vals, 1-alpha)))
}                       


#objective function (depends also on sample y) - minimize quantile

objective.fun <- function(partial){
 partial.transf <- 2/pi * atan(partial) * margin
 vals <- 
 .Call('acrt_testvals', PACKAGE = 'acrt', y, Coefpremult, partial.transf, 
       X, Wmat, Bmat, R, n, dim(y)[2], q, cores, Eicker)
 return(-quantile(vals, 1-alpha))
}

#create storage space for results

fs.results <- matrix(NA, nrow = M1, ncol = 3*(ar.order.max + 1)+1)

#generate starting values
show("Randomized search for starting values initiated.")
y <- matrix(rnorm(N0 * n), ncol = N0, nrow = n)

if( ar.order.max > 2 ) {
 order.starting <- c(2, 1:ceiling(ar.order.max/5) * 5)
 order.starting <- order.starting[order.starting <= ar.order.max]
 if((ar.order.max %% 5) != 0){
 order.starting <- c(order.starting, ar.order.max)
 }
} else {
order.starting <- ar.order.max
}

for(ar.start in order.starting){
for(i in 1:Mp){
tpara <- c(gen.start(ar.start), rep(0, ar.order.max - ar.start))
val <- objective.fun(tan(pi/2 * tpara))
M <- rbind(fs.results[, 1:(ar.order.max + 1)], c(tpara, val))
M <- M[order(M[,ar.order.max + 1]),]
fs.results[, 1:(ar.order.max + 1)] <- M[1:M1,]
}
}

#initiate first stage optimizations
show("First stage optimizations initiated.")
for(i in 1:M1){
y <- matrix(rnorm(N1 * n), ncol = N1, nrow = n)
startval <- tan(pi/2 * fs.results[i,1:ar.order.max])
max.quantile <- optim(startval, objective.fun, method = opt.method.1,
                      control = control.1)
fs.results[i,(ar.order.max+2):(2*ar.order.max + 1)] <- 
                                                  2/pi * atan(max.quantile$par)
fs.results[i,2*ar.order.max+2] <- max.quantile$value
}

#initiate second stage optimizations

index.top <- order(fs.results[,2*ar.order.max+2])[1:M2]
show("Second stage optimizations initiated.")
for(i in index.top){
y <- matrix(rnorm(N2 * n), ncol = N2, nrow = n)
startval <- tan(pi/2 * fs.results[i,(ar.order.max+2):(2*ar.order.max + 1)])
max.quantile <- optim(startval, objective.fun, method = opt.method.2,
                      control = control.2)
fs.results[i,(2*ar.order.max+3):(3*ar.order.max + 2)] <- 
                                                  2/pi * atan(max.quantile$par)
fs.results[i,3*ar.order.max+3] <- max.quantile$value
fs.results[i,3*ar.order.max+4] <- max.quantile$convergence
}

#return starting values, corresponding optimizers, and maximal quantiles

return(list(
"starting.parameters" = fs.results[,1:ar.order.max],
"starting.quantiles"     = - fs.results[,ar.order.max+1],
"first.stage.parameters" = fs.results[,
          (ar.order.max+2):(2*ar.order.max + 1)]*margin,
"first.stage.quantiles" = - fs.results[,2*ar.order.max+2],
"second.stage.parameters" = fs.results[index.top,
          (2*ar.order.max+3):(3*ar.order.max + 2)]*margin, 
"second.stage.quantiles" = - fs.results[index.top,3*ar.order.max+3],
"convergence" = fs.results[index.top,3*ar.order.max+4],
"critical.value" = max(- fs.results[index.top, 3*ar.order.max+3])
            )
      )
} 