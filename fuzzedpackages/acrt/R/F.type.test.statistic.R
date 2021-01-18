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

F.type.test.statistic <- function(y, R, r, X, bandwidth, ker, Eicker = FALSE, 
                                  cores = 1){
                        
#transform y to a matrix in case it is a vector

if(is.vector(y) == TRUE){
y <- as.matrix(y, nrow = length(y))
} 
                        
#input checks

check.X.R.order(X, R, 0)
check.y(y, X)
check.r(r, R)
check.bandwidth(bandwidth)
check.ker(ker)
check.Eicker(Eicker)
check.cores(cores)

#computation of input to c++ functions

qrX <- qr(X)
umat <- qr.resid(qrX, y) 
Rbmat <- R %*% qr.coef(qrX, y) - matrix(r, byrow = F, nrow = dim(R)[1],
                                                      ncol = dim(y)[2])
Wmat <- wm(dim(X)[1], bandwidth, ker)
Bmat <- Bfactor.matrix(qrX, R)

#call c++ functions 

if(Eicker){
test.val <-
.Call('acrt_ctestE', PACKAGE = 'acrt', umat, Rbmat, Wmat, Bmat, cores)
} else {
test.val <- 
.Call('acrt_ctest', PACKAGE = 'acrt', umat, Rbmat, Wmat, Bmat, cores)
}
return(list("test.val" = test.val))
}      
                                                                      