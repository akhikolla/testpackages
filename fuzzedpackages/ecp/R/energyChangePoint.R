 getWithin <- function (alpha_, X_) 
.Call("getWithin", alpha_, X_, PACKAGE = "ecp")

 getBetween <- function (alpha_, X_, Y_) 
.Call("getBetween", alpha_, X_, Y_, PACKAGE = "ecp")

 splitPointC <- function (s_, e_, D_, min_size_) 
.Call("splitPointC", s_, e_, D_, min_size_, PACKAGE = "ecp")

getBounds <- function(n_, lvl_, eps_)
.Call("getBounds", n_, lvl_, eps_, PACKAGE = "ecp")


eFastC <- function(Z_, K_, minsize_, alpha_, verbose_ )
.Call("eFastC", Z_, K_, minsize_, alpha_, verbose_, PACKAGE = "ecp")

eFastC_delta <- function(Z_, K_, delta_, alpha_, verbose_ )
.Call("eFastC_delta", Z_, K_, delta_, alpha_, verbose_, PACKAGE = "ecp")

ksFastC <- function(Z_, K_, minsize_, verbose_ )
.Call("ksFastC", Z_, K_, minsize_, verbose_, PACKAGE = "ecp")

ksFastC_delta <- function(Z_, K_, minsize_, verbose_ )
.Call("ksFastC_delta", Z_, K_, minsize_, verbose_, PACKAGE = "ecp")


srcGetV <- function(K_)
.Call("srcGetV", K_, PACKAGE = "ecp")

srcGetBandwidth <- function(X_, rows_)
.Call("srcGetBandwidth", X_, rows_, PACKAGE = "ecp")

srcKcpa <- function(II_, V_, H_)
.Call("srcKcpa", II_, V_, H_, PACKAGE = "ecp")