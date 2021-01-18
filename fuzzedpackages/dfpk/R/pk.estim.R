#' @useDynLib dfpk, .registration = TRUE
#' @export
pk.estim <-
function(par,t,dose,conc){
    ka <- par[1]
    CL <- par[2]
    V <- par[3]
    s<-dose*ka/V/(ka-CL/V)* (exp(-CL/V*t)-exp(-ka*t))
    if (any(is.nan(s)) ) cat("ka=", ka, "\n V=",V,"\n CL=",CL, "\n")
    s[s==0] = rep(2^(-1074),length(s[s==0] ))
    sum( (log(s)-log(conc))^2 )
}
