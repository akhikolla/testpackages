#' Numerically stable summation of logs
#' @keywords internal
log_sum_exp <- function(x){
    if(length(x) == 1)
        return(x)
    xmax <- which.max(x)
    if(is.finite(x[xmax]))
        return(log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax])
    ## in case the maximum is not finite, we let R figure out what is
    ## the correct result (usually -Inf or +Inf)
    return(log(sum(exp(x))))
} 
