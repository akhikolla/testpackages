`[.clifford` <- function(C, index, ...){
    if(is.clifford(index)){
        stop("cannot extract a clifford; try A[terms(B)]")
    } else if(is.list(index)){
        dots <- index
    } else {
        dots <- c(list(index),list(...))
    }
    clifford(dots,getcoeffs(C,dots))
}  

`[<-.clifford` <- function(C, index, ..., value){
    
    if(missing(index)){ # C[] <- value
        dots <- list(...)
        if(is.clifford(value)){
            return(as.clifford(c_overwrite(
                terms(C),coeffs(C),
                terms(value),coeffs(value),
                maxyterm(C,value)
            )))
        } else { # value a scalar
            stopifnot(length(value) == 1)
            return(clifford(terms(C),value + numeric(length(coeffs(C)))))
        }
    } else {  # index supplied, dots interpreted as more terms
        dots <- c(list(index),list(...))
        if(value==0){
            jj <- clifford(dots,1)
            return(as.clifford(c_overwrite(
                terms(C),coeffs(C),
                terms(jj),coeffs(jj),
                maxyterm(C,jj)
            ))-jj)
        } else { # value != 0
            jj <- clifford(dots,value) # sic; this is legit!
            return(as.clifford(c_overwrite(
                terms(C),coeffs(C),
                terms(jj),coeffs(jj),
                maxyterm(C,jj)
            )))
        }
    }
}

