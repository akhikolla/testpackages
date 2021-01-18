"Ops.freealg" <- function (e1, e2 = NULL) 
{
    oddfunc <- function(...){stop("odd---neither argument has class free?")}
    unary <- nargs() == 1
    lclass <- nchar(.Method[1]) > 0
    rclass <- !unary && (nchar(.Method[2]) > 0)
    
    if (!is.element(.Generic, c("+", "-", "*", "^", "=="))){
        stop("operator '", .Generic, "' is not implemented for free algebra objects")
    }

    if(unary){
        if (.Generic == "+") {
            return(e1)
        } else if (.Generic == "-") {
            return(free_negative(e1))
        } else {
            stop("Unary operator '", .Generic, "' is not implemented for free algebra objects")
        }
    }
    
    if (.Generic == "*") {
        if (lclass && rclass) {
            return(free_times_free(e1, e2))
        } else if (lclass) {
            return(free_times_scalar(e1, e2))
        } else if (rclass) {
            return(free_times_scalar(e2, e1))
        } else {
            oddfunc()
        }
    } else if (.Generic == "+") {
        if (lclass && rclass) {
            return(free_plus_free(e1, e2))  # S1+S2
        } else if (lclass) {
          return(free_plus_numeric(e1, e2)) # S+x
        } else if (rclass) {
          return(free_plus_numeric(e2, e1)) # x+S
        } else {
          oddfunc()
        }
    } else if (.Generic == "-") {
      if (lclass && rclass) {
        return(free_plus_free(e1, free_negative(e2)))  # S1-S2
      } else if (lclass) {
        return(free_plus_numeric(e1, -e2))                # S-x
      } else if (rclass) {
        return(free_plus_numeric(free_negative(e2), e1)) # x-S
      } else {
        oddfunc()
      }
    } else if (.Generic == "^") {
      if(lclass && !rclass){
        return(free_power_scalar(e1,e2)) # S^n
      } else {
        stop("Generic '^' not implemented in this case")
      }
    } else if (.Generic == "==") {
      if(lclass && rclass){
        return(free_eq_free(e1,e2))
      } else {
        stop("Generic '==' only compares two freealg objects with one another")
      }          
    } else if (.Generic == "!=") {
      if(lclass && rclass){
        return(!free_eq_free(e1,e2))
      } else {
        stop("Generic '!=' only compares two free objects with one another")
      }
    } else if (.Generic == "/") {
      if(lclass && !rclass){
        return(free_times_scalar(e1,1/e2))
      } else {
        stop("Generic '/' not supported for freealg objects")
      }
    }
}

`free_negative` <- function(S){
    if(is.zero(S)){
        return(S)
    } else {
        return(freealg(words(S), -coeffs(S)))
    }
}

`free_times_free` <- function(S1,S2){
  if(is.zero(S1) || is.zero(S2)){
    return(constant(0))
  } else {
      jj <- lowlevel_free_prod(
          words1=S1[[1]],coeffs1=S1[[2]],
          words2=S2[[1]],coeffs2=S2[[2]]
      )
      return(freealg(jj[[1]],jj[[2]]))
  }
}

`free_times_scalar` <- function(S,x){
freealg(S[[1]],x*S[[2]])
}

`free_plus_free` <- function(e1,e2){
  if(is.zero(e1)){
        return(e2)
    } else if(is.zero(e2)){
        return(e1)
    } else {
        jj <- lowlevel_free_sum(
          words1=e1[[1]],coeffs1=e1[[2]],
          words2=e2[[1]],coeffs2=e2[[2]]
        )
        return(freealg(jj[[1]],jj[[2]]))
    }
}

`free_plus_numeric` <- function(S,x){
    free_plus_free(S,numeric_to_free(x))
}

free_power_scalar <- function(S,n){
  stopifnot(n==round(n))
  if(n<0){
    stop("negative powers not implemented")
  } else if(n==0){
    return(as.freealg(1))
  } else {
      jj <- lowlevel_free_power(S[[1]],S[[2]],n)
      return(freealg(jj[[1]],jj[[2]]))
  }
}

`free_eq_free` <- function(e1,e2){
  is.zero(e1-e2)  # nontrivial; S1 and S2 might have different orders
}
