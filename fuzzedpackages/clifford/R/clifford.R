`clifford` <- function(terms,coeffs=1){
    if(length(coeffs)==1){coeffs <- coeffs+numeric(length(terms))}
    stopifnot(is_ok_clifford(terms,coeffs))
    m <- mymax(c(terms,recursive=TRUE))
    out <- c_identity(terms,coeffs,m)
    class(out) <- "clifford"  # this is the only place class clifford is set
    return(out)
}

`terms` <- function(x){ x[[1]] }  # accessor methods start
`coeffs` <- function(x){ x[[2]] }

`getcoeffs` <- function(C,B){ # accessor methods end
    c_getcoeffs(
        L = terms(C),
        c = coeffs(C),
        m = maxyterm(C),
        B = B)
}

`const` <- function(C,drop=TRUE){
  out <- getcoeffs(C,list(numeric(0)))
  if(drop){
    return(out)
  } else {
    return(as.clifford(out))
  }
}


`is.1vector` <- function(x){all(grades(x)==1)}
`is.basisblade` <- function(x){ (nterms(x)==1) || is.scalar(x) }
`is.blade` <- function(x){stop("Not implemented: factorization is hard")}
`coeffs<-` <- function(x,value){UseMethod("coeffs<-")}
`coeffs<-.clifford` <- function(x,value){
    stopifnot(length(value) == 1)
    return(clifford(terms(x),value + 0*coeffs(x)))
}

`const<-` <- function(x,value){UseMethod("const<-")}
`const<-.clifford` <- function(x,value){
    stopifnot(length(value) == 1)
    x <- x-const(x)
    return(x+value)
}



`mymax` <- function(s){
    if(length(s)==0){
        return(0)
    } else {
        return(suppressWarnings(max(s)))
    }
}

`is_ok_clifford` <- function(terms,coeffs){
    stopifnot(is.list(terms))
    term_elements <- c(terms,recursive = TRUE)
    stopifnot(length(terms) == length(coeffs))

    if(!is.null(term_elements)){
      stopifnot(all(term_elements > 0))
    }
    
    term_elements_increase <- c(lapply(terms,diff),recursive=TRUE)
    stopifnot(all(term_elements_increase > 0))

    return(TRUE)
}

`as.clifford` <- function(x){
    if(inherits(x,"clifford")){
        return(x)
    } else if(is.list(x)){
        return(clifford(x[[1]],x[[2]]))
    } else if(is.numeric(x)){
        return(numeric_to_clifford(x))
    } else if(is.null(x)){
        return(clifford(list(numeric(0)),0))
    } else {
        stop("not recognised")
    }
}

`numeric_to_clifford` <- function(x){
  if(length(x)==1){
    return(as.scalar(x))
  } else {
    return(as.1vector(x))
  }
}

`is.clifford` <- function(x){inherits(x,"clifford")}

`nterms` <- function(x){length(coeffs(x))}
`is.zero` <- function(C){nterms(C)==0}
`nbits` <- function(x){max(c(terms(x),recursive=TRUE))}
`grades` <- function(x){unlist(lapply(terms(x),length))}
`is.scalar` <- function(C){
  if(is.zero(C)){
    return(TRUE)
  } else {
    return((length(terms(C))==1) && (length(terms(C)[[1]])==0))
  }
}

`scalar` <- function(x=1){clifford(list(numeric(0)),x)}
`as.scalar` <- `scalar`
`as.1vector` <- function(x){clifford(as.list(seq_along(x)),x)}
`pseudoscalar` <- function(n,x=1){clifford(list(seq_len(n)),x)}
`as.pseudoscalar` <- `pseudoscalar`
`is.pseudoscalar` <- function(C){
    if(is.zero(C)){return(TRUE)}
    return(
        is.clifford(C)         &&
        (length(terms(C))==1) &&
        all(terms(C)[[1]] == seq_along(terms(C)[[1]]))
    )
}

`antivector` <- function(v,n=length(v)){
    stopifnot(n>=length(v))
    clifford(sapply(seq_along(v),function(i){seq_len(n)[-i]},simplify=FALSE),v)
}

`is.antivector` <- function(C, include.pseudoscalar=FALSE){

  if(!is.homog(C)){return(FALSE)}
  if(include.pseudoscalar && is.pseudoscalar(C)){return(TRUE)}
  return(grades(C)[1] == maxyterm(C)-1)
}

`basis` <- function(n,x=1){clifford(list(n),x)}
`e` <- basis

`rcliff` <- function(n=9,d=6,g=4,include.fewer=TRUE){
  if(include.fewer){
    f <- function(...){sample(g,1)}
  } else {
    f <- function(...){g}
  }
  clifford(replicate(n,sort(sample(d,f())),simplify=FALSE),sample(n))
} 

`rblade` <- function(d=9, g=4){
  Reduce(`%^%`,sapply(seq_len(g),function(...){as.1vector(rnorm(d))},simplify=FALSE))
}

`rev.clifford` <- function(x){
  f <- function(u){ifelse(length(u)%%4 %in% 0:1, 1,-1)}
  clifford(
      terms(x),
      coeffs(x) * unlist(lapply(terms(x),f))
  )
}

`Conj.clifford` <- function(z){
  z <- rev(z)
  clifford(
      terms(z),
      coeffs(z) * ifelse(gradesminus(z)%%2==0,1,-1)
  )
}

`print.clifford` <- function(x,...){
  cat("Element of a Clifford algebra, equal to\n")

  out <- ""
  for(i in seq_along(terms(x))){
    co <- coeffs(x)[i]
    if(co>0){
      pm <- " + " # pm = plus or minus
    } else {
      pm <- " - "
    }
    co <- capture.output(cat(abs(co)))
    jj <- catterm(terms(x)[[i]])
    out <- paste(out, pm, co, jj, sep="")
  }
  if(is.zero(x)){
      out <- "the zero clifford element (0)"
  } else if(is.scalar(x)){
      out <- paste("scalar (",capture.output(cat(coeffs(x))),")")
  }
  cat(paste(strwrap(out, getOption("width")), collapse="\n"))
  cat("\n")
  return(x)
}

`drop` <- function(C){
    if(is.zero(C)){
        return(0)
    } else if(is.scalar(C)){
        return(const(C))
    } else {
        return(C)
    }
}

`grade` <- function(C,n,drop=TRUE){
  C <- as.clifford(C)
  out <- as.clifford(c_grade(terms(C),coeffs(C),maxyterm(C),n))
  if(drop){out <- drop(out)}
  return(out)
}

`is.homog` <- function(C){
    C <- as.clifford(C)
    if(is.zero(C)){return (TRUE)}
    if(is.scalar(C)){ return(TRUE) }
    g <- grades(C)
    if(min(g) == max(g)){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

`scalprod` <- function(C1,C2=rev(C1),drop=TRUE){grade(C1*C2,0,drop=drop)}
`eucprod` <- function(C1,C2=C1,drop=TRUE){grade(C1*Conj(C2),0,drop=drop)}

`Mod.clifford` <- function(z){sqrt(eucprod(z))}

`is.even` <- function(C){all(grades(C)%%2==0)}
`is.odd`  <- function(C){all(grades(C)%%2==1)}

`evenpart` <- function(C){
    wanted <- which(unlist(lapply(terms(C),length))%%2==0)
    clifford(terms(C)[wanted],coeffs=coeffs(C)[wanted])
}

`oddpart` <- function(C){
    wanted <- which(unlist(lapply(terms(C),length))%%2==1)
    clifford(terms(C)[wanted],coeffs=coeffs(C)[wanted])
}

`allcliff` <- function(n){
    clifford(apply(expand.grid(rep(list(0:1),n))>0,1,which),1)
}

`zap` <- function(x,drop=TRUE,digits=getOption("digits")){
    out <- clifford(terms(x),base::zapsmall(coeffs(x),digits=digits))
    if(drop){out <- drop(out)}
    return(out)
}

`catterm` <- function(a){
  if(length(a)==0){
    return("")
  } else {
    if(isTRUE(getOption("separate"))){
      return(paste("e",a,collapse=" ",sep=""))
    } else {
      jj <- getOption("basissep")
      if(is.null(jj)){jj <- ""}
      return(paste("e_",paste(a,collapse=jj),sep=""))
    }
  }
}

`as.character.clifford` <- function(x,...){
  out <- ""
  for(i in seq_along(terms(x))){
    co <- coeffs(x)[i]
    if(co>0){
      pm <- " + " # pm = plus or minus
    } else {
      pm <- " - "
    }
    co <- abs(co)
    jj <- catterm(terms(x)[[i]])
    out <- paste(out, pm, co, jj, sep="")
  }

  if(is.zero(x)){
    out <- "0 "
  } else if(is.scalar(x)){
    out <- as.character(coeffs(x))
  } 
  return(out)
}

`gradesplus` <- function(x){
    sig <- signature()
    if(sig==0){
        return(grades(x))
    } else if(sig>0){
        return(unlist(lapply(terms(x),function(o){sum(o <= sig)})))
    } else if(sig<0){
        return(grades(x)*NA)
    } else {
        stop("this cannot happen")
    }
}

`gradesminus` <- function(x){ grades(x) - gradesplus(x) }

`dual` <- function(C,n){ C*clifford_inverse(pseudoscalar(n)) }

`neg` <- function(C,n){clifford(terms(C),coeffs(C)*ifelse(grades(C) %in% n,-1,1))}

`first_n_last` <- function(x){
  n <- nterms(x)
  paste(
      as.character(clifford(list(x[[1]][[1]]),x[[2]][1])), " ...",
      as.character(clifford(list(x[[1]][[n]]),x[[2]][n]))
  )
}

`summary.clifford` <- function(object, ...){
  out <- list(
      first_n_last = first_n_last(object),
      nterms       = nterms(object),
      magnitude    = eucprod(object)
  )
  class(out) <- "summary.clifford"
  return(out)
}

"print.summary.clifford" <- function(x, ...){
  cat("Element of a Clifford algebra \n")  
  cat("Typical terms: ", x[[1]],"\n") 
  cat("Number of terms:", x[[2]],"\n") 
  cat("Magnitude:", x[[3]],"\n") 
}

`as.vector.clifford` <- function(x, mode="any"){
    tx <- terms(x)
    stopifnot(all(lapply(tx,length)==1))
    tx <- unlist(tx)
    out <- rep(0,max(tx))
    out[tx] <- coeffs(x)
    return(out)
}

