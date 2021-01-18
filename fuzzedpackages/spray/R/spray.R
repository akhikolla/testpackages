setClass("spray",
         representation = representation(index="numeric",value="numeric"),
         prototype      = list(index=numeric(),value=numeric()),
         )

`spraymaker` <- function(L,addrepeats=FALSE,arity=ncol(L[[1]])){    # formal; this is the *only* way to create a spray object; 
    stopifnot(is_valid_spray(L))
    if(is.empty(L)){
      M <- matrix(0,0,arity)
      x <- numeric(0)
      out <- list(M,x)
    } else {
        M <- L[[1]]
        x <- L[[2]]
        if(!addrepeats & anyDuplicated(M)){ 
            stop("repeated indices; yet argument 'addrepeats' is FALSE")
        }
        out <- spray_maker(M,x)   # see RcppExports.R; strips out zeros and puts index rows in (some inscrutable) order
    }

    class(out) <- "spray"  # this is the *only time class<-() is called
    return(out)
}

`is_valid_spray` <- function(L){  # L = list(M,x)
    stopifnot(is.list(L))
    stopifnot(length(L)==2)
    if(is.empty(L)){
        return(TRUE)
    } else {
        stopifnot(is.matrix(L[[1]]))
        stopifnot(is.vector(L[[2]]))
        stopifnot(nrow(L[[1]])==length(L[[2]]))
        return(TRUE)
    }
}

`is.spray` <- function(S){inherits(S,"spray")}

`is.empty` <- function(L){
  if( is.null(L[[1]]) & is.null(L[[2]])){
    return(TRUE)
  } else if(  identical(nrow(L[[1]]),0L) & (length(L[[2]])==0)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
 
`is.zero` <- is.empty

`spray` <- function(M,x,addrepeats=FALSE){
    if(is.vector(M)){
        M <- rbind(M)
    } else if(inherits(M,"partition")){
        M <- unclass(t(M)) # saves typing t() in 'spray(parts(4))' which drives me nuts
    } else if(is.spray(M)){
        if(length(x)>1){
            stop("spray(M,x) not defined if length(x)>1")
        }
        M <- index(M)
    } else {
        M <- as.matrix(M) # e.g., coerces expand.grid(...) to a matrix
    }
    if(missing(x)){x <- rep(1,nrow(M))}
    if(length(x)==1){x <- rep(x,nrow(M))}
    return(spraymaker(list(M,x),addrepeats=addrepeats))
}

`index` <- function(S){S[[1]]}    # these two functions are the only
`value` <- function(S){S[[2]]}    # 'accessor' functions in the package

`value<-` <- function(S,value){
    stopifnot(length(value)==1)
    spray(index(S),value)
}

`as.spray` <- function(arg1, arg2, addrepeats=FALSE, offbyone=FALSE){  # tries quite hard to coerce things to a spray
    if(is.spray(arg1)){
        return(arg1)  
    } else if(is.list(arg1)){
        return(spraymaker(arg1,addrepeats=addrepeats))
    } else if(is.matrix(arg1) & !missing(arg2)){
        return(spraymaker(list(arg1,arg2),addrepeats=addrepeats))
    } else if(is.array(arg1)){
        ind <- which(arg1!=0,arr.ind=TRUE)
        if(offbyone){ ind <- ind-1 }
        return(spray(ind, arg1[arg1!=0]))
    } else {
        stop("as.spray() needs a list or an array")
    }  
}

`arity` <- function(S){
    ncol(index(S))
}

setGeneric("dim")
setGeneric("deriv")

`dim.spray` <- function(x){
    ind <- index(x)
    if(any(ind<0)){
        stop("dim() does not make sense for sprays with negative indices")
    } else {
        return(apply(ind,2,max))
    }
}

`as.array.spray` <- function(x, offbyone=FALSE, compact=FALSE, ...){
    ind <- index(x)
    dS <- dim(x) 
    if(compact){
        jj <- apply(ind,2,min)
        dS <- dS-jj+1
        ind <- 1+sweep(ind,2,jj)
    } 
    if(offbyone) {
      ind <- ind+1
      dS <- dS + 1
    }

    if(any(ind==0)){  # ind<0 detected by dS <- dim(x) above
        stop("There are zero index elements")
    }
    
    out <- array(0,dS)
    out[ind] <- value(x)
    dimnames(out) <- names(ind)
    return(out)
}


`spray_missing_accessor` <- function(S,dots){

 ## This function is a placeholder for a discussion of why idiom such
 ## as S[1,,p,1:3] is problematic (note the adjacent commas).
 

 return(0)
##  browser()



## Not implemented.

## the problem is illustrated by: removing the return() statement from
## above and uncomment the browser() line.  Then change the first
## line of the body of `[.spray`() to something like


##         spray_missing_accessor(S,dots=match.call())



## Then load the package and type
##
##
## R>      p <- 1:3
## R>      S <- spray(diag(1:6)) ;
## R>      S[1,,p,1:3]
##
##
##
##  On typing this, you enter the browser;
## 'dots' is a 6-element list.  The first, dots[[1]], is the function name `[.spray`
##
## and the second, dots[[2]], is just S.   This is the *symbol* S, not the variable.  You
## can evaluate the symbol by typing eval(dots[[2]]) which shows the value of S.

## The third element, dots[[3]], is just a numeric "1".  The next element, dots[[4]], 
## corresponding to the gap between the two consecutive commas, is some sort of ghostly
## not NULL, not a list, nothing that I can make sense of.  For
## example, 'dput(dots[[4]])' gives an empty line; 'dput(dots[4])'
## gives "()" .... previous versions of R gave things like
## structure(list(), .Names = ""); but note that
## structure(list(), .Names = "") will not parse on the commandline.

## So, at least prior to R-3.0.0,  one cannot use constructions such as
## 'identical(dots[2], structure(list(), .Names = ""))'.


## Now dots[[5]] and indeed dput(dots[[5]]) both return the string
## `p`.  This can be evaluated with something like eval(dots[[5]]) but
## you need to be careful about the environment in which p is to be
## evaluated.  This is why Jonty had to resort to peculiar tricks with
## eval.parent() in Oarray.


## Also note that dots[[4]] and dput(dots[[4]]) show "1:3"; note that
##  class(dots[[6]])is "call", not "numeric" as one might expect.

## ... quite apart from all this, there is another reason why I have
## not implemented things like S[1,,1:3].  Such constructions simply
## do not fit nicely into the STL map methods in spray_ops.cpp.  At
## least I can't see any nice way.  You'd need some sort of partial
## matching system.

## I also need to fix OArray.


}

`[.spray` <- function(S, ...,drop=FALSE){

    dots <- list(...)
    first <- dots[[1]]
    if(is.matrix(first)){
        M <- first
    } else if(is.spray(first)){
        M <- index(first)
    } else {
        M <- as.matrix(expand.grid(dots))
    }

    if(ncol(M) != arity(S)){
      stop("incorrect number of dimensions specified")
    }
    jj <- spray_accessor(index(S),value(S),M) # returns a NumericVector
    if(drop){
      return(jj)
    } else {
      return(spray(M,jj))
    }
}  

`[<-.spray` <- function(S, index, ..., value){

    if(missing(index)){ # S[] <- something
        if(is.spray(value)){
            return(
                spraymaker(spray_overwrite(
                    spray::index(S    ),spray::value(S    ),
                    spray::index(value),spray::value(value))))
            } else {
                return(spray(spray::index(S),value))
            }
    }

    if(is.matrix(index)){
        M <- index
    } else if(is.spray(index)){
        M <- spray::index(index)
    } else {
        M <- as.matrix(expand.grid(c(list(index), list(...))))
    }
    if(ncol(M) != arity(S)){
        stop("incorrect number of dimensions specified")
    }

    if(length(value)==1){value <- rep(value,nrow(M))}
    stopifnot(length(value)==nrow(M))
    return(spraymaker(spray_setter(spray::index(S),spray::value(S),M,value)))
}

`as.function.spray` <- function(x,...){
    function(X){
        if(!is.matrix(X)){X  <- rbind(X)}
        stopifnot(ncol(X)==arity(x))
        jj <- matrix(1,nrow(X),nrow(index(x)))
        for(i in seq_len(arity(x))){ # iterate over index columns
            jj <- jj * outer(X[,i],index(x)[,i],"^")
        }
        return(rowSums(sweep(jj,2,value(x),"*")))
    }
}

`constant` <- function(S,drop=FALSE){ # returns 'the constant (term) of S'
  S[t(rep(0,arity(S))),drop=drop]
}

`constant<-` <- function(S,value){ "[<-"(S,t(rep(0,arity(S))),value=value)}

`zero` <- function(d){spray(matrix(0,0,d),numeric(0))}
`product` <- function(power){spray(rbind(power))}
`homog` <- function(d,power=1){spray(partitions::compositions(power,d))}
`linear` <- function(x,power=1){spray(diag(power,length(x)),x)}
`lone` <- function(n,d=n){
  jj <- rep(0,d)
  jj[n] <- 1
  return(product(jj))
}

`xyz` <- function(d){product(rep(1,d))}

`one` <- function(d){
  if(is.spray(d)){d <- arity(d)}
  return(spray(t(rep(0,d))))
}

`ooom` <- function(S,n){
  out <- one(S)
  for(i in seq_len(n)){
    out <- 1  + S*out
  }
  return(out)
}

`print_spray_matrixform` <- function(S){
    if(is.empty(S)){
        cat(paste('empty sparse array with ', arity(S), ' columns\n',sep=""))
    } else {
        jj <-
            data.frame(index(S),symbol= " = ", val=round(value(S),getOption("digits")))
        colnames(jj) <- c(rep(" ",arity(S)+1),'val')
        print(jj,row.names=FALSE)
    }
    return(invisible(S))
}

`print_spray_polyform` <- function(S){

    vector_of_3vars <- c("x","y","z")
    multiply_symbol1 <- "*"
    multiply_symbol2 <- "*"

    if(is.empty(S)){
    cat(paste('the NULL multinomial of arity ', arity(S), '\n',sep=""))
    return(S)
  }

  if(is.null(getOption("sprayvars"))){
    if(arity(S) <= 3){
      variables <- vector_of_3vars
    } else {
      variables <- paste("x",seq_len(arity(S)),sep="")
    }
  } else {
    variables <- getOption("sprayvars")
    }
  
  if(is.empty(constant(S))){
    string <- ""
  } else {
    string <- constant(S,drop=TRUE)
  }
  
  ind <- index(S)
  val <- value(S)

    for(i in seq_len(nrow(ind))){ 
        coeff <- val[i]
        coeff_abs <- abs(coeff)
        coeff_sign <- sign(coeff)
        if(all(ind[i,]==0)){next}  # constant term
        if(coeff_abs == 1){
            term <- ""
        } else {
            term <- paste(round(coeff_abs, getOption("digits")), "",multiply_symbol1, sep="")
        }
        
        if(coeff_sign == 1){
            term <- paste(" +", term, sep="")
        } else {
            term <- paste(" -", term, sep="")
        }

        jj <- ind[i,]
        vars <- variables[which(jj !=0)]
        powers <- jj[jj!=0]
        lv <- length(vars)
        if(lv==1){   # just one variable
            if(powers==1){
                term <- paste(term, vars, sep="") 
            } else {
                term <- paste(term, vars, "^", powers, sep="")
            }
        } else {   ## length(vars) > 1
            for(i in seq_len(lv-1)){
                if(powers[i] == 1){
                    term <- paste(term, vars[i], multiply_symbol2,sep="")
                } else {
                    term <- paste(term, vars[i],"^",powers[i],multiply_symbol2,sep="") 
                }
            }
            
            if(powers[lv] == 1){
                term <- paste(term, vars[lv],sep="") 
            } else {
                term <- paste(term, vars[lv],"^",powers[lv],sep="") 
            }
        }
        string <- paste(string, term, sep="")

      }  # row of index loop closes
    
  string <- paste(string,'\n',sep="")
  for(i in strwrap(string)){
    cat(noquote(i))
    cat("\n")
  }
  return(invisible(S))
}
  
`print.spray` <- function(x, ...){
  if(isTRUE(getOption("polyform",default=FALSE))){
    out <- print_spray_polyform(x)
  } else {
    out <- print_spray_matrixform(x)
  }
  return(out)
}

`asum` <- function(S, dims,drop=TRUE, ...) {UseMethod("asum")}

#      M <- matrix(sample(0:2,5*29,rep=T),ncol=5)

`process_dimensions` <- function(S,dims){
  if(is.logical(dims)){ dims <- which(dims)  }
  stopifnot(all(dims <= arity(S)))
  stopifnot(all(dims >0))
  stopifnot(all(dims == round(dims)))
  return(dims)
}

`asum.spray` <- function(S, dims, drop=TRUE, ...){   # off-by-one dealt with in C++
  dims <- process_dimensions(S,dims)
  jj <- spraymaker(spray_asum_include(index(S),value(S),dims))
  if(drop){
    return(spray(index(jj)[,-dims,drop=FALSE],value(jj),addrepeats=TRUE))
  } else {
    return(jj)
  }
}

`asum_inverted` <- function(S, dims){   # off-by-one dealt with in C++
  dims <- process_dimensions(S,dims)
  return(spraymaker(spray_asum_exclude(index(S),value(S),dims)))
}

`subs` <- function(S,dims,x){

  dims <- process_dimensions(S,dims)

  a <- index(S)[,-dims,drop=FALSE]
  b <- index(S)[, dims,drop=FALSE]
  
  jj <- apply(sweep(b,2,x,function(x,y){y^x}),1,prod)* value(S)
  spray(a, jj, addrepeats=TRUE)

}

`deriv.spray` <- function(expr, i, derivative=1, ...){
  S <- as.spray(expr)
  orders <- rep(0,arity(S))
  orders[i] <- derivative
  return(aderiv(S,orders))
}

`aderiv` <- function(S,orders){
  orders <- round(orders)
  stopifnot(length(orders) == arity(S))
  stopifnot(all(orders >= 0))
  return(spraymaker(spray_deriv(index(S),value(S),orders)))
}

`pmax` <- function(...){ UseMethod("pmax") }
`pmax.default` <- function(..., na.rm=FALSE){ base::pmax(..., na.rm=FALSE) }

`pmin` <- function(...){ UseMethod("pmin") }
`pmin.default` <- function(..., na.rm=FALSE){ base::pmin(..., na.rm=FALSE) }

`pmax.spray` <- function(x, ...) {
    if(nargs()<3){
        return(maxpair_spray(x, ...))
    } else {
        return(maxpair_spray(x, maxpair_spray(...)))
    }
}

`pmin.spray` <- function(x, ...) {
    if(nargs()<3){
        return(minpair_spray(x, ...))
    } else {
        return(minpair_spray(x, minpair_spray(...)))
    }
}

`maxpair_spray` <- function(S1,S2){
    if(missing(S2)){ return(S1) }
    if(!is.spray(S2)){   # S2 interpreted as a scalar
        if(length(S2)>1){
            stop("pmax(S,x) not defined if length(x)>1")
        } else if (S2>0){
            stop("pmax(S,x) not defined if x>0")
        } else {
            return(spray(index(S1),pmax(value(S1),S2)))
        }
    } else {
        return(spraymaker(spray_pmax(index(S1),value(S1),index(S2),value(S2))))
    }
}

`minpair_spray` <- function(S1,S2){
    if(missing(S2)){ return(S1) }
    if(!is.spray(S2)){
        if(length(S2)>1){
            stop("pmin(S,x) not defined if length(x)>1")
        } else if (S2<0){
            stop("pmin(S,x) not defined if x<0")
        } else {
            return(spray(index(S1),pmin(value(S1),S2)))
        }
    } else {
        return(spraymaker(spray_pmin(index(S1),value(S1),index(S2),value(S2))))
    }
}

`rspray` <- function(n,vals=1,arity=3,powers=0:2){
    return(spray(matrix(sample(powers,n*arity,replace=TRUE),ncol=arity),addrepeats=TRUE,vals))
}

`knight` <- function(d=2){
  n <- d * (d - 1)
  out <- matrix(0, n, d)
  out[cbind(rep(seq_len(n), each=2), c(t(which(diag(d)==0, arr.ind=TRUE))))] <- seq_len(2)
  spray(rbind(out, -out, `[<-`(out, out==1, -1),`[<-`(out, out==2, -2)))
}

`king` <- function(d=2){
  out <- spray(expand.grid(replicate(d,list(c(-1L,0L,1L)))))
  return(out-1)
}

`nterms` <- function(S){ nrow(index(S)) }

setGeneric("zapsmall")
setMethod("zapsmall","spray",function(x,digits){
    zap(x,digits=digits)
})

setMethod("zapsmall","ANY",function(x,digits){
    base::zapsmall(x,digits=digits)
})

`zap` <- function(x, digits=getOption("digits")){
  spray(index(x),base::zapsmall(value(x),digits=digits))
}
