`freealg` <- function(words,coeffs){ # formal
  stopifnot(is_ok_free(words,coeffs))
  out <- lowlevel_simplify(words,coeffs)  # simplify() is defined in
                                 # RcppExports.R; it returns a list

  class(out) <- "freealg"   # this is the only time class() is set to "freealg"
  return(out)
}

`words` <- function(x){x[[1]]}
`coeffs` <- function(x){x[[2]]}  # accessor methods end here

`coeffs<-` <- function(x,value){UseMethod("coeffs<-")}
`coeffs<-.freealg` <- function(x,value){
    if(length(value) != 1){
        stop('order of coefficients not defined.  Idiom "coeffs(x) <- value" is meaningful only if value is unchanged on reordering, here we require "value" to have length 1') 
    }
    x[[2]][] <- value
    return(freealg(x[[1]],x[[2]]))
}

`as.freealg` <- function(x,...){
  if(is.freealg(x)){
    return(x)
  } else if(is.list(x)){
    return(freealg(x[[1]],x[[2]]))
  } else if(is.numeric(x) &&(length(x)==1)){
    return(numeric_to_free(x))
  } else if(is.numeric(x) &&(length(x) > 1)){
    return(vector_to_free(x))
  } else if(is.character(x)){
    return(natural_char_to_freealg(x))
  } else {
    stop("not recognised")
  }
}

`numeric_to_free` <- function(x){
  stopifnot(length(x)==1)
  freealg(list(numeric(0)),x)
  }

`is.zero` <- function(x){  length(coeffs(x))==0 }

`is.constant` <- function(x){
  jj <- words(x)
  (length(jj)==1) & identical(jj[[1]],integer(0))
}

"constant" <- function(x){UseMethod("constant")}
"constant<-" <- function(x, value){UseMethod("constant<-")}

`constant.numeric` <- function(x){numeric_to_free(x)}

`constant.freealg` <- function(x){
  wanted <- sapply(words(x),function(x){length(x)==0})
  if(any(wanted)){
    out <- coeffs(x)[wanted]
  } else {
    out <- 0
  }
  return(out)
}

`constant<-.freealg` <- function(x,value){
  wanted <- sapply(words(x),function(x){length(x)==0})
  if(any(wanted)){
    co <- coeffs(x)
    co[wanted] <- value
    w <- words(x)
    } else {
      co <- c(coeffs(x),value)
      w <- c(words(x),list(numeric(0)))
    }
  freealg(w,co)
}

`is.freealg` <- function(x){inherits(x,"freealg")}

`is_ok_free` <- function(words,coeffs){
    if( (length(words)==0) && (length(coeffs)==0)){
        return(TRUE)  # zero element
    }

    if(identical(words,list()) && length(coeffs)==1){
      return(TRUE)
      }
    stopifnot(unlist(lapply(words,is.numeric)))
    stopifnot(is.numeric(coeffs))

    stopifnot(length(words)==length(coeffs))
    return(TRUE)
}

`rfalg` <- function(n=7, distinct=3, maxsize=4, include.negative=FALSE){
  distinct <- seq_len(distinct)
  if(include.negative){distinct <- c(distinct,-distinct)}
  freealg(replicate(n,sample(distinct,min(1+rgeom(1,1/maxsize),maxsize),replace=TRUE),simplify=FALSE), seq_len(n))
}

`print.freealg` <- function(x,...){
  cat("free algebra element algebraically equal to\n")

  if(isTRUE(getOption("usecaret"))){
      symbols <- c(letters,paste(letters,"^-1",sep=""))
  } else {
      symbols <- c(letters,LETTERS)
  }
  n <- 26
  
  out <- ""
  for(i in seq_along(words(x))){
    co <- coeffs(x)[i]
    if(co>0){
      pm <- " + " # pm = plus or minus
    } else {
      pm <- " - "
      co <- abs(co)
    }
    jj <- words(x)[i][[1]]
    if(length(jj)>0){mulsym <- "*"} else {mulsym <- ""}
    if(any(jj<0)){jj[jj<0] <- n-jj[jj<0]}
    jj <- symbols[jj]
    jj <- paste(jj,collapse="")

    out <- paste(out, pm, co, mulsym, jj, sep="")
  }
  if(is.zero(x)){out <- "0"}
  cat(out)
  cat("\n")
  return(x)
}

`vector_to_free` <- function(v,coeffs){
  if(missing(coeffs)){coeffs <- rep(1,length(v))}
  freealg(as.list(v),coeffs)
}

`string_to_freealg` <- function(string){
  string <- gsub("^\\+","",string)  # strip initial "+"
  string <- gsub("\\*","",string)   # strip all "*"
  minus <- length(grep("^-",string))>0
  if(minus){
    sign <- (-1)
    string <- gsub("^-","",string)
  } else {
    sign <- +1
  }

  if(length(grep("[0-9]",string))>0){
    coeff <- as.numeric(gsub("[a-z]|[A-Z]","",string))
    string <- gsub("[0-9]","",string)
  } else {
    coeff <- 1
  }
  out <- match(strsplit(string,"")[[1]], c(letters,LETTERS))
  if(any(out>26)){out[out>26] <- 26-out[out>26]}
  freealg(list(out),coeffs=sign*coeff)
}

`char_to_freealg` <- function(ch){ Reduce(`+`,lapply(ch,string_to_freealg))  }

`natural_char_to_freealg` <- function(string){
  string <- paste(string, collapse = " ")
  string <- gsub(" ","",string)  # strip spaces
  string <- gsub("\\+"," +",string) # 'A+B" -> "A +B"
  string <- gsub("\\-", " -",string) # "A-B" -> "A -B"
  char_to_freealg(strsplit(string," ")[[1]]) }
