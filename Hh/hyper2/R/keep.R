         
`tidy` <- function(H){
  wanted <- sort(unique(c(brackets(H),recursive=TRUE)))  # numeric

  o <- order(wanted)  # == seq_along(wanted)

  bracketout <- list()
  powerout <- NULL
  for(i in seq_along(H)){
    b <- brackets(H)[[i]]  # numeric (not necessarily sorted)
    if(any(b %in% wanted)){
      bracketout <- c(bracketout, list(which(apply(outer(b,wanted,"=="),2,any)))) # easily the ugliest line in any of my code, anywhere
      powerout <- c(powerout, powers(H)[i])
    }
  }
  if(identical(pnames(H),NA)){
    pout <- NA
  } else {
    pout <- pnames(H)[wanted]
  }
  return(hyper2(bracketout,powerout,pout))
}

`keep` <- function(H, keep, tidy=TRUE){
  p <- pnames(H)    # might be NA
  if(is.character(keep)){
    stopifnot(all(keep %in% p))
    keep <- which(p %in% keep) # 'keep' now numeric
  } else {
    jj <- seq_len(size(H))
    stopifnot(all(keep %in% jj))
    keep <- which(jj %in% keep) # 'keep' now numeric
  }

  bracketout <- list()
  powerout <- NULL
  for(i in seq_along(H)){
    b <- brackets(H)[[i]]
    jj <- b[b %in% keep]   # the meat
    if(length(jj)>0){
      bracketout <- c(bracketout, list(jj))
      powerout <- c(powerout, powers(H)[i])
    }
  }
  out <-hyper2(L=bracketout,d=powerout,pnames=p)
  if(tidy){out <- tidy(out)}
  return(out)
}

`discard` <- function(H, discard, tidy=TRUE){
  p <- pnames(H)
  if(is.character(discard)){
    stopifnot(all(discard %in% p))
    keep <- which(!(p %in% discard))
  } else {
    jj <- seq_len(size(H))
    stopifnot(all(discard %in% jj))
    keep <- which(!(jj %in% discard))
  }
  return(keep(H,keep,tidy=tidy))  # the meat
}
