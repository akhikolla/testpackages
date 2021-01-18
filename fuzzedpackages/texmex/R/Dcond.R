Dcond <-
function(x,a,b,c,d,zi,zk)
  {
    res <- a+(x^(b-1))*zi - (x^(d-1))*zk
    return(res)
  }
