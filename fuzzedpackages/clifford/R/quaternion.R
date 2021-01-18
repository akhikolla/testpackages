`clifford_to_quaternion` <- function(C){
  stopifnot(identical(sort(grades(C)),c(0L,2L,2L,2L)))
  out <- matrix(c(const(C), -getcoeffs(C,list(c(1,2),c(1,3),c(2,3)) )))
  class(out) <- c("onion","quaternion")
  return(out)
}

`quaternion_to_clifford` <- function(Q){
  Q <- as.numeric(Q)
  stopifnot(length(Q)==4)
  clifford(list(numeric(0),c(1,2),c(1,3),c(2,3)),c(Q[1],-Q[2:4]))
}
