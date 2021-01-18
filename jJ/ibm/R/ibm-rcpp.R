local.int.cpp <- function(posx, posy, R1, R2){
  .Call("local_int_cpp", posx, posy, R1, R2, PACKAGE = "ibm" )
}





