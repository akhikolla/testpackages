
"checkbc"<-function(BreedNames, bc){
  if(is.null(bc)){stop("Breed composition 'bc' is not defined.\n")}
  
  BreedNames <- unique(BreedNames)
  bc <- c(bc)
  storage.mode(bc)<-"double"
  
  if(!is.vector(bc)){             stop("Argument 'bc' is not a vector.\n")}
  if(any(is.na(bc))){             stop("Vector 'bc' contains NA.\n")}
  if(any(bc< -1e-07)){            stop("Vector 'bc' contains negative values.\n")}
  if(abs(sum(bc)-1) > 1e-07){     stop("The values in vector 'bc' do not sum up to 1.\n")}
  if(is.null(names(bc))){         stop("Vector 'bc' needs comonent names.\n")}
  if(!all(names(bc) %in% BreedNames)){stop("Some breeds in vector 'bc' do not appear in column Breed of phen.\n")}

  bcomp <- setNames(rep(0, length(BreedNames)), BreedNames)
  bcomp[names(bc)] <- bc
  bcomp[bcomp<0.001]   <-0.001
  bcomp            <- bcomp/sum(bcomp)
  bcomp
}
