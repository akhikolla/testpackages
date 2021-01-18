augment.input <- function(input){
  
  if(is.list(input)){
    S = length(input)
  }else{
    stop("Input is not a list!")
  }
  
  input.union = list()
  for(t in 1:S){
    input.union[[t]] = do.call(rbind, input[t:S])
  }

  
  # get the missing inputs at each level
  input.miss = list()
  for(t in 1:(S-1)){
    input.miss[[t]] = (match.input(input.union[[t+1]], input[[t]])$A)
  }
  
  # only combine missing and observed inputs with repetition
  for(t in 1:(S-1)){
    input.union[[t]] = rbind(input[[t]], input.miss[[t]])
  }
  input.union[[S]] = input[[S]]

  return(list(union=input.union, miss=input.miss))
  
}

## for univariate models
create.w.new <- function(t, input, input.miss, y, ym){
  n = nrow(input.miss[[t]])
  w = matrix(NA, n)
  indtemp = ismember(input.miss[[t]], input.miss[[t-1]])
  w[indtemp$IIA] = ym[indtemp$IA]

  if(length(indtemp$IIA)<n){
    indtemp = ismember(input.miss[[t]], input)
    w[indtemp$IIA] = y[indtemp$IA]    
  }


  return(w)
}

## for univiariate models
create.w.univ <- function(t, input, input.miss, y, ym){

  n = nrow(input[[t]])
  w = matrix(NA, n, 1)
  indtemp = ismember(input.miss, input[[t]])
  w[indtemp$IA, 1] = ym[indtemp$IIA]

  # if(length(indtemp$IA)<n){
  #   input.temp = input[[t]][indtemp$NIB, , drop=FALSE]
  #   ID = ismember(input.temp, input[[t-1]])$IA 
  #   w[indtemp$NIB, 1] = y[ID]
  # }

  if(length(indtemp$IA)<n){
    input.temp = input[[t]][indtemp$NIA, , drop=FALSE]
    ID = ismember(input.temp, input[[t-1]])$IA 
    w[indtemp$NIA, 1] = y[ID]
  }

  return(w)
}



## This routine creates W's at missing inputs for multivaraite models
create.w <- function(t, input, input.miss, y, ym){

	n = nrow(input[[t]])
	N = ncol(y)
	w = matrix(NA, n, N)
	indtemp = ismember(input.miss, input[[t]])
	w[indtemp$IA, ] = ym[indtemp$IIA, ]

	# if(length(indtemp$IA)<n){
	# 	input.temp = input[[t]][indtemp$NIB, , drop=FALSE]
	# 	ID = ismember(input.temp, input[[t-1]])$IA 
	# 	w[indtemp$NIB, 1] = y[ID]
	# }

	if(length(indtemp$IA)<n){
		input.temp = input[[t]][indtemp$NIA, , drop=FALSE]
		ID = ismember(input.temp, input[[t-1]])$IA 
		w[indtemp$NIA, ] = y[ID, ]
	}

	return(w)
}


## This routine creates W's at new inputs for multivariate models
create.w.pred <- function(t, input, input.miss, y, ym){
  n = nrow(input.miss[[t]])
  N = ncol(y)
  w = matrix(NA, n, N)
  indtemp = ismember(input.miss[[t]], input.miss[[t-1]])
  w[indtemp$IIA, ] = ym[indtemp$IA, ]

  if(length(indtemp$IIA)<n){
    indtemp = ismember(input.miss[[t]], input)
    w[indtemp$IIA, ] = y[indtemp$IA, ]    
  }


  return(w)
}


ismember <- function(matA, matB){
	
	if(!is.matrix(matA)){
		message("\n The first argument is not a matrix!\n")
		matA = as.matrix(matA)
	}

	if(!is.matrix(matB)){
		message("\n The second argument is not a matrix!\n")
		matA = as.matrix(matB)
	}

	nB = nrow(matB)

	dfA = data.frame(t(matA))
	dfB = data.frame(t(matB))

	# find positions in matB such that all rows of matA are in matB
	ind = match(dfA, dfB)   
    indA = ind[!is.na(ind)] # matA == matB[indA, ]
    indtemp = 1:nB 
    NIA = indtemp[-indA] # matB[NIA, ] is not in matA 

    # find positions in matA such that all rows of matA are in matB
    IIA = which(!is.na(ind)) # matA[IIA, ] == matB[indA, ]
    NIB = which(is.na(ind))

	return(list(IA=indA, IIA=IIA, NIA=NIA, NIB=NIB))
	# matA[IIA, ] = matB[IA, ]
	# matB[NIA, ] does not belong to any rows of matA
	# matA[NIB, ] does not belog to any rows of matB
}




match.input = function(input1, input2){

  if(!is.matrix(input1)){
    message("\n The first argument is not a matrix!\n")
    input1 = as.matrix(input1)
  }

  if(!is.matrix(input2)){
    message("\n The second argument is not a matrix!\n")
    input2 = as.matrix(input2)
  }

  n1 = dim(input1)[1]
  n2 = dim(input2)[1]

  dfA = data.frame(t(input1))
  dfB = data.frame(t(input2))

  ind = match(dfA, dfB)
  indB = ind[!is.na(ind)] 

  indtemp = 1:n2

  indA = which(!is.na(ind))
  NIB = which(is.na(ind))

  A = input1[NIB, ,drop=FALSE]
  
  if(length(indA)==0){
    indA = NULL
  }

  if(length(indB)==0){
    indB = NULL
  }
  
  return(list(IA=indA, A=A, IB=indB))
  # IA: positions of elements from input1 that are also in input2
  # A: elements from input1 that are not in input2
  # IB: positions of elements from input2 that are also in input1
}





############################################################################
#' @title Compute continous rank probability score for normal distributions
#' @description This function compute the continous rank probability score for
#' normal distributions. It is mainly used to evaluate the validility of 
#' predictive distributions.
#' @param x a vector of true values (held-out data)
#' @param mu a vector of predictive means
#' @param sig a vector of predictive standard deviations
#'
#' @author Pulong Ma <mpulong@gmail.com>
#' @export

CRPS <- function(x, mu, sig){

  xo = (x-mu)/sig
  crps = sig*(xo*(2*pnorm(xo)-1) + 2*pnorm(xo) - 1/sqrt(pi))
  
  return(crps)
}



