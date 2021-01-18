#' @importFrom PK auc
#' @useDynLib dfpk, .registration = TRUE
#' @export
AUC.estim <-
function(t,conc,dose, method = 2){
	if(method == "1"){
		out_pk= optim(c(1,5,50), pk.estim, t=t, dose=dose, conc=conc, method = "L-BFGS-B", lower = c(0.1,0.2,1))
		dose/out_pk$par[2]
	}else if(method == "2"){
		a <- auc(conc, t, design="complete")
		a$est
	}else{
		stop("Error in AUCmethod: Invalid method for the AUC estimation. See the R documentation of the function for the possible methods")
}
}
