#' Next dose determination of a phase I clinical trial.
#'
#' @param model A character string to specify the selected dose-finding model. See for details \code{\link{dtox}}, \code{\link{pkcov}}, \code{\link{pkcrm}}, \code{\link{pktox}}, \code{\link{pkpop}}, \code{\link{pklogit}}.
#' @param y A binary vector of the toxicity outcomes from previous patients; 1 indicates a toxicity, 0 otherwise.
#' @param AUCs A vector with the computed AUC values of each patient for pktox, pkcrm, pklogit and pkpop; defaults to NULL. 
#' @param doses A vector with the doses panel.
#' @param x A vector with the dose level assigned to the patients.
#' @param theta The toxicity threshold.
#' @param options A list with the Stan model's options.
#' @param prob The threshold of the posterior probability of toxicity for the stopping rule in the selected model; defaults to 0.9. See for details \code{\link{dtox}}, \code{\link{pkcov}}, \code{\link{pkcrm}}, \code{\link{pktox}}, \code{\link{pkpop}}, \code{\link{pklogit}}.
#' @param betapriors A vector with the value for the prior distribution of the regression parameters in the model; defaults to NULL.
#' @param thetaL A second threshold of AUC in the \code{\link{pkcrm}} model; defaults to theta in the PKCRM model and NULL for the models \code{\link{dtox}}, \code{\link{pkcov}}, \code{\link{pktox}}, \code{\link{pkpop}} and \code{\link{pklogit}}.
#' @param p0 The skeleton of CRM for \code{\link{pkcrm}}; defaults to NULL.
#' @param L The AUC threshold to be set before starting the trial for \code{\link{pkcrm}}; defaults to NULL.
#' @param deltaAUC A vector of the difference between computed individual AUC and the AUC of the population at the same dose level (defined as an average); argument for \code{\link{pkcov}}; defaults to NULL.
#' @param CI A logical constant indicating the estimated 95\% credible interval; defaults to TRUE. 
#'
#' @description
#' nextDose is used to perform parameter estimation at each step during a dose-finding trial. Determines the next or recommended dose level in a phase I clinical trial.
#'
#' @return  An object of class "dose" is returned, consisting of determination of the next recommended dose and estimations. Objects generated 
#' by nextDose contain at least the following components:
#'
#' \item{N}{The total number of enrolled patients.}
#' \item{y}{A binary vector of toxicity outcomes from previous patients; 1 indicates a toxicity, 0 otherwise.}
#' \item{AUCs}{A vector with the computed AUC values of each patient.}
#' \item{doses}{A vector with the doses panel.}
#' \item{x}{A vector with the dose level assigned to the patients.}
#' \item{theta}{The tocixity threshold.}
#' \item{options}{List with the Stan model's options.}
#' \item{newDose}{The next recommended dose (RD) level; equals to 0 if the trial has stopped, according to the stopping rules.}
#' \item{pstim}{The mean values of the estimated probabilities of toxicity.}
#' \item{pstimQ1}{The 1st quartile of the estimated probabilities of toxicity if CI = TRUE, otherwise is NULL.}
#' \item{pstimQ3}{The 3rd quartile of the estimated probabilities of toxicity if CI = TRUE, otherwise is NULL.}
#' \item{parameters}{The estimated model's parameters.}
#' \item{model}{A character string to specify the selected dose-finding model. See for details \code{\link{dtox}}, \code{\link{pkcov}}, \code{\link{pkcrm}}, \code{\link{pktox}}, \code{\link{pkpop}}, \code{\link{pklogit}}.}
#'
#' @author Artemis Toumazi \email{artemis.toumazi@@gmail.com}, Moreno Ursino \email{moreno.ursino@@inserm.fr}, Sarah Zohar \email{sarah.zohar@@inserm.fr}
#'
#' @references Ursino, M., et al, (2017) Dose-finding methods for Phase I clinical trials using pharmacokinetics in small populations, Biometrical Journal, <doi:10.1002/bimj.201600084>.
#'
#' Toumazi, A., et al, (2018) dfpk: An R-package for Bayesian dose-finding designs using pharmacokinetics (PK) for phase I clinical trials, Computer Methods and Programs in Biomedicine, <doi:10.1016/j.cmpb.2018.01.023>.
#' 
#' @examples
#'   \dontrun{
#' doses <- c(12.59972,34.65492,44.69007,60.80685,83.68946,100.37111)
#' theta <- 0.2 
#' options <- list(nchains = 4, niter = 4000, nadapt = 0.9)     
#' AUCs <- c(1.208339,  5.506040,  6.879835,  3.307928,  3.642430, 
#'           10.271291,  3.885522,  3.086622,  2.537158,  5.525917,  
#'           8.522176,  4.642741, 11.048531, 10.246976,  5.226807)
#' x <- c(1, 2, 3, 4, 5, 6, 4, 4, 4, 5, 5, 4, 4, 5, 5)
#' y <- c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0)
#' nextD <- nextDose(model = "pktox", y=y, AUCs=AUCs, doses=doses, 
#'                   x=x, theta=theta, options=options)
#' }
#'
#' @seealso \code{\link{nsim}}
#'
#' @import ggplot2
#' @import rstan
#' @useDynLib dfpk, .registration = TRUE
#' @export
nextDose <- function(model, y, AUCs, doses, x, theta, options, prob = 0.9, betapriors = NULL, 
					 thetaL=NULL, p0 = NULL, L = NULL, deltaAUC = NULL, CI = TRUE){
	model1 = NULL
	eval(parse(text = paste("model1 =", model, sep="")))
	N <- length(x)
	if (model == "pktox" & is.null(betapriors)){betapriors = c(10, 10000, 20, 10)
	}else if(model == "pkcrm" & is.null(betapriors)){betapriors = c(10, 10000)
	}else if (model == "pkpop" & is.null(betapriors)){betapriors = c(10, 10000, 10, 5)
	}else if (model == "dtox" & is.null(betapriors)){betapriors = c(0, 16.71, 0, 6.43)
	}else if(model == "pkcov" & is.null(betapriors)){betapriors = c(-14.76, 0, 3.23 + 5)
	}else if (model == "pklogit" & is.null(betapriors)){betapriors = c(10, 10000, 20, 10)}
	m <- model1(y, AUCs, d = doses, x, theta, prob = prob, betapriors = betapriors, thetaL=NULL, options = options, p0 = p0, L = L, deltaAUC = deltaAUC, CI = CI)
	MTD <- m$newDose
	if(is.na(MTD) == TRUE){
		MTD = 0
	}

	pstim <- m$pstim
	if(CI == TRUE && length(m$p_sum) != 0){
		pstim_Q1 <- m$p_sum[,2]
    	pstim_Q3 <- m$p_sum[,5]
	}else{
		pstim_Q1 <- NULL
    	pstim_Q3 <- NULL
	}

	parameters <- m$parameters
	new("dose", N = N, y = y, AUCs = AUCs, doses = doses, x = x, theta = theta, options = options, 
		newDose = MTD, pstim = pstim, pstimQ1 = pstim_Q1, pstimQ3 = pstim_Q3, parameters = parameters, 
		model = model)
}
