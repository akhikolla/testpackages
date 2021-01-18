setClassUnion("ClassNewDose", c("numeric", "logical", "NULL"))

#' An S4 class to represent a dosefinding results.
#'
#' @slot pid The patient's ID provided in the study.
#' @slot N  The total sample size per trial.
#' @slot time The sampling time points.
#' @slot doses A vector with the doses panel.
#' @slot conc The estimated concentration values for each patient at each dose.
#' @slot p0  The skeleton of CRM for \code{\link{pkcrm}}; defaults to NULL. 
#' @slot L  The AUC threshold to be set before starting the trial for \code{\link{pkcrm}}; defaults to NULL.
#' @slot nchains The number of chains for the Stan model.
#' @slot niter The number of iterations for the Stan model.
#' @slot nadapt  The number of warmup iterations for the Stan model.
#' @slot newDose  The next maximum tolerated dose (MTD) if TR=1 otherwise the percentage of MTD selection for each dose level after all trials starting from dose 0; equals to 0 if the trial has stopped before the end, according to the stopping rules.
#' @slot MTD A vector containing the next maximum tolerated doses (MTD) of each trial (TR); equals to 0 if the trial has stopped before the end, according to the stopping rules.
#' @slot MtD The final next maximum tolerated (MTD) dose after all the trials.
#' @slot theta  The toxicity target.
#' @slot doseLevels A vector of dose levels assigned to patients in the trial.
#' @slot toxicity The estimated toxicity outcome.
#' @slot AUCs A vector with the computed AUC values of each patient.
#' @slot TR The total number of trials to be simulated.
#' @slot preal The prior toxicity probabilities.
#' @slot pstim The estimated mean probabilities of toxicity.
#' @slot pstimQ1 The 1st quartile of estimated probability of toxicity.
#' @slot pstimQ3 The 3rd quartile of estimated probability of toxicity.
#' @slot model A character string to specify the selected dose-finding model. See for details \code{\link{dtox}}, \code{\link{pkcov}}, \code{\link{pkcrm}}, \code{\link{pktox}}, \code{\link{pkpop}}, \code{\link{pklogit}}..
#' @slot seed The seed of the random number generator that is used at the beginning of each trial.
#' @import methods
#' @useDynLib dfpk, .registration = TRUE
#' @export
setClass("dosefinding", slots = list(pid="numeric", N ="numeric", time="numeric", doses = "numeric", conc="numeric", 
        p0 = "numeric", L = "numeric",  nchains = "numeric", niter = "numeric", nadapt = "numeric", newDose = "ClassNewDose", 
        MTD = "ClassNewDose", MtD = "numeric", theta = "numeric", doseLevels="matrix", toxicity= "matrix", AUCs="matrix", TR="numeric", 
        preal = "numeric", pstim = "list", pstimQ1 = "list", pstimQ3 = "list", model = "character", seed= "matrix"))


#' An S4 class to represent a simulated scenarios.
#'
#' @slot PKparameters Subject's pharmacokinetic's (PK) parameters from the population distributions defined by the population mean.
#' @slot nPK The length of the time points.
#' @slot time The sampling time points.
#' @slot idtr The id number of the corresponding simulated dataset.
#' @slot N  The total sample size per trial.
#' @slot doses A vector with the doses panel. 
#' @slot preal The prior toxicity probabilities.
#' @slot limitTox The toxicity threshold. 
#' @slot omegaIIV  The inter-individual variability for the clearance and the volume of distribution.
#' @slot omegaAlpha The patient's sensitivity parameter.
#' @slot conc The concentration computed at the PK population values.
#' @slot concPred The concentration values with proportional errors for each patient at each dose. 
#' @slot tox The toxicity outcome.
#' @slot tab A summary matrix containing the sampling time points at the first row followed by concPred, parameters and alphaAUC. It used by the simulation function nsim.
#' @slot parameters The simulated PK parameters of each patient.
#' @slot alphaAUC A vector with the computed AUC values of each patient. 
#' @import methods
#' @useDynLib dfpk, .registration = TRUE
#' @export
setClass("scen", slots = list(PKparameters="numeric", nPK="numeric", time="numeric", idtr="numeric",
        N = "numeric", doses="numeric", preal = "numeric", limitTox="numeric", omegaIIV="numeric", 
        omegaAlpha="numeric", conc="matrix", concPred="numeric",
        tox="matrix", tab="matrix", parameters = "matrix", alphaAUC="numeric"))


#' An S4 class to perform parameter estimation at each step during a dose-finding trial.
#'
#' @slot N The total number of enrolled patients.
#' @slot y A binary vector of toxicity outcomes from previous patients; 1 indicates a toxicity, 0 otherwise.
#' @slot AUCs A vector with the computed AUC values of each patient.
#' @slot doses A vector with the doses panel.
#' @slot x A vector with the dose level assigned to the patients.
#' @slot theta The toxicity threshold.
#' @slot options a list of Stan model's options.
#' @slot newDose The next recommended dose (RD) level; equals to 0 if the trial has stopped, according to the stopping rules.
#' @slot pstim The estimated mean probabilities of toxicity.
#' @slot pstimQ1 The 1st quartile of estimated probability of toxicity.
#' @slot pstimQ3 The 3rd quartile of estimated probability of toxicity.
#' @slot parameters The Stan model's estimated parameters.
#' @slot model A character string to specify the selected dose-finding model. See for details \code{\link{dtox}}, \code{\link{pkcov}}, \code{\link{pkcrm}}, \code{\link{pktox}}, \code{\link{pkpop}}, \code{\link{pklogit}}.
#' @import methods
#' @useDynLib dfpk, .registration = TRUE
#' @export
setClass("dose", slots = list(N = "numeric", y = "numeric", AUCs = "numeric", doses ="numeric", x = "numeric", 
        theta = "numeric", options = "list", newDose="ClassNewDose", pstim="numeric", pstimQ1="ClassNewDose", pstimQ3="ClassNewDose", 
        parameters="numeric", model = "character"))


setGeneric("show")
#' @export 
setMethod(f = "show", signature = "dosefinding", definition = function(object)
    {
        cat("Today: ", date(), "\n") 
        cat("\n","A. Data Summary (", object@model,"model)", "\n")
        cat("Number of simulations:", object@TR, "\n")
        cat("Total number of patients in the trial:", object@N, "\n")
        cat("The time sampling:", round(object@time, digits = 3), "\n")
        n <- object@N
        cat("Levels of Doses:", round(object@doses, digits=3), "\n")
        cat("Concentration of the drug:", round(object@conc, digits = 3), "\n")
        if(object@model == "pkcrm"){
        	cat("Skeleton of CRM:", object@p0, "\n")
        	cat("Threshold set before starting the trial:", object@L, "\n")
        }
        cat("\n","B. STAN Model's Options \n")
        cat("The Stan model runs with", object@nchains, "MCMC chains ")
        cat("which each chain has", object@niter, "iterations ")
        cat("and", object@nadapt, "warmup iterations \n")
        #cat("The seed of the random number generator that is used at the beginning of each trial is given below: \n")
        #s <- matrix(NA, nrow = 1, ncol = object@TR + 1)
        #colnames(s) <- c("Initial", paste("Trial ", 1:object@TR, sep = ""))
        #rownames(s) <- "Seed"
        #s[1, ] <- c(object@seed[1]-1, object@seed)
        #print(s)
        if(object@TR == "1"){
            cat("\n","C. Dose-Finding Results: \n")
            cat("PID", "\t", "Level", "\t", "Toxicity", "\t", "AUCs", "\n")
            for (i in 1:n){
                cat(i,"\t", object@doseLevels[i],"\t", object@toxicity[i],"\t", "\t", round(object@AUCs[i], digits=3) ,"\n")
            }
            cat("\nThe prior probabilities of toxicity:", round(object@preal, digits=4), "\n")
            cat("Next recommended dose level:", object@newDose, "\n")
            cat("Recommendation is based on a target toxicity probability of:",object@theta, "\n")
        }else{
            cat("\n","C. Dose-Finding Results: \n")
            doselevels <- as.vector(object@doseLevels)
            t <- matrix(NA, nrow=4, ncol=length(object@doses)+1)
            rownames(t) <- c("Dose", "Truth Probabilities", "Dose-Allocation (%)", "Selected % MTD")
            colnames(t) <- rep("", length(object@doses)+1)
            stop <- as.character("STOP")
            dd <- paste("", 1:length(object@doses), sep = "")
            t[1, ] <- c(stop, dd)
            t[2, ] <- c("NA", round(object@preal, digits=3))
            for(i in 1:length(object@doses)){
              n_levels = length(which(doselevels == i))
              t[3, i+1] <- round(n_levels/length(doselevels), digits=2)
            }
            zeroDose <- length(which(doselevels == "NA"))
            if(zeroDose == "0"){
              t[3,1] <- "NA"
            }else{
              t[3,1] <- zeroDose / length(doselevels)
            }
            t[4, ] <- round(object@newDose, digits=2)
            print(t, quote = FALSE)
            cat("Recommendation is based on a target toxicity probability of:",object@theta, "\n")
            #cat("Seed of each trial:",object@seed, "\n")
        }
    }
)


setGeneric("show")
#' @export
setMethod(f = "show",
          signature = "scen",
          definition = function(object){
               cat("Today: ", date(), "\n")
               cat("\n","Scenarios Settings", "(TR:", object@idtr, ")", "\n","\n")
               cat("Total number of patients in the trial:", "\t", object@N, "\n")
               cat("The subject's PK parameters (ka, CL, V):", object@PKparameters)
               cat(" with a standard deviation of CL and V equals to:", object@omegaIIV, "\n")
               cat("The doses of the drug:", round(object@doses, digits=3), "\n")
               cat("The real probabilities of toxicity:", round(object@preal,digits=3), "\n")
               cat("Time after the drug administration (hours):", "\n", round(object@time, digits = 3), "\n")
               cat("Threshold on the toxicity:", object@limitTox, "\n")
               cat("\n", "Simulated Scenarios \n")
               cat("\n", "Amount of drug in a given volume of plasma is (i.e. concentration):", "\n")
               print(round(object@conc, digits = 3))
               cat("\n")
               cat("\n", "AUC with the sensitivity parameter", "\n")
               print(round(object@alphaAUC, digits = 3))
               cat("\n")
               cat("\n", "with standard deviation omega for parameter alpha=", object@omegaAlpha, "\n")
               cat("\n", "Toxicity (0 indicates no toxicity and 1 toxicity) :","\n")
               print(object@tox)
               cat("\n")
               cat("The PK parameter's estimations for each patient are:", "\n")
               print(round(object@parameters, digits=3))
               cat("\n","NB. pid = Patient's ID number \n")
        }
)


setGeneric("show")
#' @export 
setMethod(f = "show",
          signature ="dose",
          definition = function(object) {
               cat("Today: ", date(), "\n")
               cat("Model:", object@model, "\n")
               cat("Total number of enrolled patients in the trial: ", object@N, "\n")
               cat("Levels of doses: ", object@doses, "\n")
               cat("The Next Recommended Dose: ", "\n")
               print(object@newDose)
               cat("\n")
               cat("Estimated probability of toxicity: ", "\n")
               print(round(object@pstim, digits=4))
               cat("\n")
               cat("Estimated model's parameters: ", "\n")
               print(object@parameters)
               cat("\n")
        }
)

###########################################
################ Plots ####################
###########################################

setGeneric("plot")
#' The graphical representation of dose-finding results. 
#' 
#' @param x a "dosefinding" object.
#' @param y the "y" argument is not used in the plot-method for "dosefinding" object.
#' @param TR The number of the selected trial that user wants to plot; defaults to 1.
#' @param ask Choose plot or not; defaults to TRUE.
#' @param CI Indicate if the "dosefinding" object includes the 95\% credible interval for the posterior dose response plot; defaults to TRUE.
#' @param \dots other arguments to the \code{\link[=graphics]{plot.default}} function can be passed here.
#'
#' @description A plot selection showing either the dose escalation allocation of the selected trial or the plot of the final posterior distributions of the probability of toxicity at each dose or the boxplot of the sampling distribution of the probability of toxicity at each dose in the end of the trial over the total number of trials. 
#' @author Artemis Toumazi \email{artemis.toumazi@@artemis.com}, Moreno Ursino \email{moreno.ursino@@inserm.fr}, Sarah Zohar \email{sarah.zohar@@inserm.fr}
#' 
#' @references Ursino, M., et al, (2017) Dose-finding methods for Phase I clinical trials using pharmacokinetics in small populations, Biometrical Journal, <doi:10.1002/bimj.201600084>.
#'
#' Toumazi, A., et al, (2018) dfpk: An R-package for Bayesian dose-finding designs using pharmacokinetics (PK) for phase I clinical trials, Computer Methods and Programs in Biomedicine, <doi:10.1016/j.cmpb.2018.01.023>.
#' 
#' @import methods
#' @import stats
#' @import graphics
#' @importFrom grDevices  rainbow
#' @useDynLib dfpk, .registration = TRUE
#' @importFrom utils menu
#' @export
setMethod(f = "plot", signature =c("dosefinding", "missing"), definition = function(x, y=NA, TR=1, ask=TRUE, CI = TRUE,...){
    if(CI == "TRUE" && length(x@pstimQ1) != 0){
      choices <- c("1: Plot trial summary", "2: Boxplot of sampling dose response", "3: Plot posterior dose response with 95% CI\n")
    }else{
      choices <- c("1: Plot trial summary", "2: Boxplot of sampling dose response\n")
    }

    if (ask == "TRUE") {
        cat("Make a plot selection (or 0 to exit)\n\n")
        for (i in 1:length(choices)) {
            cat(choices[i], "\n")
        }
        pick <- readline("Selection: ")
    } else {
        pick <- 2
    }
    num_choices <- c(0:length(choices))
    while(!(pick %in% num_choices)) {
        cat("Invalid choice. Enter an item from the menu, or 0 to exit")
        pick <- readline("Selection: ")
    }
    if (pick == 0) {
        return(invisible(x));
    } else if (pick == 1) {
        par(las=1)
        n <- x@N                    
        nontox <- which(x@toxicity[TR,] == "0")
        notNa <- which(is.na(x@doseLevels[TR,]) == "FALSE")
        if(x@MtD == 0) warning("Plot not completed! The trial has stopped according to the stopping rules! \n \n", call. = FALSE)
        plot(x@pid[nontox], x@doseLevels[TR,nontox], pch="O", ylim=c(1,max(x@doseLevels[TR,notNa])), xlim=c(1,n), 
        	 xlab="Patient number", ylab="Dose level", ...)
        points((1:length(x@toxicity[TR,]))[-nontox],x@doseLevels[TR,-nontox], pch="X")
        mtext("Each point represents a patient", line=2)
        mtext("A circle indicates no toxicity, a cross toxicity", line=0.5)
    } else if (pick == 3) {
        par(las=1)
        n <- x@N
        ndoses <- length(x@doses)
        # PropTox <- matrix(NA, ncol = 6, nrow = ndoses)
        # for(i in 1: ndoses){
        #    PropTox[i,] <-  rbind(summary(x@pstim_post[i,]))
        # }
        if (x@MtD == 0) stop("Unable to plot! The trial stopped based on the stopping rules \n \n", call. = FALSE)
        if (length(x@pstimQ1) == 0 && length(x@pstimQ3) == 0) stop("Error in plot.window(...) : need finite 'CI' values. \n In addition: Warning messages: \n The nsim and plot function must use the same logical value for CI. \n \n", call. = FALSE)
        plot(1:ndoses, x@pstim[[TR]][1:ndoses,n],type="l",xlab="Dose level",ylab="Probability of toxicity", ylim=c(0,max(x@pstim[[TR]][1:ndoses,n]) + 0.15))
        points(1:ndoses,x@pstim[[TR]][1:ndoses,n], pch="X")
        lines(1:ndoses,x@preal, lty=2)
        points(1:ndoses,x@preal, pch="O")
        abline(h=x@theta, lwd=2, lty=4, col = "red")
        lines(1:ndoses,x@pstimQ1[[TR]][1:ndoses,n], lty=3, col = "blue")
        lines(1:ndoses,x@pstimQ3[[TR]][1:ndoses,n], lty=3, col = "blue")
        # lines(1:ndoses, PropTox[,2], lwd=2, lty=3, col = "orange")
        # lines(1:ndoses, PropTox[,5], lwd=2, lty=3, col = "orange") 
        mtext("Prior (dashed) and updated (solid) dose-toxicity curves", line=2)
        mtext("95% CI (dotted) of the updated dose-toxicity curve", line=0.5)
    } else if (pick == 2){
        par(las=1)
        ndoses <- length(x@doses)
        if (x@MtD == 0) stop("Unable to plot! The trial stopped based on the stopping rules \n \n", call. = FALSE)
        PropTox <- matrix(NA, ncol = 6, nrow = ndoses)
        for(i in 1: ndoses){
            #PropTox[i,] <-  rbind(summary(x@pstim[i,]))
            PropTox[i,] <- rbind(summary(unlist(lapply(x@pstim, `[`,i,))))
        }
        d <- c(1:ndoses)
        boxplot(PropTox~d, xlab="Dose level", ylab="Probability of toxicity", ylim=c(0,max(PropTox) + 0.15))
        abline(h=x@theta, lwd=2, lty=4, col = "red")
        mtext("Boxplot dose response", line=1)
        # if (x@newDose == "NA") mtext("(Note: The trial stopped based on the stopping rules)", line=0)
    }
})

#' The graphical representation of the drug's concentration in the plasma at time t after the drug administration. 
#'  
#' @param x a "scen" object or a list of the selected trial from a "scen" object.
#' @param y the "y" argument is not used in the plot-method for "scen" object.
#' @param col the color argument to the \code{\link[=graphics]{plot.default}} function.
#' @param xlab the label of x-axis.
#' @param ylab the label of y-axis.
#' @param main the title of the graph.
#' @param \dots other arguments to the \code{\link[=graphics]{plot.default}} function can be passed here.
#'
#' @author Artemis Toumazi \email{artemis.toumazi@@gmail.com}, Moreno Ursino \email{moreno.ursino@@inserm.fr}, Sarah Zohar \email{sarah.zohar@@inserm.fr}
#'
#' @references Ursino, M., et al, (2017) Dose-finding methods for Phase I clinical trials using pharmacokinetics in small populations, Biometrical Journal, <doi:10.1002/bimj.201600084>.
#'
#' Toumazi, A., et al, (2018) dfpk: An R-package for Bayesian dose-finding designs using pharmacokinetics (PK) for phase I clinical trials, Computer Methods and Programs in Biomedicine, <doi:10.1016/j.cmpb.2018.01.023>.
#' 
#' @import methods
#' @import stats
#' @import graphics
#' @importFrom grDevices  rainbow
#' @useDynLib dfpk, .registration = TRUE
#' @export
setMethod(f = "plot", signature =c("scen", "missing"), definition = function(x, y=NA, col = rainbow(length(x@doses)), xlab="Time (hours)", 
          ylab="Concentration (mg/L)", main="Pharmacokinetics: Concentration vs Time", ...){
    plot(x@time, x@conc[,1], type="l", col=col[1], xlab=xlab, ylab=ylab, main=main, ylim=c(0,max(x@conc)), ...)
    for(i in 2:length(x@doses)){
        lines(x@time, x@conc[,i], col=col[i], lty=i, ...)
    }
}
)


#' The graphical representation of dose escalation for each patient in the trial. 
#'
#' @param x a "dose" object.
#' @param y the "y" argument is not used in the plot-method for "dose" object.
#' @param ask Choose plot or not; defaults to TRUE.
#' @param CI Indicate if the "dose" object includes the 95\% credible interval for the posterior dose response plot; defaults to TRUE.
#' @param \dots other arguments to the \code{\link[=graphics]{plot.default}} function can be passed here.
#'
#' @author Artemis Toumazi \email{artemis.toumazi@@gmail.com}, Moreno Ursino \email{moreno.ursino@@inserm.fr}, Sarah Zohar \email{sarah.zohar@@inserm.fr}
#'
#' @references Ursino, M., et al, (2017) Dose-finding methods for Phase I clinical trials using pharmacokinetics in small populations, Biometrical Journal, <doi:10.1002/bimj.201600084>.
#'
#' Toumazi, A., et al, (2018) dfpk: An R-package for Bayesian dose-finding designs using pharmacokinetics (PK) for phase I clinical trials, Computer Methods and Programs in Biomedicine, <doi:10.1016/j.cmpb.2018.01.023>.
#' 
#' @import methods
#' @import stats
#' @import graphics 
#' @importFrom grDevices  rainbow
#' @useDynLib dfpk, .registration = TRUE
#' @export
setMethod(f = "plot", signature =c("dose", "missing"), definition = function(x, y=NA, ask=TRUE, CI = TRUE, ...){
    if(CI == "TRUE" && length(x@pstimQ1) != 0){
      choices <- c("1: Plot trial summary", "2: Plot posterior dose response with 95% CI\n")
    }else{
      choices <- c("1: Plot trial summary\n")
    }
    if (ask == "TRUE") {
        cat("Make a plot selection (or 0 to exit)\n\n")
        for (i in 1:length(choices)) {
            cat(choices[i], "\n")
        }
        pick <- readline("Selection: ")
    } else {
        pick <- 1
    }
    num_choices <- c(0:length(choices))
    while(!(pick %in% num_choices)) {
        cat("Invalid choice. Enter an item from the menu, or 0 to exit")
        pick <- readline("Selection: ")
    }
    if (pick == 0) {
        return(invisible(x));
    } else if (pick == 1) {
        par(las=1)
		n <- x@N
		pid <- 1:n              
		nontox <- which(x@y == "0")
		plot(pid[nontox], x@x[nontox], pch="O", ylim=c(1,max(x@x)), xlim=c(1,n+1), xlab="Patient number", ylab="Dose level",...)
		points((1:length(x@y))[-nontox],x@x[-nontox], pch="X")
		mtext("Each point represents a patient", line=2)
		mtext("A circle indicates no toxicity, a cross toxicity", line = 0.5)
    } else {
        par(las=1)
        n <- x@N
        ndoses <- length(x@doses)
        plot(1:ndoses, x@pstim, type="l", xlab="Dose level", ylab="Probability of toxicity", ylim=c(0,max(x@pstim) + 0.15))
        points(1:ndoses,x@pstim, pch="X")
        abline(h=x@theta, lwd=2, lty=4, col = "red")
        lines(1:ndoses,x@pstimQ1, lty=3, col = "blue")
        lines(1:ndoses,x@pstimQ3, lty=3, col = "blue")
        mtext("Updated (solid) dose-toxicity curves with the", line=1.5)
        mtext("95% CI (dotted) of the updated dose-toxicity curve", line=0.5)
    }
})
