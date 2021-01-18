## Cpp code: Generate Rcpp code for mparse
## (based on idea in SimInf_mparse in SimInf package)
Rcpp_mparse <- function(transitions, matchCrit, obsProcess, addVars, stopCrit, tspan, incidence, afterTstar, runFromR) {
    
    ## read source code
    Rcpp_code <- readLines(system.file("", "simFunction.R", package = "SimBIID"))
    
    ## set compilation type
    if(runFromR) {
        compType <- "Rcpp_object <- Rcpp::cppFunction"
        ## set output type
        if(is.null(matchCrit)) {
            if(!tspan){
                compType <- paste0(compType, "('NumericVector ")
            } else {
                if(!runFromR) {
                    compType <- paste0(compType, "('NumericMatrix ")
                } else {
                    compType <- paste0(compType, "('List ")
                }
            }
        } else {
            compType <- paste0(compType, "('NumericVector ")
        }
    } else {
        compType <- "Rcpp_object <- RcppXPtrUtils::cppXPtr('SEXP "
    }
    
    ## set matching critera
    if(!is.null(matchCrit)) {
        Rcpp_code[1] <- paste0(compType, "simFunction(NumericVector gdata, double tstart, double tstop, IntegerVector u, IntegerVector counts")
    } else {
        Rcpp_code[1] <- paste0(compType, "simFunction(NumericVector gdata, double tstart, double tstop, IntegerVector u")
    }
    if(tspan){
        Rcpp_code[1] <- paste0(Rcpp_code[1], ", NumericVector tspan")
    }
    
    ## add additional variables to parser
    if(!is.null(addVars)) {
        Rcpp_code[1] <- paste(Rcpp_code[1], addVars)
    }
    Rcpp_code[1] <- paste0(Rcpp_code[1], ") { ")
    
    ## extract rate markers
    ratelines <- sort(c(grep("RATELINES", Rcpp_code), grep("MATCHCRIT", Rcpp_code), 
                        grep("TSPAN", Rcpp_code), grep("AFTER_TSTAR", Rcpp_code),
                        grep("RETURNCRIT", Rcpp_code)))
    if(length(ratelines) != 10){
        stop("Something wrong with simFunction code file")
    }
    
    ## number of transitions
    nrates <- length(transitions)
    currline <- ratelines[1]
    lines <- paste0("    NumericVector rates(", nrates, ");")
    Rcpp_code <- c(Rcpp_code[1:(currline - 1)], lines, Rcpp_code[(currline + 1):length(Rcpp_code)])
    currline <- currline + length(lines) - 1
    ratelines <- ratelines[-1] + length(lines) - 1
    
    ## update rates
    upRates <- c("", "    // update rates")
    for(i in 1:length(transitions)) {
        temp <- transitions[[i]]$propensity
        temp <- paste0("    rates[", i - 1, "] = ", temp, ";")
        upRates <- c(upRates, temp)
    }
    upRates <- c(upRates, "    totrate = sum(rates);")
    
    ## reset incidence if using particle filter
    if(incidence & !is.null(matchCrit)) {
        resetInc <- c("", paste0(paste(rep(" ", 4), collapse = ""), "// reset incidence values if new time point"))
        resetInc <- c(resetInc, paste0(paste(rep(" ", 4), collapse = ""), "if(tstart > 0.0) {"))
        resetInc <- c(resetInc, paste0(paste(rep(" ", 8), collapse = ""), "for(i = (int) (u.size() / 2); i < u.size(); i++) {"))
        resetInc <- c(resetInc, paste0(paste(rep(" ", 12), collapse = ""), "u[i] = 0;"))
        resetInc <- c(resetInc, paste0(paste(rep(" ", 8), collapse = ""), "}"))
        resetInc <- c(resetInc, paste0(paste(rep(" ", 4), collapse = ""), "}"))
        nupRates <- length(upRates)
        upRates <- c(upRates, resetInc)
    }
    
    Rcpp_code <- c(Rcpp_code[1:currline], upRates, Rcpp_code[(currline + 1):length(Rcpp_code)])
    ratelines <- ratelines + length(upRates)
    
    if(incidence & !is.null(matchCrit)) {
        upRates <- upRates[1:nupRates]
    }
    
    ## update output size
    if(is.null(matchCrit)){
        ## initialise tspan outputs
        if(tspan){
            if(!is.data.frame(obsProcess)){
                outsize <- paste0(paste(rep(" ", 4), collapse = ""), "NumericMatrix out(tspan.size() + 1, u.size() + 2);")
            } else {
                outsize <- paste0(paste(rep(" ", 4), collapse = ""), "NumericMatrix out(tspan.size() + 1, u.size() + ", nrow(obsProcess) + 2, ");")
            }
        } else {
            if(!is.data.frame(obsProcess)){
                outsize <- paste0(paste(rep(" ", 4), collapse = ""), "NumericVector out(u.size() + 2);")
            } else {
                outsize <- paste0(paste(rep(" ", 4), collapse = ""), "NumericVector out(u.size() + ", nrow(obsProcess) + 2, ");")
            }
        }
    } else {
        outsize <- paste0(paste(rep(" ", 4), collapse = ""), "NumericVector out(u.size() + 1);")
    }
    currline <- ratelines[1]
    Rcpp_code <- c(Rcpp_code[1:(currline - 1)], outsize, Rcpp_code[(currline + 1):length(Rcpp_code)])
    ratelines <- ratelines[-1]
    
    ## update tspan
    if(tspan) {
        tempnSpace <- 4
        tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        
        ## check tspan
        checkTspan <- paste0("")
        checkTspan <- c(checkTspan, paste0(tempSpace, "// check tspan"))
        checkTspan <- c(checkTspan, paste0(tempSpace, "for(i = 1; i < tspan.size(); i++){"))
        tempSpace <- paste0(rep(" ", tempnSpace + 4), collapse = "")
        checkTspan <- c(checkTspan, paste0(tempSpace, "if(tspan[i] <= tspan[i - 1]) {"))
        tempSpace <- paste0(rep(" ", tempnSpace + 8), collapse = "")
        checkTspan <- c(checkTspan, paste0(tempSpace, "stop(\"tspan not ordered\");"))
        tempSpace <- paste0(rep(" ", tempnSpace + 4), collapse = "")
        checkTspan <- c(checkTspan, paste0(tempSpace, "}"))
        checkTspan <- c(checkTspan, paste0(tempSpace, "if(tstart >= tspan[0]){"))
        tempSpace <- paste0(rep(" ", tempnSpace + 8), collapse = "")
        checkTspan <- c(checkTspan, paste0(tempSpace, "stop(\"tstart >= tspan[0]\");"))
        tempSpace <- paste0(rep(" ", tempnSpace + 4), collapse = "")
        checkTspan <- c(checkTspan, paste0(tempSpace, "}"))
        checkTspan <- c(checkTspan, paste0(tempSpace, "if(tstop < tspan[tspan.size() - 1]){"))
        tempSpace <- paste0(rep(" ", tempnSpace + 8), collapse = "")
        checkTspan <- c(checkTspan, paste0(tempSpace, "stop(\"tstop < tspan[n]\");"))
        tempSpace <- paste0(rep(" ", tempnSpace + 4), collapse = "")
        checkTspan <- c(checkTspan, paste0(tempSpace, "}"))
        tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        checkTspan <- c(checkTspan, paste0(tempSpace, "}"))
        
        ## update tspan
        upTspan <- paste0(tempSpace, "while(tspan[k] < t && k < tspan.size()) {")
        tempSpace <- paste0(rep(" ", tempnSpace + 4), collapse = "")
        upTspan <- c(upTspan, paste0(tempSpace, "out(k, 0) = NA_REAL;"))
        upTspan <- c(upTspan, paste0(tempSpace, "out(k, 1) = tspan[k];"))
        upTspan <- c(upTspan, paste0(tempSpace, "for(i = 0; i < u.size(); i++) {"))
        upTspan <- c(upTspan, paste0(tempSpace, "    out(k, i + 2) = (double) u[i];"))
        upTspan <- c(upTspan, paste0(tempSpace, "}"))
        if(incidence) {
            upTspan <- c(upTspan, paste0(tempSpace, "for(i = u.size() / 2; i < u.size(); i++) {"))
            upTspan <- c(upTspan, paste0(tempSpace, "    u[i] = 0;"))
            upTspan <- c(upTspan, paste0(tempSpace, "}"))
        }
        if(is.data.frame(obsProcess)){
            for(i in 1:nrow(obsProcess)){
                upTspan <- c(upTspan, paste0(tempSpace, "out(k, u.size() + ", i - 1 + 2, ") = ", obsProcess$compiled[i]))
            }
        }
        upTspan <- c(upTspan, paste0(tempSpace, "k++;"))
        tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        upTspan <- c(upTspan, paste0(tempSpace, "}"))
        
        currline <- ratelines[1]
        Rcpp_code <- c(
            Rcpp_code[1:(currline - 1)], 
            checkTspan,
            "",
            paste0(tempSpace, "// update tspan"), 
            paste0(tempSpace, "int k = 0;"), 
            upTspan, 
            Rcpp_code[(currline + 1):length(Rcpp_code)]
        )
        ratelines <- ratelines[-1] + length(checkTspan) + length(upTspan) + 2
    } else {
        Rcpp_code <- Rcpp_code[-ratelines[1]]
        ratelines <- ratelines[-1] - 1
    }
    
    ## update states
    tempnSpace <- 12
    tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
    if(length(transitions) == 1) {
        ## remove u_tmp and cumrate
        temp <- grep("u_tmp = 0.0, ", Rcpp_code)
        Rcpp_code[temp] <- gsub("u_tmp = 0.0, ", "", Rcpp_code[temp])
        temp <- grep("cumrate = 0.0;", Rcpp_code)
        Rcpp_code[temp] <- gsub(", cumrate = 0.0", "", Rcpp_code[temp])
        
        upStates <- NULL
        ## update negative states
        temp <- which(transitions[[1]]$S < 0)
        if(length(temp) > 1){
            stop("Cannot update more than one state for each transition")
        }
        if(length(temp) > 0) {
            upStates <- c(upStates, paste0(tempSpace, "u[", temp - 1, "]--;"))
        }
        ## update positive states
        temp <- which(transitions[[1]]$S > 0)
        if(length(temp) > 1){
            stop("Cannot update more than one state for each transition")
        }
        if(length(temp) > 0) {
            upStates <- c(upStates, paste0(tempSpace, "u[", temp - 1, "]++;"))
            ## add incidence curve if required
            if(incidence) {
                upStates <- c(upStates, paste0(tempSpace, "u[", length(transitions) + temp, "]++;"))
            }
        }
        tempnSpace <- tempnSpace - 4
        tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
    } else {
        upStates <- paste0(tempSpace, "cumrate = rates[0];")
        for(i in 1:length(transitions)) {
            if(i < length(transitions)) {
                if(i > 1) {
                    upStates <- c(upStates, paste0(tempSpace, "cumrate += rates[", i - 1, "];"))
                }
                upStates <- c(upStates, paste0(tempSpace, "if(u_tmp < cumrate) {"))
                tempnSpace <- tempnSpace + 4
                tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
            }
            ## update negative states
            temp <- which(transitions[[i]]$S < 0)
            if(length(temp) > 1){
                stop("Cannot update more than one state for each transition")
            }
            if(length(temp) > 0) {
                upStates <- c(upStates, paste0(tempSpace, "u[", temp - 1, "]--;"))
            }
            ## update positive states
            temp <- which(transitions[[i]]$S > 0)
            if(length(temp) > 1){
                stop("Cannot update more than one state for each transition")
            }
            if(length(temp) > 0) {
                upStates <- c(upStates, paste0(tempSpace, "u[", temp - 1, "]++;"))
                ## add incidence curve if required
                if(incidence) {
                    upStates <- c(upStates, paste0(tempSpace, "u[", length(transitions) + temp, "]++;"))
                }
            }
            tempnSpace <- tempnSpace - 4
            tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
            if(i < length(transitions)) {
                upStates <- c(upStates, paste0(tempSpace, "} else {"))
                tempnSpace <- tempnSpace + 4
                tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
            }
        }
        for(i in 1:(length(transitions) - 1)) {
            upStates <- c(upStates, paste0(tempSpace, "}"))
            tempnSpace <- tempnSpace - 4
            tempSpace <- paste0(rep(" ", tempnSpace), collapse = "")
        }
    }
    
    ## update after_tstar
    if(!is.null(afterTstar)){
        afterTstar <- strsplit(afterTstar, "\n")[[1]]
        afterTstar <- sapply(afterTstar, function(x) paste0(tempSpace, x))
        currline <- ratelines[1]
        ratelines <- ratelines[-1]
        Rcpp_code <- c(Rcpp_code[1:(currline - 1)], afterTstar, Rcpp_code[(currline + 1):length(Rcpp_code)])
        ratelines <- ratelines + length(afterTstar) - 1
    } else {
        Rcpp_code <- Rcpp_code[-ratelines[1]]
        ratelines <- ratelines[-1] - 1
    }
    
    ## update tspan
    if(tspan) {
        currline <- ratelines[1]
        upTspan[1] <- "    while(tspan[k] < tstar && k < tspan.size()) {"
        upTspan <- paste(paste(rep(" ", 7), collapse = ""), upTspan)
        Rcpp_code <- c(
            Rcpp_code[1:(currline - 1)], 
            "",
            paste0(paste(rep(" ", 12), collapse = ""), "// update tspan"),
            upTspan, 
            Rcpp_code[(currline + 1):length(Rcpp_code)]
        )
        ratelines <- ratelines[-1] + length(upTspan) + 1
    } else {
        Rcpp_code <- Rcpp_code[-ratelines[1]]
        ratelines <- ratelines[-1] - 1
    }
    
    ## update rates
    upRates <- sapply(as.list(upRates), function(x) {
        paste0("        ", x)
    })
    ## if only one event, then remove random event draw
    if(length(transitions) == 1) {
        currline <- ratelines[1]
        Rcpp_code <- c(Rcpp_code[1:(currline - 2)], Rcpp_code[currline:length(Rcpp_code)])
        ratelines <- ratelines - 1
    }
    currline <- ratelines[1]
    ratelines <- ratelines[-1]
    Rcpp_code <- c(Rcpp_code[1:(currline - 1)], upStates, upRates, Rcpp_code[(currline + 1):length(Rcpp_code)])
    ratelines <- ratelines + length(upStates) + length(upRates) - 1
    
    ## update after_tstar
    if(!is.null(afterTstar)){
        afterTstar <- sapply(afterTstar, function(x) paste0("    ", x))
        currline <- ratelines[1]
        ratelines <- ratelines[-1]
        Rcpp_code <- c(Rcpp_code[1:(currline - 1)], afterTstar, Rcpp_code[(currline + 1):length(Rcpp_code)])
        ratelines <- ratelines + length(afterTstar) - 1
    } else {
        Rcpp_code <- Rcpp_code[-ratelines[1]]
        ratelines <- ratelines[-1] - 1
    }
    
    ## add stopping criteria
    if(!is.null(stopCrit)) {
        if(runFromR & tspan) {
            retCrit <- "    List out1(2);"
            retCrit <- c(retCrit, "    out1[0] = out(tspan.size(), _);")
            retCrit <- c(retCrit, "    out1[1] = out(Range(0, tspan.size() - 1), Range(1, out.ncol() - 1));")
            retCrit <- c(retCrit, "    return out1;")
            if(!is.null(stopCrit)){
                tempind <- grep("return out;", stopCrit)
                temp <- gsub("return out;", "", stopCrit[tempind])
                tretCrit <- gsub("    ", "", retCrit)
                tretCrit <- paste0(temp, tretCrit)
                stopCrit <- c(stopCrit[c(1:(tempind - 1))], tretCrit, stopCrit[(tempind + 1):length(stopCrit)])
            }
        }
        currline <- ratelines[1]
        Rcpp_code <- c(Rcpp_code[1:(currline - 1)], "", stopCrit, Rcpp_code[(currline + 1):length(Rcpp_code)])
        ratelines <- ratelines[-1] + length(stopCrit)
    } else {
        Rcpp_code <- Rcpp_code[-ratelines[1]]
        ratelines <- ratelines[-1] - 1
    }
    if(is.null(matchCrit)) {
        tempSpace <- paste(rep(" ", 4), collapse = "")
        if(!tspan){
            matchCrit <- paste0(tempSpace, "// set output\n",
                                tempSpace, "out[0] = (totrate == 0.0 ? 1:0);\n",
                                tempSpace, "out[1] = t;\n",
                                tempSpace, "out[Range(2, u.size() + 1)] = as<NumericVector>(u);")
            if(is.data.frame(obsProcess)){
                for(i in 1:nrow(obsProcess)){
                    matchCrit <- paste0(matchCrit, "\n", tempSpace, "out[u.size() + ", i - 1 + 2, "] = ", obsProcess$compiled[i])
                }
            }
        } else {
            matchCrit <- paste0(tempSpace, "// set output\n",
                                tempSpace, "while(k < tspan.size()) {\n",
                                tempSpace, "    out(k, 0) = NA_REAL;\n",
                                tempSpace, "    out(k, 1) = tspan[k];\n",
                                tempSpace, "    for(i = 0; i < u.size(); i++) {\n",
                                tempSpace, "        out(k, i + 2) = (double) u[i];\n",
                                tempSpace, "    }")
            if(incidence) {
                matchCrit <- paste0(matchCrit, "\n",
                                tempSpace, "    for(i = u.size() / 2; i < u.size(); i++) {\n",
                                tempSpace, "        u[i] = 0;\n",
                                tempSpace, "    }")
            }
            if(is.data.frame(obsProcess)){
                for(i in 1:nrow(obsProcess)){
                    matchCrit <- paste0(matchCrit, "\n", tempSpace, "    out(k, u.size() + ", i - 1 + 2, ") = ", obsProcess$compiled[i])
                }
            }
            matchCrit <- paste0(matchCrit, "\n", tempSpace, "    k++;\n",
                                tempSpace, "}\n")
            matchCrit <- paste0(matchCrit,
                                tempSpace, "out(tspan.size(), 0) = (totrate == 0.0 ? 1:0);\n",
                                tempSpace, "out(tspan.size(), 1) = t;\n",
                                tempSpace, "for(i = 0; i < u.size(); i++) {\n",
                                tempSpace, tempSpace, "out(tspan.size(), i + 2) = (double) u[i];\n",
                                tempSpace, "}")
            if(is.data.frame(obsProcess)){
                for(i in 1:nrow(obsProcess)){
                    matchCrit <- paste0(matchCrit, "\n", tempSpace, 
                                        "out(tspan.size(), u.size() + ", i - 1 + 2, ") = (out(tspan.size(), 1) == out(tspan.size() - 1, 1) ? out(tspan.size() - 1, u.size() + ", i - 1 + 2, "):", gsub(";", "", obsProcess$compiled[i]), ");")
                }
            }
        }
    }
    currline <- ratelines[1]
    ratelines <- ratelines[-1]
    Rcpp_code <- c(Rcpp_code[1:(currline - 1)], matchCrit, Rcpp_code[(currline + 1):length(Rcpp_code)])
    ratelines <- ratelines + length(matchCrit) - 1
    
    ## set return criteria
    if(runFromR & tspan) {
        retCrit <- "    List out1(2);"
        retCrit <- c(retCrit, "    out1[0] = out(tspan.size(), _);")
        retCrit <- c(retCrit, "    out1[1] = out(Range(0, tspan.size() - 1), Range(1, out.ncol() - 1));")
        retCrit <- c(retCrit, "    return out1;")
    } else {
        retCrit <- "    return out;"
    }
    currline <- ratelines[1]
    Rcpp_code <- c(Rcpp_code[1:(currline - 1)], retCrit, Rcpp_code[(currline + 1):length(Rcpp_code)])
    
    ## return parsed code
    Rcpp_code
}

