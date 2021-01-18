## returns a matrix of simulated states
bootStates <- function(dataset, func, pars, u, npart = 50, ...) {
    
    ## 'func' must be SimBIID_model object
    if(class(func) != "SimBIID_model"){
        stop("'func' must be 'SimBIID_model' object")
    }
    
    ## check parameters
    checkInput(pars, c("numeric", "matrix"), ncol = length(func$pars))
    if(!identical(colnames(pars), func$pars)){
        stop("'colnames(pars)' does not match model parameters")
    }
    
    ## check data set
    checkInput(dataset, "data.frame")
    if(colnames(dataset)[1] != "t"){
        stop("First column of 'data' must be 't'")
    }
    if(ncol(dataset) < 2) {
        stop("Must have at least one count column in 'data'")
    }
    checkInput(dataset$t, "numeric")
    for(j in 2:ncol(dataset)) {
        checkInput(dataset[, j, drop = TRUE], "numeric", int = TRUE)
    }
    if(!identical(func$obsProcess$dataNames, colnames(dataset)[-1])){
        stop("'data' columns do not match 'obsProcess$datanames'")
    }
    dataset <- as.matrix(dataset)
    
    ## check number of particles
    checkInput(npart, c("numeric", "vector"), 1, int = TRUE, gt = 0)
    
    ## checks on function / input
    if(func$tspan){
        warning("'SimBIID_model' object will have 'tspan' set to F for particle filter\n")
    }
    if(is.null(func$obsProcess[1])){
        stop("'SimBIID_model' must have non-NULL 'obsProcess'")
    }
    if(!is.null(func$addVars[1])){
        stop("'SimBIID_model' can't have non-NULL 'addVars'")
    }
    if(!is.null(func$stopCrit[1])){
        stop("'SimBIID_model' can't have non-NULL 'stopCrit'")
    }
    
    ## check u
    checkInput(u, c("vector", "numeric"), int = TRUE, gte = 0)
    checkInput(sum(u), "numeric", int = TRUE, gt = 1)
    checkInput(u, length = length(func$compartments))
    if(!identical(names(u), func$compartments)) {
        stop("'names(u)' does not match 'func$compartments'")
    }
    
    if(func$incidence) {
        func$compartments <- func$compartments[-grep("_inc$", func$compartments)]
    }
    ## generate model
    func <- mparseRcpp(
        transitions = func$transitions,
        compartments = func$compartments,
        pars = func$pars,
        obsProcess = func$obsProcess,
        addVars = NULL,
        stopCrit = NULL,
        tspan = FALSE,
        incidence = func$incidence,
        afterTstar = NULL,
        PF = TRUE,
        runFromR = FALSE
    )
    
    ## compile model
    compfunc <- compileRcpp(func)
    
    ## run function
    output <- bootstrapPartFilterState(npart, pars, dataset, u, compfunc)
    pars <- pars[!is.na(output) & !map_lgl(output, is.null), ]
    output <- output[!is.na(output)]
    output <- output[!map_lgl(output, is.null)]
    if(length(output) == 0) {
        stop("No particles with non-zero weight at the end of the fitted period")
    }
    output <- lapply(1:length(output), function(i, x, time){
        x <- x[[i]]
        x <- cbind(rep(i, nrow(x)), time, x)
        x
    }, x = output, time = dataset[, 1])
    output <- do.call("rbind", output)
    colnames(output) <- c("rep", "t", names(u))
    output <- as.data.frame(output)
    
    ## add observation process sampling if required
    if(is.data.frame(func$obsProcess)){
        obsProcess <- func$obsProcess
        ntime <- nrow(output[output$rep == 1, ])
        for(i in 1:nrow(obsProcess)){
            if(obsProcess$dist[i] != "pois"){
                obsProcess$compiled[i] <- paste0("r", obsProcess$dist[i], 
                                                 "(", ntime, ", ", obsProcess$p1[i], ", ", obsProcess$p2[i], ")")
            } else {
                obsProcess$compiled[i] <- paste0("r", obsProcess$dist[i], 
                                                 "(", ntime, ", ", obsProcess$p1[i], ")")
            }
        }
        for(j in 1:nrow(obsProcess)){
            eval(parse(text = paste0("output$", obsProcess$dataNames, " <- NA")))
        }
        for(i in 1:nrow(pars)){
            for(j in 1:nrow(obsProcess)){
                temp <- with(c(output[output$rep == i, ], pars[i, ]), eval(parse(text = obsProcess$compiled[j])))
                eval(parse(text = paste0("output$", obsProcess$dataNames, "[output$rep == ", i, "] <- temp")))
            }
        }
    }
    ## return output
    list(pars = pars, output = output)
}

