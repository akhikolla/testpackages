test_that("ABCRef works", {
    
    ## set up SIR simulation model
    transitions <- c(
        "S -> beta * S * I -> I",
        "I -> gamma * I -> R"
    )
    compartments <- c("S", "I", "R")
    pars <- c("beta", "gamma")
    model <- mparseRcpp(
        transitions = transitions,
        compartments = compartments,
        pars = pars
    )
    model <- compileRcpp(model)

    ## generate function to run simulators
    ## and produce final epidemic size and time
    ## summary statistics
    simRef <- function(pars, model) {
        ## run model over a 100 day period with
        ## one initial infective in a population
        ## of 120 individuals
        sims <- model(pars, 0, 20, c(119, 1, 0))

        ## return vector of summary statistics
        c(finaltime = sims[2], finalsize = sims[5])
    }

    ## set priors
    priors <- data.frame(
        parnames = c("beta", "gamma"),
        dist = rep("gamma", 2),
        stringsAsFactors = FALSE
    )
    priors$p1 <- c(10, 10)
    priors$p2 <- c(10^4, 10^2)

    ## produce reference table by sampling from priors
    ## (add additional arguments to 'func' at the end)
    ## first test should fail as no sumNames
    expect_error(ABCRef(
        npart = 100,
        priors = priors,
        func = simRef,
        model = model
    ))
    set.seed(50)
    refTable <- ABCRef(
        npart = 100,
        priors = priors,
        func = simRef,
        sumNames = c("finaltime", "finalsize"),
        model = model
    )
    ## file to save results
    tmp <- "ABCRef"
    ## the first run always succeeds, but warns
    expect_known_output(refTable, tmp, print = TRUE)
    
    ## produce reference table by adding in parameters
    ## as matrix (first test should fail since pars and priors
    ## can't be set together)
    expect_error(ABCRef(
        pars = as.matrix(refTable[, 1:2]),
        npart = 100,
        priors = priors,
        func = simRef,
        sumNames = c("finaltime", "finalsize"),
        model = model
    ))
    set.seed(50)
    refTable <- ABCRef(
        pars = as.matrix(refTable[, 1:2]),
        npart = 100,
        func = simRef,
        sumNames = c("finaltime", "finalsize"),
        model = model
    )
    ## file to save results
    tmp <- "ABCRef1"
    ## the first run always succeeds, but warns
    expect_known_output(refTable, tmp, print = TRUE)
    
    ## check for reproducibility in serial
    set.seed(50)
    refTable1 <- ABCRef(
        pars = as.matrix(refTable[, 1:2]),
        npart = 100,
        func = simRef,
        sumNames = c("finaltime", "finalsize"),
        model = model
    )
    
    expect_identical(refTable, refTable1)
    
    ## check for reproducibility in parallel
    if(.Platform$OS.type != "windows" & requireNamespace("parallel", quietly = TRUE)) {
        set.seed(50)
        refTable <- ABCRef(
            pars = as.matrix(refTable[, 1:2]),
            npart = 100,
            func = simRef,
            sumNames = c("finaltime", "finalsize"),
            model = model,
            parallel = TRUE,
            mc.cores = 2
        )
        set.seed(50)
        refTable1 <- ABCRef(
            pars = as.matrix(refTable[, 1:2]),
            npart = 100,
            func = simRef,
            sumNames = c("finaltime", "finalsize"),
            model = model,
            parallel = TRUE,
            mc.cores = 2
        )
        
        expect_identical(refTable, refTable1)
        
        ## check for changing seeds in parallel
        set.seed(50)
        refTable <- ABCRef(
            pars = as.matrix(refTable[, 1:2]),
            npart = 100,
            func = simRef,
            sumNames = c("finaltime", "finalsize"),
            model = model,
            parallel = TRUE,
            mc.cores = 2
        )
        refTable1 <- ABCRef(
            pars = as.matrix(refTable[, 1:2]),
            npart = 100,
            func = simRef,
            sumNames = c("finaltime", "finalsize"),
            model = model,
            parallel = TRUE,
            mc.cores = 2
        )
        
        expect_false(identical(refTable, refTable1))
    }
    
})
