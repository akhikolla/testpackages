test_that("mparse/run works", {
    
    ## test that single transition model compiles
    
    ## set up transition rates
    transitions <- "S -> beta * S * I -> I"
    
    ## set compartment and parameter names
    compartments <- c("S", "I")
    pars <- "beta"
    
    ## parse model
    model <- mparseRcpp(
        transitions = transitions, 
        compartments = compartments,
        pars = pars
    )
    expect_error(compileRcpp(model), NA)
    
    ## parse model
    transitions <- c(
        "S -> beta * S * I -> I",
        "I -> gamma * I -> R"
    )
    compartments <- c("S", "I", "R")
    pars <- c("beta", "gamma")
    
    ## check standard parsed model
    model <- mparseRcpp(
        transitions = transitions,
        compartments = compartments,
        pars = pars
    )
    ## file to save results
    tmp <- "mparse"
    ## the first run always succeeds, but warns
    expect_known_output(model, tmp, print = TRUE)
    
    ## set seed
    set.seed(50)
    ## compile and run model
    sims <- run(
        model = model,
        pars = c(beta = 0.001, gamma = 0.1),
        tstart = 0,
        tstop = 20,
        u = c(S = 119, I = 1, R = 0)
    )
    ## file to save results
    tmp <- "run"
    ## the first run always succeeds, but warns
    expect_known_output(sims, tmp, print = TRUE)
    
    ## check model with incidence
    model <- mparseRcpp(
        transitions = transitions,
        compartments = compartments,
        pars = pars,
        incidence = TRUE
    )
    ## file to save results
    tmp <- "mparseinc"
    ## the first run always succeeds, but warns
    expect_known_output(model, tmp, print = TRUE)
    
    ## set seed
    set.seed(50)
    ## compile and run model
    sims <- run(
        model = model,
        pars = c(beta = 0.001, gamma = 0.1),
        tstart = 0,
        tstop = 20,
        u = c(S = 119, I = 1, R = 0, S_inc = 0, I_inc = 1, R_inc = 0)
    )
    ## file to save results
    tmp <- "runinc"
    ## the first run always succeeds, but warns
    expect_known_output(sims, tmp, print = TRUE)
    
    ## check model with tspan
    model <- mparseRcpp(
        transitions = transitions,
        compartments = compartments,
        pars = pars,
        tspan = TRUE
    )
    ## file to save results
    tmp <- "mparsetspan"
    ## the first run always succeeds, but warns
    expect_known_output(model, tmp, print = TRUE)
    
    ## set seed
    set.seed(50)
    ## compile and run model
    sims <- run(
        model = model,
        pars = c(beta = 0.001, gamma = 0.1),
        tstart = 0,
        tstop = 20,
        u = c(S = 119, I = 1, R = 0),
        tspan = 1:20
    )
    ## file to save results
    tmp <- "runtspan"
    ## the first run always succeeds, but warns
    expect_known_output(sims, tmp, print = TRUE)
    
    ## check model with incidence and tspan
    model <- mparseRcpp(
        transitions = transitions,
        compartments = compartments,
        pars = pars,
        incidence = TRUE,
        tspan = TRUE
    )
    ## file to save results
    tmp <- "mparseinctspan"
    ## the first run always succeeds, but warns
    expect_known_output(model, tmp, print = TRUE)
    
    ## set seed
    set.seed(50)
    ## compile and run model
    sims <- run(
        model = model,
        pars = c(beta = 0.001, gamma = 0.1),
        tstart = 0,
        tstop = 20,
        u = c(S = 119, I = 1, R = 0, S_inc = 0, I_inc = 1, R_inc = 0),
        tspan = 1:20
    )
    ## file to save results
    tmp <- "runinctspan"
    ## the first run always succeeds, but warns
    expect_known_output(sims, tmp, print = TRUE)
    
    ## CHECK PARALLELISATION REPRODUCIBILITY
    
    ## check serial runs with set seed are identical
    
    ## set seed
    set.seed(50)
    ## compile and run model
    sims <- run(
        model = model,
        pars = c(beta = 0.001, gamma = 0.1),
        tstart = 0,
        tstop = 20,
        u = c(S = 119, I = 1, R = 0, S_inc = 0, I_inc = 1, R_inc = 0),
        tspan = 1:20
    )
    
    ## set seed
    set.seed(50)
    ## compile and run model
    sims1 <- run(
        model = model,
        pars = c(beta = 0.001, gamma = 0.1),
        tstart = 0,
        tstop = 20,
        u = c(S = 119, I = 1, R = 0, S_inc = 0, I_inc = 1, R_inc = 0),
        tspan = 1:20
    )
    
    ## check same output
    expect_identical(sims, sims1)
    
    ## check reproducibility in parallel for the same mc.cores
    
    if(.Platform$OS.type != "windows" & requireNamespace("parallel", quietly = TRUE)) {
        ## set seed
        set.seed(50)
        ## compile and run model
        sims <- run(
            model = model,
            pars = c(beta = 0.001, gamma = 0.1),
            tstart = 0,
            tstop = 20,
            u = c(S = 119, I = 1, R = 0, S_inc = 0, I_inc = 1, R_inc = 0),
            tspan = 1:20,
            nrep = 2,
            parallel = TRUE,
            mc.cores = 2
        )
        
        ## set seed
        set.seed(50)
        ## compile and run model
        sims1 <- run(
            model = model,
            pars = c(beta = 0.001, gamma = 0.1),
            tstart = 0,
            tstop = 20,
            u = c(S = 119, I = 1, R = 0, S_inc = 0, I_inc = 1, R_inc = 0),
            tspan = 1:20,
            nrep = 2,
            parallel = TRUE,
            mc.cores = 2
        )
        
        ## check same output
        expect_identical(sims, sims1)
        
        ## check result when not setting seed second time
        
        ## set seed
        set.seed(50)
        ## compile and run model
        sims <- run(
            model = model,
            pars = c(beta = 0.001, gamma = 0.1),
            tstart = 0,
            tstop = 20,
            u = c(S = 119, I = 1, R = 0, S_inc = 0, I_inc = 1, R_inc = 0),
            tspan = 1:20,
            nrep = 2,
            parallel = TRUE,
            mc.cores = 2
        )
        
        ## compile and run model
        sims1 <- run(
            model = model,
            pars = c(beta = 0.001, gamma = 0.1),
            tstart = 0,
            tstop = 20,
            u = c(S = 119, I = 1, R = 0, S_inc = 0, I_inc = 1, R_inc = 0),
            tspan = 1:20,
            nrep = 2,
            parallel = TRUE,
            mc.cores = 2
        )
        
        ## check same output
        expect_false(identical(sims, sims1))
    }
    
})
