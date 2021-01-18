test_that("PMCMC works", {
    
    ## set up data to pass to PMCMC
    flu_dat <- data.frame(
        t = 1:14,
        Robs = c(3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4)
    )

    ## set up observation process
    obs <- data.frame(
        dataNames = "Robs",
        dist = "pois",
        p1 = "R + 1e-5",
        p2 = NA,
        stringsAsFactors = FALSE
    )

    ## set up model (no need to specify tspan
    ## argument as it is set in PMCMC())
    transitions <- c(
        "S -> beta * S * I / (S + I + R + R1) -> I",
        "I -> gamma * I -> R",
        "R -> gamma1 * R -> R1"
    )
    compartments <- c("S", "I", "R", "R1")
    pars <- c("beta", "gamma", "gamma1")
    model <- mparseRcpp(
        transitions = transitions,
        compartments = compartments,
        pars = pars,
        obsProcess = obs
    )

    ## set priors
    priors <- data.frame(
        parnames = c("beta", "gamma", "gamma1"),
        dist = rep("unif", 3),
        stringsAsFactors = FALSE)
    priors$p1 <- c(0, 0, 0)
    priors$p2 <- c(5, 5, 5)

    ## define initial states
    iniStates <- c(S = 762, I = 1, R = 0, R1 = 0)
    
    set.seed(50)

    ## run PMCMC algorithm for first three days of data
    post <- PMCMC(
        x = flu_dat[1:3, ],
        priors = priors,
        func = model,
        u = iniStates,
        npart = 75,
        niter = 5000,
        nprintsum = 1000
    )
    ## file to save results
    tmp <- "PMCMC"
    ## the first run always succeeds, but warns
    expect_known_output(post, tmp, print = TRUE)

    ## run predictions forward in time
    post_pred <- predict(
        window(post, start = 2000, thin = 3),
        tspan = 4:14
    )
    ## file to save results
    tmp <- "PMCMCpred"
    ## the first run always succeeds, but warns
    expect_known_output(post_pred, tmp, print = TRUE)
    
})
