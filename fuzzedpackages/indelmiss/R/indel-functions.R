plottree <- function(x, toilabel = TRUE, colors = NULL, ...) {
  if (is.null(colors)) 
    colors <- grDevices::rainbow(length(x$bg))
  colo <- NULL
  for (i in 1:nrow(x$tree$edge)) {
    tmp <- which(unlist(lapply(lapply(x$bg, FUN = function(X) c(x$tree$edge[i, 1], x$tree$edge[i, 2]) %in% X), FUN = "all"), use.names = FALSE))
    if ((length(tmp) > 1)) {
      colo[i] <- tmp[which.min(unlist(lapply(x$bg, FUN = length)))]
    } else{
      colo[i] <- tmp
    }
  }
  if (toilabel) 
    x$tree$tip.label[x$toi] <- paste(x$tree$tip.label[x$toi], "+", sep = "")
  ape::plot.phylo(x$tree, edge.color = colors[colo], ...)
}

plot.indelmiss <- function(x, model = NULL, ci = TRUE, cil = 95, ...) {
  if (is.null(model)) 
    stop("model option is required")
  plotrates(x, model, ci, cil, ...)
  if (model %in% c("M2", "M4")) 
    plotp(x, model, ci, cil, ...)
}

plotp <- function(x, model, ci = TRUE, cil = 95, ...) {
  if (!model %in% c("M2", "M4")) 
    stop("This plot is reserved for models where p is estimated.")
  if (model %in% c("M2", "M4")) {
    pest <- x$results[[model]]$parsep$p
    pest_se <- x$results[[model]]$parsep$se$p
  }
  cil = qnorm(1 - (1 - cil/100)/2)
  plot(x = x$toi, y = pest, ylab = "Estimates", xlab = "Taxa", pch = 16, ylim = c(min(pest - cil * pest_se), max(pest + cil * pest_se)), 
       xaxt = "n")
  axis(1, at = x$toi, labels = x$toi)
  if (ci == TRUE) {
    segments(x$toi, pest - cil * pest_se, x$toi, pest + cil * pest_se)
    epsilon = 0.05
    segments(x$toi - epsilon, pest - cil * pest_se, x$toi + epsilon, pest - cil * pest_se)
    segments(x$toi - epsilon, pest + cil * pest_se, x$toi + epsilon, pest + cil * pest_se)
  }
}

plotrates <- function(x, model, ci = TRUE, cil = 95, ...) {
  rest <- x$results[[model]]$parsep$rates
  rest_se <- x$results[[model]]$parsep$se$rates
  xax <- 1:(2 * length(x$bg))
  xaxlab <- paste(c("mu", "nu"), rep(1:length(x$bg), each = 2, len = length(x$bg) * 2), sep = "")
  if (model %in% c("M1", "M2")) {
    rest <- rest[1, ]
    rest_se <- rest_se[1, ]
    xax <- 1:length(x$bg)
    xaxlab <- paste(c("mu"), 1:length(x$bg), sep = "")
  }
  cil = qnorm(1 - (1 - cil/100)/2)
  plot(x = xax, y = rest, ylab = "Estimates", xlab = "Clade", pch = 16, ylim = c(min(rest - cil * rest_se), max(rest + cil * rest_se)), 
       xaxt = "n")
  axis(1, at = xax, labels = xaxlab)
  if (ci == TRUE) {
    segments(xax, rest - cil * rest_se, xax, rest + cil * rest_se)
    epsilon = 0.05
    segments(xax - epsilon, rest - cil * rest_se, xax + epsilon, rest - cil * rest_se)
    segments(xax - epsilon, rest + cil * rest_se, xax + epsilon, rest + cil * rest_se)
  }
}

print.indelmiss <- function(x, ...) {
  cat("Call:\n")
  print.default(x$call)
  cat("\n", x$taxa, "taxa with", x$phyl, "gene families and", (length(x$w) - 1), "different phyletic gene patterns.\n")
  cat("-----------------------------------\n")
  cat("Groups of nodes with the same rates:\n")
  print.default(x$bg)
  cat("-----------------------------------\n")
  for (i in x$modelnames) {
    cat(i, "\n")
    if (x$results[[i]]$convergence == 0) {
      print.default(x$results[[i]]$parsep)
      if (i %in% c("M2", "M4")) {
        cat("Number of genes estimated as missing corresponding to the missing data proportions is:\n ")
        print(round(x$results[[i]]$estmissgenes))
      }
      if(x$optmethod=="nlminb") cat("\nLoglikelihood for model", i, ":", -x$results[[i]]$objective, "\n")
      if(x$optmethod=="optim") cat("\nLoglikelihood for model", i, ":", -x$results[[i]]$value, "\n")
      cat("AIC           for model", i, ":", x$results[[i]]$AIC, "\n")
      cat("BIC           for model", i, ":", x$results[[i]]$BIC, "\n")
      cat("-----------------------------------\n")
    } else {
      cat("Convergence not achieved.\n")
      cat("Increase number of iterations and/or function calls. ")
      cat("See control arguments for the optimization method used.\n")
      cat("-----------------------------------\n")
    }
  }
  if (!all(x$conv == 0)) 
    cat("Not all models converged. \n")
  cat("Time taken:", x$time[3], "seconds.\n")
}
checkbglist <- function(bg, tree){
  all_branches <- TRUE
  for (i in 1:length(tree$edge.length)) {
    if (!any(unlist(lapply(lapply(bg, FUN = function(j) tree$edge[i,] %in% j),FUN="all")))) all_branches <- FALSE
  }
  return(all_branches)
}
indelrates <- function(verbose = FALSE, usertree = NULL, userphyl = NULL, matchtipstodata = FALSE, datasource = "user", seed = 1, taxa = 5, brlensh = c(1, 4), mu = 1, 
                       nu = 1, phyl = 5000, pmiss = 0, toi = 1, bgtype = "listofnodes", bg = NULL, zerocorrection = TRUE, rootprob = "stationary", rpvec = NULL, optmethod = "nlminb", init = 0.9, lowlim = 0.001, uplim = 100,
                       numhessian = TRUE, modelnames = c("M1", "M2", "M3", "M4"), ...) {
  #indelrate is M1; indelandp is M2; insanddelrate is M3; insanddelandp is M4.
  ptm <- proc.time()
  bgtype <- match.arg(bgtype, c("listofnodes", "ancestornodes"))
  modelnames <- match.arg(modelnames, c("M1", "M2", "M3", "M4"), several.ok = TRUE)
  if (!is.null(rootprob)) 
    rootprob <- match.arg(rootprob, c("equal", "stationary", "maxlik", "user"))
  if (datasource == "user" & is.null(usertree)) 
    stop("usertree option is required.")
  if (any(pmiss < 0 | pmiss > 1)) stop("pmiss is the proportion of gene families that are (or are treated as) in the genome but are not observed. Should be between 0 and 1.")
  if ((is.null(userphyl) | !(is.matrix(userphyl) | is.data.frame(userphyl))) & datasource == "user") 
    stop("userphyl option is required.")
  if (bgtype == "listofnodes" & !is.null(bg) & !is.list(bg)) 
    stop("bg should be a vector with listofnodes argument.")
  if (bgtype == "listofnodes" & !is.null(bg) & is.list(bg)) {
    if (!checkbglist(bg, usertree)) stop("Some branch is missing in bg.")
  }
  if (bgtype == "ancestornodes" & !is.null(bg) & !is.vector(bg)) 
    stop("bg should be a list with ancestornodes argument.")
  if (is.null(rootprob)) 
    stop("rootprob option is required.")
  if (datasource == "user") {
    if (length(table(as.matrix(userphyl))) != 2) {
      stop("Number of discrete states in the data provided does not equal 2. The discrete states allowed are 0 or 1.")
    }
  }
  
  if (datasource == "user") {
    if ((nrow(userphyl) < ncol(userphyl))) {
      warning("Looks like there are more taxa than gene families. The rows of the matrix given as the value of the userphyl argument should represent the gene families and the columns the taxa.")
    }
  }
  if (is.null(datasource)) 
    stop("datasource option is required")
  if (datasource == "user") {
    if (is.data.frame(userphyl)) {
      userphyl <- as.matrix(userphyl)
    }
  }
  
  #   if (is.null(userphyl) & datasource == "user") {
  #     stop("userphyl option is required")
  #}
  #if(is.numeric(toi) & !length(toi)==length(nmiss))
  #stop("Length of vector of the number of present genes to be removed should be the same as the length of the vector of taxa of interest.")
  set.seed(seed)
  # singlebranch <- logical(length = 1)
  
  #############Libraries###################
  libload <- function(funcval) {
    if (requireNamespace(paste(funcval), quietly = TRUE) == FALSE) 
      stop("Please install package dependencies before continuing.") 
    suppressMessages(loadNamespace(funcval))
  }
  libload("ape")
  libload("numDeriv")
  libload("Rcpp")
  #########Transition rate matrix and substitution rate matrix#######
  TPM_taxa <- function(rates, ad, ti) {
    #Insertion rate is nu; deletion rate is mu
    j <- unlist(lapply(lapply(bg, FUN = function(X) c(ad[1], ad[2]) %in% X), FUN = "all"))
    if ((sum(j) > 1)) j[j][-which.min(unlist(lapply(bg, FUN = length)))] <- FALSE
    mu <- rates[1, j] #
    nu <- rates[2, j] #
    return(matrix(c(mu + nu * exp(-(mu + nu) * ti), nu - nu * exp(-(mu + nu) * ti), mu - mu * exp(-(mu + nu) * ti), nu + mu * exp(-(mu + 
                                                                                                                                      nu) * ti)) * 1/(mu + nu), 2, 2, byrow = TRUE))
  }
  ###########Data Generation/ Import#########
  # lowlim = 0.01
  # uplim = 100
  alphabet <- c(0, 1)
  al <- length(alphabet)
  if (datasource == "user") {
    if (ape::is.binary.phylo(usertree) == FALSE | ape::is.rooted(usertree) == FALSE) {
      usertree <- ape::multi2di(usertree)
      cat("Tree either not binary or not rooted.\n")
      cat("Transformed using multi2di from ape.", "\n")
      cat("Labels of nodes of interest might change. Please check the tree in the fitted object.\n")
    }
    tree1 <- usertree
    tree1 <- ape::reorder.phylo(tree1, order = "postorder")
    if (matchtipstodata) {
      fin <- userphyl[,pmatch(usertree$tip.label, colnames(userphyl))] 
      datab <- fin
    } else {
      datab <- userphyl
    }
    # datab <- t(userphyl)
    phyl <- nrow(datab)
    nooftaxa <- length(tree1$tip.label)
    mu <- NULL
    nu <- NULL
    brlensh <- NULL
  }
  if (datasource == "simulation") {
    libload("phangorn")
    nooftaxa <- taxa
    tree1 <- ape::rtree(nooftaxa, br = rbeta(n = (nooftaxa - 2), shape1 = brlensh[1], shape2 = brlensh[2]))
    tree1 <- ape::reorder.phylo(tree1, order = "postorder")
    data <- phangorn::simSeq(tree1, l = phyl, type = "USER", levels = alphabet, bf = c(mu/(mu + nu), nu/(mu + nu)), Q = mu)
    datab <- matrix(as.numeric(as.character(data)), nrow = nooftaxa)
    datab <- t(datab)
  }
  if (bgtype == "listofnodes" & is.null(bg)) {
    bg <- list(c(1:(nooftaxa * 2 - 1)))
  } else {# if (bgtype == \listofnodes\" & !is.null(bg)) "
    bg <- bg
  }
  if (bgtype == "ancestornodes") {
    uncl <- sort(bg, decreasing = TRUE) #unique clades
    for (i in 1:length(uncl)) {
      if (!is.na(uncl[i])) {
        try <- phangorn::Ancestors(tree1, uncl[i], type = "all")
        res <- match(try, uncl)
        uncl[na.omit(res)] <- NA
      }
    }
    bg <- uncl <- na.omit(uncl)
    ########
    if (is.null(bg)) 
      stop("bg not specified.")
    bgo <- bg
    bg <- list()
    for (i in 1:length(bgo)) {
      bg[[i]] <- c(bgo[i], phangorn::Descendants(tree1, bgo[i], type = "all"))
    }
    bg[[length(bg) + 1]] <- c(setdiff(1:(nooftaxa * 2 - 1), c(unlist(bg))), bgo)
  }
  ptogenes <- function(model) {
    presentgenes <- NULL
    est_missing <- NULL
    spect <- 1:length(toi)
    for (i in 1:length(toi)) {
      presentgenes <- sum((databp_red)[, toi[i]] * w)
      p <- res$parsep$p[spect[i]]
      est_missing[i] <- presentgenes/(1 - p) - presentgenes
    }
    return(est_missing)
  }
  ##############Parameters and Variables################
  nointnodes <- nooftaxa - 1
  tips <- 1:nooftaxa
  len_tips <- length(tips)
  csp <- length(bg) #how many different clade-specific parameters
  #if(simplify==TRUE) 
  datab[which(datab > 1)] <- 1
  if (is.character(toi)) {
    if (toi == "all") {
      toi = 1:nooftaxa
    }
  } else {
    toi <- toi
  } #tip label of taxa of interest
  if (any(pmiss > 0)) {
    for (e in 1:length(toi)) {
      datab[sample(which(datab[, toi[e] ] == 1), pmiss[e] * sum(datab[, toi[e] ] == 1), replace = FALSE), toi[e] ] <- 0
      #randomly setting some "P" to "A"
    }
  }
  
  check.equal <- function(x, y) {
    isTRUE(all.equal(y, x, check.attributes = FALSE))
  }
  zeroentr <- which(apply(datab, 1, check.equal, y = rep(0, nooftaxa)))
  ifelse(length(zeroentr) > 0, databp <- datab[-zeroentr, ], databp <- datab)
  # databp <- datab
  
  w <- summary(as.factor(apply(databp, 1, paste, collapse = " ")), maxsum = nrow(databp))
  b <- attr(w, "names")
  ww <- unname(cbind(b, w))
  databp_red <- matrix(na.omit(as.numeric(unlist(strsplit(ww[,1],split = " ")))), ncol = nooftaxa, byrow = T)
  databp_red <- rbind(databp_red, rep(0, nooftaxa))
  if (datasource == "user") 
    colnames(databp_red) <- colnames(userphyl)
  w <- c(w, 1) #Count for no gene correction...
  logll <- NULL
  results <- list()
  nodelist <- c(setdiff(tree1$edge[, 2], tips), length(tips) + 1)
  F <- cbind(tree1$edge, tree1$edge.length)
  rates <- matrix(NA, nrow = 2, ncol = csp)
  #rates: parameters to be estimated
  ###############Log-likelihood function###############
  bamsp <- function(previtval, model) { #branch and model specific assignment
    le_prev <- length(previtval)
    if (model == "M1") {
      if (rootprob == "maxlik") {
        for (i in 1:csp) rates[, i] <- (previtval[-le_prev])[i]
        return(list(rates = rates, iroot = previtval[le_prev]))
      } else {
        for (i in 1:csp) rates[, i] <- previtval[i]
        return(list(rates = rates))
      }
    } else if (model == "M2") {
      if (rootprob == "maxlik") {
        for (i in 1:csp) {
          rates[, i] <- previtval[1:(le_prev - length(toi) - 1)][i]
        }
        p <- previtval[(le_prev - length(toi)):(le_prev - 1)]
        return(list(rates = rates, p = p, iroot = previtval[le_prev]))
      } else {
        for (i in 1:csp) rates[, i] <- previtval[1:(le_prev - length(toi))][i]
        p <- previtval[(le_prev - length(toi) + 1):le_prev]
        return(list(rates = rates, p = p))
      }
    } else if (model == "M3") {
      if (rootprob == "maxlik") {
        rates[] <- previtval[-le_prev] #by default should enter by column
        return(list(rates = rates, iroot = previtval[le_prev]))
      } else {
        rates[] <- previtval #by default should enter by column
        return(list(rates = rates))
      }
    } else if (model == "M4") {      
      if (rootprob == "maxlik") {
        rates[] <- previtval[1:(le_prev - length(toi) - 1)]
        #rates, p, genefamilypresencerootprobability
        p <- previtval[(le_prev - length(toi)):(le_prev - 1)]
        return(list(rates = rates, p = p, iroot = previtval[le_prev]))
      } else {
        rates[] <- previtval[1:(length(previtval) - length(toi))]
        p <- previtval[(length(previtval) - length(toi) + 1):length(previtval)]
        return(list(rates = rates, p = p))
      }
    }
  }
  # Lixi <- matrix(0, nooftaxa + nointnodes, al)
  alseq <- 1:length(alphabet)
  update_Lixi_init <- function(miss_curr){
    p <- miss_curr
    #miss_curr is the toi-ordered vector of current values for the missing data proportion.
    for (i in 1:length(w)) {
      if (length(toi) > 1) {
        Lixi_init[[i]][toi, ][databp_red[i, ][toi] == 1, ] <- cbind(0, 1 - p[databp_red[i, ][toi] == 1])
        #prob. that we observe gene given that it is truly absent or present
        Lixi_init[[i]][toi, ][databp_red[i, ][toi] == 0, ] <- cbind(1, p[databp_red[i, ][toi] == 0])
        #prob. that we don't observe gene given that it is truly absent or present
      }
      if (length(toi) == 1) {
        if (databp_red[i, ][toi] == 1) 
          Lixi_init[[i]][toi, ] <- c(0, 1 - p)
        else if (databp_red[i, ][toi] == 0) 
          Lixi_init[[i]][toi, ] <- c(1, p)
      }
    }
    return(Lixi_init)
  }
  tmp <- numeric(al)
  ll <- function(model, pm, rootp, Lixi_in_ll) {
    tmp <- allpatt_loopC(nodelist, al, tree1$edge[, 1], tree1$edge[, 2], pm, Lixi_in_ll, len_tips)
    logll <- log(rowSums(sweep(tmp, MARGIN = 2, rootp, `*`)))#log(sum(rootp * tmp))
    return(logll)
  }
  #################Objective function################
  totalll <- function(previtval, model, Lixi_in_totll) {
    ptbrun <- bamsp(previtval, model)
    if (model == "M2" | model == "M4") {
      Lixi_in_totll <- update_Lixi_init(ptbrun$p)
    }
    ptb_loc <- ptbrun$rates
    rootp <- rootpfun(rootprob, currpar = ptbrun)
    pm <- lapply(1:nrow(F), function(i) TPM_taxa(ptb_loc, ad = F[i, c(1:2)], ti = F[i, 3]))
    phy <- w * ll(model, pm, rootp, Lixi_in_ll = Lixi_in_totll)
    if (zerocorrection == TRUE) 
      return(-(sum(phy[-length(phy)]) - nrow(databp) * log(1 - exp(phy[length(phy)]))))
    if (zerocorrection == FALSE) 
      return(-(sum(phy[-length(phy)])))
    #-ve of the value because we are minimizing by default.
  }
  rootpfun <- function(rootprob, currpar) {
    if (rootprob == "equal") {
      return(rep(1/al, al))
    } else if (rootprob == "maxlik") {
      irootprob <- currpar$iroot
      return(c(irootprob, 1 - irootprob))
    } else if (rootprob == "user") {
      return(rpvec)
    } else {
      rates_rp <- currpar$rates
      mptb2 <- mean.default(rates_rp[2, ])
      mptb1 <- mean.default(rates_rp[1, ])
      return(c(mptb1/(mptb1 + mptb2), mptb2/(mptb1 + mptb2)))
    }
  }
  indelinit <- init
  missingguess <- rep(0, length(toi))
  for (e in 1:length(toi)) {
    if (toi[e] < nooftaxa) 
      missingguess[e] <- abs(1 - summary(as.factor(databp[, toi[e] ]))[2]/summary(as.factor(databp[, (toi[e] + 1) ]))[2])
    else if (toi[e] == nooftaxa) 
      missingguess[e] <- abs(1 - summary(as.factor(databp[, toi[e] ]))[2]/summary(as.factor(databp[, (toi[e] - 1) ]))[2])
  }
  missingguess[missingguess >= 1] <- 0
  ################Estimation################
  options(digits = 7)
  if (rootprob == "maxlik") {
    modelop <- list(M1 = list(start = c(rep(indelinit, csp), 0.5), df = 1 * csp + 1, pb = 1, lower = c(rep(lowlim, csp), 0.01), upper = c(rep(uplim, csp), 0.99), model = "M1"), M2 = list(start = c(rep(indelinit, csp), missingguess, 0.5), df = 1 * csp + 1 * length(toi) + 1, pb = 1, lower = c(rep(lowlim, csp), rep(0, length(toi)), 0.01), upper = c(rep(uplim, csp), rep(1, length(toi)), 0.99), model = "M2"), M3 = list(start = c(rep(indelinit, csp), rep(indelinit, csp), 0.5), df = 2 * csp + 1, pb = 2, lower = c(rep(c(lowlim, lowlim), csp), 0.01), upper = c(rep(c(uplim, uplim), csp), 0.99), model = "M3"), M4 = list(start = c(rep(indelinit, csp), rep(indelinit, csp), missingguess, 0.5), df = 2 * csp + 1 * length(toi) +  1, pb = 2, lower = c(c(rep(c(lowlim, lowlim), csp)), rep(0, length(toi)), 0.01), upper = c(c(rep(c(uplim, uplim), csp)), rep(1, length(toi)), 0.99), model = "M4"))
  } else {
    modelop <- list(M1 = list(start = c(rep(indelinit, csp)), df = 1 * csp, pb = 1, lower = c(rep(lowlim, csp)), upper = c(rep(uplim, csp)), model = "M1"), M2 = list(start = c(rep(indelinit, csp), missingguess), df = 1 * csp + 1 * length(toi), pb = 1, lower = c(rep(lowlim, csp), rep(0, length(toi))), upper = c(rep(uplim, csp), rep(1, length(toi))), model = "M2"), M3 = list(start = c(rep(indelinit, csp), rep(indelinit, csp)), df = 2 * csp, pb = 2, lower = c(rep(c(lowlim, lowlim), csp)), upper = c(rep(c(uplim, uplim), csp)), model = "M3"), M4 = list(start = c(rep(indelinit, csp), rep(indelinit, csp), missingguess), df = 2 * csp + 1 * length(toi), pb = 2, lower = c(c(rep(c(lowlim, lowlim), csp)), rep(0, length(toi))), upper = c(c(rep(c(uplim, uplim), csp)), rep(1, length(toi))), model = "M4"))
  }
  parvec <- NULL
  sevec <- NULL
  convcheck <- NULL
  
  for (i in modelnames) {
    Lixi_init <- rep(list(matrix(data = 0, nrow = nooftaxa + nointnodes, ncol = al)), length = length(w))
    for (k in 1:length(w)) {
      for (u in alseq) {
        Lixi_init[[k]][which(databp_red[k, ] == alphabet[u]), u] <- 1
      }
    }
    if (optmethod == "nlminb") {
      res <- nlminb(start = modelop[[i]]$start, objective = totalll, model = i, lower = modelop[[i]]$lower, upper = modelop[[i]]$upper, Lixi_in_totll = Lixi_init, ...)
      if (numhessian) 
        res$hessian <- numDeriv::hessian(totalll, model = i, res$par, Lixi_in_totll = Lixi_init)
      res$parsep <- bamsp(res$par, i)
      rownames(res$parsep$rates) <- c("mu", "nu")
    }
    if (optmethod == "optim") {
      res <- optim(par = modelop[[i]]$start, fn = totalll, model = i, lower = modelop[[i]]$lower, upper = modelop[[i]]$upper, Lixi_in_totll = Lixi_init, hessian = TRUE, 
                   method = "L-BFGS-B", ...)
      res$parsep <- bamsp(res$par, i)
      rownames(res$parsep$rates) <- c("mu", "nu")
    }
    res$df <- modelop[[i]]$df
    convcheck <- c(convcheck, res$convergence)
    if (verbose == TRUE) 
      cat("---------------------------", "\n")
    if (res$convergence == 0) {
      estimate <- res$par
      parvec <- c(parvec, res$par[1:modelop[[i]]$pb])
      if (numhessian) {
        res$se <- try(sqrt(diag(solve(res$hessian))), silent = TRUE)
        res$parsep$se <- try(bamsp(res$se, i), silent = TRUE)
        rownames(res$parsep$se$rates) <- try(c("mu", "nu"), silent = TRUE)
        sevec <- try(c(sevec, res$se), silent = TRUE)
      }
      if (optmethod == "nlminb") 
        llhat <- -round(res$objective, 3)
      if (optmethod == "optim") 
        llhat <- -round(res$value, 3)
      res$AIC <- 2 * llhat - 2 * res$df
      res$BIC <- 2 * llhat - log(phyl) * res$df
      if (i %in% c("M2", "M4")) 
        res$estmissgenes <- ptogenes(i)
    }
    results[[i]] <- res
    if (verbose == TRUE) {
      if (res$convergence == 0) {
        cat("Estimate from model", i, ":", estimate, "\n")
        if (numhessian) 
          cat("SEs for model parameters:", res$se, "\n")
        cat("Estimated log-likelihood   :", llhat, "\n")
        cat("---------------------------", "\n")
      } else {
        cat("Convergence not achieved.\n")
        cat("Increase number of iterations and/or function calls. ")
        cat("See control arguments for the optimization method used.\n")
      }
    }
  }
  if (lowlim %in% parvec || uplim %in% parvec) {
    cat("Estimated parameters on interval bounds.")
  }
  if (!all(is.finite(sevec)) & numhessian) {
    cat("Something is not right with the standard errors.")
    cat("Check Hessian matrix estimate.\n")
    cat("Consider calculating bootstrap errors (make sure to use numhessian=FALSE).\n")
  }
  #  ifelse(length(toi) > 1, removed <- nmiss/colSums(databp[, toi ] == 1), removed <- nmiss/sum(databp[, toi ] == 1))
  # removed <- pmiss
  
  timetaken <- proc.time() - ptm
  val <- list(call = match.call(), conv = convcheck, time = timetaken, bgtype = bgtype, bg = bg, results = results, tree = tree1, verbose = verbose, 
              data_red = databp_red, w = w, mu = mu, nu = nu, datasource = datasource, seed = seed, pmiss = pmiss, toi = toi, zerocorrection = zerocorrection, 
              optmethod = optmethod, brlensh = brlensh, taxa = nooftaxa, phyl = nrow(databp), rootprob = rootprob, modelnames = modelnames)
  class(val) <- "indelmiss"
  return(invisible(val))
}