

# helper function to create Xmeans for smooth plotting
create_xmeans <- function(xvars_static, Xorig, nvals, xvars_np, xvars_npargs,
                          xvars_bivariate, xvars_basis,
                          npdegree, has_intercept, betanms) {

  X1b <- NULL
  if (length(xvars_bivariate) > 0) {
    X1b <- as.matrix(Xorig[, unlist(xvars_bivariate)])
    X1b <- matrix(colMeans(X1b), ncol=ncol(X1b), nrow=nvals, byrow=T)
    colnames(X1b) <- unlist(xvars_bivariate)
  }

  X1 <- NULL
  if (length(xvars_static) > 0) {
    X1 <- as.matrix(Xorig[, xvars_static])
    X1 <- matrix(colMeans(X1), ncol=ncol(X1), nrow=nvals, byrow=T)
    colnames(X1) <- xvars_static
  }

  # tps / trunc.poly
  X2 <- NULL
  X3 <- NULL
  
  if (length(xvars_np) > 0) {
    # tps
    xvars_np_tps <- unlist(xvars_npargs[xvars_basis == "tps"])
    
    X2 <- as.matrix(Xorig[, xvars_np_tps])
    X2 <- matrix(colMeans(X2), ncol=ncol(X2), nrow=nvals, byrow=T)
    colnames(X2) <- xvars_np_tps

    # trunc.poly
    xvars_np_truncpoly <- xvars_np[xvars_basis == "trunc.poly"]
    xvars_np_truncpoly_degrees <- npdegree[xvars_basis == "trunc.poly"]
    X3mx1 <- as.matrix(Xorig[, xvars_np_truncpoly])
    X3mx1 <- matrix(colMeans(X3mx1), ncol=ncol(X3mx1), nrow=nvals, byrow=T)
    X3df <- as.data.frame(X3mx1)
    colnames(X3df) <- xvars_np_truncpoly

    getX <- mapply(np, x1=as.list(X3df), degree=as.list(xvars_np_truncpoly_degrees),
                   MoreArgs = list(basis="trunc.poly"), SIMPLIFY=FALSE)
    X3 <- lapply(getX, function(xx) xx$X)
    X3 <- do.call(cbind, X3)

    # remove trailing
    np_nms <- unlist(xvars_npargs)[xvars_basis == "trunc.poly"]

    if (!is.null(X3)) {
      cnms <- mapply(function(nms, maxdegree) {
        paste0(nms, 1:maxdegree)
      }, nms = as.list(np_nms), maxdegre=as.list(xvars_np_truncpoly_degrees), SIMPLIFY=FALSE)

      colnames(X3) <- unlist(cnms)
    }

  }
  
  # remove dup columns in xmeans
  if (!is.null(X1b) & !is.null(X2)) {
    cnX1b <- colnames(X1b)
    cnX2 <- colnames(X2)
    if(all(cnX1b %in% cnX2)) {
      X1b <- NULL
    }
  }
  

  Xmeans <- cbind(X1, X1b, X2, X3)
  if (has_intercept) {
    Xmeans <- cbind(1, Xmeans)
    colnames(Xmeans)[1] <- betanms[1]
  }

  Xmeans
}

# create a single X dataset for smooth plotting
create_single_smooth_data <- function(xnm, Xorig, Xmeans, xvars_static, xvars_np,
                                      xvars_npargs, xvars_basis,
                                      xvars_bivariate,
                                      npdegree, random_intercept, nvals, knots, zvars,
                                      betanms, unms, betavals, uvals, linkname,
                                      applylink, rg, probs, names_y, simresults,
                                      interval = "simulation", multresponse, mcmcres) {

  
  # classify type of variable
  if (length(xnm) == 2) {
    type <- "bivariate"
  } else if (xnm %in% xvars_static) {
    type <- "static"
  } else {
    which_np_var <- which(xvars_np %in% xnm)
    type <- xvars_basis[which_np_var]
    degree <- npdegree[which_np_var]
    truexnm <- unlist(xvars_npargs)[which_np_var]
  }
  
  xplotvals <- Xorig[, xnm]
  xplotvals <- seq(min(xplotvals), max(xplotvals), length.out=nvals)
  Xplot <- Xmeans
  yplotvals <- 0
  
  # process X data based on type
  if (type %in% c("static", "tps")) {
    Xplot[, xnm] <- xplotvals
  } else if (type == "trunc.poly") {
    npow <- 1:degree
    Xtpoly <- matrix(xplotvals^(rep(npow, each=length(xplotvals))),ncol=degree,byrow=FALSE)
    colnames(Xtpoly) <- paste0(truexnm, 1:degree)
    Xplot[, which(colnames(Xplot) %in%  colnames(Xtpoly))] <- Xtpoly
  }
  
  Zplot <- NULL
  if (length(xvars_np) > 0) {
    
    univariate <- sapply(xvars_npargs, length) == 1
    knots_univ <- knots[univariate]
    xvars_basis_univ <- xvars_basis[univariate]
    npdegree_univ <- npdegree[univariate]
    
    Xnpdat <- as.data.frame(Xplot[, xvars_np])
    
    npres <- mapply(np,
                    x1 = as.list(Xnpdat),
                    knots = knots_univ,
                    basis = as.list(xvars_basis_univ),
                    degree = as.list(npdegree_univ), SIMPLIFY=FALSE)
    Zplot <- do.call(cbind, lapply(npres, function(xx) xx$Z))
  }
  
  Zplot2 <- NULL

  if (length(xvars_bivariate) > 0) {
    which_bivariate <- sapply(xvars_npargs, length) == 2

    knots_bivariate <- knots[which_bivariate]

    Xbivdatlst <- lapply(xvars_bivariate, function(xx) {
      Xorig[, xx]
    })

    npres <- mapply(create_bivariate_smooth,
                    X = Xbivdatlst,
                    knots = knots_bivariate,
                    MoreArgs=list(nvals=nvals), SIMPLIFY=F)

    Zplot2 <- do.call(cbind, lapply(npres, function(xx) xx$Z))
    Xplotbiv <- do.call(cbind, lapply(npres, function(xx) xx$X))

    xplotvals <- Xplotbiv[, 2]
    yplotvals <- Xplotbiv[, 3]

    which_bivariate <- sapply(xvars_npargs, length) == 2
    xvars_npargs_biv <- xvars_npargs[which_bivariate]
    
    bivcols <- which(colnames(Xplot) %in% unlist(xvars_npargs_biv))
    
    
    if (type == 'bivariate') {
      Xplot <- matrix(rep(colMeans(Xplot), nrow(Xplotbiv)), byrow=TRUE,
                      nrow=nrow(Xplotbiv))
      Xplot[, bivcols] <- cbind(xplotvals, yplotvals)
    } else {
      
      Xplot[, bivcols] <- matrix(colMeans(cbind(xplotvals, yplotvals)), 
                       nrow=nrow(Xplot), ncol=length(bivcols), byrow=TRUE)
      Zplot2 <- matrix(rep(colMeans(Zplot2), nvals), byrow = TRUE, nrow=nvals)
    }

  }
  
  # process all Z
  Z1 <- NULL
  if (random_intercept & type != 'bivariate') {
    Z1 <- matrix(0, nrow=nvals, ncol=zvars[1])
  } else if (random_intercept & type == 'bivariate') {
    Z1 <- matrix(0, nrow=nrow(Zplot2), ncol=zvars[1])
  } 
  
  if (type == "bivariate" & !is.null(Zplot)) {

    zvals <- colMeans(Zplot)
    Zplot <- matrix(rep(zvals, nrow(Zplot2)), nrow=nrow(Zplot2), 
                    byrow=TRUE)
  }

  # include random intercepts if needed
  # Z1: random intercept Z
  # Zplot: univariate smoothing Z
  # Zplot2: bivariate Z
  Zplot <- cbind(Z1, Zplot, Zplot2)

  # fitted and lower/upper for y
  Cgrid <- cbind(Xplot, Zplot)

  if (applylink) {
    linkinfo <- make.link(linkname)
  } else {
    linkinfo <- make.link("identity")
  }
  
  invlink <- linkinfo$linkinv

  if (interval == "simulation") {
    simvals <- getSamples(simresults, nsamp=nrow(Cgrid), results=mcmcres)

    #------------------------------------------------------------------
    #
    # multresponse
    #
    #------------------------------------------------------------------

    if (multresponse) {
      num_y <- length(names_y)
      nms_simvals <- colnames(simvals)

      var_nums <- sapply(nms_simvals, get_multresponse_ynums)

      var_nums_uniq <- sort(unique(var_nums))
      
      # split simvals matrix (X,Z)
      simvals_lst <- lapply(var_nums_uniq,
                            function(xx, dat=simvals) {
        cnums <- which(var_nums == xx)

        # subset for response
        dat_subs <- dat[, cnums]
        
        # get corresponding y
        yfitsim_subs <- Cgrid %*% t(dat_subs)
        yfit_subs <- invlink(apply(t(yfitsim_subs), 2, quantile, probs=0.50))
        ylower_subs <- invlink(apply(t(yfitsim_subs), 2, quantile, probs=probs[1]))
        yupper_subs <- invlink(apply(t(yfitsim_subs), 2, quantile, probs=probs[2]))

        list(dat_subs = dat_subs,
             yfitsim_subs = yfitsim_subs,
             yfit_subs = yfit_subs,
             ylower_subs = ylower_subs,
             yupper_subs = yupper_subs)

      })

      # get simvals, yfit values for plotting
      yfitsim_lst <- lapply(simvals_lst, function(xx) xx$yfitsim_subs)
      yfitsim <- do.call(rbind, yfitsim_lst)

      yfit <- sapply(simvals_lst, function(xx) xx$yfit_subs)
      ylower <- sapply(simvals_lst, function(xx) xx$ylower_subs)
      yupper <- sapply(simvals_lst, function(xx) xx$yupper_subs)

    } else {
      yfitsim <- Cgrid %*% t(simvals)
      yfit <- invlink(apply(t(yfitsim), 2, quantile, probs=0.50))
      ylower <- invlink(apply(t(yfitsim), 2, quantile, probs=probs[1]))
      yupper <- invlink(apply(t(yfitsim), 2, quantile, probs=probs[2]))
    }



  } else {
    stop("invalid interval parameter")
  }

  colnames(Cgrid) <- c(betanms, unms)

  if (type == "bivariate")
    xnm <- paste(rev(xnm), collapse=' ~ ')

  retdf <- data.frame(y = yfit,
                    ylower = ylower,
                    yupper = yupper,
                    Cgrid,
                    xplotvals = xplotvals,
                    yplotvals = yplotvals,
                    grouping = xnm,
                    type = type)

  if (!multresponse)
    colnames(retdf)[1] <- names_y

  return(retdf)

}

# process bivariate
create_bivariate_smooth <-function(X, knots, nvals) {


  x1 <- X[, 1]; x2 <- X[, 2]
  num_knots <- nrow(knots)

  xvals <- seq(min(x1), max(x1), length=nvals)
  yvals <- seq(min(x2), max(x2), length=nvals)

  # create design matrices
  Xnew <- expand.grid(int = 1,
                      xvals=xvals,
                      yvals=yvals)
  Xfit <- as.matrix(Xnew)

  biv <- create_bivariate_design(X1 = Xfit[, 2],
                                 X2 = Xfit[, 3],
                                 knots = knots)

  # random effects design matrix
  retval <- list(X = Xfit,
                 Z = biv$Z)
  return(retval)

}

# plot predictions
ypred_plot <- function(object, lower.prob=0.025, upper.prob=0.975, ...) {
  
  ypredlst <- object@pp
  ynms <- object@model@names_y
  
  # create plot object
  highlow <- lapply(ypredlst, function(xx) {
    apply(xx, 2, quantile, probs=c(lower.prob, upper.prob))
  })
  
  ymed <- lapply(ypredlst, function(xx) {
    apply(xx, 2, median)
  })
  
  p1all <- mapply(function(interval, pointest, yname) {
    x <- 1:ncol(interval)
    pdata <- data.frame(x=x, 
                        pointest=as.numeric(pointest), 
                        lower=interval[1, ], 
                        upper = interval[2, ])
    
    p1 <- ggplot2::ggplot(pdata, 
                          mapping=ggplot2::aes(x=x, y=pointest))
    
    p1 <- p1 + ggplot2::xlab("") + ggplot2::ylab(yname)
    
    p1 <- p1 + ggplot2::geom_ribbon(data=pdata, 
                                    mapping=ggplot2::aes_string(ymin="lower",
                                                 ymax="upper"),
                                    alpha=0.2)
    p1 <- p1 + ggplot2::geom_line(data=pdata,
                                  colour="blue",
                                  mapping=ggplot2::aes(x=x, y=pointest))
    
    p1 <- p1 + ggplot2::theme_bw()
  }, interval=highlow, pointest=ymed, yname=ynms, SIMPLIFY=FALSE)
  
  n <- length(p1all)
  nCol <- floor(sqrt(n))
  do.call(gridExtra::grid.arrange, c(p1all, ncol=nCol))

  
}


