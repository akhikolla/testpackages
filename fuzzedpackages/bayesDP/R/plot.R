#' @title bdpnormal Object Plot
#' @description \code{plot} method for class \code{bdpnormal}.
#' @param x object of class \code{bdpnormal}. The result of a call to the
#'   \code{\link{bdpnormal}} function.
#' @param type character. Optional. Select plot type to print.
#'   Supports the following: "discount" gives the discount function;
#'   "posteriors" gives the posterior plots of the event rates; and
#'   "density" gives the augmented posterior density plot(s).  Leave NULL to
#'   display all plots in sequence.
#'
#' @param print logical. Optional. Produce a plot (\code{print = TRUE}; default) or
#'   return a ggplot object (\code{print = FALSE}). When \code{print = FALSE},
#'   it is possible to use \code{ggplot2} syntax to modify the plot appearance.
#'
#' @details The \code{plot} method for \code{bdpnormal} returns up to three
#'   plots: (1) posterior density curves; (2) posterior density of the effect
#'   of interest; and (3) the discount function. A call to \code{plot} that
#'   omits the \code{type} input returns one plot at a time and prompts the
#'   user to click the plot or press return to proceed to the next plot.
#'   Otherwise, the user can specify a plot \code{type} to display the
#'   requested plot.
#'
#' @import methods
#' @importFrom utils head
#' @importFrom ggplot2 aes_string ggtitle ylim guides guide_legend theme element_blank
#' @importFrom graphics par
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("plot", signature(x = "bdpnormal"), function(x, type=NULL, print=TRUE){
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  discount_function   <- x$args1$discount_function
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c
  N_t                 <- x$args1$N_t
  N_c                 <- x$args1$N_c
  arm2                <- x$args1$arm2
  method              <- x$args1$method

  information_sources <- NULL

  ### Create data frames for plotting via ggplot2
  D1 <- D2 <- D3 <- D5 <- D6 <- NULL

  if (arm2){
    dens1 <- density(posterior_control$posterior_mu,
                     adjust=0.5)
    D1 <- data.frame(information_sources="Posterior",
                     group="Control",
                     x=dens1$x,
                     y=dens1$y)

    if(!is.null(posterior_control$posterior_flat_mu)){
      dens2 <- density(posterior_control$posterior_flat_mu,adjust=0.5)
      D2 <- data.frame(information_sources="Current Data",
                       group="Control",
                       x=dens2$x,
                       y=dens2$y)
    }

    if(!is.null(posterior_control$prior_mu)){
      dens3 <- density(posterior_control$prior_mu,adjust=0.5)
      D3 <- data.frame(information_sources="Historical Data",
                       group="Control",
                       x=dens3$x,
                       y=dens3$y)
    }
  }


  dens4 <- density(posterior_treatment$posterior_mu,adjust=0.5)
  D4 <- data.frame(information_sources="Posterior",
                   group="Treatment",
                   x=dens4$x,
                   y=dens4$y)


  if(!is.null(posterior_treatment$posterior_flat_mu)){
    dens5 <- density(posterior_treatment$posterior_flat_mu,adjust=0.5)
    D5 <- data.frame(information_sources="Current Data",
                     group="Treatment",
                     x=dens5$x,
                     y=dens5$y)
  }

  if(!is.null(posterior_treatment$prior_mu)){
    dens6 <- density(posterior_treatment$prior_mu,adjust=0.5)
    D6 <- data.frame(information_sources="Historical Data",
                     group="Treatment",
                     x=dens6$x,
                     y=dens6$y)
  }

  D <- rbind(D1,D2,D3,D4,D5,D6)

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior","Current Data","Historical Data")))
  D$group <- factor(D$group, levels=c("Treatment", "Control"))

  ##############################################################################
  ### Posterior Type Plots
  ##############################################################################
  post_typeplot <- ggplot(D,aes_string(x="x",y="y")) +
    geom_line(size=2,
              aes_string(colour="information_sources",
                         lty="information_sources")) +
    facet_wrap(~group, ncol=1,scales='free') +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Posterior Type Plot") +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())


  ##############################################################################
  ### Density Plots
  ##############################################################################
  densityplot <- ggplot(subset(D,information_sources=="Posterior"),
                        aes_string(x="x",y="y")) +
    geom_line(size=2,aes_string(colour="group")) +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Density Plot") +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())


  ##############################################################################
  ### Discount function plot
  ### - Only makes sense to plot if Current/historical treatment are present or
  ###   both current/historical control are present
  ##############################################################################
  p_hat <- seq(0,1,length.out=100)

  discountfun_plot <- NULL

  if(!is.null(N0_t) & !is.null(N_t)){

    ### Create discount function based on distribution type
    if(discount_function == "weibull"){
      discount_function_treatment <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[1],
                                              scale=x$args1$weibull_scale[1])
    } else if(discount_function == "scaledweibull"){
      max_p <- pweibull(1,
                        shape=x$args1$weibull_shape[1],
                        scale=x$args1$weibull_scale[1])

      discount_function_treatment <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[1],
                                              scale=x$args1$weibull_scale[1])/max_p
    } else if(discount_function == "identity"){
      discount_function_treatment <- p_hat
    }

    D1 <- data.frame(group = "Treatment",
                     y     = discount_function_treatment,
                     x     = seq(0,1,length.out=100))
    D2 <- data.frame(group = c("Treatment"),
                     p_hat = median(posterior_treatment$p_hat))
    D3 <- data.frame(group = c("Treatment"),
                     p_hat = median(posterior_treatment$alpha_discount))

    discountfun_plot <- ggplot() +
      geom_line(data=D1,aes_string(y="y",x="x",color="group"),size=1) +
      geom_vline(data=D2, aes_string(xintercept="p_hat",color="group"),lty=2) +
      geom_hline(data=D3, aes_string(yintercept="p_hat",color="group"),lty=2)
  }


  if(arm2){
    if(!is.null(N0_c) & !is.null(N_c)){
      if(is.null(discountfun_plot)){
        discountfun_plot <- ggplot()
      }

      ### Create discount function based on distribution type
      if(discount_function == "weibull"){
        discount_function_control <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[2],
                                              scale=x$args1$weibull_scale[2])
      } else if(discount_function == "scaledweibull"){
        max_p <- pweibull(1,
                          shape=x$args1$weibull_shape[2],
                          scale=x$args1$weibull_scale[2])

        discount_function_control <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[2],
                                              scale=x$args1$weibull_scale[2])/max_p
      } else if(discount_function == "identity"){
        discount_function_control <- p_hat
      }

      D4 <- data.frame(group = "Control",
                       y     = discount_function_control,
                       x     = seq(0,1,length.out=100))
      D5 <- data.frame(group="Control", p_hat=median(posterior_control$p_hat))
      D6 <- data.frame(group="Control", p_hat=median(posterior_control$alpha_discount))

      discountfun_plot  <- discountfun_plot +
        geom_line(data=D4,aes_string(y="y",x="x",color="group"),size=1) +
        geom_vline(data=D5,aes_string(xintercept="p_hat",color="group"),lty=2) +
        geom_hline(data=D6,aes_string(yintercept="p_hat",color="group"),lty=2)
    }
  }


  if(!is.null(discountfun_plot)){
    discountfun_plot <- discountfun_plot +
      facet_wrap(~group, ncol=1) +
      theme_bw() +
      ylab("Alpha Discount Value") +
      xlab("Stochastic comparison (Current vs Historical Data)") +
      ggtitle("Discount Function") +
      ylim(0,1) +
      guides(fill=guide_legend(title=NULL)) +
      theme(legend.title=element_blank())
  }


  if(print){
    if(is.null(type)){
      op <- par(ask=TRUE)
      plot(post_typeplot)
      plot(densityplot)
      if(!is.null(discountfun_plot)){
        plot(discountfun_plot)
      }
      par(op)
    } else if (type=="discount"){
      if(!is.null(discountfun_plot)){
        plot(discountfun_plot)
      }
    } else if (type=="posteriors"){
      plot(post_typeplot)
    } else if (type=="density"){
      plot(densityplot)
    }
  } else{
    if(is.null(type)){
      out <- list(posteriors = post_typeplot,
                  density = densityplot)
      if(!is.null(discountfun_plot)){
        out <- c(out, discount = list(discountfun_plot))
        out
      } else{
        out
      }
    } else if (type=="discount"){
      if(!is.null(discountfun_plot)){
        discountfun_plot
      }
    } else if (type=="posteriors"){
      post_typeplot
    } else if (type=="density"){
      densityplot
    }
  }
})


#' @title bdpbinomial Object Plot
#' @description \code{plot} method for class \code{bdpbinomial}.
#' @param x object of class \code{bdpbinomial}. The result of a call to the
#'   \code{\link{bdpbinomial}} function.
#' @param type character. Optional. Select plot type to print.
#'   Supports the following: "discount" gives the discount function;
#'   "posteriors" gives the posterior plots of the event rates; and
#'   "density" gives the augmented posterior density plot(s).  Leave NULL to
#'   display all plots in sequence.
#'
#' @param print logical. Optional. Produce a plot (\code{print = TRUE}; default) or
#'   return a ggplot object (\code{print = FALSE}). When \code{print = FALSE},
#'   it is possible to use \code{ggplot2} syntax to modify the plot appearance.
#'
#' @details The \code{plot} method for \code{bdpbinomial} returns up to three
#'   plots: (1) posterior density curves; (2) posterior density of the effect
#'   of interest; and (3) the discount function. A call to \code{plot} that
#'   omits the \code{type} input returns one plot at a time and prompts the
#'   user to click the plot or press return to proceed to the next plot.
#'   Otherwise, the user can specify a plot \code{type} to display the
#'   requested plot.
#'
#' @import methods
#' @importFrom utils head
#' @importFrom ggplot2 aes_string ggtitle ylim guides guide_legend theme element_blank
#' @importFrom graphics par
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("plot", signature(x = "bdpbinomial"), function(x, type=NULL, print=TRUE){
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  discount_function   <- x$args1$discount_function
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c
  N_t                 <- x$args1$N_t
  N_c                 <- x$args1$N_c
  arm2                <- x$args1$arm2
  method              <- x$args1$method

  information_sources <- NULL

  ### Create data frames for plotting via ggplot2
  D1 <- D2 <- D3 <- D5 <- D6 <- NULL

  if (arm2){
    dens1 <- density(posterior_control$posterior,
                     adjust=0.5)
    D1 <- data.frame(information_sources="Posterior",
                     group="Control",
                     x=dens1$x,
                     y=dens1$y)

    if(!is.null(posterior_control$posterior_flat)){
      dens2 <- density(posterior_control$posterior_flat,adjust=0.5)
      D2 <- data.frame(information_sources="Current Data",
                       group="Control",
                       x=dens2$x,
                       y=dens2$y)
    }

    if(!is.null(posterior_control$prior)){
      dens3 <- density(posterior_control$prior,adjust=0.5)
      D3 <- data.frame(information_sources="Historical Data",
                       group="Control",
                       x=dens3$x,
                       y=dens3$y)
    }
  }


  dens4 <- density(posterior_treatment$posterior,adjust=0.5)
  D4 <- data.frame(information_sources="Posterior",
                   group="Treatment",
                   x=dens4$x,
                   y=dens4$y)


  if(!is.null(posterior_treatment$posterior_flat)){
    dens5 <- density(posterior_treatment$posterior_flat,adjust=0.5)
    D5 <- data.frame(information_sources="Current Data",
                     group="Treatment",
                     x=dens5$x,
                     y=dens5$y)
  }

  if(!is.null(posterior_treatment$prior)){
    dens6 <- density(posterior_treatment$prior,adjust=0.5)
    D6 <- data.frame(information_sources="Historical Data",
                     group="Treatment",
                     x=dens6$x,
                     y=dens6$y)
  }

  D <- rbind(D1,D2,D3,D4,D5,D6)

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior","Current Data","Historical Data")))
  D$group <- factor(D$group, levels=c("Treatment", "Control"))

  ##############################################################################
  ### Posterior Type Plots
  ##############################################################################
  post_typeplot <- ggplot(D,aes_string(x="x",y="y")) +
    geom_line(size=2,aes_string(color="information_sources",lty="information_sources")) +
    theme_bw() +
    facet_wrap(~group, ncol=1,scales='free') +
    ylab("Density (PDF)") +
    xlab("Values") +
    ggtitle("Posterior Type Plot") +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())


  ##############################################################################
  ### Density Plots
  ##############################################################################
  densityplot <- ggplot(subset(D,information_sources=="Posterior"),
                        aes_string(x="x",y="y")) +
    geom_line(size=2,aes_string(color="group")) +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Density Plot") +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())


  ##############################################################################
  ### Discount function plot
  ### - Only makes sense to plot if Current/historical treatment are present or
  ###   both current/historical control are present
  ##############################################################################
  p_hat <- seq(0,1,length.out=100)

  discountfun_plot <- NULL

  if(!is.null(N0_t) & !is.null(N_t)){
    ### Create discount function based on distribution type
    if(discount_function == "weibull"){
      discount_function_treatment <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[1],
                                              scale=x$args1$weibull_scale[1])
    } else if(discount_function == "scaledweibull"){
      max_p <- pweibull(1,
                        shape=x$args1$weibull_shape[1],
                        scale=x$args1$weibull_scale[1])

      discount_function_treatment <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[1],
                                              scale=x$args1$weibull_scale[1])/max_p
    } else if(discount_function == "identity"){
      discount_function_treatment <- p_hat
    }

    D1 <- data.frame(group = "Treatment",
                     y     = discount_function_treatment,
                     x     = seq(0,1,length.out=100))
    D2 <- data.frame(group="Treatment", p_hat=median(posterior_treatment$p_hat))
    D3 <- data.frame(group="Treatment", p_hat=median(posterior_treatment$alpha_discount))

    discountfun_plot <- ggplot() +
      geom_line(data=D1,aes_string(y="y",x="x",color="group"),size=1) +
      geom_vline(data=D2, aes_string(xintercept="p_hat",color="group"),lty=2) +
      geom_hline(data=D3, aes_string(yintercept="p_hat",color="group"),lty=2)
  }


  if(arm2){
    if(!is.null(N0_c) & !is.null(N_c)){

      if(is.null(discountfun_plot)){
        discountfun_plot <- ggplot()
      }

      ### Create discount function based on distribution type
      if(discount_function == "weibull"){
        discount_function_control <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[2],
                                              scale=x$args1$weibull_scale[2])
      } else if(discount_function == "scaledweibull"){
        max_p <- pweibull(1,
                          shape=x$args1$weibull_shape[2],
                          scale=x$args1$weibull_scale[2])

        discount_function_control <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[2],
                                              scale=x$args1$weibull_scale[2])/max_p
      } else if(discount_function == "identity"){
        discount_function_control <- p_hat
      }

      D4 <- data.frame(group = "Control",
                       y     = discount_function_control,
                       x     = seq(0,1,length.out=100))
      D5 <- data.frame(group="Control",p_hat=median(posterior_control$p_hat))
      D6 <- data.frame(group="Control",p_hat=median(posterior_control$alpha_discount))

      discountfun_plot  <- discountfun_plot +
        geom_line(data=D4,aes_string(y="y",x="x",color="group"),size=1) +
        geom_vline(data=D5,aes_string(xintercept="p_hat",color="group"),lty=2) +
        geom_hline(data=D6,aes_string(yintercept="p_hat",color="group"),lty=2)
    }
  }


  if(!is.null(discountfun_plot)){
    discountfun_plot <- discountfun_plot +
      facet_wrap(~group, ncol=1) +
      theme_bw() +
      ylab("Alpha Discount Value") +
      xlab("Stochastic comparison (Current vs Historical Data)") +
      ggtitle("Discount Function") +
      ylim(0,1) +
      guides(fill=guide_legend(title=NULL)) +
      theme(legend.title=element_blank())
  }


  if(print){
    if(is.null(type)){
      op <- par(ask=TRUE)
      plot(post_typeplot)
      plot(densityplot)
      if(!is.null(discountfun_plot)){
        plot(discountfun_plot)
      }
      par(op)
    } else if (type=="discount"){
      if(!is.null(discountfun_plot)){
        plot(discountfun_plot)
      }
    } else if (type=="posteriors"){
      plot(post_typeplot)
    } else if (type=="density"){
      plot(densityplot)
    }
  } else{
    if(is.null(type)){
      out <- list(posteriors = post_typeplot,
                  density = densityplot)
      if(!is.null(discountfun_plot)){
        out <- c(out, discount = list(discountfun_plot))
        out
      } else{
        out
      }
    } else if (type=="discount"){
      if(!is.null(discountfun_plot)){
        discountfun_plot
      }
    } else if (type=="posteriors"){
      post_typeplot
    } else if (type=="density"){
      densityplot
    }
  }


})


#' @title bdpsurvival Object Plot
#' @description \code{plot} method for class \code{bdpsurvival}.
#' @param x object of class \code{bdpsurvival}. The result of a call to the
#'   \code{\link{bdpsurvival}} function.
#' @param type character. Optional. Select plot type to print.
#'   Supports the following: "discount" gives the discount function and
#'   "survival" gives the survival curves.  Leave NULL to
#'   display all plots in sequence.
#'
#' @param print logical. Optional. Produce a plot (\code{print = TRUE}; default) or
#'   return a ggplot object (\code{print = FALSE}). When \code{print = FALSE},
#'   it is possible to use \code{ggplot2} syntax to modify the plot appearance.
#'
#' @details The \code{plot} method for \code{bdpsurvival} returns up to two
#'   plots: (1) posterior survival curves and (2) the discount function. A
#'   call to \code{plot} that omits the \code{type} input returns one plot
#'   at a time and prompts the user to click the plot or press return to
#'   proceed to the next plot. Otherwise, the user can specify a plot
#'   \code{type} to display the requested plot.
#'
#' @import methods
#' @importFrom utils head
#' @importFrom ggplot2 aes_string ggtitle ylim guides guide_legend theme element_blank
#' @importFrom graphics par
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("plot", signature(x = "bdpsurvival"), function(x, type=NULL, print=TRUE){
  args1               <- x$args1
  discount_function   <- x$args1$discount_function
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  data                <- args1$data
  breaks              <- args1$breaks
  arm2                <- args1$arm2
  method              <- args1$method

  ##############################################################################
  ### Survival curve(s)
  ### - Only computed for one-arm trial
  ##############################################################################
  ### Organize data for current treatment
  time_t  <- sort(unique(args1$S_t$time))
  survival_times_posterior_flat <- lapply(time_t, ppexp,
    posterior_treatment$posterior_flat_hazard, cuts=c(0,breaks))
  survival_median_posterior_flat <- 1-sapply(survival_times_posterior_flat, median)

  D1 <- data.frame(source  = "Current Data",
                   group = "Treatment",
                   x      = time_t,
                   y      = survival_median_posterior_flat)

  ### Organize data for historical treatment
  if(!is.null(args1$S0_t)){
    time0_t <- sort(unique(args1$S0_t$time))
    survival_times_prior <- lapply(time0_t, ppexp,
      posterior_treatment$prior_hazard, cuts=c(0,breaks))
    survival_median_prior <- 1-sapply(survival_times_prior, median)

    D2 <- data.frame(source = "Historical Data",
                     group  = "Treatment",
                     x      = time0_t,
                     y      = survival_median_prior)
  } else{
    D2 <- NULL
  }

  ### Organize data for treatment posterior
  survival_times_posterior  <- lapply(time_t, ppexp,posterior_treatment$posterior_hazard,cuts=c(0,breaks))
  survival_median_posterior <- 1-sapply(survival_times_posterior, median)

  D3 <- data.frame(source = "Posterior",
                   group  = "Treatment",
                   x      = time_t,
                   y      = survival_median_posterior)

  ### Organize data for current control
  if(!is.null(args1$S_c)){
    time_c               <- sort(unique(args1$S_c$time))
    survival_times_posterior_flat <- lapply(time_c, ppexp,
                                           posterior_control$posterior_flat_hazard, cuts=c(0,breaks))
    survival_median_posterior_flat <- 1-sapply(survival_times_posterior_flat, median)

    D4 <- data.frame(source = "Current Data",
                     group  = "Control",
                     x      = time_c,
                     y      = survival_median_posterior_flat)
  } else{
    D4 <- NULL
  }

  ### Organize data for historical control
  if(!is.null(args1$S0_c)){
    time0_c               <- sort(unique(args1$S0_c$time))
    survival_times_prior <- lapply(time0_c, ppexp,
                                   posterior_control$prior_hazard, cuts=c(0,breaks))
    survival_median_prior <- 1-sapply(survival_times_prior, median)

    D5 <- data.frame(source = "Historical Data",
                     group  = "Control",
                     x      = time0_c,
                     y      = survival_median_prior)
  } else{
    D5 <- NULL
  }

  ### Organize data for control posterior
  if(!is.null(args1$S_c) & !is.null(args1$S0_c)){
    survival_times_posterior  <- lapply(time_c, ppexp,posterior_control$posterior_hazard,cuts=c(0,breaks))
    survival_median_posterior <- 1-sapply(survival_times_posterior, median)

    D6 <- data.frame(source = "Posterior",
                     group  = "Control",
                     x      = time_c,
                     y      = survival_median_posterior)
  } else if(is.null(args1$S_c) & !is.null(args1$S0_c)){
    survival_times_posterior  <- lapply(time0_c, ppexp,posterior_control$posterior_hazard,cuts=c(0,breaks))
    survival_median_posterior <- 1-sapply(survival_times_posterior, median)

    D6 <- data.frame(source = "Posterior",
                     group  = "Control",
                     x      = time0_c,
                     y      = survival_median_posterior)
  } else if(!is.null(args1$S_c) & is.null(args1$S0_c)){
    survival_times_posterior  <- lapply(time_c, ppexp,posterior_control$posterior_hazard,cuts=c(0,breaks))
    survival_median_posterior <- 1-sapply(survival_times_posterior, median)

    D6 <- data.frame(source = "Posterior",
                     group  = "Control",
                     x      = time_c,
                     y      = survival_median_posterior)
  } else{
    D6 <- NULL
  }

  D       <- rbind(D4,D5,D6, D1,D2,D3)
  D$group <- factor(D$group, levels=c("Treatment", "Control"))

  ### Plot survival curve
  survival_curves <- ggplot(D, aes_string(x="x",y="y")) +
    geom_line(size=1.4, aes_string(color="source", lty="source")) +
    facet_wrap(~group, ncol=1,scales='free') +
    ylab("Survival probability") +
    xlab("Time") +
    theme_bw() +
    ggtitle("Survival Curve(s)") +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())


  ##############################################################################
  ### Discount function plot
  ### - Only makes sense to plot if Current/historical treatment are present or
  ###   both current/historical control are present
  ##############################################################################
  p_hat <- seq(0,1,length.out=100)

  discountfun_plot <- NULL

  if(!is.null(args1$S_t) & !is.null(args1$S0_t)){

    ### Create discount function based on distribution type
    if(discount_function == "weibull"){
      discount_function_treatment <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[1],
                                              scale=x$args1$weibull_scale[1])
    } else if(discount_function == "scaledweibull"){
      max_p <- pweibull(1,
                        shape=x$args1$weibull_shape[1],
                        scale=x$args1$weibull_scale[1])

      discount_function_treatment <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[1],
                                              scale=x$args1$weibull_scale[1])/max_p
    } else if(discount_function == "identity"){
      discount_function_treatment <- p_hat
    }


    D1 <- data.frame(group = "Treatment",
                     y     = discount_function_treatment,
                     x     = seq(0,1,length.out=100))
    D2 <- data.frame(group="Treatment", p_hat=c(median(posterior_treatment$p_hat)))
    D3 <- data.frame(group="Treatment", p_hat=c(median(posterior_treatment$alpha_discount)))

    discountfun_plot <- ggplot() +
      geom_line(data=D1,aes_string(y="y",x="x",color="group"),size=1) +
      geom_vline(data=D2, aes_string(xintercept="p_hat",color="group"),lty=2) +
      geom_hline(data=D3, aes_string(yintercept="p_hat",color="group"),lty=2)
  }


  if(arm2){
    if(!is.null(args1$S_c) & !is.null(args1$S0_c)){
      if(is.null(discountfun_plot)){
        discountfun_plot <- ggplot()
      }

      ### Create discount function based on distribution type
      if(discount_function == "weibull"){
        discount_function_control <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[2],
                                              scale=x$args1$weibull_scale[2])
      } else if(discount_function == "scaledweibull"){
        max_p <- pweibull(1,
                          shape=x$args1$weibull_shape[2],
                          scale=x$args1$weibull_scale[2])

        discount_function_control <- pweibull(p_hat,
                                              shape=x$args1$weibull_shape[2],
                                              scale=x$args1$weibull_scale[2])/max_p
      } else if(discount_function == "identity"){
        discount_function_control <- p_hat
      }

      D4 <- data.frame(group = "Control",
                       y     = discount_function_control,
                       x     = seq(0,1,length.out=100))
      D5 <- data.frame(group="Control", p_hat=c(median(posterior_control$p_hat)))
      D6 <- data.frame(group="Control", p_hat=c(median(posterior_control$alpha_discount)))

      discountfun_plot  <- discountfun_plot +
        geom_line(data=D4,aes_string(y="y",x="x",color="group"),size=1) +
        geom_vline(data=D5,aes_string(xintercept="p_hat",color="group"),lty=2) +
        geom_hline(data=D6,aes_string(yintercept="p_hat",color="group"),lty=2)
    }
  }


  if(!is.null(discountfun_plot)){
    discountfun_plot <- discountfun_plot +
      facet_wrap(~group, ncol=1) +
      theme_bw() +
      ylab("Alpha Discount Value") +
      xlab("Stochastic comparison (Current vs Historical Data)") +
      ggtitle("Discount Function") +
      ylim(0,1) +
      guides(fill=guide_legend(title=NULL)) +
      theme(legend.title=element_blank())
  }


  if(print){
    if(is.null(type)){
      op <- par(ask=TRUE)
      plot(survival_curves)
      if(!is.null(args1$S0_t) | (!is.null(args1$S_c) & !is.null(args1$S0_c))){
        plot(discountfun_plot)
      }
      par(op)
    } else if(type=="discount"){
      if(!is.null(args1$S0_t) | (!is.null(args1$S_c) & !is.null(args1$S0_c))){
        plot(discountfun_plot)
      }
    } else if(type=="survival"){
      plot(survival_curves)
    }
  } else{
    if(is.null(type)){
      out <- list(survival = survival_curves)
      if(!is.null(args1$S0_t) | (!is.null(args1$S_c) & !is.null(args1$S0_c))){
        out <- c(out, discount = list(discountfun_plot))
        out
      } else{
        out
      }
    } else if(type=="discount"){
      if(!is.null(args1$S0_t) | (!is.null(args1$S_c) & !is.null(args1$S0_c))){
        discountfun_plot
      }
    } else if(type=="survival"){
      survival_curves
    }
  }
})
