# FUNCTION FOR LOOKING AT POSSIBLE VARIABLE STRATA

#' Stratification function
#'
#' This function assists in evaluating whether the supplied variable is useful
#' for stratification when fitting a cox proportional hazards model.
#' @param x variable that may be used for stratification, can be categorical or
#'   continuous.
#' @param survival.time vector of survival time corresponding to input vector x.
#' @param split specifies how to split a continuous variable. Default is median
#'   value.
#' @return Generates a plot and table. Table displays the quartiles of the
#'   groups of x. A boxplot is also generated to display the distributions of
#'   the groups in x visually.
#' @keywords gradient boosting
#' @export
#' @examples
#' data <- simulate_survival_cox(true_beta=c(1,1,1,1,1,0,0,0,0,0))
#' strata.boosting(data$strata_idx, data$time)
#' 
strata.boosting <- function(x, survival.time, split="median"){
  data <- data.frame(x, survival.time)
  
  give.n <- function(x){
    return(c(y = mean(x), label = length(x)))
  }
  
  if(is.vector(x) & !is.list(x)){ # check that x is a vector - just one variable
    if(is.factor(x)){
      p <- ggplot2::ggplot(data, ggplot2::aes(x=x, y=survival.time)) + ggplot2::geom_boxplot() + ggplot2::theme_bw() +
        ggplot2::stat_summary(fun.data = give.n, geom = "text") + ggplot2::ylab("Survival Time") + ggplot2::xlab(NULL)
      print(p)
      
      # print summary table
      plyr::ddply(data, plyr::.(as.factor(x)), plyr::summarise,  Min=stats::quantile(survival.time, 0), Q1=stats::quantile(survival.time, 0.25),
                  Median=stats::quantile(survival.time, 0.5), Q3=stats::quantile(survival.time, 0.75), Max=stats::quantile(survival.time, 1))
    }
    else if(length(unique(x))<10){ # limit is 10 strata
      p <- ggplot2::ggplot(data, ggplot2::aes(x=as.factor(x), y=survival.time)) + ggplot2::geom_boxplot() + ggplot2::theme_bw()+
        ggplot2::stat_summary(fun.data = give.n, geom = "text") + ggplot2::ylab("Survival Time") + ggplot2::xlab(NULL)
      print(p)
      
      # print summary table of quartiles
      plyr::ddply(data, plyr::.(as.factor(x)), plyr::summarise,  Min=stats::quantile(survival.time, 0), Q1=stats::quantile(survival.time, 0.25),
                  Median=stats::quantile(survival.time, 0.5), Q3=stats::quantile(survival.time, 0.75), Max=stats::quantile(survival.time, 1))
      
    }
    else{ # continuous case
      if(split=="median"){
        med <- stats::median(data$x)
        x.split <- rep(0,length(x))
        x.split[which(x > med)] <- 1
        p <- ggplot2::ggplot(data, ggplot2::aes(x=as.factor(x.split), y=survival.time)) + ggplot2::geom_boxplot() + ggplot2::theme_bw()+
          ggplot2::stat_summary(fun.data = give.n, geom = "text") + ggplot2::ylab("Survival Time") + ggplot2::xlab(NULL) +
          ggplot2::scale_x_discrete(labels=c("0" = "< Median", "1" = "> Median"))
        print(p)
        plyr::ddply(data, plyr::.(as.factor(x.split)), plyr::summarise,  Min=stats::quantile(survival.time, 0), Q1=stats::quantile(survival.time, 0.25),
                    Median=stats::quantile(survival.time, 0.5), Q3=stats::quantile(survival.time, 0.75), Max=stats::quantile(survival.time, 1))
      }
    }
  }
  # suppose that x is a vector of variables - may want a combination
  
}