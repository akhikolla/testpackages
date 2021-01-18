#' @title plot_Moran.I
#' @description Implements Global Moran I test to evaluate spatial autocorrelation in a units' risk propensity in the data.  Positive values indicate spatial clustering of similar values.
#'
#' @param data data.
#' @param var_duration variable that measures duration until censoring or failure.
#' @param var_id ID's unique identifier.
#' @param var_time variable that measures time.
#' @param n number of observation per id.
#' @param t value of the confidence interval.
#'
#'
#' @return A ggplot object
#'
#' @examples
#'
#' library(BayesSPsurv)
#' dataw <- spduration::add_duration(data = BayesSPsurv::Walter_2015_JCR,
#'                                   y = "renewed_war",
#'                                   unitID = "ccode",
#'                                   tID = "year",
#'                                   freq = "year",
#'                                   ongoing = FALSE)
#'
#'
#' plot_Moran.I(data = dataw ,
#'             var_duration = "duration",
#'             var_id = "ccode",
#'             var_time = "year",
#'             n = 12)
#'
#'
#' @export


plot_Moran.I <- function(data,
                         var_duration = character(),
                         var_id = character(),
                         var_time = character(),
                         n = 1,
                         t = 1.645){

    wdata <- data_plots(data = data, var_id = var_id, var_time = var_time, n = n)
    w2   <- wdata[[2]]
    mats <- wdata[[1]]
    t_ <- t
    morans <- as.numeric()
    uci    <- as.numeric()
    lci    <- as.numeric()
    for(i in 1:length(w2)){
        m         <- ape::Moran.I(w2[[i]][, var_duration], mats[[i]], scaled = FALSE, na.rm = FALSE, alternative = "two.sided")
        morans[i] <- m$observed
        uci[i]    <- m$observed+(t_*((m$sd)/sqrt(nrow(mats[[i]]))))
        lci[i]    <- m$observed-(t_*((m$sd)/sqrt(nrow(mats[[i]]))))
    }
    year      <- as.numeric(names(mats))
    morandata <- as.data.frame(cbind(year,morans,uci,lci))

    ggplot2::ggplot(morandata, ggplot2::aes(x = year, y = morans)) +
        ggplot2::geom_point(size = 2)  +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::xlab("Year") +
        ggplot2::ylab("Observed Moran's I") +
        ggplot2::geom_errorbar(ggplot2::aes(ymax = uci, ymin = lci))+
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 10), axis.title = ggplot2::element_text(size = 15, face = "bold"))


}


#' @title plot_JointCount
#' @description Uses Joint Count tests to assess spatial clustering or dispersion of categorical variables in the data. Negative values indicate positive spatial clustering.
#'
#' @param data data.
#' @param var_cured binary indicator of immunity.
#' @param var_id ID's unique identifier.
#' @param var_time variable that measures time.
#' @param n number of observation per id.
#' @param t value of the confidence interval.
#'
#' @return A ggplot object
#'
#' @examples
#'
#' library(BayesSPsurv)
#' dataw  <- spduration::add_duration(data = BayesSPsurv::Walter_2015_JCR,
#'                                    y = "renewed_war",
#'                                    unitID = "ccode",
#'                                    tID = "year",
#'                                    freq = "year",
#'                                    ongoing = FALSE)
#'
#'
#' plot_JointCount(data = dataw,
#'                var_cured = "cured",
#'                var_id = "ccode",
#'                var_time = "year",
#'                n = 12)
#'
#'
#' @export

plot_JointCount <- function(data,
                           var_cured = character(),
                           var_id = character(),
                           var_time = character(),
                           n = 1,
                           t = 1.645){


    wdata <- data_plots(data = data, var_id = var_id, var_time = var_time, n = n)
    w2   <- wdata[[2]]
    mats <- wdata[[1]]
    t_ <- t
    jcs  <-as.numeric()
    juci <-as.numeric()
    jlci <-as.numeric()
    exp  <-as.numeric()
    diff <-as.numeric()

    for(i in 1:length(w2)){
        m       <- JointCount(w2[[i]][, var_cured], mats[[i]])
        jcs[i]  <- m$Obs
        exp[i]  <- m$Exp
        d       <- m$Obs-m$Exp
        diff[i] <- d
        juci[i] <- d+(t_*((m$sd)/sqrt(nrow(mats[[i]]))))
        jlci[i] <- d-(t_*((m$sd)/sqrt(nrow(mats[[i]]))))
    }
    year      <- as.numeric(names(mats))
    jcdata    <- as.data.frame(cbind(year,jcs,juci,jlci))

    ggplot2::ggplot(jcdata, ggplot2::aes(x = year, y = diff)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::xlab("Year") +
        ggplot2::ylab("Observed Joint Counts") +
        ggplot2::geom_errorbar(ggplot2::aes(ymax = juci, ymin = jlci)) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 10), axis.title = ggplot2::element_text(size=15, face = "bold"))
}


