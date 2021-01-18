
#' Air pollution data, separately for summer and winter months
#'
#' Air pollution data from Leeds (U.K.) city centre, collected from 1994 to
#' 1998. The \code{summer} data set corresponds to the months of April to July
#' inclusive. The \code{winter} data set corresponds to the months of November
#' to February inclusive. Some outliers have been removed, as discussed by
#' Heffernan and Tawn, 2004.
#'
#'
#' @rdname airPollution
#' @name summer and winter data
#' @aliases summer winter
#' @docType data
#' @format Data frames with 578 (summer) and 532 (winter) observations on the
#' following 5 variables.
#' \describe{
#' \item{O3}{Daily maximum ozone in
#' parts per billion.}
#' \item{NO2}{Daily maximum NO2 in parts per
#' billion.}
#' \item{NO}{Daily maximum NO in parts per billion.}
#' \item{SO2}{Daily maximum SO2 in parts per billion.}
#' \item{PM10}{Daily maximum PM10 in micrograms/metre^3}
#' }
#' @references J. E. Heffernan and J. A. Tawn, A conditional approach for
#' multivariate extreme values, Journal of the Royal Statistical society B, 66,
#' 497 -- 546, 2004
#' @source Provided as online supplementary material to Heffernan and Tawn,
#' 2004:
#'
#' http://www.blackwellpublishing.com/rss/Readmefiles/heffernan.htm
#' @keywords datasets
#' @examples
#'
#' data(summer)
#' data(winter)
#'
NULL





#' Liver related laboratory data
#'
#' Liver related laboratory data from a randomized, blind, parallel group
#' clinical trial with 4 doses of a drug.
#'
#' Dose A is the lowest dose, dose, B the next, C the next, and D the highest
#' dose. The baseline values were taken prior to any treatment being received,
#' and the clinical trial had a single post-baseline visit.
#'
#' @name liver
#' @docType data
#' @usage data(liver)
#' @format A data frame with 606 observations on the following 9 variables.
#' \describe{
#'
#' \item{ALP.B}{Alkaline phosphatase at baseline. A numeric vector.}
#' \item{ALT.B}{Alanine aminotransferase at baseline. A numeric vector.}
#' \item{AST.B}{Aspartate aminotransferase at baseline. A numeric vector.}
#' \item{TBL.B}{Total bilirubin at baseline. A numeric vector.}
#' \item{ALP.M}{Alkaline phosphatase after treatment. A numeric vector.}
#' \item{ALT.M}{Alanine aminotransferase after treatment. A numeric vector.}
#' \item{AST.M}{Aspartate aminotransferase after treatment. A numeric vector.}
#' \item{TBL.M}{Total bilirubin after treatment. A numeric vector.}
#' \item{dose}{The treatment group (i.e. dose group). A factor with levels \code{A} \code{B} \code{C} \code{D}}
#' }
#' @source AstraZeneca data on file.
#' @keywords datasets
NULL


#' Rain, wavesurge, portpirie and nidd datasets.
#'
#' Rainfall, wave-surge, Port Pirie and River Nidd data sets.
#'
#' The rain, wave-surge and Port Pirie datasets are used by Coles and appear in
#' the \code{ismev} package. The River Nidd data appear in the \code{evir}
#' package.
#'
#' @name rain, wavesurge and portpirie
#' @rdname rain
#' @aliases rain wavesurge portpirie nidd
#' @docType data
#' @format The format of the rain data is: num [1:17531] 0 2.3 1.3 6.9 4.6 0 1
#' 1.5 1.8 1.8 ...
#'
#' The wave-surge data is bivariate and is used for testing functions in
#' \code{texmex}.
#'
#' The Port Pirie data has two columns: 'Year' and 'SeaLevel'.
#'
#' The River Nidd data represents 154 measurements of the level of the River
#' Nidd at Hunsingore Weir (Yorkshire, UK) between 1934 and 1969. Each
#' measurement breaches the threshold of $65 m^3/2$. Various authors have
#' analysed this dataset, as described by Papastathopoulos and Tawn~\cite{egp},
#' there being some apparent difficulty in identifying a threshold above which
#' GPD models are suitable.
#' @references S. Coles, An Introduction to Statistical Modeling of Extreme
#' Values, Springer, 2001
#'
#' I. Papastathopoulos and J. A. Tawn, Extended Generalised Pareto Models for
#' Tail Estimation, Journal of Statistical Planning and Inference, 143, 134 --
#' 143, 2011
#' @source Copied from the \code{ismev} package and the \code{evir} package
#' @keywords datasets
NULL



#' Extreme value modelling
#'
#' Extreme values modelling, including the conditional multivariate approach of
#' Heffernan and Tawn (2004).
#'
#' The
#' package was originally called `texmex' for Threshold EXceedances and
#' Multivariate EXtremes. However, it is no longer the case that only threshold
#' excess models are implemented, so the `tex' bit doesn't make sense. So, the
#' package is called `texmex' because it used to be called `texmex'.
#'
#' \code{\link{evm}}: Fit extreme value distributions to data, possibly with
#' covariates. Use maximum likelihood estimation, maximum penalized likelihood
#' estimation, simulate from the posterior distribution or run a parametric
#' bootstrap. Extreme value families include the generalized Pareto
#' distribution (\code{gpd}) and generalized extreme value (\code{gev})
#' distribution.
#'
#' \code{\link{mex}}: Fit multiple, independent generalized Pareto models to
#' the the upper tails of the columns of a data set, and estimate the
#' conditional dependence structure between the columns using the method of
#' Heffernan and Tawn.
#'
#' \code{\link{bootmex}}: Bootstrap estimation for parameters in generalized
#' Pareto models and in the dependence structure.
#'
#' \code{\link{declust}}: Estimation of extremal index and subsequent
#' declustering of dependent sequences using the intervals estimator of Ferro
#' and Segers.
#'
#' @name texmex-package
#' @aliases texmex-package texmex
#' @docType package
#' @author Harry Southworth, Janet E. Heffernan, Paul D. Metcalfe
#'
#' Maintainer: Harry Southworth <harry.southworth@@gmail.com>
#'
#' URL: https://github.com/harrysouthworth/texmex
#' @references J. E. Heffernan and J. A. Tawn, A conditional approach for
#' multivariate extreme values, Journal of the Royal Statistical society B, 66,
#' 497 -- 546, 2004.
#'
#' C.A.T Ferro and J. Segers, Inference for Clusters of Extreme Values, Journal
#' of the Royal Statistical society B, 65, 545 -- 556, 2003.
#' @keywords models multivariate package
#' @examples
#'
#' # Analyse the winter data used by Heffernan and Tawn
#' mymex <- mex(winter, mqu = .7, penalty="none", dqu=.7, which = "NO")
#' plot(mymex)
#' # Only do 10 replicates to keep CRAN checks happy. Do many more in any
#' # real application
#' myboot <- bootmex(mymex, R=10)
#' plot(myboot)
#' mypred <- predict(myboot,  pqu=.95)
#' summary(mypred , probs = c( .025, .5, .975 ))
#'
#' # Analyse the liver data included in the package
#' library(MASS) # For the rlm function
#'
#' liver <- liver[liver$ALP.M > 1,] # Get rid of outlier
#' liver$ndose <- as.numeric(liver$dose)
#'
#' alt <- resid(rlm(log(ALT.M) ~ log(ALT.B) + ndose, data=liver, method="MM"))
#' ast <- resid(rlm(log(AST.M) ~ log(AST.B) + ndose, data=liver, method="MM"))
#' alp <- resid(rlm(log(ALP.M) ~ log(ALP.B) + ndose, data=liver, method="MM"))
#' tbl <- resid(rlm(log(TBL.M) ~ log(TBL.B) + ndose, data=liver, method="MM"))
#'
#' r <- data.frame(alt=alt, ast=ast, alp=alp, tbl=tbl)
#'
#' Amex <- mex(r[liver$dose == "A",], mqu=.7)
#' Bmex <- mex(r[liver$dose == "B",], mqu=.7)
#' Cmex <- mex(r[liver$dose == "C",], mqu=.7)
#' Dmex <- mex(r[liver$dose == "D",], mqu=.7)
#'
#' par(mfcol=c(3,3))
#' plot(Amex)
#'
#' plot(Dmex, col="blue")
#'
#' ## Take a closer look at the marginal behaviour of ALT
#' \donttest{
#' r$ndose <- liver$ndose
#'
#' altmod1 <- evm(alt, qu=.7, phi = ~ ndose, xi = ~ ndose, data=r)
#' altmod2 <- evm(alt, qu=.7, phi = ~ ndose, data=r)
#' altmod3 <- evm(alt, qu=.7, xi = ~ ndose, data=r)
#' altmod4 <- evm(alt, qu=.7, data=r)
#'
#' # Prefer model 3, with term for xi on basis of AIC
#'
#' balt3 <- evm(alt, qu=.7, xi = ~ ndose, data=r, method="simulate")
#' par(mfrow=c(3,3))
#' plot(balt3)
#'
#' # use longer burn-in and also thin the output
#'
#' balt3 <- thinAndBurn(balt3,burn=1000,thin=5)
#' plot(balt3)
#'
#' # Get some simulated values for dose D
#'
#' DParam <- predict(balt3,type="lp",newdata=data.frame(ndose=4),all=TRUE)$obj$link[[1]]
#'
#' simD <- rgpd(nrow(DParam), sigma=exp(DParam[,"phi"]), xi=DParam[,"xi"], u=quantile(alt, .7))
#'
#' # These are simulated residuals. Get some baselines and transform all
#' # to raw scale
#'
#' b <- sample(log(liver$ALT.M), size=nrow(balt3$param), replace=TRUE)
#' res <- exp(b + simD)
#'
#' # estimate quantiles on raw scale
#' quantile(res, prob=c(.5, .75, .9, .95, .99))
#'
#' # estimate proportion exceeding 3*upper limit of normal mean(res >
#' # 36 * 3) # 36 is the upper limit of normal for ALT
#' }
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom grDevices dev.interactive nclass.Sturges terrain.colors
#'     axisTicks
#' @importFrom graphics abline axTicks axis box filled.contour hist
#'     lines matplot mtext pairs panel.smooth par plot points polygon
#'     rug segments title
#' @importFrom stats AIC D acf coef cov2cor density formula lag lowess
#'     mahalanobis median model.frame model.matrix model.response
#'     na.omit optim pexp ppoints predict qexp qnorm qlogis quantile rcauchy
#'     rlogis resid rexp runif sd simulate spline uniroot update var bw.nrd
#'     rnorm
#' @importFrom utils tail
#' @importFrom ggplot2 ggplot scale_x_continuous scale_x_discrete
#'     scale_y_continuous ggtitle aes geom_area geom_segment
#'     stat_density geom_point position_jitter stat_smooth theme
#'     geom_line geom_polygon element_text geom_rug geom_abline
#'     geom_rect geom_histogram geom_text geom_hline geom_vline
#'     geom_density labs coord_cartesian geom_smooth
#' @importFrom Rcpp evalCpp
#' @useDynLib texmex, .registration=TRUE
NULL



