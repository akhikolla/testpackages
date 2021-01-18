#' Elastic Functional Data Analysis
#'
#' A library for functional data analysis using the square root
#' velocity framework which performs pair-wise and group-wise
#' alignment as well as modeling using functional component
#' analysis
#'
#' @name fdasrvf
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @references J. D. Tucker, W. Wu, and A. Srivastava, ``Phase-Amplitude Separation of Proteomics Data Using Extended Fisher-Rao Metric," Electronic Journal of Statistics, Vol 8, no. 2. pp 1724-1733, 2014.
#' @references J. D. Tucker, W. Wu, and A. Srivastava, ``Analysis of signals under compositional noise With applications to SONAR data," IEEE Journal of Oceanic Engineering, Vol 29, no. 2. pp 318-330, Apr 2014.
#' @references Tucker, J. D. 2014, Functional Component Analysis and Regression using Elastic Methods. Ph.D. Thesis, Florida State University.
#' @references Robinson, D. T. 2012, Function Data Analysis and Partial Shape Matching in the Square Root Velocity Framework. Ph.D. Thesis, Florida State University.
#' @references Huang, W. 2014, Optimization Algorithms on Riemannian Manifolds with Applications. Ph.D. Thesis, Florida State University.
#' @references Cheng, W., Dryden, I. L., and Huang, X. (2016). Bayesian registration of functions and curves. Bayesian Analysis, 11(2), 447-475.
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @references Cheng, W., Dryden, I. L., and Huang, X. (2016). Bayesian registration of functions and curves. Bayesian Analysis, 11(2), 447-475.
#' @references W. Xie, S. Kurtek, K. Bharath, and Y. Sun, A geometric approach to visualization of variability in functional data, Journal of American Statistical Association 112 (2017), pp. 979-993.
#' @references Lu, Y., R. Herbei, and S. Kurtek, 2017: Bayesian registration of functions with a Gaussian process prior. Journal of Computational and Graphical Statistics, 26, no. 4, 894–904.
#' @references Lee, S. and S. Jung, 2017: Combined analysis of amplitude and phase variations in functional data. arXiv:1603.01775 [stat.ME], 1–21.
#' @references J. D. Tucker, J. R. Lewis, and A. Srivastava, “Elastic Functional Principal Component Regression,” Statistical Analysis and Data Mining, vol. 12, no. 2, pp. 101-115, 2019.
#' @references J. D. Tucker, J. R. Lewis, C. King, and S. Kurtek, “A Geometric Approach for Computing Tolerance Bounds for Elastic Functional Data,” Journal of Applied Statistics, 10.1080/02664763.2019.1645818, 2019.
#' @references T. Harris, J. D. Tucker, B. Li, and L. Shand, "Elastic depths for detecting shape anomalies in functional data," Technometrics, 10.1080/00401706.2020.1811156, 2020.
#' @docType package
#' @useDynLib fdasrvf, .registration=TRUE
#' @import foreach mvtnorm matrixcalc splines parallel doParallel Rcpp fields
#' @importFrom graphics layout legend matplot plot title lines image plot.new
#' @importFrom grDevices rainbow
#' @importFrom coda traceplot mcmc
#' @importFrom viridisLite viridis
#' @importFrom stats approx cov optim predict quantile rnorm runif sd smooth.spline var spline median rgamma optimize
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom tolerance mvtol.region
#' @importFrom testthat test_check
#' @aliases fdasrvf fdasrvf-package
NULL
#' Simulated two Gaussian Dataset
#'
#' A functional dataset where the individual functions are given by:
#' \eqn{y_i(t) = z_{i,1} e^{-(t-1.5)^2/2} + z_{i,2}e^{-(t+1.5)^2/2}}, \eqn{t \in [-3, 3], ~i=1,2,\dots, 21},
#' where \eqn{z_{i,1}} and \eqn{z_{i,2}} are \emph{i.i.d.} normal with mean one and standard deviation
#' 0.25. Each of these functions is then warped according to: \eqn{\gamma_i(t) = 6({e^{a_i(t+3)/6} -1 \over e^{a_i} - 1}) - 3}
#' if  \eqn{a_i \neq 0}, otherwise \eqn{\gamma_i = \gamma_{id}} (\eqn{gamma_{id}(t) = t})
#' is the identity warping). The variables are as follows: f containing the
#' 21 functions of 101 samples and time which describes the sampling
#'
#'
#' @docType data
#' @keywords datasets
#' @name simu_data
#' @usage data("simu_data")
#' @format A list which contains f and time
NULL
#' Distributed Gaussian Peak Dataset
#'
#' A functional dataset where the individual functions are given by a Gaussian
#' peak with locations along the \eqn{x}-axis. The variables are as follows:
#' f containing the 29 functions of 101 samples and time which describes the
#' sampling
#'
#'
#' @docType data
#' @keywords datasets
#' @name toy_data
#' @usage data("toy_data")
#' @format A list which contains f and time
NULL
#' Berkley Growth Velocity Dataset
#'
#' Combination of both boys and girls growth velocity from the Berkley Dataset
#'
#'
#' @docType data
#' @keywords datasets
#' @name growth_vel
#' @usage data("growth_vel")
#' @format A list which contains f and time
NULL
#' Aligned Distributed Gaussian Peak Dataset
#'
#' A functional dataset where the individual functions are given by a Gaussian
#' peak with locations along the \eqn{x}-axis. The variables are as follows:
#' f containing the 29 functions of 101 samples and time which describes the
#' sampling which as been aligned
#'
#'
#' @docType data
#' @keywords datasets
#' @name toy_warp
#' @usage data("toy_warp")
#' @format A list which contains the outputs of the time_warping function
NULL
#' Aligned Simulated two Gaussian Dataset
#'
#' A functional dataset where the individual functions are given by:
#' \eqn{y_i(t) = z_{i,1} e^{-(t-1.5)^2/2} + z_{i,2}e^{-(t+1.5)^2/2}}, \eqn{t \in [-3, 3], ~i=1,2,\dots, 21},
#' where \eqn{z_{i,1}} and \eqn{z_{i,2}} are \emph{i.i.d.} normal with mean one and standard deviation
#' 0.25. Each of these functions is then warped according to: \eqn{\gamma_i(t) = 6({e^{a_i(t+3)/6} -1 \over e^{a_i} - 1}) - 3}
#' if  \eqn{a_i \neq 0}, otherwise \eqn{\gamma_i = \gamma_{id}} (\eqn{gamma_{id}(t) = t})
#' is the identity warping). The variables are as follows: f containing the
#' 21 functions of 101 samples and time which describes the sampling which has been aligned
#'
#'
#' @docType data
#' @keywords datasets
#' @name simu_warp
#' @usage data("simu_warp")
#' @format A list which contains the outputs of the time_warping function
NULL
#' Aligned Simulated two Gaussian Dataset using Median
#'
#' A functional dataset where the individual functions are given by:
#' \eqn{y_i(t) = z_{i,1} e^{-(t-1.5)^2/2} + z_{i,2}e^{-(t+1.5)^2/2}}, \eqn{t \in [-3, 3], ~i=1,2,\dots, 21},
#' where \eqn{z_{i,1}} and \eqn{z_{i,2}} are \emph{i.i.d.} normal with mean one and standard deviation
#' 0.25. Each of these functions is then warped according to: \eqn{\gamma_i(t) = 6({e^{a_i(t+3)/6} -1 \over e^{a_i} - 1}) - 3}
#' if  \eqn{a_i \neq 0}, otherwise \eqn{\gamma_i = \gamma_{id}} (\eqn{gamma_{id}(t) = t})
#' is the identity warping). The variables are as follows: f containing the
#' 21 functions of 101 samples and time which describes the sampling which has been aligned
#'
#'
#' @docType data
#' @keywords datasets
#' @name simu_warp_median
#' @usage data("simu_warp_median")
#' @format A list which contains the outputs of the time_warping function finding the median
NULL
#' MPEG7 Curve Dataset
#'
#' Contains the MPEG7 curve data set which is 20 curves in 65 classes. The array
#' is structured with dimension (2,100,65,20)
#'
#'
#' @docType data
#' @keywords datasets
#' @name beta
#' @usage data("mpeg7")
#' @format an array of shape (2,100,65,20)
NULL
#' Example Image Data set
#'
#' Contains two simulated images for registration
#'
#'
#' @docType data
#' @keywords datasets
#' @name im
#' @usage data("image")
#' @format a list containing two images of dimension  (64,64)
NULL
