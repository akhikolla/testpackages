# Package description ####

#' multivariance: Measuring Multivariate Dependence Using Distance Multivariance
# Multivariance: detecting and measuring multivariate dependence
#'
#' The multivariance package provides basic functions to calculate distance multivariance and related quantities. To test independence use \code{\link{multivariance.test}}, it provides an interface (via its arguments) to all the tests based on distance (m-/total-)multivariance. The package offers also several other functions related to distance multivariance, e.g. a detection and visualization of dependence structures \code{\link{dependence.structure}}. See below for details on the full content of the package.
#'
# It also includes a function to perform a dependence analysis.
#'
#' Distance multivariance is a multivariate dependence measure, which can be used to detect dependencies between an arbitrary number of random vectors each of which can have a distinct dimension. The necessary functions are implemented in this package, and examples are given. For the theoretic background we refer to the papers [1,2,3,4,5,6]. Paper [3] includes a summary of the first two. It is the recommended starting point for users with an applied interest. Paper [4] is concerned with new (faster) p-value estimates for the independence tests, [5] introduces the copula versions of distance multivariance, [6] discusses the quantification of dependence using distance multicorrelations.
#'
#' The (current) code is speed improved in comparison to the former releases. Certainly there is still room for improvement and development. Questions, comments and remarks are welcome: \email{bjoern.boettcher@@tu-dresden.de}
#'
# Users on Windows machines might get a considerable speed up using MRO instead of the standard R release - since this is particularly faster for large matrix operations.
#'
#' For infos on the latest changes and/or updates to the package use \code{news(package="multivariance")}.
#'
#' To cite this package use the standard citation for R packages, i.e., the output of \code{citation("multivariance")}.
#'
#' @section Multivariance:
#'
#'  \code{\link{multivariance}} computes the distance multivariance
#'
#'  \code{\link{total.multivariance}} computes the total distance multivariance
#'
#'  \code{\link{m.multivariance}} computes the m-multivariance (introduced in [3])
#'
#'  It might be convenient to compute these simultaneously using \code{\link{multivariances.all}}.
#'
#'  \code{\link{copula.multivariance}} computes the copula versions of the above (introduced in [5])
#'
#'  \code{\link{multicorrelation}} computes the multicorrelations (discussed specifically in [6])
#'
#'
#' @section Functions to use and interpret multivariance:
#'
#'  \code{\link{rejection.level}} computes a (conservative) rejection level for a given significance level. This can be used for a conservative interpretation of distance multivariance. The counterpart is \code{\link{multivariance.pvalue}}, which computes a conservative p-value for a given distance multivariance. Both methods are distribution-free.
#'
#'  \code{\link{resample.rejection.level}} and \code{\link{resample.pvalue}} are the distribution dependent versions of the above. They are approximately sharp, but computational more expensive. Any resampling is done by \code{\link{resample.multivariance}}.
#'
#'  Using the methods developed in [4] approximate p-value estimates are provided by \code{\link{pearson.pvalue}}. This method is much faster than the resampling method.
#'
#'  \code{\link{multivariance.test}} provides the corresponding tests of independence. The former provides output as common for tests in R.
#'
#' \code{\link{cdm}} and \code{\link{cdms}} compute the doubly centered distance matrix and matrices, respectively. These can be used to speed up repeated computations of distance multivariance.
#'
#' In [4] various methods to estimate the moments of the test statistic under H0 were developed, these are (implicitly) implemented in this package only for the moments used in \code{\link{pearson.pvalue}}. Further and explicit functions can be added upon request. Please feel free to contact the author.
#'
#' \code{\link{emp.transf}} computes the Monte Carlo empirical transform of the data. This data yields the copula version of distance multivariance. Hereto note, that values become randomized due to the "Monte Carlo empirical transform", i.e., the copula versions yield in a finite sample setting not identical values for repeated runs.
#'
#' For planing of large projects or studies it might be convenient to estimate the computation time of multivariance via \code{\link{multivariance.timing}}.
#'
#' @section Dependence structures:
#'
#'  \code{\link{dependence.structure}} performs the dependence structure detection algorithm as described in [3].
#'
#'  \code{\link{find.cluster}} is the basic building block of \code{\link{dependence.structure}}. It is recommended to use \code{\link{dependence.structure}}.
#'
#' @section Examples:
#'
#' \code{\link{coins}} and \code{\link{tetrahedron}} generate samples of pairwise independent random variables, with dependence of higher order.
#'
#' \code{\link{dep_struct_iterated_13_100}}, \code{\link{dep_struct_ring_15_100}}, \code{\link{dep_struct_several_26_100}} and \code{\link{dep_struct_star_9_100}} are example data sets for the dependence structure detection. These might also serve as benchmark examples.
#'
#' \code{\link{anscombe.extended}} provides an extension of Anscombe's Quartett. It illustrates that a large value of Pearson's correlation can occur for very different dependencies and that this is not a small-sample problem. These dependencies are at least partly differentiated by values of distance multicorrelation.
#'
#' @references
#' [1] B. Böttcher, M. Keller-Ressel, R.L. Schilling, Detecting independence of random vectors: generalized distance covariance and Gaussian covariance. Modern Stochastics: Theory and Applications, Vol. 5, No. 3(2018) 353-383. \url{https://www.vmsta.org/journal/VMSTA/article/127/info}
#'
#' [2] B. Böttcher, M. Keller-Ressel, R.L. Schilling, Distance multivariance: New dependence measures for random vectors. The Annals of Statistics, Vol. 47, No. 5 (2019) 2757-2789. \url{https://projecteuclid.org/euclid.aos/1564797863}
#'
#' [3] B. Böttcher, Dependence and Dependence Structures: Estimation and Visualization using the Unifying Concept of Distance Multivariance. Open Statistics, Vol. 1, No. 1 (2020) 1-46. \url{https://doi.org/10.1515/stat-2020-0001}
#'
#' [4] G. Berschneider, B. Böttcher, On complex Gaussian random fields, Gaussian quadratic forms and sample distance multivariance. Preprint. \url{https://arxiv.org/abs/1808.07280}
#'
#' [5] B. Böttcher, Copula versions of distance multivariance and dHSIC via the distributional transform -- a general approach to construct invariant dependence measures. Statistics, (2020) 1-18. \url{https://doi.org/10.1080/02331888.2020.1748029}
#'
#' [6] B. Böttcher, Notes on the interpretation of dependence measures -- Pearson's correlation, distance correlation, distance multicorrelations and their copula versions. Preprint. \url{https://arxiv.org/abs/2004.07649}
#'
#' @docType package
#' @name multivariance-package
NULL
#' @useDynLib multivariance
NULL
#' @importFrom Rcpp sourceCpp
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("Welcome to 'multivariance' version ",utils::packageVersion("multivariance"),".

For usage hints and details on the theoretic backgound see
help(\"multivariance-package\")
and the references given therein, starting with:
https://doi.org/10.1515/stat-2020-0001

			To suppress this message use:
			suppressPackageStartupMessages(library(multivariance))"))
}


################# Multivariance ###########


#' rejection level for the test statistic
#'
#' Under independence the probability for the normalized and Nscaled (squared) multivariance to be above this level is less than \code{alpha}. The same holds for the normalized, Nscaled and Escaled (squared) total multivariance and m-multivariance.
#'
#' @param alpha level of significance
#' @details
#' This is based on a distribution-free approach. The value might be very conservative. This is the counterpart to \code{\link{multivariance.pvalue}}. For a less conservative approach see \code{\link{resample.rejection.level}}.
#'
#' The estimate is only valid for \code{alpha} smaller than 0.215.
#'
#' @examples
#' rejection.level(0.05) #the rejection level, for comparison with the following values
#' total.multivariance(matrix(rnorm(100*3),ncol = 3)) #independent sample
#' total.multivariance(coins(100)) #dependent sample which is 2-independent
#'
#' # and the p values are (to compare with alpha)
#' multivariance.pvalue(total.multivariance(matrix(rnorm(100*3),ncol = 3))) #independent sample
#' multivariance.pvalue(total.multivariance(coins(100))) #dependent sample which is 2-independent
#'
#' \dontrun{
#' # visualization of the rejection level
#' curve(rejection.level(x),xlim = c(0.001,0.215),xlab = "alpha")
#' }
#'
#' @export
rejection.level = function(alpha) {
  if (any(alpha > 0.215)) warning("alpha too large. Only valid for alpha smaller than 0.215! \n")
  return((stats::qnorm(1-alpha/2)^2))
  # identical with qchisq(1-alpha,1)
}

#' transform multivariance to p-value
#'
#' Computes a conservative p-value for the hypothesis of independence for a given multivariance / m-multivariance / total multivariance.
#'
#' @param x value of a normalized \code{\link{multivariance}} scaled by the sample size (i.e., computed with \code{Nscale = TRUE})
#'
#' @details
#' This is based on a distribution-free approach. The p-value is conservative, i.e. it might be much smaller. This is the counterpart to \code{\link{rejection.level}}. For a less conservative approach see \code{\link{resample.pvalue}} or \code{\link{pearson.pvalue}}.
#'
#' p-values larger than 0.215 might be incorrect, since the distribution-free estimate on which the computation is based only holds up to 0.215.
#'
#' @references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
multivariance.pvalue = function(x) {
  if (any(x < 0,na.rm = TRUE)) print(paste("Negative multivariance = ",x[which(x<0)]))
  2-2*stats::pnorm(sqrt(x))
}


#' computes a doubly centered distance matrix
#'
#' computes the doubly centered distance matrix
#'
#' @param x matrix, each row of the matrix is treated as one sample
#' @param normalize logical, indicates if the matrix should be normalized
#' @param psi if it is \code{NULL}, the euclidean distance will be used. In the case of \code{isotropic = TRUE}: a real valued negative definite function of one variable (accepting vectors as arguments; returning a vector of the same length). In the case of \code{isotropic = FALSE}: a real valued function of two variables (or vectors) to compute the distance of two samples based on a continuous negative definite function.
#' @param isotropic logical, indicates if psi of the Euclidean distance matrix should be computed, i.e., if an isotropic distance should be used.
#' @param p numeric, if it is a value between 1 and 2 then the Minkowski distance with parameter p is used.
#' @param external.dm.fun here one can supply an external function, which computes the distance matrix given \code{x}.
#'
#' @details
#' The doubly centered distance matrices are required for the computation of (total / m-) multivariance.
#'
#'If \code{normalize = TRUE} then the value of multivariance is comparable and meaningful. It can be compared to the \code{\link{rejection.level}} or its p-value \code{\link{multivariance.pvalue}} can be computed.
#'
#' More details: If \code{normalize = TRUE} the matrix is scaled such that the multivariance based on it, times the sample size, has in the limit - in the case of independence - the distribution of an L^2 norm of a Gaussian process with known expectation.
#'
#' As default the Euclidean distance is used. The parameters \code{psi}, \code{p}, \code{isotropic} and \code{external.dm.fun} can be used to select a different distance. In particular, \code{external.dm.fun} can be used to provide any function which calculates a distance matrix for the rows of a given matrix.
#'
#' @references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' x = coins(100)
#' cdm(x) # fast euclidean distances
#' cdm(x,psi = function(x,y) sqrt(sum((x-y)^2))) # this is identical to the previous (but slower)
#'
#' # the function cdm does the following three lines in a faster way
#' N = nrow(x)
#' C = diag(N) - matrix(1/N,nrow = N,ncol = N)
#' A = - C %*% as.matrix(stats::dist(x,method="euclidean")) %*% C #'
#' all(abs(A- cdm(x,normalize = FALSE)) < 10^(-12))
#'
#' @export
cdm = function(x, normalize = TRUE, psi = NULL, p = NULL, isotropic = FALSE, external.dm.fun = NULL) {
  if (!is.matrix(x)) x = as.matrix(x)
  if (is.null(psi) & is.null(p) & is.null(external.dm.fun)) {
    #dm = dist.to.matrix(stats::dist(x,method="euclidean"))
    #DEVELOPING NOTE: here as.matrix was slow, dist.to.matrix is faster. Instead one could just use the vector....
    # even faster for the euclidean case is fastdist defined via Rcpp
    #dm = fastdist(as.matrix(x))

    return(fastEuclideanCdm(x,normalize))

  } else {
    if (!is.null(p)) {
      if ((p<1) || (p>2)) warning("p is not in [1,2]. \n")
      dm = dist.to.matrix(stats::dist(x,method="minkowski", p = p))
    } else { # case: psi is given
      if (!is.null(external.dm.fun)) {
        dm = external.dm.fun(as.matrix(x))
      } else {
        if (isotropic) {
          #dm = psi(dist.to.matrix(stats::dist(x,method="euclidean")))
          dm = psi(fastdist(as.matrix(x)))
        } else {
          x = as.matrix(x)
          n = nrow(x)
          d = ncol(x)
          dm = matrix(apply(cbind(x[rep(1:n,n),],x[rep(1:n,each = n),]), #create all combinations
                            1, # apply to each row
                            function(y) psi(y[1:d], y[(d+1):(2*d)])),nrow = n)
          # DEVELOPING NOTE: could double the speed if only the upper triangular matrix is computed, using the idea of dist.to.matrix
        }
      }
    }
  }

  return(doubleCenterSymMat(dm,normalize))

  #alternative (even slower) implementations:
  # colm = colMeans(dm)
  # m = mean(colm)  # equals mean(dm)
  #
  # if (m == 0) warning("It seems that one variable is constant. Constants are always independent. \n")
  # if (normalize && (m != 0)) {
  #   return((-dm + outer(colm,colm, FUN ="+") - m)/ m)
  # } else {
  #   return(-dm + outer(colm,colm, FUN ="+") - m)
  # }

  #alternative (even slower) implementations:
  #cdm1 = sweep(dm,1,colm)
  #cdm2 = -sweep(cdm1,2,rowMeans(dm)) - m
  #cdm2 = -(x - rep(colm, ncol(dm)) - rep(rowMeans(dm),each = ncol(dm))) - m # for quadratic matrix

}

#' computes the doubly centered distance matrices
#' @param x matrix, each row is a sample
#' @param vec vector which indicates which columns are treated as one sample
#' @param membership depreciated. Now use \code{vec}.
#' @param ... these are passed to \code{\link{cdm}}
#'
#' @return It returns a list of distance matrices.
#'
#' @export
cdms = function(x,vec = 1:ncol(x),membership = NULL,...) {
  if (!is.null(membership)) {
    vec = membership
    warning("Use 'vec' instead of 'membership' as argument to 'cdms'. 'membership' is depreciated. \n")
  }
  if (anyNA(vec)) vec = 1:ncol(x)
  n = max(vec)
  # N = nrow(x)
  # array.cdm = array(,dim = c(N,N,n))
  # for (i in 1:n) array.cdm[,,i] = cdm(x[,(vec == i)],...)
  # return(array.cdm)

  return(lapply(1:n, function(i) cdm(x[,(vec == i),drop = FALSE],...)))
}

#' double centering of a matrix
#'
#' # changed default after 2.0.0
#' @keywords internal
double.center = function(dm,normalize = TRUE) {
  if (is.list(dm)) { # double center a list of matrices
    return(lapply(dm, function(l) double.center(l,normalize) ))
    # n = lenght(tm)
    # array.cdm = array(,dim=dim(dm))
    # for (i in 1:n) array.cdm[,,i] = double.center(dm[,,i],normalize)
    # return(array.cdm)
  } else { # double center a matrix
    colm = colMeans(dm)
    m = mean(colm)  # equals mean(dm)

    if (m == 0) warning("It seems that one variable is constant. Constants are always independent. \n")
    if (normalize && (m != 0)) {
      return((-dm + outer(colm,colm, FUN ="+") - m)/ m)
    } else {
      return(-dm + outer(colm,colm, FUN ="+") - m)
    }
  }
}

#' distance multivariance
#'
#' Computes the distance multivariance, either for given data or a given list of doubly centered distance matrices.
#'
#' @param x either a data matrix or a list of doubly centered distance matrices
#' @param vec if x is a matrix, then this indicates which columns are treated together as one sample; if x is a list, these are the indexes for which the multivariance is calculated. The default is all columns and all indexes, respectively.
#' @param Nscale if \code{TRUE} the multivariance is scaled up by the sample size (and thus it is exactly as required for the test of independence)
#' @param squared if \code{FALSE} it returns the actual multivariance, otherwise the squared multivariance (less computation)
#' @param ... these are passed to \code{\link{cdms}} (which is only invoked if \code{x} is a matrix)
#' @param correlation depreciated, please use the function \code{\link{multicorrelation}} instead.
#'
#' @details
#'
#' If \code{x} is a matrix and \code{vec} is not given, then each column is treated as a separate sample. Otherwise \code{vec} has to have as many elements as \code{x} has columns and values starting from 1 up to the number of 'variables', e.g. if \code{x} is an \code{N} by 5 matrix and \code{vec = c(1,2,1,3,1)} then the multivariance of the 1-dimensional variables represented by column 2 and 4 and the 3-dimensional variable represented by the columns 1,3,5 is computed.
#'
#' As default it computes the normalized Nscaled squared multivariance, for a multivariance without normalization the argument \code{normalize = FALSE} has to be passed to \code{cdms}.
#'
#'
#' \code{correlation = TRUE} yields values between 0 and 1. These can be interpreted similarly to classical correlations, see also \code{\link{multicorrelation}}.
#'
#' As a rough guide to interpret the value of distance multivariance note:
#' \itemize{
#' \item If the random variables are not (n-1)-independent, large values indicate dependence, but small values are meaningless. Thus in this case use \code{\link{total.multivariance}}.
#' \item If the random variables are (n-1)-independent and \code{Nscale = TRUE}, values close to 1 and smaller indicate independence, larger values indicate dependence. In fact, in the case of independence the test statistic is a Gaussian quadratic form with expectation 1 and samples of it can be generated by \code{\link{resample.multivariance}}.
#' \item If the random variables are (n-1)-independent and \code{Nscale = FALSE}, small values (close to 0) indicate independence, larger values indicate dependence.
#' }
#'
#' Finally note, that due to numerical (in)precision the value of multivariance might become negative. In these cases it is set to 0. A warning is issued, if the value is negative and further than the usual (used by \code{\link[base]{all.equal}}) tolerance away from 0.
#'
#' @references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' multivariance(matrix(rnorm(100*3),ncol = 3)) #independent sample
#' multivariance(coins(100)) #dependent sample which is 2-independent
#'
#' x = matrix(rnorm(100*2),ncol = 2)
#' x = cbind(x,x[,2])
#' multivariance(x) #dependent sample which is not 2-independent (thus small values are meaningless!)
#' multivariance(x[,1:2]) #these are independent
#' multivariance(x[,2:3]) #these are dependent
#'
#' multivariance(x[,2:3],correlation = TRUE)
#'
#' @export
multivariance = function(x,vec = NA,Nscale = TRUE,correlation = FALSE, squared = TRUE, ...) {

  if (correlation) warning("The option 'correlation' is depreciated. Please use the function 'multicorrelation' instead. \n")
  if (is.list(x)) {
    if (is.data.frame(x)) stop("Input is a data.frame. Please provide a matrix or a list of doubly centered distance matrices.")
  } else { # if the input is a matrix, the distance matrices are computed
    if (is.array(x) & (length(dim(x))>2)) stop("Please provide a list instead of an array. Changed since version 2.0.0.\n")

    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }

  if (anyNA(vec)) vec = 1:length(x)

  if (anyNA(x)) stop("Provided x contains NA. \n")

  #if (length(vec) > dim(x)[2]) warning("More data columns than rows.")
  #Aprod = x[,,vec[1]]
  #for (i in 2:length(vec)) Aprod = Aprod * x[,,vec[i]]
  #result = mean(Aprod)

  result = mean(Reduce("*", x[vec]))

  if (Nscale && !correlation) result = result *nrow(x[[1]])

  if (correlation) {
    n = length(vec)
    norm.n = function(x) mean(abs(x)^n)^(1/n)
    #Anorm = mean(abs(x[,,vec[1]]^n))^(1/n)
    #for (i in 2:length(vec)) Anorm = Anorm * mean(abs(x[,,vec[i]]^n))^(1/n)
    #result = result / Anorm

    result = result / Reduce("*",lapply(x[vec], norm.n))

  }
  # DEVELOPING NOTE: The following is much slower .... It might be faster, if we store the matrices only as vectors with the upper triangular as elements.
  #ut = upper.tri(Aprod)
  #diat = diag(x[vec[1],,])
  #test = x[vec[1],,][ut]
  #for (i in 2:length(vec)) {
  #  diat = diat * diag(x[vec[i],,])
  #  test = test * x[vec[i],,][ut]
  #}
  #erg = sum(diat,2*test)/ncol(x)^2

  if (result < 0) {
    if (!isTRUE(all.equal(result,0))) warning(paste("Value of multivariance was negative (",result,"). This is usually due to numerical (in)precision. It was set to 0. \n"))
    result = 0
  }

  if (squared | (Nscale && !correlation)) { return(result)
  } else { return(sqrt(result))}


  # return(mean(apply(x[vec,,],c(2,3),prod))) #this vector version is also much much slower!!!
  # DEVELOPING NOTE: mean is slow due to error correction. sum()/length() is faster, but not as precise.
}

#' total distance multivariance
#'
#' computes the total distance multivariance
#'
#' @inheritParams multivariance
#' @param lambda a scaling parameter >0. Each k-tuple multivariance gets weight \code{lambda^(n-k)}.
#' @param Escale if \code{TRUE} then it is scaled by the number of multivariances which are theoretically summed up (in the case of independence this yields for normalized distance matrices an estimator with expectation 1)
#'
#' @details
#' Total distance multivariance is per definition the scaled sum of certain distance multivariances, and it characterize dependence.
#'
#'  As a rough guide to interpret the value of total distance multivariance note:
#' \itemize{
#' \item Large values indicate dependence.
#' \item For \code{Nscale = TRUE} values close to 1 and smaller indicate independence, larger values indicate dependence. In fact, in the case of independence the test statistic is a Gaussian quadratic form with expectation 1 and samples of it can be generated by \code{\link{resample.multivariance}}.
#' \item For \code{Nscale = FALSE} small values (close to 0) indicate independence, larger values indicate dependence.
#' }

#'
#' Finally note, that due to numerical (in)precision the value of total multivariance might become negative. In these cases it is set to 0. A warning is issued, if the value is negative and further than the usual (used by \code{\link[base]{all.equal}}) tolerance away from 0.
#'
#'@references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' x = matrix(rnorm(100*3),ncol = 3)
#' total.multivariance(x) #for an independent sample
#' # the value coincides with
#' (multivariance(x[,c(1,2)],Nscale = TRUE) + multivariance(x[,c(1,3)],Nscale = TRUE)+
#'  multivariance(x[,c(2,3)],Nscale = TRUE) + multivariance(x,Nscale = TRUE))/4
#'
#' total.multivariance(coins(100)) #value for a dependent sample which is 2-independent
#'
#' @export
total.multivariance = function(x,vec = NA,lambda = 1, Nscale = TRUE,Escale = TRUE,squared = TRUE,...) {
  if (is.list(x)) {
    if (is.data.frame(x)) stop("Input is a data.frame. Please provide a matrix or a list of doubly centered distance matrices.")
  } else { # if the input is a matrix, the distance matrices are computed
    if (is.array(x) & (length(dim(x))>2)) stop("Please provide a list instead of an array. Changed since version 2.0.0.")

    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) vec = 1:length(x)

  if (anyNA(x)) stop("provided x contains NA")

  n = length(vec)

  result = mean(Reduce("*", lapply(x[vec],function(y) lambda + y))) - lambda^n

  #Aprod = lambda + x[,,vec[1]]
  #for (i in 2:length(vec)) Aprod = Aprod * (lambda + x[,,vec[i]])
  #result = mean(Aprod)-lambda^(length(vec))

  if (result < 0) {
    if (!isTRUE(all.equal(result,0))) warning(paste("Value of total multivariance was negative (",result,"). This is usually due to numerical (in)precision. It was set to 0. \n"))
    result = 0
  }

  if (Nscale) result = result * nrow(x[[1]])
  if (Escale) result = result/((1+lambda)^n - n*lambda^(n-1) - lambda^n)

  if (squared | Nscale) { return(result)
  } else { return(sqrt(result))}
}

#' m distance multivariance
#'
#' Computes m distance multivariance.
#'
#' @details
#'
#' m-distance multivariance is per definition the scaled sum of certain distance multivariances, and it characterize m-dependence.
#'
#'  As a rough guide to interpret the value of total distance multivariance note:
#' \itemize{
#' \item Large values indicate dependence.
#' \item If the random variables are (m-1)-independent and \code{Nscale = TRUE}, values close to 1 and smaller indicate m-independence, larger values indicate dependence. In fact, in the case of independence the test statistic is a Gaussian quadratic form with expectation 1 and samples of it can be generated by \code{\link{resample.multivariance}}.
#' \item If the random variables are (m-1)-independent and \code{Nscale = FALSE}, small values (close to 0) indicate m-independence, larger values indicate dependence.
#' }
#'
#' Since random variables are always 1-independent, the case \code{m=2} characterizes pairwise independence.
#'
#' Finally note, that due to numerical (in)precision the value of m-multivariance might become negative. In these cases it is set to 0. A warning is issued, if the value is negative and further than the usual (used by \code{\link[base]{all.equal}}) tolerance away from 0.
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @inheritParams multivariance
#' @param m \code{=2} or \code{3} the m-multivariance will be computed.
#' @param Escale if \code{TRUE} then it is scaled by the number of multivariances which are theoretically summed up (in the case of independence this yields for normalized distance matrices an estimator with expectation 1)
#'
#'
#'
#' @examples
#' x = matrix(rnorm(3*30),ncol = 3)
#'
#' # the following values are identical
#' m.multivariance(x,m =2)
#' 1/choose(3,2)*(multivariance(x[,c(1,2)]) +
#'                multivariance(x[,c(1,3)]) +
#'                multivariance(x[,c(2,3)]))
#'
#' # the following values are identical
#' m.multivariance(x,m=3)
#' multivariance(x)
#'
#' # the following values are identical
#' 1/4*(3*(m.multivariance(x,m=2)) + m.multivariance(x,m=3))
#' total.multivariance(x, Nscale = TRUE)
#' 1/4*(multivariance(x[,c(1,2)], Nscale = TRUE) +
#'      multivariance(x[,c(1,3)], Nscale = TRUE) +
#'      multivariance(x[,c(2,3)], Nscale = TRUE) + multivariance(x, Nscale = TRUE))
#'
#' @export
m.multivariance = function(x, vec= NA, m = 2, Nscale = TRUE, Escale = TRUE, squared = TRUE,...) {
  if (is.list(x)) {
    if (is.data.frame(x)) stop("Input is a data.frame. Please provide a matrix or a list of doubly centered distance matrices.")
  } else { # if the input is a matrix, the distance matrices are computed
    if (is.array(x) & (length(dim(x))>2)) stop("Please provide a list instead of an array. Changed since version 2.0.0.")

    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) vec = 1:length(x)

  if (anyNA(x)) stop("provided x contains NA")

  n = length(vec)

  xvec = x[vec]
  vec = 1:n

  #k = 2
  if (m == 2) {
    # Asum = x[,,vec[1]]
    # A2sum = Asum^2 #x[vec[1],,]^2
    # for (i in 2:length(vec)) {
    #   tempFactor = x[,,vec[i]]
    #   Asum = Asum + tempFactor
    #   A2sum = A2sum + tempFactor^2
    # }

    Asum = Reduce("+", xvec)
    A2sum = Reduce("+", lapply(xvec,function(y) y^2))

    result = mean(Asum^2 - A2sum)/2
  }

  #k more general
  #k = 3
  if (m == 3) {
    if (n < 3) {
      result = NA
    } else {
      # Asum = x[,,vec[1]]
      # A2sum = Asum^2 #x[vec[1],,]^2
      # A3sum = A2sum * Asum #x[vec[1],,]^3
      # for (i in 2:length(vec)) {
      #   tempFactor = x[,,vec[i]]
      #   Asum = Asum + tempFactor
      #   summand = tempFactor^2
      #   A2sum = A2sum + summand #x[vec[i],,]^2
      #   A3sum = A3sum + summand *tempFactor #x[vec[i],,]^3
      # }

      Asum = Reduce("+", xvec)
      A2sum = Reduce("+", lapply(xvec,function(y) y^2))
      A3sum = Reduce("+", lapply(xvec,function(y) y^2*y))

      result = mean(Asum^2*Asum- 3* Asum *A2sum + 2 * A3sum)/ 6
      #    result = mean(Asum^3 - choose(3,2)* Asum *A2sum + 2 * A3sum)/ factorial(3)
    }
  }

  if (m > 3) {
    warning("m > 3, not implemented. \n")
    result = NA
  }

  if (is.na(result)) {
    return(result)
  } else {
    if (result < 0) {
      if (!isTRUE(all.equal(result,0))) warning(paste("Value of m-multivariance was negative (",result,"). This is usually due to numerical (in)precision. It was set to 0. \n"))
      result = 0
    }

    if (Nscale) result = result *nrow(Asum)
    if (Escale) result = result/(choose(n,m))

    if (squared | Nscale) { return(result)
    } else { return(sqrt(result))}
  }
}

#' simultaneous computation of multivariance and total/ 2-/ 3-multivariance
#'
#' Computes simultaneously multivariance, total multivariance, 2-multivariance and 3-multivariance.
#'
#' @inheritParams multivariance
#'
#' @seealso \code{\link{multivariance}}, \code{\link{total.multivariance}}, \code{\link{m.multivariance}}
#'
#' @details
#' The computation is faster than the separate computations.
#'
#' @return Returns a vector with multivariance, total.multivariance, 2-multivariance and 3-multivariance
#'
#' @examples
#' x = coins(100,k = 3)
#' multivariances.all(x)
#' # yields the same as:
#' multivariance(x)
#' total.multivariance(x)
#' m.multivariance(x,m=2)
#' m.multivariance(x,m=3)
#'
#'
#' @export
#'
multivariances.all = function(x, vec= NA, Nscale = TRUE, squared = TRUE,...) {
  if (is.list(x)) {
    if (is.data.frame(x)) stop("Input is a data.frame. Please provide a matrix or a list of doubly centered distance matrices.")
  } else { # if the input is a matrix, the distance matrices are computed
    if (is.array(x) & (length(dim(x))>2)) stop("Please provide a list instead of an array. Changed since version 2.0.0.")

    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) vec = 1:length(x)

  if (anyNA(x)) stop("provided x contains NA")

  n = length(vec)

  xvec = x[vec]
  vec = 1:n

  # Asum = x[,,vec[1]]
  # A2sum = Asum^2 # = x[vec[1],,]^2
  # A3sum = A2sum * Asum # = x[vec[1],,]^3
  # Aprod = Asum
  # Aplusprod = 1 + Asum # = x[vec[1],,]
  # for (i in 2:n) {
  #   tempFactor = x[,,vec[i]]
  #   Asum = Asum + tempFactor
  #   Aprod = Aprod * tempFactor
  #   Aplusprod = Aplusprod * (1 + tempFactor)
  #   summand = tempFactor^2
  #   A2sum = A2sum + summand #x[vec[i],,]^2
  #   A3sum = A3sum + summand *tempFactor #x[vec[i],,]^3
  # }
  #
  Aprod = Reduce("*", xvec)
  Aplusprod = (Reduce("*", lapply(xvec,function(y) 1 + y)))
  Asum = Reduce("+", xvec)
  A2sum = Reduce("+", lapply(xvec,function(y) y^2))

  m = mean(Aprod)
  if (n==2) {
    # this and the following formular are different
    # within numerical tolerance, to avoid confusing
    # differences we use the same formular in this case
    mt = m
    m2 = m
  } else {
    mt = (mean(Aplusprod)-1)/(2^n - n - 1)
    m2 = mean(Asum^2 - A2sum)/(n *(n-1))
  }
  if (n > 2) {
    if (n == 3) {
      # this and the following formular are different
      # within numerical tolerance, to avoid confusing
      # differences we use the same formular in this case
      m3 = m
    } else {
      A3sum = Reduce("+", lapply(xvec,function(y) y^2*y))
      m3 = mean(Asum^2*Asum - 3* Asum *A2sum + 2 * A3sum)/ (n *(n-1)*(n-2))
    }
  } else {
    m3 = NA
  }
  result = c(multi = m,total = mt, m.multi.2 = m2,m.multi.3 = m3)

  neg.res = result<0
  if (any(neg.res,na.rm = TRUE)) {
    if (!isTRUE(all.equal(result[neg.res],rep(0,sum(neg.res)),check.names = FALSE))) {
      warning(paste0("Value of ",c("","total-","2-","3-")[neg.res],"multivariance was negative (",result[neg.res],"). This is usually due to numerical (in)precision. It was set to 0. \n"))
    }
    result[neg.res] = 0
  }

  if (Nscale) result = result *nrow(Asum)

  if (squared | Nscale) { return(result)
  } else { return(sqrt(result))}
}


#' distance multicorrelation
#'
#' Computes various types of sample distance multicorrelation as defined and discussed in [3,4,6].
#'
#' @inheritParams multivariance
#' @param type default: "total.lower.upper", for details and other options see below
#' @param multicorrelation.type one of \code{"normalized","unnormalized"}
#' @param estimator.type one of \code{"biased","bias.corrected"}
#'
#' @details
#'
#' There exist many variants of distance multicorrelation as discussed in [6] -- and only in specific cases a direct comparison of the values is meaningful.
#'
#' The implemented options are:
#' \itemize{
#' \item \code{total.upper.lower normalized bias.corrected}: default; bounded by 1; fast; population limit characterizes independence by 0
#' \item \code{pairwise normalized bias.corrected}: bounded by 1; fast; population limit characterizes pairwise independence by 0
#' \item \code{total.upper normalized biased}: biased versions of the above
#' \item \code{total.lower normalized biased}
#' \item \code{pairwise normalized biased}
#' \item \code{multi normalized biased}: population limit characterizes only in case of lower independence the independence of all variables by 0
#' \item \code{m.multi.3 normalized biased}: population limit characterizes only in case of pairwise independence the 3-independence of all variables by 0
#' \item \code{pairwise unnormalized biased} population limit characterizes pairwise independence by 0 and relation by similarity transforms by 1
#' \item \code{multi unnormalized biased}: population limit characterizes only in case of lower independence the independence of all variables by 0 and relation by similarity transforms by 1
#' \item \code{m.multi.3 unnormalized biased}: population limit characterizes only in case of pairwise independence the 3-independence of all variables by 0 and relation by similarity transforms by 1
#' }
#'
#'
#' Further details:
#'
#' The \code{"bias.corrected"} versions require a data matrix, since they compute bias corrected centered distance matricies.
#'
#' For \code{"multi"} the unnormalized and normalized version coincide if an even number of variables is considered. They usually differ if an odd number of variables is considered. If all variables are related by similarity transforms the unnormalized \code{"unnormalized"} multicorrelations are 1.
#'
#' For \code{"pairwise"} an alias is \code{"m.multi.2"}.
#'
#' For total multicorrelation there is currently only a feasible empirical estimator for a lower or upper bound. These are upper and lower bounds for in the population setting. When using bias corrected estimators these are in general no proper bounds, but their range can be used as values for comparisons.
#'
#' @return
#'
#' Value of the multicorrelation(s).
#'
#' @references
#' For the theoretic background see the references [2,3,6] given on the main help page of this package: \link{multivariance-package}.

#' @examples
#' y = rnorm(100)
#' x = cbind(y,y*2,(y-2)/3,y+1,y*5) # all variables are related by similarity transforms
#'
#' # compute all types of correlations for x:
#' for (ty in c("total.lower","total.upper","pairwise","m.multi.3","multi"))
#'  for (mty in c("normalized"))
#'   print(paste(format(multicorrelation(
#'   x,type=ty,multicorrelation.type = mty,estimator.type = "biased")
#'   ,digits=3,nsmall = 3,width = 7),mty,ty,"correlation - biased estimate"))
#'
#' for (ty in c("total.upper.lower","pairwise"))
#'  for (mty in c("normalized"))
#'   print(paste(format(multicorrelation(
#'   x,type=ty,multicorrelation.type = mty,estimator.type = "bias.corrected")
#'   ,digits=3,nsmall = 3,width = 7),mty,ty,"correlation - bias corrected estimate"))
#'
#' for (ty in c("m.multi.2","m.multi.3","multi"))
#'  for (mty in c("unnormalized"))
#'   print(paste(format(multicorrelation(
#'   x,type=ty,multicorrelation.type = mty,estimator.type = "biased")
#'   ,digits=3,nsmall = 3,width = 7),mty,ty,"correlation - biased estimate"))
#'
#' @export
multicorrelation = function(x, vec = 1:ncol(x), type = "total.upper.lower", multicorrelation.type = "normalized", estimator.type = "bias.corrected",squared = TRUE, ...) {

  if (estimator.type == "bias.corrected") {
    # implemented options:
    # normalized: total.upper.lower, total.upper pairwise
    # all
    if (!(multicorrelation.type == "normalized")) {
      stop("For the bias corrected estimator currently only 'normalized' with type 'total.upper.lower' and 'pairwise' are implemented.")
    }
    if (!is.matrix(x)) stop("'x' must be a matrix (for the current options).")
    return(multicorrelation.bias.corrected(x = x,vec = vec,type = type,squared = squared,...))
  }

  if (!(estimator.type == "biased")) {
    stop(paste0("Unknown 'estimator.type': ",estimator.type))
  }

  if (is.list(x)) {
    if (is.data.frame(x)) stop("Input is a data.frame. Please provide a matrix or a list of doubly centered distance matrices.")
  } else { # if the input is a matrix, the distance matrices are computed
    if (is.array(x) & (length(dim(x))>2)) stop("Please provide a list instead of an array. Changed since version 2.0.0.")

    if (anyNA(vec)) vec = 1:ncol(x)
    x = cdms(x,vec,normalize = FALSE,...)
    vec = 1:max(vec)
  }
  if (anyNA(vec)) vec = 1:length(x)

  if (anyNA(x)) stop("provided x contains NA")

  n = max(vec)

  if (type == "total.lower") {
    type = "total"
  }
  if (type == "total.upper") {
    n = 2
    type = "total"
  }

  if (type == "pairwise") type = "m.multi.2"
  if (type == "m.multi.2") n = 2
  if (type == "m.multi.3") n = 3

  #list.cdms = cdms(x,normalize = FALSE)
  switch(multicorrelation.type,
    normalized = {
      list.norms = lapply(x, function(x) (mean(abs(x)^n))^(1/n)) # for lower bound
      list.norms[list.norms == 0] = 1 # for constant random variables no scaling is required, since the matrix is 0 anyway
    },
    unnormalized = {
      list.norms = lapply(x, function(x) (mean(x^n))^(1/n))
      if (any(list.norms == 0)) {
        check.val = lapply(x[list.norms == 0], function(x) (mean(x^2)))
        if (any(check.val != 0))  {# non constant with normalizing constant 0
          warning("Division by 0 due to a normalizing constants of 0. This can happen e.g. for an even number of samples of bernoulli distributed random variables.")
          return(NaN)
        } else {
          list.norms[list.norms == 0] = 1 # for constant random variables no scaling is required, since the matrix is 0 anyway
        }
      }
      if (type == "total") stop("The bound is only valid for 'normalized' total multivariance.")
    },
    {stop(paste("unkown multicorrelation.type:",type))}
  )
  list.cdms = lapply(1:length(x), function(i) x[[i]]/list.norms[[i]])

  switch(type,
    multi = {
      return(c( multicorrelation = multivariance(list.cdms,Nscale = FALSE,squared = squared)))
    },
    total = {
      if (n == 2) {
        return( c(total.multicorrelation.upper.bound = total.multivariance(list.cdms,Nscale = FALSE,squared = squared)))
      } else {
      return( c(total.multicorrelation.lower.bound = total.multivariance(list.cdms,Nscale = FALSE,squared = squared)))
      }
    },
    m.multi.2 = {
      return(c(multicorrelation.2 = m.multivariance(list.cdms,Nscale = FALSE,squared = squared)))
    },
    m.multi.3 = {
      return(c(multicorrelation.3 = m.multivariance(list.cdms,m = 3,Nscale = FALSE,squared = squared)))
    },
    #all = { # this does not work, since the norms are different for the cases
    #  val = c(multivariances.all(list.cdms,m = 3,Nscale = FALSE,squared = squared))
    #  names(val) = paste0(c("Mcor","tMcorlb","M2cor","M3cor"),ifelse(multicorrelation.type == "unnormalized","-unnormalized","-normalized"))
#      return(val)
#    },
    {stop(paste("unkown type:",type))}
  )
}

#' @rdname multicorrelation
#' @export
Mcor <- multicorrelation


#' test for independence
#'
#' Depreciated. Use \code{\link{multivariance.test}} instead. It provides all options and returns test result in a standard R format.
#'
#' This computes a test of independence for the columns of a sample matrix (required for the resampling test) or for given doubly centered distance matrices (only possible for the distribution-free test).
#'
#' @inheritParams multivariance
#' @param alpha significance level
#' @param type one of \code{"pearson_approx","distribution_free","resample"}
#' @param verbose logical, if TRUE meaningful text output is generated.
#'
#' @return Returns \code{TRUE} if the hypothesis of independence is NOT rejected, otherwise \code{FALSE}.
#' @details The \code{"pearson_approx"} and \code{"resample"} are approximately sharp. The latter is based on a resampling approach and thus much slower. The \code{"distribution_free"} test might be very conservative.
#' The doubly centered distance matrices can be prepared by \code{\link{cdms}}. But note that for the test based on Pearson's approximation and for the resampling test, the data matrix has to be given.
#'
#' @references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' independence.test(coins(100)) #dependent sample which is 2-independent
#' independence.test(coins(100),type = "resample") #dependent sample which is 2-independent
#'
#' independence.test(coins(100)[,2:3]) # independent sample
#' independence.test(coins(100)[,2:3],type = "resample") # independent sample
#'
#' independence.test(coins(10),type = "resample") #dependent sample which is 2-independent
#' independence.test(coins(10)[,2:3],type = "resample") #dependent sample which is 2-independent
#'
#' @export
independence.test = function(x,vec = 1:ncol(x),alpha = 0.05,type = "distribution_free",verbose = TRUE,...) {
  tm = total.multivariance(x,vec,...)

  warning("The function 'independence.test' is depreciated, use 'multivariance.test(...)$p.value < alpha' instead to compute the test result.")

  switch(type,
         distribution_free = {
           R = rejection.level(alpha)
           outtext = paste("\nDistribution free test (classical): The value of the test statistic is",tm,"and values above",R,"are rejected.\n")
           result = tm>R
         },
         resample = {
           p.value = resample.pvalue(tm,x=x,vec=vec,times = 300,type="total",...)
           outtext = paste("\nResampling test: The value of the test statistic is",tm,"and (by resampling) its p-value is",p.value,"\n")
           result = p.value<alpha
         },
         pearson_approx = {
           p.value = pearson.pvalue(x=x,vec=vec,type="total",...)
           outtext = paste("\nTest using Pearson's approximation: The value of the test statistic is",tm,"and its p-value is",p.value,"\n")
           result = p.value<alpha
         }
  )

  if (verbose) {
    cat(outtext)
    if (result) {cat("The hypothesis of independence is rejected.\n")
    }
    else {cat("The hypothesis of independence is NOT rejected.\n")
    }
  }
  invisible(result)
}

#' independence tests based on (total-/2-/3-) multivariance
#'
#' This performs the (specified by \code{type} and \code{p.value.type}) independence test for the columns of a sample matrix.
#'
#' @inheritParams cdms
#' @param type one of \code{"independence"}, \code{"pairwise independence"}, \code{"multi"}, \code{"total"}, \code{"m.multi.2"}, \code{"m.multi.3"}
#' @param p.value.type one of \code{"pearson_approx"}, \code{"distribution_free"}, \code{"resample"}, \code{"pearson_unif"}
#' @param verbose logical, if TRUE meaningful text output is generated.
#'
#' @return  A list with class "\code{htest}" containing the following components:
#' \describe{
#'   \item{\code{statistic}}{the value of the test statistic,}
#'   \item{\code{p.value}}{the p-value of the test statistic,}
#'   \item{\code{method}}{a character string indicating the type of test performed,}
#'   \item{\code{data.name}}{a character string giving the name(s) of the data.}
#' }
#' @details
#' For the use of \code{vec} see the examples below and the more detailed explanation of this argument for \code{\link{multivariance}}.
#'
#' The types \code{"independence"} and \code{"total"} are identical: an independence test is performed.
#'
#' Also the types \code{"pairwise independence"} and \code{"m.multi.2"} are identical:  a test of pairwise independence is performed.
#'
#' The type \code{"m.multi.3"}, performs a test for 3-independence, assuming pairwise independence. The type  \code{"multi"} performs a test for n-independence, assuming (n-1)-independence.
#'
#' There are several ways (determined by \code{p.value.type}) to estimate the p-value: The \code{"pearson_approx"} and \code{"resample"} are approximately sharp. The latter is based on a resampling approach and thus much slower. The \code{"distribution_free"} test might be very conservative, its p-value estimates are only valid for p-values lower than 0.215 - values above should be interpreted as "values larger than 0.215". Finally, \code{"pearson_unif"} uses fixed parameters in Pearson's estimate, it is only applicable for univariate uniformly distributed marginals
#'
#' All tests are performed using the standard euclidean distance. Other distances can be supplied via the \code{...}, see \code{\link{cdm}} for the accepted arguments.
#'
#' @references
#' For the theoretic background see the references given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' # an independence test
#' multivariance.test(dep_struct_several_26_100,p.value.type = "distribution_free") # conservative
#' multivariance.test(dep_struct_several_26_100,p.value.type = "resample") #sharp but slow
#' multivariance.test(dep_struct_several_26_100,p.value.type = "pearson_approx") #
#'
#' # as an example, all tests for one data set:
#' coins100 = coins(100)
#' for (ty in c("total","m.multi.2","m.multi.3","multi"))
#'  for (pvt in c("distribution_free","resample","pearson_approx"))
#'   print(multivariance.test(coins100,type=ty,p.value.type = pvt))
#'
#' # using the vec argument:
#' x = matrix(rnorm(50*6),ncol = 10) # a 50x6 data matrix
#' vec = c(1,2,3,4,5,6) # each column is treated as one variable
#' multivariance.test(x,vec,p.value.type = "distribution_free") # is the same as the default
#'
#' vec = c(1,2,2,1,3,1)
#' # column 1,4,6 are treated as one variable
#' # column 2,3 are treated as one variable
#' # column 5 is treated as one variable
#' multivariance.test(x,vec,p.value.type = "distribution_free")
#'
#' @export
multivariance.test = function(x,vec = 1:ncol(x),type = "total",p.value.type = "pearson_approx",verbose = TRUE,...) {

  if (is.data.frame(x)) stop("Input is a data.frame. Please provide a matrix or a list of doubly centered distance matrices.")

  if (p.value.type == "conservative") p.value.type = "distribution_free" #undocumented

  data.name = deparse(substitute(x))

  if ((p.value.type %in% c("pearson_approx","pearson_unif")) & (is.list(x))) {
    if (is.null(x$list.cdm)) {
      stop("'pearson_approx' requires the data. A list of doubly centered distance matrices is not sufficient.")
    } else {
      x.cdms.mu.bcd = x
      x = x$list.cdm
    }
  }

  if (type == "independence") type = "total"
  if (type == "pairwise independence") type = "m.multi.2"

  switch(type,
    multi = {
      method = "multivariance: test of n-independence (assuming (n-1)-independence)"
      fun = function(dat = x) multivariance(dat,vec = vec,...)
      statistic = c(multivariance = NA)
    },
    total = {
      method = "total multivariance: test of independence"
      fun = function(dat = x) total.multivariance(dat,vec = vec,...)
      statistic = c(total.multivariance = NA)
    },
    m.multi.2 = {
      method = "2-multivariance: test of pairwise independence"
      fun = function(dat = x) m.multivariance(dat,vec = vec,...)
      statistic = c(m.multivariance = NA)
    },
    m.multi.3 = {
      method = "3-multivariance: test of 3-independence (assuming pairwise independence)"
      fun = function(dat = x) m.multivariance(dat,vec = vec,m = 3,...)
      statistic = c(m.multivariance = NA)
    },
    {stop(paste("unkown type:",type))}
  )

  #method = "total multivariance: test of independence"
  #m = total.multivariance(x,vec,...)
  #statistic = c(total.multivariance = m)

  #statistic = NULL
  p.value = NULL
  #method = NULL
  parameter = NULL
  alternative = NULL

  switch(p.value.type,
    distribution_free = {
      method = paste0(method, "; distribution-free p-value (conservative for p-values in (0,0.215))")
      statistic[1] = fun()
      p.value = multivariance.pvalue(statistic)
    },
    resample = {
      method = paste0(method, "; resampling test")
      res = resample.multivariance(x=x,vec=vec,type=type,...)
      statistic[1] = res$original
      p.value = res$p.value

    },
    pearson_approx = {
      method = paste0(method, "; p-value via Pearson's approximation (approx. sharp)")
      # statistic[1] = fun()
      if (is.matrix(x)) {
        dots <- list(...) # we filter the argument lambda which might be used for total.multivariance, we keep everything else to preserve other error messages.
        x.cdms.mu.bcd = do.call('cdms.mu.bcd', c(list(x = x, vec = vec), dots[!(names(dots) %in% "lambda")]))
        # = cdms.mu.bcd(x,vec,...) # problem by extra argments in ... e.g. "lambda"
      p.value = pearson.pvalue(x=x.cdms.mu.bcd,vec=vec,type=type,...)
      statistic[1] = fun(x.cdms.mu.bcd$list.cdm)

      } else {
        statistic[1] = fun()
        p.value = pearson.pvalue(x=x.cdms.mu.bcd,vec=vec,type=type,...)
      }
    },
    pearson_unif = {
      method = paste0(method, "; p-value via Pearson's approximation (approx. sharp) using theoretic values based on the Euclidean distance. Only valid for 1-dim uniform marginals.")
      statistic[1] = fun()
      if (is.matrix(x)) {
        p.value = pearson.pvalue.unif(x=x,vec=vec,type=type,...)
      } else {
        n = length(x) # x is now list.cdm

        unif.cmb = list(list.cdm = x,
          mu = matrix(rep(c(1/3,2/45,8/945),n),ncol=n),
          bcd = matrix(rep(c(1/6,7/60,1/9),n),ncol=n),
          mean =rep(1/3,n))

        p.value = pearson.pvalue(x=unif.cmb,vec=vec,type=type,...)
      }
    },
    {stop(paste("unkown p.value.type:",p.value.type))}
  )

  result = list(statistic = statistic, p.value = p.value, method = method, data.name = data.name) #, parameter = parameter, alternative = alternative
  class(result) = "htest"

  return(result)
}

################ Bias corrected #########

#' distance matrix
#'
#' # currently only used for the bias corrected multicorrelations
#' # it might be used globally to remove redundancies
#'
#' @keywords internal
dm = function(x, psi = NULL, p = NULL, isotropic = FALSE, external.dm.fun = NULL) {
  if (!is.matrix(x)) x = as.matrix(x)
  if (is.null(psi) & is.null(p) & is.null(external.dm.fun)) {
    #dm = dist.to.matrix(stats::dist(x,method="euclidean"))
    #DEVELOPING NOTE: here as.matrix was slow, dist.to.matrix is faster. Instead one could just use the vector....
    # even faster for the euclidean case is fastdist defined via Rcpp

    return(fastdist(x))

  } else {
    if (!is.null(p)) {
      if ((p<1) || (p>2)) warning("p is not in [1,2]. \n")
      dm = dist.to.matrix(stats::dist(x,method="minkowski", p = p))
    } else { # case: psi is given
      if (!is.null(external.dm.fun)) {
        dm = external.dm.fun(x)
      } else {
        if (isotropic) {
          #dm = psi(dist.to.matrix(stats::dist(x,method="euclidean")))
          dm = psi(fastdist(x))
        } else {
          n = nrow(x)
          d = ncol(x)
          dm = matrix(apply(cbind(x[rep(1:n,n),],x[rep(1:n,each = n),]), #create all combinations
            1, # apply to each row
            function(y) psi(y[1:d], y[(d+1):(2*d)])),nrow = n)
          # DEVELOPING NOTE: could double the speed if only the upper triangular matrix is computed, using the idea of dist.to.matrix
        }
      }
    }
  }
}


#' list of distance matrices
#'
#' # currently only used for the bias corrected multicorrelations
#' # it might be used globally, to remove redundancies
#'
#' @keywords internal
dms = function(x,vec = 1:ncol(x),...) {
  return(lapply(1:max(vec),function(i) dm(x[,vec == i,drop = FALSE],...)))
}

#' bias corrected total multicorrelations
#'
#' @keywords internal
multicorrelation.bias.corrected = function(x,
  vec = 1:ncol(x), squared = FALSE, type = "all",...)
{
  # convention: 0/0 = 0

  if (type == "total.upper.lower")
    return(total.multicorrelation.bias.corrected.upper.lower(x,vec,squared,...))

  if ((type == "m.multi.2") | (type == "pairwise"))
    return(pairwise.multicorrelation.bias.corrected(x,vec,squared,...))

  if (type == "total.upper")
    return(total.multicorrelation.bias.corrected.upper(x,vec,squared,...))


  ## below slow for 'all'
  ## undocumented

  if (!(type == "all")) stop(paste0("unknown type: ",type))


  N = nrow(x)
  n = max(vec)
  const.lambda = 1 # for lambda total multicorrelation

  #dms.list = lapply(1:n,function(i) fastdist(as.matrix(x[,vec == i])))

  dms.list = dms(x,vec,...)

  cdms.list = lapply(dms.list,doubleCenterBiasCorrected)
  #    double.centerUB(dms.list,normalize = FALSE)

  values = c(unnormalized = NA,upper = NA, normalized = NA, lower = NA, pairwise = NA)

  for (typ in c("lower","upper"))  {
    switch(typ,
      lower = {
        fun = function(x) abs(x)^n # for lower
        fun.inv = function(x) x^(1/n) # for lower
      },
      upper = {
        fun = function(x) x^2 # for upper and pairwise
        fun.inv = function(x) sqrt(x)
      })

    # by the following roots a lot of precision is lost
    norms.list = lapply(cdms.list, function(x) (N/((N-3))*mean(fun(x))))

    too.small = norms.list < .Machine$double.eps
    if (isTRUE(any(c(too.small)))) {
      warning("some norming constant is approximately 0")
      norms.list[too.small] = 0
    }

    #if ((n==2) & isTRUE(any(list.norms.upper == 0)))
    #  return(c(0,0,0))

    cdms.normed.list = lapply(1:length(cdms.list), function(i) cdms.list[[i]]/fun.inv(norms.list[[i]]))


    global.factor = ((1+const.lambda)^n - n*const.lambda^(n-1) - const.lambda^n)

    values[typ] = N/((N-3))*(mean(Reduce("*", lapply(cdms.normed.list,function(y) const.lambda + y))) - const.lambda^n)/global.factor

    if (typ == "upper") {
      # pairwise multicorrelation / 2-multicorrelation (needs cdms.normed.list used for upper)
      Asum = Reduce("+", cdms.normed.list)
      A2sum = Reduce("+", lapply(cdms.normed.list,function(y) y^2))

      values["pairwise"] = N/((N-3))*mean(Asum^2 - A2sum)/(2*(choose(n,2)))
    }
  }

  ### full normalized and unnormalized

  if (n > 15) warning(paste0("Normalized and unnormalize multicorreltions become very slow for large number of variables. For the current dataset in these need to compute ",2^n-n-1," terms. In contrast 'upper.lower' requires only 1 of such terms."))

  for (typ in c("normalized","unnormalized")) {

    switch(typ,
      normalized = {
        fun = function(x) abs(x) # for normalizied
      },
      unnormalized = {
        fun = function(x) x # for unnormalized
      })

    vals = NULL
    for (k in 2:n) { # for each k-tuple
      norms.list = lapply(cdms.list, function(x) (N/((N-3))*mean(fun(x)^k))^(1/k))

      if (typ == "normalized") {
        norms.list = lapply(norms.list,function(x) ifelse(x < .Machine$double.eps^(1/k),0,x))
      }

      cdms.normed.list = lapply(1:length(cdms.list), function(i) cdms.list[[i]]/(norms.list[[i]]))

      S = utils::combn(n,k) # all possible k-tuples

      vals = c(vals,apply(S,2,function(s) (N/((N-3))*(mean(Reduce("*", cdms.normed.list[s]))))))
    }

    values[typ] = mean(vals)
  }

  if (!squared) values = signed.sqrt(values)

  # names(values) = paste0("Mcor.",names(values))
  return(values)
}

#' # included for speed. it is faster than upper.lower
#'
#' @keywords internal
total.multicorrelation.bias.corrected.upper = function(x,vec = 1:ncol(x),squared = FALSE,...){
  n = max(vec)
  global.factor = 2^n-n-1
  N = nrow(x)

  dms.list = dms(x,vec,...)

  cdms.list = lapply(dms.list,function(mat) doubleCenterBiasCorrectedUpper(mat))

  ret = c(total.upper = ((N-1)/(N-3)*(mean(Reduce("*", lapply(cdms.list
    ,function(y) 1 + y))) - 1)/global.factor))

  if (abs(ret) < .Machine$double.eps) ret[1] = 0
  if (!squared) ret = signed.sqrt(ret)

  return(ret)
}

#'
#'
#' @keywords internal
total.multicorrelation.bias.corrected.upper.lower = function(x,vec = 1:ncol(x),squared = FALSE,...){
  n = max(vec)
  const.lambda = 1

  # global.factor = ((1+const.lambda)^n - n*const.lambda^(n-1) - const.lambda^n) # for lambda total mutlicorrelation
  global.factor = 2^n - n - 1

  N = nrow(x)
  #dms.list = lapply(1:n,function(i) fastdist(as.matrix(x[,vec == i])))
  dms.list = dms(x,vec,...)

  cdms.ul.list = lapply(dms.list
    ,function(mat) doubleCenterBiasCorrectedUpperLower(mat,n))

  ret = c(total.upper = ((N-1)/(N-3)*(mean(Reduce("*", lapply(
    cdms.ul.list,function(y) 1 + y$out/y$upper))) - 1)/global.factor),
    total.lower = ((N-1)/(N-3)*(mean(Reduce("*", lapply(
      cdms.ul.list,function(y) 1 + y$out/y$lower))) - 1)/global.factor)
  )

  if (isTRUE(abs(ret[1]) < .Machine$double.eps)) ret[1] = 0
  if (isTRUE(abs(ret[2]) < .Machine$double.eps)) ret[2] = 0
  if (!squared) ret = signed.sqrt(ret)

  return(ret)


  #  For lambda total multicorrelation:
  #  return(
  #   c(upper = signed.sqrt((N-1)/(N-3)*(mean(Reduce("*", lapply(  cdms.ul.list,function(y) const.lambda + y$out/y$upper))) - const.lambda^n)/global.factor),
  #   lower = signed.sqrt((N-1)/(N-3)*(mean(Reduce("*", lapply(
  #   cdms.ul.list,function(y) const.lambda + y$out/y$lower))) - const.lambda^n)/global.factor)
  # ))
}

#' pairwise multicorrelation
#'
#' @keywords internal
pairwise.multicorrelation.bias.corrected = function(x,vec = 1:ncol(x),squared = FALSE,...){
  n = max(vec)
  N = nrow(x)
  #dms.list = lapply(1:n,function(i) fastdist(as.matrix(x[,vec == i])))
  dms.list = dms(x,vec,...)

  cdms.list = lapply(dms.list
    ,function(mat) doubleCenterBiasCorrectedUpper(mat))

  Asum = Reduce("+", cdms.list)
  A2sum = Reduce("+", lapply(cdms.list,function(y) y^2))

  ret = c(pairwise = (N-1)/((N-3))*mean(Asum^2 - A2sum)/(2*(choose(n,2))))

  if (abs(ret) < .Machine$double.eps) {
    #warning(paste0(ret," rounded to 0."))
    ret[1] = 0
  }
  if (!squared) ret = signed.sqrt(ret)

  return(ret)

}




################# Example data ##########

#' dependence example: tetrahedron sampling
#'
#' This function creates samples of a tetrahedron-dice colored r, g, b and rgb. Each sample indicates if for the thrown dice the colors r, g and b are contained on the bottom side of the dice.
#'
#' @param N number of samples
#' @return It returns the samples of the events r, g and b as rows of a \code{N} by 3 matrix (the first column corresponds to r, the second to g,...). TRUE indicates that this color is on the bottom side of the dice. The columns are dependent but 2-independent.
#' @examples
#' tetrahedron(10)
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
tetrahedron = function(N = 1000) {
  # rolls the tetrahedron with sides r,g,b,rgb
  # this is an explicit implementation, the distribution is the same
  # as for 'coins(N,3,type = "even")'
  side = sample.int(4,N,replace=TRUE)
  x = (side == 1)|(side == 4)
  y = (side == 2)|(side == 4)
  z = (side == 3)|(side == 4)
  return(unname(cbind(x,y,z)))
}

#' dependence example: k-independent coin sampling
#'
#' This function creates samples which are dependent but k-independent.
#' @param N number of samples
#' @param k each k-tuple will be independent
#' @param type one of \code{"even"} or \code{"odd"}
#' @return It returns the samples as rows of an \code{N} by \code{k+1} matrix. The columns are dependent but k-independent.
#'
#' @details Throw \code{k} independent fair coins. Now consider
#' the k+1 events: The first shows head, the second shows head,... the \code{k}-th shows head,
#' there is an \code{even} (or \code{odd} as selected via \code{type}) number of heads. Each row
#' contains the state of these k+1 events.
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' coins(200,4)
#'
#' @export
coins = function(N = 1000, k = 2, type = "even") {
  k.coins = matrix(sample.int(2,k*N,replace = TRUE)-1,ncol = k)
  switch(type,
         even = { d = ((rowSums(as.matrix(k.coins)) %% 2) == 0) },
         odd =  { d = ((rowSums(as.matrix(k.coins)) %% 2) == 1) }
  )
  return(unname(cbind(k.coins,d)))
}


# Resampling ######

#' resample the columns of a matrix
#' @param x matrix
#' @param vec vector, indicates which columns belong together
#' @param replace boolean, sampling with or without replacement
#' @param incl.first boolean, if \code{TRUE} also the first component is resampled
#'
#' @return Returns a matrix with the same dimensions as \code{x}. The columns are resampled from the original columns. The resampling is done with replacement (\code{replace = TRUE}) or without (\code{replace = FALSE}). Columns which belong together (indicated by vec) are resampled identically, i.e., all values in rows of these are kept together.
#'
#' @examples
#' sample.cols(matrix(1:15,nrow = 5),vec = c(1,1,2))
#'
#' @export
sample.cols = function(x,vec = 1:ncol(x),replace = TRUE,incl.first = TRUE) {
  if (anyNA(vec)) vec = 1:ncol(x)
  N = nrow(x)
  xnew = x
  for (i in (2-incl.first):max(vec)) {
    neworder = sample.int(N,replace = replace)
    xnew[,vec == i] = x[neworder,vec == i]
    #print(neworder)
  }
  return(xnew)
}

#' resamples doubly centered distance matrices
#' @param list.cdm a list of doubly centered distance matrices
#' @inheritParams sample.cols
#'
#' @return Returns a list of doubly centered distance matrices, each matrix corresponds to the resampled columns of the corresponding sample, using resampling with replacement (bootstrap) or without replacement (permutations).
#'
#' @export
sample.cdms = function(list.cdm,replace = FALSE, incl.first = FALSE) {
  N = nrow(list.cdm[[1]])
  # n = dim(array.cdm)[3]
  #
  # for (i in 2:n) { # for i == 1 the resampling is not necessary.
  #   neworder = sample.int(N)
  #   array.cdm[,,i] = array.cdm[neworder,neworder,i]
  # }
  # return(array.cdm)
  #
  #return(lapply(list.cdm, function(y) {neworder = sample.int(N, replace = replace) y[neworder,neworder]})) # resamples also the first
  return(
  lapply(1:length(list.cdm),
    function(i) {
      if (i > (!incl.first)) {
        neworder = sample.int(N,replace = replace)
        #print(neworder)
        return(list.cdm[[i]][neworder,neworder])
      } else { # no resampling for the first
        return(list.cdm[[i]])
      }
    })
  )
}

#' resampling (total /m-) multivariance
#'
#' The distribution of the test statistic under the hypothesis of independence is required for the independence tests. This function generates approximate samples of this distribution either by sampling without replacement (permutations) or by sampling with replacement (bootstrap).
#'
#' @details
#' The resampling is done by sampling from the original data either without replacement (\code{"permutation"}) or with replacement (\code{"bootstrap"}). Using resampling without replacement is (much) faster (due to special identities which only hold in this case).
#'
#' For convenience also the actual (total /m-) multivariance is computed and its p-value.
#'
#' @param x matrix, the rows should be iid samples
#' @param vec vector, which indicates which columns of \code{x} are treated together as one sample
#' @param times integer, number of samples to generate
#' @param type one of \code{"multi","total","m.multi.2","m.multi.3","all"}
#' @param resample.type one of \code{"permutation", "bootstrap"}. The samples are generated without replacement (permutations) or with replacement (bootstrap).
#' @param ... is passed to \code{\link{cdms}, \link{multivariance}, \link{total.multivariance}, \link{m.multivariance}}, respectively.
#'
#' @return A list with elements
#' \describe{
#'   \item{\code{resampled}}{the (total/m-)multivariances of the resampled data,}
#'   \item{\code{original}}{the (total/m-)multivariance of the original data,}
#'   \item{\code{p.value}}{the p-value of the original data, computed using the resampled data}
#' }
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' re.m = resample.multivariance(matrix(rnorm(30*2),nrow = 30),
#'                         type= "multi",times = 300)$resampled
#' curve(ecdf(re.m)(x), xlim = c(0,4),main = "empirical distribution of the test statistic under H_0")
#' @export
resample.multivariance = function(x,vec = 1:ncol(x),times = 300,type = "multi",resample.type = "permutation",...) {

  res.vec = vec # the vec argument used for the resampling methods; some methods use the doubly centered distance matrices, and there vec has to be 1:max(vec)
  switch(resample.type,
         #distinct.permutation = {resample = function() matrix(x[derangements.without.fixpoint(N,n, distinctcols,vec)],ncol = n)},
       permutation = {
         if (is.list(x)) {
           list.cdm = x[vec]
           res.vec = 1:length(vec)
         } else {
           dots <- list(...) # we filter the argument lambda which might be used for total.multivariance, we keep everything else to preserve other error messages.
           #argnames <- names(formals(cdm))
           #list.cdm = do.call('cdms', c(list(x = x, vec = vec), dots[names(dots) %in% argnames]))
           list.cdm = do.call('cdms', c(list(x = x, vec = vec), dots[!(names(dots) %in% "lambda")]))

           # list.cdm = cdms(x,vec,...) #!! TODO: note there is the argument lambda for total.multivariance this should be excluded here!
           res.vec = 1:max(vec) # all distance matrices shall be used.
         }
         resample = function() sample.cdms(list.cdm, replace = FALSE)}, # resampling of the doubly centered distance matrices - this is faster than permutation.orig
    bootstrap = {
         warning("Note that bootstrap resampling is (for multivariance) much slower than permutation resampling. Because certain simplifying identities fail to hold. \n")
         res.vec = 1:max(vec) # all distance matrices shall be used.
         resample = function() cdms(sample.cols(x,vec,replace =TRUE))
      },
    permutation.orig = {resample = function() sample.cols(x,vec,replace =FALSE)}, # resampling of the original data
    {stop(paste("unkown resample.type:",resample.type))}
  )

  switch(type,
         multi = {fun = function (x,v) multivariance(x,vec = v,...)},
         total = {fun = function (x,v) total.multivariance(x,vec = v,...)},
         m.multi.2 = {fun = function (x,v) m.multivariance(x,vec = v,...)},
         m.multi.3 = {fun = function (x,v) m.multivariance(x,vec = v,m = 3,...)},
         all = {fun = function (x,v) multivariances.all(x,vec = v,...)}#doAll = TRUE}
    ,
    {stop(paste("unkown type:",type))}
  )

  results = matrix(,nrow = times, ncol = 1 + (type == "all")*3)
  for (i in 1:times) {# we use a for loop instead of replicate to prevent trouble with '...'
    results[i,] = fun(resample(),res.vec)
  }

  multi = fun(x,vec)
  p.value = rowSums(t(results) >= multi)/times
  names(p.value) = names(multi)

  invisible(list(resampled = results,
                 original = multi,
                 p.value = p.value))

}


#' rejection level via resampling
#'
#' Uses the resample method to sample from the test statistic under the hypothesis of independence. The alpha quantile of these samples is returned.
#'
#' @param alpha numeric, the significance value
#' @param ... passed to \code{\link{resample.multivariance}}. Required is the data matrix \code{x}.
#'
#'@references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' resample.rejection.level(0.05,matrix(rnorm(30*2),nrow = 30))
#' resample.rejection.level(0.05,matrix(rnorm(30*3),nrow = 30),vec = c(1,1,2))
#'
#' @export
resample.rejection.level = function(alpha = 0.05,...) {
  samples = resample.multivariance(...)$resampled
  stats::quantile(samples,probs= 1-alpha)
}

#' p-value via resampling
#'
#' Use a resampling method to generate samples of the test statistic under the hypothesis of independence. Based on these the p.value of a given value of a test statistic is computed.
#'
#' @return It returns 1 minus the value of the empirical distribution function of the resampling samples evaluated at the given value.

#' @param value numeric, the value of (total-/m-)multivariance for which the p-value shall be computed
#' @param ... passed to \code{\link{resample.multivariance}}. Required is the data matrix \code{x}.
#'
#' @details This function is useful if a p-value of a test statistic shall be computed based on the resampling values of the test statistic of a different sample. For the p-value based on the same sample \code{\link{resample.multivariance}(...)$p.value} is sufficient.
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#' x = coins(100)
#' resample.pvalue(multivariance(x),x=x,times = 300)
#' resample.pvalue(multivariances.all(x),x=x,times = 300,type = "all")
#'
#' @export
resample.pvalue = function(value,...){
  samples = resample.multivariance(...)$resampled
  #sum(samples >= value) / length(samples)
  #slower:  1-stats::ecdf(samples)(value)
  result = rowSums(t(samples) >= value)/ncol(t(samples))
  names(result) = names(value)
  return(result)
}

##### Moments ####

#' given the distance matrix the unbiased estimate for mu3 is computed
#' @keywords internal
mu3.unbiased = function(B,b2ob = sum(tcrossprod(B)*B)) {
  #B2 = tcrossprod(B) # B%*%B ; note that the matrices are symmetric; this seems to be a faster implementation.
  #B3 = tcrossprod(B,B2) # B2%*%B

  cb = colSums(B)
  cbob = colSums(B^2)

  b = sum(cb) #sum(B)
  b2 = sum(cb^2) #sum(B2)
  b3 = cb%*%B%*%cb #sum(B3)

  bob = sum(cbob) #sum(B^2)
  bobob = sum(B^2*B)
  lboblb = sum(cbob * cb) #sum(tcrossprod(B,B^2))

  csbo3 = sum(cb^2*cb)

  N = ncol(B)

  facN = function(k) 1/(prod(N:(N-k)))

  f = facN(3)*(b3-b2ob-2*lboblb+bobob)

  ei = facN(2)*b2ob

  y = facN(4)*(-4*b3 - 4*bobob - 2*csbo3 + 2*b2ob + 10*lboblb + b*b2 - b*bob)

  u = facN(5)*(-48*lboblb - 8*b2ob + 16*bobob + 16*csbo3 + 24*b3 - 12*b*b2 + 6*b*bob + b^2*b)

  mu3 = - ei + 3*f - 3*y + u

  return(mu3)
}


#' given the sample of a single variable the doubly centered distance matrix, mu (the limit moments) and bcd (the terms for the finite sample moments) are computed
#'
#' The normalization should be postponed to the moment calculation.
#'
# NOTE: speedup might be possible by incorporating mu3 and some matrix-norm-identities
#' @keywords internal
cdm.mu.bcd = function(x, normalize = FALSE, psi = NULL, p = NULL, isotropic = FALSE, unbiased.moments = TRUE, external.dm.fun = NULL) {
  if (normalize) stop("normalized not implemented")

  ##### lines copied from cmd
  if (!is.matrix(x)) x = as.matrix(x)
  if (is.null(psi) & is.null(p) & is.null(external.dm.fun)) {
    #dm = dist.to.matrix(stats::dist(x,method="euclidean"))
    #DEVELOPING NOTE: here as.matrix was slow, dist.to.matrix is faster. Instead one could just use the vector....
    # even faster for the euclidean case is fastdist defined via Rcpp

    dm = fastdist(x)
  } else {
    if (!is.null(p)) {
      if ((p<1) || (p>2)) warning("p is not in [1,2]. \n")
      dm = dist.to.matrix(stats::dist(x,method="minkowski", p = p))
    } else { # case: psi is given
      if (!is.null(external.dm.fun)) {
        dm = external.dm.fun(as.matrix(x))
      } else {
        if (isotropic) {
          #dm = psi(dist.to.matrix(stats::dist(x,method="euclidean")))
          dm = psi(fastdist(as.matrix(x)))
        } else {
          x = as.matrix(x)
          n = nrow(x)
          d = ncol(x)
          dm = matrix(apply(cbind(x[rep(1:n,n),],x[rep(1:n,each = n),]), #create all combinations
                            1, # apply to each row
                            function(y) psi(y[1:d], y[(d+1):(2*d)])),nrow = n)
          # DEVELOPING NOTE: could double the speed if only the upper triangular matrix is computed, using the idea of dist.to.matrix
        }
      }
    }
  }
  ### end copy of cdm


  colm = colMeans(dm)
  m = mean(colm)  # equals mean(dm)

  # if (unbiased & normalize) warning("Unbiased normalized cdm, m2, m3 not implemented. \n")

  if (m == 0) warning("It seems that one variable is constant. Constants are always independent. \n")

  #cdm = (-dm + outer(colm,colm, FUN ="+") - m)
  cdm = doubleCenterSymMat(dm,normalize)

  ###### for mu (biased)
  B = dm
  #if (normalize) B = B/mean(B) # normalize

  B2 = tcrossprod(B) #mymm(B,B) # B%*%B note that the matrices are symmetric; this seems to be a faster implementation.
  #B3 = tcrossprod(B,B2) #mymm(B,B2) #B2%*%B

  N = nrow(B)

  cb = colSums(B)
  mB = 1/N^2 * sum(cb) #mean(B)
  mB2 = 1/N^2 * sum(cb^2) #mean(B2)
  mBsq = mean(B^2)

  mB2oB = mean(B2*B)

  ###### for bcd (and mu unbiased)
  if (unbiased.moments) {
    #    bcd = c(N^2/(N*(N-1))*mBsq,N^2/(N*(N-1)*(N-2)) * (mB2-mBsq),N^2/(N*(N-1)*(N-2)*(N-3))*(N^2*mB^2+2*mBsq-4*mB2)) #b,c,d
    bcd = c(N/((N-1))*mBsq, N/((N-1)*(N-2)) * (mB2-mBsq), N/((N-1)*(N-2)*(N-3))*(N^2*mB^2+2*mBsq-4*mB2)) #b,c,d

  #  if (normalize) {
  #
  #  } else {
      mB = N/(N-1)*mB #!! since we are in the case without normalization
   # }
    m2 = sum(bcd * c(1,-2,1))
    m3 = mu3.unbiased(B,b2ob = mB2oB*N^2)

    mu = c(mB,m2,m3) # unbiased estimators for limit moments

  } else {

    bcd = c(mBsq,1/N*mB2,mB^2) # biased estimators for b,c,d

    mB3 = 1/N^2 * cb%*%B%*%cb #mean(B3)

    m2 = mBsq-2/N*mB2+mB^2
    m3 = -1/N*mB2oB + 3/N^2*mB3-3/N*mB2*mB + mB^2*mB

    mu = c(mB,m2,m3) # biased estimators for the limit moments

  }

  return(list(cdm = cdm, mu = mu, bcd = bcd, mean = m))
}

#' computes the doubly centered distance matrices, mus and bcds
#' @param x matrix, each row is a sample
#' @param vec vector which indicates which columns are treated as one sample
#' @param membership depreciated. Now use \code{vec}.
#' @param ... these are passed to \code{\link{cdm}}
#'
#' @return  A list containing the following components:
#' \describe{
#'   \item{\code{list.cdm}}{list of the doubly centered distance matrices - these are always normalized if 'cdm.normalize = TRUE'!!!,}
#'   \item{\code{mu}}{matrix with the limit moments in a column for each variable,}
#'   \item{\code{bcd}}{matrix with b, c, d (which are required for the computation of the finite sample moments) in columns for each variable,}
#'   \item{\code{mean}}{vector with the mean of each distance matrix.}
#' }
#'
#' @keywords internal
cdms.mu.bcd = function(x,vec = 1:ncol(x),membership = NULL,cdm.normalize = TRUE,...) {
  if (!is.null(membership)) {
    vec = membership
    warning("Use 'vec' instead of 'membership' as argument to 'cdms'. 'membership' is depreciated. \n")
  }
  if (anyNA(vec)) vec = 1:ncol(x)
  n = max(vec)
  N = nrow(x)
  #array.cdm = array(,dim = c(N,N,n))
  list.cdm = list()
  mu = matrix(ncol = n,nrow = 3)
  bcd = matrix(ncol = n,nrow = 3)
  m = numeric(n)
  for (i in 1:n) {
    res = cdm.mu.bcd(x[,(vec == i),drop = FALSE],...)
    #array.cdm[,,i] = res$cdm
    if((res$mean ==0)|(!cdm.normalize)) {
      list.cdm[[i]] = res$cdm
    } else {
      list.cdm[[i]] = res$cdm/res$mean # returns always the normalized matrices!
    }
    mu[,i] = res$mu
    bcd[,i] = res$bcd
    m[i] = res$mean
  }
  return(list(list.cdm = list.cdm, mu = mu, bcd = bcd, mean = m))
}



coef.7.cases = structure(c(0, 0, 0, 0, 0, 0, 0, 0, -6, -2, 8, -1, -4, 4, 1,
                           0, 11, 2, -12, 1, 4, -6, 0, -1, -6, 0, 4, 0, 0, 2, 0, 2, 1, 0,                           0, 0, 0, 0, 0, 0), .Dim = c(8L, 5L))

coef.array = structure(c(0, 0, 0, 0, 0, 6, -24, 18, -1, 2, -4, 3, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, -24, 18, -1, -2, 12, -9, 0,
                         -2, 4, -2, 0, 1, -2, 1, 0, 0, 0, 0, 0, 6, -24, 18, -1, 0, 4,
                         -3, 0, -1, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, -24, 18, -1,
                         -2, 12, -9, 2, 0, 0, -2, -1, 0, 0, 1, 0, 0, 0, 0, 0, 6, -24,
                         18, -1, -4, 20, -15, 1, 0, -4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         6, -24, 18, -1, 0, 4, -3, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 6, -24, 18, -1, -10, 44, -33, 2, 4, -24, 18, -1, 0, 4, -3), .Dim = c(4L,
                                                                                                 5L, 7L))

#' Computes the explicit coefficients for the finite sample variance for a sample of size N
#' @keywords internal
N.coefficients = function(N) {
  freq.o.c = as.vector(coef.7.cases %*% N^(0:4))
  case.coef = t(apply(coef.array,3,function(x) x %*%(N^{0:4})))
  dimnames(case.coef)[[2]]=c("a","b","c","d")
  #  names(freq.o.c)=c("dd","bb","cc")
  return(list(freq.o.c = freq.o.c, case.coef = case.coef))
}

#' functions which are required for the calculation of the finite sample expectation and variance for m-multivariance and total multivariance
#' @keywords internal
d2 = function(a,b) sum(a)*sum(b)-sum(a*b)
d3 = function(a,b,c) sum(c)*d2(a,b) - d2(a*c,b) - d2(a,b*c)
d4 = function(a,b,c,d) sum(d)*d3(a,b,c)- d3(a*d,b,c)- d3(a,b*d,c)- d3(a,b,c*d)

G2 = function(a,b,c) {
  d2(c,c)/2 + d3(a,b,c) + d4(a,a,b,b)/4
}

d5 = function(a,b,c,d,e) sum(e)*d4(a,b,c,d) - d4(e*a,b,c,d) - d4(a,e*b,c,d) -d4(a,b,e*c,d) -d4(a,b,c,e*d)
d6 = function(a,b,c,d,e,f) sum(f)*d5(a,b,c,d,e) - d5(f*a,b,c,d,e) - d5(a,f*b,c,d,e) -d5(a,b,f*c,d,e) -d5(a,b,c,f*d,e) - d5(a,b,c,d,f*e)

G3 = function(a,b,c) {
  d3(c,c,c)/6 + d4(a,b,c,c)/2 + d5(a,a,b,b,c)/4 + d6(a,a,a,b,b,b)/36
}

Gt = function(a,b,c) {
  return(
    prod(a+b+c+1) - prod(b+1) *(1 + sum((a+c)/(b+1))) - prod(a+1) *(1 + sum((b+c)/(a+1)))
    + 1 + sum(a)*sum(b) - sum(a*b) + sum(a+b+c)
  )
}

#' This is the function GC which is required for the computation of the finite sample variance for m and total multivariance
#' @keywords internal
sums.of.products = function(a,b,c, type = "multi") {
  switch(type,
         multi = { temp = prod(c)},
         m.multi.2 = { temp = G2(a,b,c)},
         m.multi.3 = { temp = G3(a,b,c)},
         total = { temp = Gt(a,b,c)},
         {stop(paste("unkown type:",type))}
  )
  return(temp)
}

#' computes the moments as required for Pearson's approximation
#'
#' @param N sample size
#' @param bcd an array with b c d
#' @param mu the limit moments
#' @param mmean the means of the distance matrices
#' @details
#' Note: It is currently only implemented for the case of normalized multivariance, i.e., the OUTPUT values correspond to normalized multivariance!
# COMMENT  But also note that the skewness is (appart from the scaling with respect to the number of summands) invariant with respect to normalization.
#'
#' @keywords internal
moments.for.pearson = function(N, bcd, mu, mmean, type = "multi") {

  switch(type,
         multi = {fun = function (x) prod(x)},
         total = {fun = function (x) prod(x+1)-sum(x)-1},
         m.multi.2 = {fun = function (x) (sum(x)^2-sum(x^2))/2},
         m.multi.3 = {fun = function (x) (sum(x)^3-3*sum(x)*sum(x^2)+2*sum(x^2*x))/6},
         {stop(paste("unkown type:",type))}
  )

  #limit.variance = 2*fun(mu[2,]/mmean^2) # variance
  #limit.skewness = 8*fun(mu[3,]/mmean^3)/limit.variance^(3/2) #skewness

  one.over.mmean = 1/mmean
  one.over.mmean[mmean == 0] = 0 # convention 0/0 = 0 (since the corresponding mu[i,] below are 0)

  limit.variance = 2*fun(mu[2,]*one.over.mmean^2) # variance
  limit.skewness = 8*fun(mu[3,]*one.over.mmean^3)/limit.variance^(3/2) #skewness

  if (is.nan(limit.skewness)) { # convention 0/0 = 0
    limit.skewness = 0
  }

  n = dim(mu)[2]

  switch(type,
         multi = {scalevec = 1},
         total = {scalevec = (2^n-n-1)},
         m.multi.2 = {scalevec = (choose(n,2))},
         m.multi.3 = {scalevec = (choose(n,3))},
  )
  scalevec2 = scalevec^2

  res = N.coefficients(N)

  muvec = mu[1,]

  bcdsums = (res$case.coef[,2:4]/N^4)%*%bcd[1:3,] #bbi+cci+ddi

  ### normalized:

  sumh2 = (res$freq.o.c[c(2,3,1)]/N^4)%*%bcd[1:3,] #(C(N,2)bi+C(N,3)ci+C(N,1)di)/N^4

  #bcdsums.normalized = t(apply(bcdsums,1,function(x) x/sumh2))

  one.over.sumh2 = 1/sumh2
  one.over.sumh2[sumh2 == 0] = 0


  bcdsums.normalized = t(apply(bcdsums,1,function(x) x*one.over.sumh2))


  one = rep(1,n)

  one[mmean == 0] = 0 # fixing 0/0

  E.no = (fun(one) + (N-1)*fun((-1/(N-1))*one))/scalevec # theoretic finite sample expectation for normalized multivariances - otherwise 'one' has to be replaced by 'muvec'

  # E2.no is the second moment (unbiased estimated)
  s.no = numeric(7)
  for (i in 1:3) s.no[i] = res$freq.o.c[i]/N^2*sums.of.products((-1/(N-1))*one,(-1/(N-1))*one,bcdsums.normalized[i,],type = type)
  for (i in c(4,7)) s.no[i] = res$freq.o.c[i]/N^2*sums.of.products(one,one,bcdsums.normalized[i,],type = type)
  for (i in c(5,6)) s.no[i] = res$freq.o.c[i]/N^2*sums.of.products(one,(-1/(N-1))*one,bcdsums.normalized[i,],type = type)
  E2.no = sum(s.no/scalevec2)

  Esq.no = (E.no)^2 # this is the same as the following (modulo numeric tolerance), since for normalized multivariance the expectation does not depend on the estimated marginals...

  # Esq.no is an unbiased estimate for the first moment squared.
#  r.no = numeric(7)
#  for (i in 1:3) r.no[i]  = res$freq.o.c[i]/N^2*sums.of.products((-1/(N-1))*one,(-1/(N-1))*one,(-1/(N-1))^2*one,type = type)
#  for (i in c(4,7)) r.no[i] = res$freq.o.c[i]/N^2*sums.of.products(one,one,one,type = type)
#  for (i in c(5,6)) r.no[i] = res$freq.o.c[i]/N^2*sums.of.products(one,(-1/(N-1))*one,(-1/(N-1))*one,type = type)
#  Esq.no = sum(r.no/scalevec2)

  Var.no = E2.no-Esq.no # finite sample variance

  return(c(E.no,Var.no,limit.skewness))
}

#' approximate distribution function of a Gaussian quadratic form
#'
#' Approximation of the of the value of the distribution function of a Gaussian quadratic form based on its first three moments.
#'
#' @param x value at which the distribution function is to be evaluated
#' @param moment vector with the mean, variance and skewness of the quadratic form
#' @param lower.tail logical, indicating of the lower or upper tail of the distribution function should be calculated
#' @param verbose logical, if \code{TRUE} a warning is issued if negative moments are sanitized to 0.
#'
#' @details This is Pearson's approximation for Gaussian quadratic forms as stated in [4] (equation (4.65) in arXiv:1808.07280v2)
#'
#' @references
#' For the theoretic background see the reference [4] given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
pearson.qf = function(x,moment,lower.tail = TRUE,verbose = FALSE) {
  if (any(is.na(moment[1:3]))) return(NA) # e.g. in the case where 3-multivariance is tested but only two variables are available

  if (any(moment[1:2] <= 0)) {
    if (verbose) warning("The mean or variance was negative, using 0 instead.")
    if (lower.tail) {
      return(as.numeric(x <= max(moment[1],0)))
    } else {
      return(as.numeric(x >= max(moment[1],0))) # for a discrete distribution the p.value (upper.tail) includes the value itself
    }
  } else {
    if (moment[3] <= 0) {
      if (verbose) warning("The skewness was negative, using 0 instead.")
      stats::pnorm((x-moment[1])/sqrt(moment[2]),lower.tail = FALSE)
    } else {
      a = 2/moment[3]
      a2 = a*a
      return(stats::pgamma((x-moment[1])/sqrt(moment[2])*a+a2, shape = a2, scale = 1, lower.tail = lower.tail))
    }
  }
}


#' fast p-value approximation
#'
#' Computes the p-value of a sample using Pearson's approximation of Gaussian quadratic forms with the estimators developed by Berschneider and Böttcher in [4].
#'
#' @param x matrix, the rows should be iid samples
#' @param vec vector, which indicates which columns of \code{x} are treated together as one sample. The default case treats each column as a separate sample.
#' @param type one of \code{"multi","total","m.multi.2","m.multi.3","all"}
#' @param ... these are passed to \code{\link{cdms}}
#'
#' @details This is the method recommended in [4], i.e., using Pearson's quadratic form estimate with the unbiased finite sample estimators for the mean and variance of normalized multivariance together with the unbiased estimator for the limit skewness.
#'
#' @references
#' For the theoretic background see the reference [4] given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
pearson.pvalue = function(x,vec = NA, type = "multi",...) {
  # undocumented: x can also be a list computed via multivariance:::cdms.mu.bcd()

  if (is.list(x)) { # assuming it is of type cdms.mu.bcd
    cmb = x

    if (any(is.na(vec))) {
      vec = 1:length(x$list.cdm)
    } else {
      cmb$list.cdm = cmb$list.cdm[vec]
      cmb$mu = cmb$mu[,vec]
      cmb$bcd = cmb$bcd[,vec]
      cmb$mean = cmb$mean[vec]
    }
    N = nrow(cmb$list.cdm[[1]])
  } else {
    if (any(is.na(vec))) vec = 1:ncol(x)

    dots <- list(...) # we filter the argument lambda which might be used for total.multivariance, we keep everything else to preserve other error messages.
    cmb = do.call('cdms.mu.bcd', c(list(x = x, vec = vec), dots[!(names(dots) %in% "lambda")]))
    # cmb = cdms.mu.bcd(x,vec, psi = psi, p = p, isotropic = isotropic,...)
    N = nrow(x)
  }

  # calculate the moments
  if (type == "all") {
    moms = sapply(c("multi","total","m.multi.2","m.multi.3"),function(ty) moments.for.pearson(N,cmb$bcd, cmb$mu, cmb$mean, type = ty))
  } else {
    moms = moments.for.pearson(N,cmb$bcd, cmb$mu, cmb$mean, type = type)
  }

  if ((!any(is.na(cmb$mu[3,]))) & (any(cmb$mu[3,]<0,na.rm = TRUE))) {
    warning(paste("Pearson's approximation: For the variable(s)",paste(which(cmb$mu[3,]<0),collapse = " and "),"the estimated limit skewness was negative. This usually indicates a (too) small sample size. It is recommended to use a resampling test (e.g. 'p.value.type='resample'') instead. \n"))
  }

  # normalization already in cdms.mu.bcd
  # normalizing.factor = rep(cmb$mean, each = nrow(x)^2)
  # normalizing.factor[normalizing.factor == 0] = 1 # prevent division by 0. if 0 then cdm == 0 anyway

  ## list.cdm = alply(cmb$array.cdm/normalizing.factor,3)

  if (type == "all") {
    m = multivariances.all(cmb$list.cdm)

    return(sapply(c("multi","total","m.multi.2","m.multi.3"),function(ty) pearson.qf(unname(m[ty]),moms[,ty],lower.tail = FALSE)))
  } else {

    switch( type,
            multi = {m = multivariance(cmb$list.cdm)},
            total = {m = total.multivariance(cmb$list.cdm)},
            m.multi.2 = {m = m.multivariance(cmb$list.cdm)},
            m.multi.3 = {m = m.multivariance(cmb$list.cdm,m = 3)},
            {stop(paste("unkown type:",type))}
    )

    return(pearson.qf(m,moms,lower.tail = FALSE))
  }
}

##### Copula Multivariance ####

#' Transform a vector of samples into a vector of samples of the uniform distribution
#' such that, if applied to multiple (dependent) sample vectors, the dependence is preserved.
#'
#' @examples
#' # x = rnorm(100)
#' # plot(emp.transf.vec(x))
#'
#' @keywords internal
emp.transf.vec = function(x,unif.samples = stats::runif(length(x))) {
  N = length(x)
  emp.p.x = rank(x,ties.method = "max")/N # value of the empirical distribution function at x
  emp.d.x = stats::ave(x, x, FUN = length)/N # empirical probability of x

  emp.t = emp.p.x - unif.samples*emp.d.x # here additional uniform samples are introduced to "smooth" jumps in the original distribution function - this is the key concept of the monte carlo distributional transform.

  return(emp.t)
}

#' Monte Carlo empirical transform
#'
#' Transforms a matrix (rows: samples, columns: variables) into a matrix of uniform samples with the same dependence structure via the Monte Carlo empirical transform.
#'
#' @param x data matrix (rows: samples, columns: variables)
# ru = runif(nrow(x)) # this would introduce dependence
# apply(x,2,function(y) emp.transf.vec(y,ru))
#'
#' @references
#' For the theoretic background see the reference [5] given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
emp.transf = function(x) {
  apply(x,2,function(y) emp.transf.vec(y))
}

#' A dependent Monte Carlo emprical transform
#'
#' @keywords internal
emp.transf.dep = function(x) {
  u = stats::runif(nrow(x))
  apply(x,2,function(y) emp.transf.vec(y,u))
}


#' copula version of distance multivariance
#'
#' Formally it is nothing but distance multivariance applied to the Monte Carlo emprical transform of the data. Hence its values vary for repeated runs.
#'
#' @inheritParams multicorrelation
#'
#' @references
#' For the theoretic background see the reference [5] given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
copula.multivariance = function(x,vec = 1:ncol(x), type = "total",...) {
  if (!is.matrix(x)) stop("The copula multivariance requires x to be a data matrix.")

  switch(type,
    multi = {fun = function (x,v) multivariance(x,vec = v,...)},
    total = {fun = function (x,v) total.multivariance(x,vec = v,...)},
    m.multi.2 = {fun = function (x,v) m.multivariance(x,vec = v,...)},
    m.multi.3 = {fun = function (x,v) m.multivariance(x,vec = v,m = 3,...)},
    all = {fun = function (x,v) multivariances.all(x,vec = v,...)},
    {stop(paste("unkown type:",type))}
  )

  return(fun(emp.transf(x),vec))
}

#' coupla versions of distance multicorrelation
#'
#' Formally it is nothing but distance multicorrelation applied to the Monte Carlo emprical transform of the data. Hence its values vary for repeated runs.
#'
#' @inheritParams multicorrelation
#' @param ... are passed to \code{\link{multicorrelation}}
#'
#' @references
#' For the theoretic background see the reference [5] given on the main help page of this package: \link{multivariance-package}.
#'
#' @seealso
#' \code{\link{multicorrelation}}
#'
#' @aliases CMcor
#' @export
copula.multicorrelation = function(x,vec = 1:ncol(x),...)  {
  if (!is.matrix(x)) stop("The copula multivariance requires x to be a data matrix.")
  multicorrelation(emp.transf(x),vec = vec,...)
}


#' @rdname copula.multicorrelation
#' @export
CMcor <- copula.multicorrelation


#' independence tests using the copula versions of distance multivariance
#'
#' Formally it is nothing but tests for distance multivariance applied to the Monte Carlo emprical transform of the data. Hence its values vary for repeated runs.
#'
#' @inheritParams multivariance.test
#'
#' @references
#' For the theoretic background see the reference [5] given on the main help page of this package: \link{multivariance-package}.
#'
#' @export
copula.multicorrelation.test = function(x,vec = 1:ncol(x),...) {
  if (!is.matrix(x)) stop("The copula multivariance requires x to be a data matrix.")

  if (all(table(vec) == 1)) {
    multivariance.test(emp.transf(x),vec = vec,p.value.type = "pearson_unif",...)
  } else {
    multivariance.test(emp.transf(x),vec = vec,p.value.type = "pearson_approx",...)
  }
}


#' compute the p-value by Pearson's approximation assuming uniform marginals and euclidean distance
#' @examples
#' multivariance:::pearson.pvalue.unif(matrix(runif(300),ncol = 3))
#'
#' \dontrun{
#' library(microbenchmark)
#' x = matrix(runif(300*3),ncol = 3)
#' microbenchmark(
#'   multivariance.test(x,p.value.type = "pearson_approx")$p.value,
#'   multivariance:::pearson.pvalue.unif(emp.transf(x))
#'   )
#' }
#'
#' @keywords internal
pearson.pvalue.unif = function(x,vec = NA,type = "total",psi = NULL, isotropic = TRUE, ...) {

  if (any(is.na(vec))) {
    vec = 1:ncol(x)
  }

  if (!(all(table(vec) == 1)))
    stop("'pearson_unif' makes only sense for one dimensional marginals, use 'pearson_approx' instead.")

  if ((!is.null(psi)) | (!isotropic) | (max(vec)<ncol(x))) {
    stop("The supplied arguments are not supported when using 'pearson_unif', use 'pearson_approx' instead.")
    #return(pearson.pvalue(x,vec = vec, type = type, psi = psi, isotropic = isotropic, ...))
  } else {
    n = max(vec)
    list.cdm = cdms(x,vec = vec,normalize = TRUE)

    unif.cmb = list(list.cdm= list.cdm,
      mu = matrix(rep(c(1/3,2/45,8/945),n),ncol=n),
      bcd = matrix(rep(c(1/6,7/60,1/9),n),ncol=n),
      mean =rep(1/3,n))

    return(pearson.pvalue(unif.cmb,type = type))
  }
}


##### Dependence structure ####



# * Data ####

# ** dep_struct_several ####
#' example dataset for \code{\link{dependence.structure}}
#'
#' It was generated by \preformatted{
#' set.seed(1348879148)
#' N = 100
#' dep_struct_several_26_100 = cbind(coins(N,2),tetrahedron(N),coins(N,4),
#'     tetrahedron(N),tetrahedron(N),coins(N,3),coins(N,3),rnorm(N))
#'save(dep_struct_several_26_100,file ="dep_struct_several_26_100.rda")
#'}
#'
#' To avoid irritation, note that the seed is just a simple integer hash value of the variable name.
#'
#' @format \code{matrix} 26 variables (columns), 100 independent samples (rows)
#'
"dep_struct_several_26_100"

# ** dep_struct_star ####
#' example dataset for \code{\link{dependence.structure}}
#'
#' It was generated by \preformatted{
#' set.seed(222454572)
#' N = 100
#' y = coins(N,2)
#' dep_struct_star_9_100 = cbind(y,y,y)
#' save(dep_struct_star_9_100,file ="dep_struct_star_9_100.rda")
#'}
#'
#' To avoid irritation, note that the seed is just a simple integer hash value of the variable name.
#'
#' @format \code{matrix} 9 variables (columns), 100 independent samples (rows)
#'
"dep_struct_star_9_100"

# ** dep_struct_iterated ####
#' example dataset for \code{\link{dependence.structure}}
#'
#' It was generated by \preformatted{
#' set.seed(532333356)
#' N = 100
#' x = matrix(sample.int(2,10*N,replace = TRUE)-1,ncol = 10)
#' for (i in c(2,5,9)) x = cbind(x,(rowSums(as.matrix(x[,1:(i-1)])) %% 2) == x[,i])
#' dep_struct_iterated_13_100 = x
#' save(dep_struct_iterated_13_100,file ="dep_struct_iterated_13_100.rda")
#'}
#'
#' To avoid irritation, note that the seed is just a simple integer hash value of the variable name.
#'
#' @format \code{matrix} 13 variables (columns), 100 independent samples (rows)
#'
"dep_struct_iterated_13_100"

# ** dep_struct_ring ####
#' example dataset for \code{\link{dependence.structure}}
#'
#' It was generated by \preformatted{
#' set.seed(436646700)
#' N = 100
#' n= 15
#' x=matrix(sample.int(2,N*n,replace = TRUE)-1,nrow =N)
#' x[,4] = rowSums(x[,1:3]) %% 2
#' x[,7] = rowSums(x[,4:6]) %% 2
#' x[,10] = rowSums(x[,7:9]) %% 2
#' x[,13] = rowSums(x[,10:12]) %% 2
#' x[,15] = rowSums(x[,c(13,14,1)]) %% 2
#' dep_struct_ring_15_100 = x
#' save(dep_struct_ring_15_100,file ="dep_struct_ring_15_100.rda")
#'}
#'
#' To avoid irritation, note that the seed is just a simple integer hash value of the variable name.
#'
#' @format \code{matrix} 15 variables (columns), 100 independent samples (rows)
#'
"dep_struct_ring_15_100"

# ** Anscombe ####
#' Extended Anscombe's Quartett
#'
#' The dataset extends 'anscombe' provided in the
#' standard R-package 'datasets'. All examples feature the same
#' correlation of 0.82, but different types of dependencies. The main aim was to extend the classical examples, which have
#' sample size 11, to larger sample sizes. This illustrates that the
#' implied problems of Pearson's correlation are not small sample
#' problems! Distance multicorrelation (which coincides in
#' this case with distance correlation) yields different values
#' for the datasets.
#'
#' Note: Anscombe's quartett features further identical parameters
#' besides Pearson's correlation. The extended set is only
#' concerned with correlation.
#'
#' @format \code{list} with elements:
#' \itemize{
#' \item\code{anscombe.extended$N11}
#'  matrix with 11 samples for 5 examples the first 4 are the
#'         classical Anscombe Quartett, the fifth is a monoton relation
#'         which also features the same correlation.
#' \item \code{anscombe.extended$N100} same as above but 100 samples
#' \item \code{anscombe.extended$N1000} same as above but 1000 samples
#' }
#'
#' @references
#' This example was introduced in the reference [6] given on the main help page of this package: \link{multivariance-package}.
#'
#' @examples
#'
#'# Code which generates plots of all included data:
#' op = par(mfrow = c(3,5),mar = c(0.5,0.5,3,0.5))
#'for (name in c("N11","N100","N1000")) {
#'  for (i in 1:5) {
#'    x = anscombe.extended[[name]][,2*i-1]
#'    y = anscombe.extended[[name]][,2*i]
#'    plot(x,y,main = paste0("cor = ",round(cor(x,y),2),
#' "\n Mcor = ",round(multicorrelation(cbind(x,y),type = "pairwise",squared = FALSE),2),
#' "\n CMcor = ",round(copula.multicorrelation(cbind(x,y),type = "pairwise",squared = FALSE),2)),
#'         axes = FALSE,xlab ="",ylab = "", cex.main=1)
#'    # for two variables 'pairwise' coincides with
#'    # both values of 'total.upper.lower'.
#'    box()
#'  }
#' }
#' par(op)
#'
"anscombe.extended"

# * detection ####

#' determines the dependence structure
#'
#' Determines the dependence structure as described in [3].
#'
#' @param x matrix, each row of the matrix is treated as one sample
#' @param vec vector, it indicates which columns are initially treated together as one sample
#' @param type the method used for the detection, one of '\code{conservative}','\code{resample}','\code{pearson_approx}' or '\code{consistent}'
#' @param structure.type either the '\code{clustered}' or the '\code{full}' structure is detected
#' @param verbose boolean, if \code{TRUE} details are printed during the detection and whenever a cluster is newly detected the (so far) detected dependence structure is plotted.
#' @param detection.aim \code{=NULL} or a list of vectors which indicate the expected detection, see below for more details
#' @param alpha numeric between 0 and 1, the significance level used for the tests
#' @param p.adjust.method a string indicating the p-value adjustment for multiple testing, see \code{\link{p.adjust.methods}}
#' @param c.factor numeric, larger than 0, a constant factor used in the case of '\code{type = "consistent"}'
#' @param stop.too.many numeric, upper limit for the number of tested tuples. A warning is issued if it is used. Use \code{stop.too.many = NULL} for no limit.
#' @param list.cdm not required, the list of doubly centered distance matrices corresponding to \code{x} speeds up the computation if given
#' @param ... these are passed to \code{\link{find.cluster}}
#'
#' @details
#' Performs the detection of the dependence structure as described in [3]. In the \code{clustered} structure variables are clustered and treated as one variable as soon as a dependence is detected, the \code{full} structure treats always each variable separately. The detection is either based on tests with significance level \code{alpha} or a \code{consistent} estimator is used. The latter yields (in the limit for increasing sample size) under very mild conditions always the correct dependence structure (but the convergence might be very slow).
#'
#' If \code{fixed.rejection.level} is not provided, the significance level \code{alpha} is used to determine which multivariances are significant using the distribution-free rejection level. As default the Holm method is used for p-value correction corresponding to multiple testing.
#'
#' The resulting graph can be simplified (pairwise dependence can be represented by edges instead of vertices) using \code{\link{clean.graph}}.
#'
#' Advanced:
#' The argument \code{detection.aim} is currently only implemented for \code{structure.type = clustered}. It can be used to check, if an expected dependence structure was detected. This might be useful for simulation studies to determine the empirical power of the detection algorithm. Hereto  \code{detection.aim} is set to a list of vectors which indicate the expected detected dependence structures (one for each run of \code{\link{find.cluster}}). The vector has as first element the \code{k} for which k-tuples are detected (for this aim the detection stops without success if no k-tuple is found), and the other elements, indicate to which clusters all present vertices belong after the detection, e.g. \code{c(3,2,2,1,2,1,1,2,1)} expects that 3-tuples are detected and in the graph are 8 vertices (including those representing the detected 3 dependencies), the order of the 2's and 1's indicate which vertices belong to which cluster. If \code{detection.aim} is provided, the vector representing the actual detection is printed, thus one can use the output with copy-paste to fix successively the expected detection aims.
#'
#' Note that a failed detection might invoke the warning:
#' \preformatted{
#' run$mem == detection.aim[[k]][-1] :
#' longer object length is not a multiple of shorter object length
#' }
#'
#' @return returns a list with elements:
#' \describe{
#'   \item{\code{multivariances}}{calculated multivariances,}
#'   \item{\code{cdms}}{calculated doubly centered distance matrices,}
#'   \item{\code{graph}}{graph representing the dependence structure,}
#'   \item{\code{detected}}{boolean, this is only included if a \code{detection.aim} is given,}
#'   \item{\code{number.of.dep.tuples}}{vector, with the number of dependent tuples for each tested order. For the full dependence structure a value of -1 indicates that all tuples of this order are already lower order dependent, a value of -2 indicates that there were more than \code{stop.too.many} tuples,}
#'   \item{\code{structure.type}}{either \code{clustered} or \code{full},}
#'   \item{\code{type}}{the type of p-value estimation or consistent estimation used,}
#'   \item{\code{total.number.of.tests}}{numeric vector, with the number of tests for each group of tests,}
#'   \item{\code{typeI.error.prob}}{estimated probability of a type I error,}
#'   \item{\code{alpha}}{significance level used if a p-value estimation procedure is used,}
#'   \item{\code{c.factor}}{factor used if a consistent estimation procedure is used,}
#'   \item{\code{parameter.range}}{significance levels (or 'c.factor' values) which yield the same detection result.}
#   \item{\code{}}{}
#' }
#'
#' @references
#' For the theoretic background see the reference [3] given on the main help page of this package: \link{multivariance-package}.
#'
#' @example inst/examples/dependence-structures.R
#' @export
#'
dependence.structure = function(x, vec = 1:ncol(x), verbose = TRUE, detection.aim = NULL, type = "conservative", structure.type = "clustered", c.factor = 2, list.cdm= NULL, alpha = 0.05, p.adjust.method = "holm",stop.too.many = NULL,...) {

  dots <- list(...) # just to preserve backwards compatibility for "fixed.rejection.level".

  fixed.rejection.level = dots[["fixed.rejection.level"]]
  if (is.null(fixed.rejection.level)) {
    #    if (!any(names(dots)=="fixed.rejection.level"))
    fixed.rejection.level = NA
  } else {
    type = "fixed.rejection.level"
    if (length(fixed.rejection.level) == 1)
      fixed.rejection.level = rep(fixed.rejection.level,max(vec))
  }


  if (verbose){ cat(paste0("\n    Dependence structure detection\n\ndata: '",deparse(substitute(x)),"' with ",nrow(x)," samples of ",max(vec)," variables\nstructure: ",structure.type,"\ndetection method: "))
    switch(type,
           consistent = {
             cat(paste0("consistent estimate -- using the factor: ",c.factor,"\n\n"))
           },
           conservative =, resample =, pearson_approx = {
             cat(paste0("test -- method for p-values: ", type," -- significance level: ",alpha,"\n\n"))
           },
           fixed.rejection.level = {cat(paste0("using the given fixed rejection level\n\n"))},
           {stop(paste0("unknown type: ",type))}
    )
  }

  switch(structure.type,
         full = return(dependence.structure.full(x,vec = vec,verbose = verbose,type = type, alpha = alpha,list.cdm = list.cdm,stop.too.many = stop.too.many,c.factor = c.factor,...)),
         clustered = {},
         {stop(paste0("unknown structure.type: ",structure.type))}
  )

  if (type == "consistent") {
    fixed.rejection.level = rep(sqrt(nrow(x))*c.factor,max(vec))
  }

  if (is.null(list.cdm)) {
    switch (type,
            consistent =, conservative =, resample =, fixed.rejection.level= {
              list.cdm = cdms(x,vec = vec) # creates the distance matrices
            },
            pearson_approx = {
              list.cdm = cdms.mu.bcd(x,vec = vec)
            },
            {stop(paste0("unknown type: ",type))}
    )

  }

  all.multivariances = numeric(0) # vector which will contain all distance multivariances which are calculated

  mem = as.numeric(1:max(vec))
  #its length is the number of vertices, its content is the number of the corresponding cluster for the current iteration!!!
  # has to be numeric, since otherwise 'identical' fails to end the loop (in the case of
  # no detected clusters in the first run)

  cluster.to.vertex = 1:max(mem) # cluster to vertex relation - gets renewed each iteration (since the names of the clusters change)

  vertex.to.cdm = 1:max(mem) # vertex to A (the doubly centered distance matrices) relation - gets appended each iteration


  previous.n.o.cdms = rep(0,max(mem)) # number of As in the previous iteration

  sig.limits = NULL

  n = max(mem) # number of clusters

  g = igraph::graph.empty(,directed=FALSE)
  g = igraph::add.vertices(g,n,label = sapply(1:max(mem),function(r) paste(colnames(x,do.NULL = FALSE,prefix = "")[vec == r],collapse = ",")),shape = "circle")

  # Loop through the tuples
  detected = TRUE
  k = 1
  number.of.tests = NULL

  while (detected) {
    if (!is.null(detection.aim)) {
      #    run = find.cluster(x,vec,list.cdm,mem,cluster.to.vertex,vertex.to.cdm,previous.n.o.cdms,all.multivariances,g,kvec = 2:detection.aim[[k]][1], verbose = verbose, sig.limits = limits, type = type, ...)
      if (length(detection.aim)<k) {
        warning("Not enough detection aims to complete detection. This is expected if this run was used to output a detection aim.\n")
        return()
      }

      run = do.call('find.cluster', c(list(
        x = x, vec = vec,list.cdm = list.cdm,
        mem = mem, cluster.to.vertex = cluster.to.vertex,
        vertex.to.cdm = vertex.to.cdm,
        previous.n.o.cdms = previous.n.o.cdms,
        all.multivariances = all.multivariances,
        g = g, kvec = 2:detection.aim[[k]][1],
        verbose = verbose, parameter.range = sig.limits,
        type = type, fixed.rejection.level = fixed.rejection.level),
        alpha = alpha, p.adjust.method = p.adjust.method,
        stop.too.many = stop.too.many,
        dots[!(names(dots) == "fixed.rejection.level")]))

      if (verbose) {
        cat("last detected structure (in detection.aim format): ")
        dput(c(run$k,run$mem))
      }
      success = all(run$mem == detection.aim[[k]][-1])
    } else {
      #     run = find.cluster(x,vec,list.cdm,mem,cluster.to.vertex,vertex.to.cdm,previous.n.o.cdms,all.multivariances,g,verbose = verbose, sig.limits =sig.limits, type = type,...)
      run = do.call('find.cluster', c(list(
        x = x, vec = vec,list.cdm = list.cdm,
        mem = mem, cluster.to.vertex = cluster.to.vertex,
        vertex.to.cdm = vertex.to.cdm,
        previous.n.o.cdms = previous.n.o.cdms,
        all.multivariances = all.multivariances,
        g = g,
        verbose = verbose, parameter.range = sig.limits,
        type = type, fixed.rejection.level = fixed.rejection.level),
        alpha = alpha, p.adjust.method = p.adjust.method,
        stop.too.many = stop.too.many,
        dots[!(names(dots) == "fixed.rejection.level")]))
    }

    detected = run$detected
    list.cdm = run$list.cdm
    mem = run$mem
    cluster.to.vertex = run$cluster.to.vertex
    vertex.to.cdm = run$vertex.to.cdm
    previous.n.o.cdms = run$previous.n.o.cdms
    all.multivariances = run$all.multivariances
    g = run$g
    sig.limits = run$parameter.range
    number.of.tests = c(number.of.tests, run$number.of.tests)

    k = k+1
    if (!is.null(detection.aim)) if (!success) break
    if (run$stopped) break

  }

  typeI.error.prob = NA

  #number.of.tests = length(run$all.multivariances)
  #groups.of.tests = nrow(sig.limits)

  groups.of.tests = sum(number.of.tests > 0)

  if (anyNA(fixed.rejection.level))
    typeI.error.prob = 1-(1-alpha)^groups.of.tests

  if (type %in% c("consistent","fixed.rejection.level")) {
    typeI.error.prob = 1- (1- multivariance.pvalue(fixed.rejection.level[1]))^sum(number.of.tests)

    sig.m = run$all.multivariances> sqrt(nrow(x))*c.factor

    mult.factor = 1
    if (type == "consistent") mult.factor = 1/sqrt(nrow(x))

    low.c = ifelse(any(!sig.m),max(run$all.multivariances[!sig.m]*mult.factor),0)
    high.c = ifelse(any(sig.m),min(run$all.multivariances[sig.m]*mult.factor),Inf)
    parameter.range = c(low.c,high.c)
  } else {
    parameter.range = cbind(sig.limits[,1],sig.limits[,2])
  }


  if (verbose) {

    cat(paste0("\ntotal number of tests: ",sum(number.of.tests),", groups of tests: ", groups.of.tests,"\n"))
   # print(number.of.tests)

    #print(sig.limits)

    switch(type,
           conservative =, resample =, pearson_approx = {
             cat(paste0("\nSame structure for any significance level in (",signif(max(sig.limits[,1]),3),",",signif(min(sig.limits[,2]),3),")\n"))

             # typeI.error.prob = 1-(1-alpha)^nrow(sig.limits)
             #cat(paste0("\nThe total number of tests was ", length(run$all.multivariances), ", these were grouped into ",nrow(sig.limits)," tests of significance level ",alpha," (within each group the p-values were adjusted by ",p.adjust.method,"'s method).", "\nThis yields (assuming independence) an approximate ", if (type == "conservative") "conservative bound for the ", "global type I error probability of ",signif(typeI.error.prob,3),".\n"))
           },
           consistent = {
             #typeI.error.prob = 1- (1- multivariance.pvalue(fixed.rejection.level[1]))^length(run$all.multivariances)
             cat(paste0("\nSame structure for any 'c.factor' in (",signif(low.c,3),",",signif(high.c,3),")\n"))
             #cat(paste0("\nApproximate upper bound for prob. of type I error: ", signif(typeI.error.prob,3),"\n"))
           },
           fixed.rejection.level = {
             #      cat(paste0("\nSame structure for rejection level in (",signif(low.c,3),",",signif(high.c,3),")\n")) #problematic, since it can be a vector with different values
           }
    )

    cat(paste0("\n",ifelse(type == "conservative","Conservative estimate","Estimate")," of the type I error probability: ",signif(typeI.error.prob,3),"\n"))

  }

  number.of.dep.tuples = c(NA,as.vector(table(factor(igraph::V(g)$level,2:max(vec)))))


  if (!is.null(detection.aim)) {
    return(invisible(list(cdms = run$list.cdm, multivariances = run$all.multivariances, graph = run$g,detected = success)))
  } else {
    return(invisible(list(cdms = run$list.cdm, multivariances = run$all.multivariances, graph = run$g,number.of.dep.tuples = number.of.dep.tuples,parameter.range = parameter.range,typeI.error.prob = typeI.error.prob,total.number.of.tests = number.of.tests,structure.type = structure.type, type = type,alpha = alpha, c.factor = c.factor)))
  }
}

#' cluster detection
#'
#' Performs the detection of dependence structures algorithm until a cluster is found. This function is the basic building block \code{\link{dependence.structure}}. Advanced users, might use it directly.
#'
#' @param x matrix with the samples
#'
#' @param vec vector, it indicates which columns are initially treated together as one sample
#' @param list.cdm list of doubly centered distance matrices
#' @param mem numeric vector, its length is the number of vertices, its content is the number of the corresponding cluster for the current iteration, i.e., vertex \code{i} belongs to cluster \code{mem[i]}
#' @param cluster.to.vertex vector, contains the cluster to vertex relations, i.e., \code{cluster.to.vertex[i]} is the index of the vertex which represents cluster \code{i}
#' @param vertex.to.cdm vector, contains the vertex to doubly centered distance matrix relations, i.e., \code{vertex.to.cdm[i]} is the index of the doubly centered distance matrix in \code{list.cdm} which corresponds to vertex \code{i}
#' @param previous.n.o.cdms vector, number of the doubly centered distance matrices in the previous iteration (it is used to ensure that previously check tuples are not checked again)
#' @param all.multivariances vector, which contains all distance multivariances which have been calculated so far. Only used to finally return all distance multivariances which have been calculated.
#' @param g dependence structure graph
#' @param alpha numeric, significance level used for the (distribution-free) tests
#' @param fixed.rejection.level vector, if not \code{NA} the \code{fixed.rejection.level[k]} is used for the k-tuples, instead of a level derived from the significance level \code{alpha}
#' @param p.adjust.method name of the method used to adjust the p-values for multiple testing, see \code{\link[stats]{p.adjust}} for all possible options.
#' @param verbose boolean, if \code{TRUE} details during the detection are printed and whenever a cluster is newly detected the (so far) detected dependence structure is plotted.
#' @param kvec vector, k-tuples are only checked for each k in \code{kvec}, i.e., for \code{kvec = 2:4} only 2,3 and 4-tuples would be check and then the algorithm stops.
#' @param parameter.range numeric matrix, which hosts the range of significance levels or '\code{c.factor}' which yield the same detected structure
#' @param type the method for the detection, one of '\code{conservative}','\code{resample}','\code{pearson_approx}' or '\code{consistent}'.
#' @param stop.too.many numeric, upper limit for the number of tested tuples. A warning is issued if it is used. Use \code{stop.too.many = NULL} for no limit.
#' @param ... are passed to \code{\link{resample.multivariance}} in the case of '\code{type = resample}'
#'
#' @details
#' For further details see \code{\link{dependence.structure}}.
#'
find.cluster = function(x,
  vec = 1:ncol(x), # which columns should be treated as one sample
  list.cdm = cdms(x,vec = vec), # creates the distance matrices
  mem = as.numeric(1:max(vec)),
  #its length is the number of vertices, its content is the number of the corresponding cluster for the current iteration!!!
  # has to be numeric, since otherwise 'identical' fails to end the loop (in the case of
  # no detected clusters in the first run)
  cluster.to.vertex = 1:max(mem), # cluster to vertex relation - gets renewed each iteration (since the names of the clusters change)
  vertex.to.cdm = 1:max(mem), # vertex to A (the doubly centered distance matrices) relation - gets appended each iteration
  previous.n.o.cdms = rep(0,max(mem)), # number of As in the iteration before. it is used to speed up the detection.
  all.multivariances = numeric(0), # vector which will contain all distance multivariances which are calculated
  g = igraph::add.vertices(igraph::graph.empty(,directed=FALSE),max(mem),label = sapply(1:max(mem),function(r) paste(colnames(x,do.NULL = FALSE,prefix = "")[vec == r],collapse = ",")),shape = "circle"), #the graph
  fixed.rejection.level = NA,
  alpha=0.05,p.adjust.method = "holm",
  verbose = TRUE,
  kvec = 2:max(mem),
  parameter.range = NULL,
  type = "conservative",
  stop.too.many = NULL,
  ...) {
  explore = FALSE # undocumented option, which would provide some more graphs during the detection
  if (verbose) graphics::plot(g)

  n = max(mem) # number of clusters
  nV = length(igraph::V(g)) #length(mem) # number of vertices at the start of the iteration

  if (type == "pearson_approx") {
    n.o.cdm = length(list.cdm$list.cdm) # number of As at the start of the iteration
  } else {
    n.o.cdm = length(list.cdm) # number of As at the start of the iteration
  }

  cluster.to.cdm = vertex.to.cdm[cluster.to.vertex] #cluster to A relation - gets renewed each iteration
  # Each cluster is represented by its 'largest' vertex

  cdm.to.vertex = NA # A to vertex relation
  for (i in 1:n.o.cdm) cdm.to.vertex[i] = which(vertex.to.cdm == i)

  number.of.tests = 0
  stopped = FALSE

  for (k in 2:min(max(kvec),max(mem))) { # look at the k-tuples of the n variables.

    tuples = utils::combn(n,k) # all k-tuples of 1,..,n

    tuples = matrix(cluster.to.cdm[tuples],ncol = k,byrow = TRUE)
    #transform tuples into the A indexes

    tuples = tuples[apply(tuples,1,max) > previous.n.o.cdms[k],]
    # to speed up, we only consider those with new A

    if (length(tuples) == k) dim(tuples) = c(1,k)
    # make sure that it is a matrix

    if (all(c(!is.null(stop.too.many),(nrow(tuples) > stop.too.many)))) {
      warning(paste0("More than ",stop.too.many," tuples. Detection stopped. The above results represent only the performed tests.\n"))
      stopped = TRUE
      break
    }

    number.of.tests[k] = nrow(tuples)

    if (verbose)
      cat(paste0(k,"-tuples x ",nrow(tuples),":"))

    switch(type,
      consistent =,
      conservative =,
      fixed.rejection.level = {
        multivariances = apply(tuples,1,function(x) multivariance(list.cdm,x,Nscale = TRUE)) #calculates all distance multivariances
        multivariances.pvalues = multivariance.pvalue(multivariances)},

      resample = {
        res = apply(tuples,1,function(x) resample.multivariance(list.cdm,x,...))
        multivariances = sapply(res,function(x) x$original)
        multivariances.pvalues = sapply(res,function(x) x$p.value)
      },

      pearson_approx = {
        multivariances = apply(tuples,1,function(x) multivariance(list.cdm$list.cdm,x,Nscale = TRUE))
        multivariances.pvalues = apply(tuples,1,function(x) pearson.pvalue(list.cdm,x))
      },

      {#default
        stop(paste("unkown type:",type))
      }
    )



    all.multivariances = c(all.multivariances,multivariances)
    # print(multivariances)
    if (explore) {
      graphics::plot(multivariances,main = paste(k,"tuple multivariances"))
      graphics::abline(h = c( rejection.level(alpha/choose(n,k)), rejection.level(alpha)), col = c("red","green"))
      readline("continue with [enter]")
    }


    adjusted.pvs = stats::p.adjust(multivariances.pvalues,method = p.adjust.method)

    critical.tuples =  which(((anyNA(fixed.rejection.level) & (adjusted.pvs < alpha))) | (multivariances > fixed.rejection.level[k]) )
    for (i in critical.tuples) {
      # for each tuple, with adjusted p value less than the significance level (or multivariance less than a prescribed fixed.rejection.level, if given) we add a vertex and the edges

      new = length(igraph::V(g))+1
      # next line adds the label to the dependency vertex
      # g = igraph::add.vertices(g,1,label = pvalue.to.starlabel(adjusted.pvs[i]), shape = "none", level=k)
      g = igraph::add.vertices(g,1,label = signif(multivariances[i],4), shape = "none", level=k)

      g = igraph::add_edges(g, as.vector(t(cbind(new,cdm.to.vertex[tuples[i,]]))), weight= NA, color = k, lty = k) # weight = adjusted.pvs[i]
      #  }
    }

    #if (type != "consistent")
    {
      sorted.adjusted.pvs = sort(adjusted.pvs)

      # print(sorted.adjusted.pvs)

      first.non.critical = match(FALSE, sorted.adjusted.pvs < alpha)
      if (is.na(first.non.critical)) {
        lower.sig = sorted.adjusted.pvs[length(sorted.adjusted.pvs)]
        upper.sig = 1
      } else {
        if (first.non.critical==1) {
          lower.sig = 0
          upper.sig = sorted.adjusted.pvs[1]

        } else {
          lower.sig = sorted.adjusted.pvs[first.non.critical-1]
          upper.sig = sorted.adjusted.pvs[first.non.critical]
        }
      }

      parameter.range = rbind(parameter.range,c(lower.sig,upper.sig))
    }

    if (verbose) {
      cat(paste0(" max. multivariance: ",format(max(multivariances),digits=3,nsmall = 3,width = 7)))
      if (anyNA(fixed.rejection.level)) {
        cat(paste0("; min. p-value (unadjusted): ",min(multivariances.pvalues),"\n",sep =""))
      } else {
        cat("\n")
      }

      #    cat(paste0("Same result for any significance level in (",lower.sig,",",upper.sig,")\n"))
    }
    #readline(paste("level",k,", press [Enter] for next level"))

    previous.n.o.cdms[k] = n.o.cdm

    #if ((anyNA(fixed.rejection.level) && (stats::p.adjust(multivariance.pvalue(max(multivariances)),method = p.adjust.method,n = length(multivariances)) < alpha)) || (!anyNA(fixed.rejection.level) && (max(multivariances) > fixed.rejection.level[k]))) {
    if (length(critical.tuples)>0) {
      #if a cluster was found exit the loop
      break
    }
  } # end loop k


  previous.n.o.cdms[(k+1):length(previous.n.o.cdms)] = 0 # ensures that for all higher tuples all combinations are checked. (Otherwise those tuples which were checked previously are skipped)

  if (verbose) graphics::plot(g)

  #print(paste("mem: ",paste(mem,collapse = " "),"clust",paste(igraph::clusters(g)$membership,collapse = " ")))
  if (identical(mem, igraph::clusters(g)$membership) || (igraph::clusters(g)$no == 1) ) {
    detected = FALSE
    if (verbose) {
      if (igraph::clusters(g)$no == 1) {cat("All vertices are in one cluster.\n")} else { cat("No new cluster detected.\n")}
    }    # will end recursion if the membership did not change or there is only one cluster
  } else {
    if (verbose) cat("New cluster detected. Not all vertices are in one cluster.\n")
    detected = TRUE
  }
  mem = igraph::clusters(g)$membership


  cluster.to.vertex = NA
  for (i in 1:max(mem)) cluster.to.vertex[i] = max(which(mem == i))
  # cluster to vertex relation. Each cluster is represented by its 'largest' vertex

  previous.n.o.cdms[2] = n.o.cdm

  for (i in cluster.to.vertex[cluster.to.vertex > nV]) {
    # for every vertex which represents a new cluster the distance matrix representing the cluster is calculated.
    #  the vertex i gets the cdm with index vertex.to.cdm[i]
    n.o.cdm = n.o.cdm + 1
    vertex.to.cdm[i] = n.o.cdm
    # print(which(mem[1:ncol(x)] == which(cluster.to.vertex== i)))
    #      array.cdm = abind::abind(array.cdm, cdm(x[,which(mem[1:ncol(x)] == which(cluster.to.vertex== i))]),along=1)
    #array.cdm = abind::abind(array.cdm, cdm(x[,which(vec %in% which(mem[1:ncol(x)] == which(cluster.to.vertex== i)))]),along=3)

    if (type == "pearson_approx") {
      new.cmb = cdm.mu.bcd(x[,which(vec %in% which(mem[1:ncol(x)] == which(cluster.to.vertex== i)))])
      new.ind = length(list.cdm$list.cdm)+1

      if (new.cmb$mean == 0) { # note that cdm.mu.bcd returns unnormalized cdm
        list.cdm$list.cdm[[new.ind]] = new.cmb$cdm
      } else {
        list.cdm$list.cdm[[new.ind]] = new.cmb$cdm/new.cmb$mean
      }

      list.cdm$mu = cbind(list.cdm$mu,new.cmb$mu)
      list.cdm$bcd = cbind(list.cdm$bcd,new.cmb$bcd)
      list.cdm$mean[new.ind] = new.cmb$mean

    } else {
      list.cdm[[length(list.cdm)+1]] = cdm(x[,which(vec %in% which(mem[1:ncol(x)] == which(cluster.to.vertex== i)))])
    }


    # which(cluster.to.vertex== i) is the number of the cluster represented by vertex i
    # which(mem ...) gives the "original" vertices which belong to the same cluster as vertex i
    # which(vec ...) finally gives the columns of x which belong to the same cluster as vertex i

  } # at the end of this for loop, n.o.cdm contains the new number of As


  invisible(list(detected = detected,list.cdm = list.cdm,mem = mem,cluster.to.vertex = cluster.to.vertex,vertex.to.cdm = vertex.to.cdm,previous.n.o.cdms = previous.n.o.cdms,all.multivariances = all.multivariances,g = g, k = k, parameter.range = parameter.range,number.of.tests = number.of.tests[-1], stopped = stopped))
  # the first element of number.of.tests is 0 for 1-tuples
}





#' functions to detect the full (without clustering) dependence structure
#' @examples
#' # multivariance:::dependence.structure.full(dep_struct_ring_15_100)
#' # dependence.structure(dep_struct_ring_15_100,structure.type = "full")
#'
#' @keywords internal
dependence.structure.full = function(x,vec = 1:ncol(x),verbose = TRUE,type = "conservative", alpha = 0.05,list.cdm = NULL,maxk = max(vec),stop.too.many = NULL, c.factor = 2,...) {

  # compute the distance matrices
  if (is.null(list.cdm)) {
    switch( type,
      conservative =, resample = , consistent = {
        list.cdm = cdms(x,vec = vec)},
      pearson_approx = {
        list.cdm = cdms.mu.bcd(x,vec = vec)},
      {stop(paste0("unknown type: ",type))}
    )
  }

  n = max(vec)
  N = nrow(x)

  #R = rejection.level(0.005)
  R = alpha
  consistent.rejection.level = sqrt(N)*c.factor
  ## generate the graph

  g = igraph::graph.empty(,directed=FALSE)
  g = igraph::add.vertices(g,n,label = sapply(1:max(vec),function(r) paste(colnames(x,do.NULL = FALSE,prefix = "")[vec == r],collapse = ",")),shape = "circle")

  if (verbose) graphics::plot(g,layout = function(g) layout_on_circles(g,n))

  number.of.dep.tuples = rep(NA,n)

  number.of.tests = rep(NA,n)

  m.values = list()
  for (k in 2:min(maxk,n)) {
    if (all(c(!is.null(stop.too.many),(choose(n,k) > stop.too.many)))) {
      warning(paste0("More than ",stop.too.many," tuples. Detection stopped. Output contains only results of tests up to order ",k-1,".\n"))
      number.of.dep.tuples[k] = -2 # as an indication that there can be no such tuple
      k = k-1
      break
    }

    if (verbose) {
      cat(paste0(k,"-tuples x ",choose(n,k),""))
    }
    tuples = t(utils::combn(n,k))

    # check which tuples are lower order dependent
    lo = apply(tuples,1,function(tup) lower.order(tup,m.values))

    number.of.tests[k] = sum(!lo)

    if (verbose & (k>2)) cat(paste0(" of which ",sum(!lo)," are ",k-1,"-independent\n"))

    if (verbose & (k == 2)) cat("\n")

    if (all(lo)) {
      if (verbose) cat(paste0("all ",k,"-tuples contain lower order dependent tuples\n"))
      number.of.dep.tuples[k] = -1 # as an indication that there can be no such tuple
      k = k-1
      break
    }

    if (type == "consistent") {
      # calculate the multivariances
      if (sum(!lo) == 1) {
        ms = multivariance(list.cdm,tuples[!lo,],...)
      } else {
        ms = apply(tuples[!lo,],1,function(tup) multivariance(list.cdm,tup,...))
      }

      sig.tuples = rep(TRUE,nrow(tuples))
      sig.tuples[!lo] = ms > consistent.rejection.level
      sig.tuples.values = rep(NA,nrow(tuples))
      sig.tuples.values[!lo] = ms

    } else {
      # calculate the p values
      if (sum(!lo) == 1) {
        pvs = multivariance.test(list.cdm,tuples[!lo,],type = "multi",p.value.type = type,...)$p.value
      } else {
        pvs = apply(tuples[!lo,],1,function(tup) multivariance.test(list.cdm,tup,type = "multi",p.value.type = type,...)$p.value)
      }

      sig.tuples.values = rep(NA,nrow(tuples))
      sig.tuples = rep(TRUE,nrow(tuples))
      #sig.tuples.values[which(!lo)] = p.adjust(pvs)
      #sig.tuples[which(!lo)] = sig.tuples.values[which(!lo)] < R
      sig.tuples.values[!lo] = stats::p.adjust(pvs)
      sig.tuples[!lo] = sig.tuples.values[!lo] < R

    }

    number.of.dep.tuples[k] = sum(sig.tuples[which(!lo)])

    if (verbose) cat(paste0(number.of.dep.tuples[k]," dependencies found\n"))

    m.values[[k]] = cbind(tuples,
      # apply(tuples,1,function(tup) multivariance(list.cdm,tup,Nscale = TRUE)>R),
      sig.tuples,
      lo,
      sig.tuples.values
    )

    for (i in 1:nrow(m.values[[k]]))
    {
      if (m.values[[k]][i,k+1] & !m.values[[k]][i,k+2]) {
        new = length(igraph::V(g))+1
        g = igraph::add.vertices(g,1,label = signif(k,4), shape = "none", level=k)
        g = igraph::add_edges(g, as.vector(t(cbind(new,m.values[[k]][i,1:k]))), weight= NA, color = k, lty = k)
      }
    }

    if (verbose) graphics::plot(g,layout = function(g) layout_on_circles(g,n))

  }

  number.of.tests = number.of.tests[!is.na(number.of.tests)] # the first is NA (for 1-tuples) and the last might be NA if it was stopped before (either because there is no independent tuple, or too.many

  groups.of.tests = sum(number.of.tests >0) #,na.rm = TRUE

  if (verbose) {
    cat(paste0("\ntotal number of tests: ",sum(number.of.tests),", groups of tests: ", groups.of.tests,"\n"))

    #print(number.of.tests)
  }


  if (type == "consistent") {
    parameter.range = matrix(NA,nrow = length(m.values),ncol = 2)

    for (k in 2:length(m.values)) {
      #largest value of the independent tuples with lower order independence
      parameter.range[k,1] = max(m.values[[k]][!m.values[[k]][,"lo"]&!m.values[[k]][,"sig.tuples"],"sig.tuples.values"],0)
      #smallest value of the dependent tuples with lower order independence
      parameter.range[k,2] = min(m.values[[k]][!m.values[[k]][,"lo"]&m.values[[k]][,"sig.tuples"],"sig.tuples.values"],Inf)
    }

    parameter.range = parameter.range/sqrt(N)

    if (verbose) cat(paste0("\nSame structure for any 'c.factor' in (",signif(max(parameter.range[,1],na.rm = TRUE),3),",",signif(min(parameter.range[,2],na.rm = TRUE),3),").\n"))

  } else {
    parameter.range = matrix(NA,nrow = length(m.values),ncol = 2)

    for (k in 2:length(m.values)) {
      #largest (p-)value of the dependent tuples with lower order independence
      parameter.range[k,1] = max(m.values[[k]][!m.values[[k]][,"lo"]&m.values[[k]][,"sig.tuples"],"sig.tuples.values"],0)
      #smallest (p-)value of the dependent tuples with lower order independence
      parameter.range[k,2]= min(m.values[[k]][!m.values[[k]][,"lo"]&!m.values[[k]][,"sig.tuples"],"sig.tuples.values"],1)
    }

    if (verbose) cat(paste0("\nSame structure for any significance level in (",signif(max(parameter.range[,1],na.rm = TRUE),3),",",signif(min(parameter.range[,2],na.rm = TRUE),3),").\n"))
  }


  if (type == "consistent") {
    typeI.error.prob = 1- (1- multivariance.pvalue(consistent.rejection.level))^sum(number.of.tests,na.rm = TRUE)
  } else {
    typeI.error.prob = 1-(1-alpha)^groups.of.tests
  }


  if (verbose) cat(paste0("\n",ifelse(type == "conservative","Conservative estimate","Estimate")," of the type I error probability: ",signif(typeI.error.prob,3),"\n"))



  g$layout = layout_on_circles(g)

  return(invisible(list(cdms = list.cdm, values = m.values, graph = g,number.of.dep.tuples = number.of.dep.tuples,parameter.range = parameter.range,typeI.error.prob = typeI.error.prob,total.number.of.tests = number.of.tests,structure.type = "full", type = type,alpha = alpha, c.factor = c.factor)))

}

#' Returns the row indices of matrix A which match with B
#' Use the fast cpp implementation 'match_rows' instead.
#' Function here just for reference.
#'
#' @examples
#' # A = t(utils::combn(10,3))
#' # B = A[sort(sample.int(nrow(A),10)),]
#' # match.rows(A,B)
#'
#' @keywords internal
match.rows = function(A,B) {
  # for ordered matrices
  # rows of B have to be a subset of rows of A
  i = 1
  res = numeric(nrow(B))
  for (k in 1:nrow(B))
  {
    while(any(A[i,] != B[k,])) {
      i = i+1
    }
    res[k] = i
  }
  return(res)
}



#' check if lower order dependencies are present for the given tuple indices
#' here 'm.values' is a list of boolean matrices. Matrix [[k]] corresponds to the k tuples. For each number of tuples, the first columns of the matrix always contain the indices of the tuples
#'
#' @keywords internal
lower.order = function(tuple,m.values) {
  k = length(tuple)
  if( k >2) {
    tuples = t(utils::combn(tuple,k-1))

    #matched.rows = match.rows(m.values[[k-1]][,1:(k-1)],tuples) #R-implementation
    matched.rows = match_rows(m.values[[k-1]][,1:(k-1)],tuples)  #cpp-implementation

    return(any(as.logical(m.values[[k-1]][matched.rows,c(k,k+1)]),na.rm = TRUE))
  } else {
    return(FALSE) # in the pairwise case all are 1 independent
  }
}

#' cleanup dependence structure graph
#'
#' Given a dependence structure graph: vertices representing the multivariances of only two vertices can be turned into an edge labeled with the label of the vertex. Moreover, only subsets of the graph can be selected.
#'
#' @param g graph, created by \code{\link{dependence.structure}}
#' @param only.level integer vector, if provided all edges and dependency nodes corresponding to dependence orders not given in 'only.level' are removed
#' @param simplify.pairs boolean, if true dependency nodes which are only connected to two variables are turned into edges
#' @param drop.label.pairs boolean, if true the labels for edges indicating pairwise dependence are removed
#' @return graph
#'
#' @details
#' Note: The option 'only.level' works only properly for a full dependence structure graph, in the case of a clustered dependence structure graph dependency nodes representing a cluster might be removed.
#'
#' @examples
#' N = 200
#' y = coins(N,2)
#' x = cbind(y,y,y)
#'
#' ds = dependence.structure(x,structure.type = "clustered")
#' plot(clean.graph(ds$graph))
#' plot(clean.graph(ds$graph,only.level = 2))
#' plot(clean.graph(ds$graph,only.level = 3)) # of limited use for a clustered graph,
#' # i.e., here the three-dependence node without edges indicates that
#' # all edges were connected to clusters
#'
#' ds = dependence.structure(x,structure.type = "full")
#' plot(clean.graph(ds$graph))
#' plot(clean.graph(ds$graph,drop.label.pairs = TRUE))
#' plot(clean.graph(ds$graph,only.level = 2))
#' plot(clean.graph(ds$graph,only.level = 2,drop.label.pairs = TRUE))
#' plot(clean.graph(ds$graph,only.level = 3))
#'
#' @export
clean.graph = function(g,only.level = NULL,simplify.pairs = TRUE, drop.label.pairs = FALSE) {
  enumerate = FALSE # option which would allow to numerate the labels of the dependency nodes - maybe useful for rearrangement of the graph

  #  g = ds$graph

  if (!is.null(only.level)) {
    to.delete = which(!(igraph::V(g)$level %in% c(only.level,NA)))
    g = igraph::delete.vertices(g,to.delete)

    if (enumerate) {
      vn = length(igraph::V(g)) #total number of vertices
      n = sum(is.na(igraph::V(g)$level)) # initial number of vertices
      #igraph::V(g)$label = numeric(vn)
      igraph::V(g)$label[!is.na(igraph::V(g)$level)] = paste0(1:(vn-n),".")
    }

  }

  if (simplify.pairs) {
    vert = which(igraph::V(g)$level == 2) # might still have more neighbors!!!
    only.two.neighbors = logical(igraph::vcount(g)) # is initialized with FALSE
    for (i in vert) {
      only.two.neighbors[i] = length(igraph::neighbors(g,igraph::V(g)[i]))==2
      if (only.two.neighbors[i]) {
        if (drop.label.pairs) {
          g = igraph::add_edges(g, igraph::neighbors(g,i), weight= NA, color = 2, lty=2)
        } else {
          g = igraph::add_edges(g, igraph::neighbors(g,i), weight= NA, color = 2, lty=2,label = igraph::V(g)[i]$label)
        }
      }
    }
    g = igraph::delete.vertices(g,only.two.neighbors)
  }

  if (!is.null(g$layout)) g$layout = layout_on_circles(g)

  return(g)
}


#' calculates the coordinates of n points on a circle of radius r
#' if only 1 inner point, then it is placed in the center
#' @keywords internal
circle.coordinates = function(n,r = 0.5,add.angle = 0) {
  if (n == 0) {
    return(NULL)
  } else {
    if (n == 1) r = 0
    #   print(paste("on circle:",n))
    angles = (1:n)*2*pi/n + add.angle
    coords = matrix(c(r*sin(angles),r*cos(angles)),ncol = 2)
    #   print(coords)
    return(coords)
  }
}

#' special igraph layout for the dependence structure visualization
#'
#' It places the variable nodes on an outer circle and the dependency nodes on an inner circle
#' @param g graph
#' @param n number of vertices on outer circle
#'
#' @details
#' This is the standard layout for the full dependence structure, since in this case there often too many nodes which make the other (usual) layout incomprehensible.
#'
#' @examples
#' N = 200
#' y = coins(N,2)
#' x = cbind(y,y,y)
#'
#' g = dependence.structure(x,structure.type = "clustered",verbose = FALSE)$graph
#' plot(g)
#' plot(g,layout = layout_on_circles(g))
#' @export
layout_on_circles = function(g,n = sum(is.na(igraph::V(g)$level))) {
  vn = length(igraph::V(g))
  return(rbind(circle.coordinates(n,1),circle.coordinates(vn-n)))
}


# Utility functions ####


#' transforms a distance matrix to a matrix
#'
#' Does for a distance matrix generated via \code{dist} the same as \code{as.matrix} only slightly faster.
#'
#' @param ds a distance matrix object, e.g. generated by \code{\link{dist}}
#'
#' @keywords internal
dist.to.matrix = function(ds) {
  N=attr(ds,"Size")
  m = matrix(0, nrow = N, ncol = N)
  m[outer(1:N,1:N,">")] = ds
  m+t(m)
}

#' checks if a matrix is doubly centered
#'
#' !works only for the biased estimators
#'
#' @param mat matrix
#'
#' @examples
#' multivariance:::is.doubly.centered(as.matrix(dist(rnorm(10))))
#' multivariance:::is.doubly.centered(multivariance:::double.center(as.matrix(dist(rnorm(10)))))
#'
#' @keywords internal
is.doubly.centered = function(mat) {
  zero = rep(0,nrow(mat))
  return(isTRUE(all.equal(rowMeans(mat),zero,check.attributes = FALSE)) &
         isTRUE(all.equal(colMeans(mat),zero,check.attributes = FALSE)))
}


#' transforms a p-value into the corresponding label
#'
#' @param pv p-value
#'
#' @keywords internal
p.value.to.star.label = function(pv) {
 # inefficient code
 lab = ifelse(pv  > 0.05  ,"ns"  ,pv)
 lab = ifelse(pv <= 0.05  ,"*"   ,lab)
 lab = ifelse(pv <= 0.01  ,"**"  ,lab)
 lab = ifelse(pv <= 0.001 ,"***" ,lab)
 lab = ifelse(pv <= 0.0001,"****",lab)
 return(lab)
}


#' estimate of the computation time
#'
#' Estimates the computation time. This is relative rough. First run with \code{determine.parameters = TRUE} (which takes a while). Then use the computed parameters to determine the computation time/or sample size.
#'
#' @param determine.parameters if \code{TRUE} then the parameters for the current computer are determined. This might take a while (3 loops to N=1000).
#' @param N number of samples. If \code{NULL} and \code{sectime} is given, then \code{N} is computed.
#' @param n number of variables
#' @param sectime desired computation time in seconds. If \code{NULL} then the required computation time is computed.
#' @param coef.cdm computation time parameter for the doubly centered distance matrices
#' @param coef.prod computation time parameter for matrix products
#' @param coef.sum computation time parameter for matrix sums
#'
#' @details When detecting the parameters, the median of the computation times is used.
#'
#' @examples
#' Ns = (1:100)*10
#' ns = 1:100
#' fulltime = outer(Ns,ns,FUN = function(N,n) multivariance.timing(N,n))
#' contour(Ns,ns,fulltime,xlab = "N",ylab = "n",
#'  main = "computation time of multivariance in secs",
#'  sub = "using default parameters -
#'  use 'determine.parameters = TRUE' to compute machine specific values")
#'
#' # Run to determine the parameters of your system:
#' # multivariance.timing(determine.parameters = TRUE)
#'
#' @export
multivariance.timing = function(N=NULL,n,sectime = NULL,coef.cdm=15.2,coef.prod=2.1,coef.sum = 1.05,determine.parameters = FALSE) {
# coef.cdm.Nssq coef.prod.Nssq  coef.sum.Nssq  machine
#   8.808793       1.455658       0.932342     i5-1035G4 CPU @ 1.10GHz, 1498 MHz, 4 cores, 8 (logical) processors

    if (!determine.parameters) {
    const = coef.cdm*n+coef.prod*(n-1)+coef.sum
    if (!is.null(N) & is.null(sectime)) {
      nanotime = const * N^2
      return(nanotime/10^9)
    }
    if (is.null(N) & !is.null(sectime)) {
      N = sqrt(sectime*10^9/const)
      return(ceiling(N))
    }
  } else {
    op = graphics::par(mfrow = c(2,2))
    Ns = (1:100)*10
    Nssq = Ns^2
    times = numeric(length(Ns))

    for (i in 1:length(Ns)) {
      N = Ns[i]
      x = stats::rnorm(N)

      mb = microbenchmark::microbenchmark(cdm(x)) # timing of cdm
      print(N)
      print(mb)
      times[i] = stats::median(mb$time)
    }
    res = stats::lm(times ~ Nssq -1)
    graphics::plot(Ns^2,times,main = paste("cdm of NxN - ",version$version.string),sub = res$coefficients)
    graphics::abline(res)
    coef.cdm = res$coefficients

    Ns = (1:100)*10
    prodtimes = numeric(length(Ns))

    mat = matrix(stats::rnorm(max(Ns)^2),nrow = max(Ns))

    for (i in 1:length(Ns)) {
      N = Ns[i]
      x = mat[1:N,1:N]

      mb = microbenchmark::microbenchmark(x*x) #timing of matrix products
      print(N)
      print(mb)
      prodtimes[i] = stats::median(mb$time)
    }
    res = stats::lm(prodtimes ~ Nssq -1)
    graphics::plot(Ns^2,prodtimes,main = paste("products of NxN - ",version$version.string),sub = res$coefficients)
    graphics::abline(res)
    coef.prod = res$coefficients

    Ns = (1:100)*10
    sumtimes = numeric(length(Ns))

    for (i in 1:length(Ns)) {
      N = Ns[i]
      x = mat[1:N,1:N]

      mb = microbenchmark::microbenchmark(sum(x)) # timing of matrix sums
      print(N)
      print(mb)
      sumtimes[i] = stats::median(mb$time)
    }
    res = stats::lm(sumtimes ~ Nssq -1)
    graphics::plot(Ns^2,sumtimes,main = paste("sum of NxN elem - ",version$version.string),sub = res$coefficients)
    graphics::abline(res)
    coef.sum = res$coefficients


    ns = 1:100

    fulltime = outer(Ns,ns,FUN = function(N,n) multivariance.timing(N,n,sectime = NULL, coef.cdm,coef.prod,coef.sum))

    graphics::contour(Ns,ns,fulltime,xlab = "N",ylab = "n", main = "computation time of multivariance in secs")

    graphics::par(op)

    print( "\nThe computed parameter are:")
    return(c( coef.cdm = coef.cdm,coef.prod = coef.prod,coef.sum = coef.sum))

  }
}

#' sign preserving square root
#' @keywords internal
signed.sqrt = function(x) {
  sign(x)*sqrt(abs(x))
}


#' Simple integer hash from text
#'
#' Used to compute a random seed for \code{set.seed} based on the example name,
#' in order to aviod arbitrary seeds like '1234'.
#'
#' @examples
#' multivariance:::simple.int.hash("dep_struct_several_26_100")
#' @keywords internal
simple.int.hash = function(x) {
  ra = charToRaw(x)
  return(as.integer(sum(strtoi(ra,16L)*2^(1:length(ra))) %% 2^31))
}

