//*--------------------------------------------------------------------------*//
//  File:               depthsCurveHfsp.cpp
//  Created by:         Pavlo Mozharovskyi
//  First released:     01.11.2018
//
//  Contains R-level Rcpp-interfaced functions for computing Tukey curve depth
//  and related functionality.
//
//  For a description of the algorithms, see:
//  Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//  Depth for curve data and applications.
//
//  Subsequent changes are listed below:
//  01.11.2018 (Pavlo Mozharovskyi): First version released.
//  21.01.2019 (Pavlo Mozharovskyi): Corrected
//    in function "depthCurveTukey" [[Rcpp::export(depth.curve.Tukey)]]
//    3.b) Subsample points of the object is necessary:
//      3.b.i) Subsample "refEmpDist", line 622
//  19.02.2019 (Pavlo Mozharovskyi): removed abind call from the example
//    to the function "distImages" [[Rcpp::export(dist.images)]]
//*--------------------------------------------------------------------------*//

// [[Rcpp::depends(RcppArmadillo)]]

#include "stdafx.h"

//' Distance for images
//'
//' Calculates distance matrix for a sample of images using the minimax metric.
//' This fucntion can be seen as a wrapper of a sequential call of
//' \code{images2curves} and \code{dist.curves}.
//'
//' @param images A 3-dimensional array with each slice (matrix in first two
//' dimensions) corresponding to an image. Each (eps-strictly) positive
//' entry is regarded as an occupied pixel (one), otherwise it is regarded
//' as an empty pixel, of an image.
//'
//' @param verbosity Level of reporting messages, the higher the more progress
//' reports are printed, set \code{0} (default) for no messages.
//'
//' @return A matrix \code{dim(images)[3] x dim(images)[3]} with each entry
//' being the distance between two images.
//'
//' @references Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//' Depth for curve data and applications.
//'
//' @examples
//' library(curveDepth)
//' # Pixel-grid filling function for an image
//' plotGridImage <- function(dims){
//'   redDims1 <- dims[1] - 1
//'   redDims2 <- dims[2] - 1
//'   for (i in 1:(dims[1] - 1)){
//'     lines(c(i / redDims1 - 0.5 / redDims1,
//'             i / redDims1 - 0.5 / redDims1),
//'           c(0 - 0.5 / redDims2, 1 + 0.5 / redDims2),
//'           lwd = 1, lty = 3, col = "lightgray")
//'     lines(c(0 - 0.5 / redDims1, 1 + 0.5 / redDims1),
//'           c(i / redDims2 - 0.5 / redDims2,
//'             i / redDims2 - 0.5 / redDims2),
//'           lwd = 1, lty = 3, col = "lightgray")
//'   }
//'   rect(0 - 0.5 / redDims1, 0 - 0.5 / redDims2,
//'        1 + 0.5 / redDims1, 1 + 0.5 / redDims2)
//' }
//' # Load two Sevens and one One, and plot them
//' data("mnistShort017")
//' # First Seven
//'   firstSevenDigit <- mnistShort017$`7`[, , 5]
//' image(as.matrix(rev(as.data.frame(firstSevenDigit))),
//'       col = gray((255:0) / 256), asp = 1,
//'       xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'       ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(firstSevenDigit)[1:2])
//' # Second Seven
//' secondSevenDigit <- mnistShort017$`7`[, , 6]
//' image(as.matrix(rev(as.data.frame(secondSevenDigit))),
//'       col = gray((255:0) / 256), asp = 1,
//'       xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'       ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(secondSevenDigit)[1:2])
//' # A One
//' aOneDigit <- mnistShort017$`1`[, , 1]
//' image(as.matrix(rev(as.data.frame(aOneDigit))),
//'       col = gray((255:0) / 256), asp = 1,
//'       xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'       ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(aOneDigit)[1:2])
//' # Caculate distances between all the images
//' threeDigits <- array(NA, dim = c(nrow(firstSevenDigit),
//'   ncol(firstSevenDigit), 3))
//' threeDigits[, , 1] <- firstSevenDigit
//' threeDigits[, , 2] <- secondSevenDigit
//' threeDigits[, , 3] <- aOneDigit
//' distMatrix <- dist.images(threeDigits)
//' # Print distance matrix
//' print(distMatrix)
// [[Rcpp::export(dist.images)]]
NumericMatrix distImages(arma::cube images, int verbosity = 0){
  // 0. Prepare C++ structures according to input-data dimensions
  int d = 2; // images' dimensions
  Curve** curves = new Curve*[images.n_slices];
  int* dims = new int[d];
  for (int i = 0; i < images.n_slices; i++){
    // Read input image
    dims[0] = images.slice(i).n_rows;
    dims[1] = images.slice(i).n_cols;
    ImageDensity imageDens(d, dims, images.slice(i).memptr());
    // Convert it to a curve
    //Rcout << "Start transforming curve " << i << std::endl;
    curves[i] = imageDens.toOrderedCurve();
    //Rcout << "Finish transforming curve " << i << std::endl;
  }
  // 1. Prepare output structure
  NumericMatrix similarities(images.n_slices, images.n_slices);
  // 2. Go through all pairs of images
  for (int i = 0; i < images.n_slices - 1; i++){
    if (verbosity > 0){
      Rcout << "Calculate distances from image " << i  << ": ";
    }
    for (int j = i + 1; j < images.n_slices; j++){
      // 2. Calculate the curve space metric
      similarities(i,j) = curves[i]->distCurve(*curves[j], false);
      similarities(j,i) = similarities(i,j);
      if (verbosity > 0){
        Rcout << j << " ";
      }
    }
    if (verbosity > 0){
      Rcout << "done for image " << i << "." << std::endl;
    }
  }
  // 3. Release memory
  delete[] dims;
  delete[] curves;
  // Return the similarity matrix
  return similarities;
}

//' Distance for curves
//'
//' Calculates distance matrix for a sample of curves using the minimax metric.
//'
//' @param curves A list where each element is a function being a list
//' containing a matrix \code{coords} (values, d columns).
//'
//' @param oneWay Whether curves should be condisered as a one-directional,
//' \code{FALSE} by default.
//'
//' @param verbosity Level of reporting messages, the higher the more progress
//' reports are printed, set \code{0} (default) for no messages.
//'
//' @return A matrix \code{length(curves) x length(curves)} with each entry
//' being the distance between two curves.
//'
//' @references Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//' Depth for curve data and applications.
//'
//' @examples
//' library(curveDepth)
//' # Pixel-grid filling function for an image
//' plotGridImage <- function(dims){
//'   redDims1 <- dims[1] - 1
//'   redDims2 <- dims[2] - 1
//'   for (i in 1:(dims[1] - 1)){
//'     lines(c(i / redDims1 - 0.5 / redDims1,
//'             i / redDims1 - 0.5 / redDims1),
//'           c(0 - 0.5 / redDims2, 1 + 0.5 / redDims2),
//'           lwd = 1, lty = 3, col = "lightgray")
//'     lines(c(0 - 0.5 / redDims1, 1 + 0.5 / redDims1),
//'           c(i / redDims2 - 0.5 / redDims2,
//'             i / redDims2 - 0.5 / redDims2),
//'           lwd = 1, lty = 3, col = "lightgray")
//'   }
//'   rect(0 - 0.5 / redDims1, 0 - 0.5 / redDims2,
//'        1 + 0.5 / redDims1, 1 + 0.5 / redDims2)
//' }
//' # Load two Sevens and one One, plot them,
//' # and transform to curves
//' data("mnistShort017")
//' # First Seven
//' firstSevenDigit <- mnistShort017$`7`[, , 5]
//' image(as.matrix(rev(as.data.frame(firstSevenDigit))),
//'       col = gray((255:0) / 256), asp = 1,
//'       xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'       ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(firstSevenDigit)[1:2])
//' firstSevenCurve <- images2curves(array(
//'   firstSevenDigit, dim = c(28, 28, 1)))[[1]]
//' # Second Seven
//' secondSevenDigit <- mnistShort017$`7`[, , 6]
//' image(as.matrix(rev(as.data.frame(secondSevenDigit))),
//'       col = gray((255:0) / 256), asp = 1,
//'       xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'       ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(secondSevenDigit)[1:2])
//' secondSevenCurve <- images2curves(array(
//'   secondSevenDigit, dim = c(28, 28, 1)))[[1]]
//' # A One
//' aOneDigit <- mnistShort017$`1`[, , 1]
//' image(as.matrix(rev(as.data.frame(aOneDigit))),
//'       col = gray((255:0) / 256), asp = 1,
//'       xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'       ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(aOneDigit)[1:2])
//' aOneCurve <- images2curves(array(
//'   aOneDigit, dim = c(28, 28, 1)))[[1]]
//' # Caculate distances between all the curves
//' distMatrix <- dist.curves(list(
//'   firstSevenCurve, secondSevenCurve, aOneCurve))
//' # Print distance matrix
//' print(distMatrix)
// [[Rcpp::export(dist.curves)]]
NumericMatrix distCurves(List curves, bool oneWay = false, int verbosity = 0){
  // 0. Read the input via Rcpp
  int n = curves.size();
  Curve* cCurves = new Curve[n];
  for (int i = 0; i < n; i++){
    cCurves[i].type = VOXCOORDS;
    NumericMatrix rcppCurveVals = as<List>(curves(i))["coords"];
    std::vector<double> stlCurveVals = Rcpp::as<std::vector<double> >
      (transpose(rcppCurveVals));
    cCurves[i].n = rcppCurveVals.nrow();
    cCurves[i].d = rcppCurveVals.ncol();
    cCurves[i].rawVals = new double[cCurves[i].n * cCurves[i].d];
    memcpy(cCurves[i].rawVals, stlCurveVals.data(),
           cCurves[i].n * cCurves[i].d * sizeof(double));
    cCurves[i].vals = asMatrix(cCurves[i].rawVals, cCurves[i].n, cCurves[i].d);
  }
  // 1. Prepare output structure
  NumericMatrix similarities(n, n);
  // 2. Go through all pairs of curves
  for (int i = 0; i < n - 1; i++){
    if (verbosity > 0){
      Rcout << "Calculate distances from curve " << i  << ": ";
    }
    for (int j = i + 1; j < n; j++){
      // 2. Calculate the curve space metric
      similarities(i,j) = cCurves[i].distCurve(cCurves[j], oneWay);
      similarities(j,i) = similarities(i,j);
      if (verbosity > 0){
        Rcout << j << " ";
      }
    }
    if (verbosity > 0){
      Rcout << "done for curve " << i << "." << std::endl;
    }
  }
  // 3. Release memory
  delete[] cCurves;
  // Return the similarity matrix
  return similarities;
}

//' Distance for curves
//'
//' Calculates distance matrix for two samples of curves using minimax metric.
//' The function can be particularly useful for parallel computation of a
//' big distance matrix.
//'
//' @param curvesRows A list where each element is a function being a list
//' containing a matrix \code{coords} (values, d columns).
//'
//' @param curvesCols A list where each element is a function being a list
//' containing a matrix \code{coords} (values, d columns).
//'
//' @param oneWay Whether curves should be condisered as a one-directional,
//' \code{FALSE} by default.
//'
//' @param verbosity Level of reporting messages, the higher the more progress
//' reports are printed, set \code{0} (default) for no messages.
//'
//' @return A matrix \code{length(curvesRows) x length(curvesCols)} with each
//' entry being the distance between two corresponding curves.
//'
//' @references Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//' Depth for curve data and applications.
//'
//' @examples
//' library(curveDepth)
//' # Pixel-grid filling function for an image
//' plotGridImage <- function(dims){
//'   redDims1 <- dims[1] - 1
//'   redDims2 <- dims[2] - 1
//'   for (i in 1:(dims[1] - 1)){
//'     lines(c(i / redDims1 - 0.5 / redDims1,
//'             i / redDims1 - 0.5 / redDims1),
//'           c(0 - 0.5 / redDims2, 1 + 0.5 / redDims2),
//'           lwd = 1, lty = 3, col = "lightgray")
//'     lines(c(0 - 0.5 / redDims1, 1 + 0.5 / redDims1),
//'           c(i / redDims2 - 0.5 / redDims2,
//'             i / redDims2 - 0.5 / redDims2),
//'           lwd = 1, lty = 3, col = "lightgray")
//'   }
//'   rect(0 - 0.5 / redDims1, 0 - 0.5 / redDims2,
//'        1 + 0.5 / redDims1, 1 + 0.5 / redDims2)
//' }
//' # Load two Sevens and one One, plot them,
//' # and transform to curves
//' data("mnistShort017")
//' # First Seven
//' firstSevenDigit <- mnistShort017$`7`[, , 5]
//' image(as.matrix(rev(as.data.frame(firstSevenDigit))),
//'       col = gray((255:0) / 256), asp = 1,
//'       xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'       ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(firstSevenDigit)[1:2])
//' firstSevenCurve <- images2curves(array(
//'   firstSevenDigit, dim = c(28, 28, 1)))[[1]]
//' # Second Seven
//' secondSevenDigit <- mnistShort017$`7`[, , 6]
//' image(as.matrix(rev(as.data.frame(secondSevenDigit))),
//'       col = gray((255:0) / 256), asp = 1,
//'       xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'       ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(secondSevenDigit)[1:2])
//' secondSevenCurve <- images2curves(array(
//'   secondSevenDigit, dim = c(28, 28, 1)))[[1]]
//' # A One
//' aOneDigit <- mnistShort017$`1`[, , 1]
//' image(as.matrix(rev(as.data.frame(aOneDigit))),
//'   col = gray((255:0) / 256), asp = 1,
//'   xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'   ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(aOneDigit)[1:2])
//' aOneCurve <- images2curves(array(
//'   aOneDigit, dim = c(28, 28, 1)))[[1]]
//' # Caculate distances between all the curves
//' distMatrix <- matrix(0, 3, 3)
//' distMatrix[3, 1:2] <- distMatrix[1:2, 3] <-
//'   dist.curves.asymm(list(
//'     firstSevenCurve, secondSevenCurve), list(aOneCurve))
//' distMatrix[2, 1] <- distMatrix[1, 2] <-
//'   dist.curves.asymm(
//'     list(firstSevenCurve), list(secondSevenCurve))
//' # Print distance matrix
//' print(distMatrix)
// [[Rcpp::export(dist.curves.asymm)]]
NumericMatrix distCurvesAsymm(List curvesRows, List curvesCols,
                              bool oneWay = false, int verbosity = 0){
  // 0. Read the input via Rcpp
  // 0.a Curves forming rows of the distance matrix
  int nRows = curvesRows.size();
  Curve* cCurvesRows = new Curve[nRows];
  for (int i = 0; i < nRows; i++){
    cCurvesRows[i].type = VOXCOORDS;
    NumericMatrix rcppCurveVals = as<List>(curvesRows(i))["coords"];
    std::vector<double> stlCurveVals = Rcpp::as<std::vector<double> >
      (transpose(rcppCurveVals));
    cCurvesRows[i].n = rcppCurveVals.nrow();
    cCurvesRows[i].d = rcppCurveVals.ncol();
    cCurvesRows[i].rawVals = new double[cCurvesRows[i].n * cCurvesRows[i].d];
    memcpy(cCurvesRows[i].rawVals, stlCurveVals.data(),
           cCurvesRows[i].n * cCurvesRows[i].d * sizeof(double));
    cCurvesRows[i].vals = asMatrix(cCurvesRows[i].rawVals,
                                   cCurvesRows[i].n, cCurvesRows[i].d);
  }
  // 0.a Curves forming columns of the distance matrix
  int nCols = curvesCols.size();
  Curve* cCurvesCols = new Curve[nCols];
  for (int i = 0; i < nCols; i++){
    cCurvesCols[i].type = VOXCOORDS;
    NumericMatrix rcppCurveVals = as<List>(curvesCols(i))["coords"];
    std::vector<double> stlCurveVals = Rcpp::as<std::vector<double> >
      (transpose(rcppCurveVals));
    cCurvesCols[i].n = rcppCurveVals.nrow();
    cCurvesCols[i].d = rcppCurveVals.ncol();
    cCurvesCols[i].rawVals = new double[cCurvesCols[i].n * cCurvesCols[i].d];
    memcpy(cCurvesCols[i].rawVals, stlCurveVals.data(),
           cCurvesCols[i].n * cCurvesCols[i].d * sizeof(double));
    cCurvesCols[i].vals = asMatrix(cCurvesCols[i].rawVals,
                                   cCurvesCols[i].n, cCurvesCols[i].d);
  }
  // 1. Prepare output structure
  NumericMatrix similarities(nRows, nCols);
  // 2. Go through all pairs of curves
  for (int i = 0; i < nRows; i++){
    if (verbosity > 0){
      Rcout << "Calculate distances from curve " << i  << ": ";
    }
    for (int j = 0; j < nCols; j++){
      // 2. Calculate the curve space metric
      similarities(i,j) = cCurvesRows[i].distCurve(cCurvesCols[j], oneWay);
      if (verbosity > 0){
        Rcout << j << " ";
      }
    }
    if (verbosity > 0){
      Rcout << "done for curve " << i << "." << std::endl;
    }
  }
  // 3. Release memory
  delete[] cCurvesRows;
  delete[] cCurvesCols;
  // Return the similarity matrix
  return similarities;
}

//' Calculate Tukey curve depth using given points
//'
//' Calculates Tukey curve depth of each curve in \code{objects} w.r.t. the
//' sample of curves in \code{data}. Calculation of partial depth of each
//' single point can be either exact or approximate. If exact, modified method
//' of Dyckerhoff and Mozharovskyi (2016) is used; if approximate,
//' approximation is performed by projections on directions - points uniformly
//' distributed on the unit hypersphere.
//'
//' @param objects A list where each element is a multivariate curve being a
//' list containing a matrix \code{coords} (values, d columns).
//'
//' @param data A list where each element is a multivariate curve being a list
//' containing a matrix \code{coords} (values, d columns). The depths are
//' computed w.r.t. this data set.
//'
//' @param nDirs Number of directions used to inspect the space, drawn from the
//' uniform distribution on the sphere.
//'
//' @param subs Whether to split each object into two disjunctive subsets (one
//' for integrating and one for estimation) when computing the depth.
//'
//' @param fracInt Portion of an object used for integrating.
//'
//' @param fracEst Portion of an object used for estimation,
//' maximum: \code{1 - fracInt}.
//'
//' @param subsamples A list indicating subsamples of points for each curve in
//' \code{objects}. Each elemnt of the list corresponds to a single curve and
//' should be given as a vector of the length equal to the number of points
//' on it, with entries indicating:
//' \itemize{
//'   \item 0 do not take the point into account at all,
//'   \item 1 use point as a reference (i.e. for integrating) and thus
//'          calculate its depth,
//'   \item 2 utilize point in depth calculation (i.e. for estimation).}
//'
//' @param exactEst Is calculation of depth for each reference point of the
//' curve exact (\code{TRUE}, by default) or approximate (\code{FALSE}).
//'
//' @param minMassObj Minimal portion of the \code{objects} distribution in the
//' halfspace to be considered when calculating depth.
//'
//' @param minMassDat minimal portion of the \code{data} distribution in the
//' halfspace to be considered when calculating depth.
//'
//' @return A vector of doubles having the same length as \code{objects}, whose
//' each entry is the depth of each element of \code{objects} w.r.t.
//' \code{data}.
//'
//' @references Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//' Depth for curve data and applications.
//'
//' Dyckerhoff, R. and Mozharovskyi P. (2016). Exact computation of the
//' halfspace depth. \emph{Computational Statistics and Data Analysis}, 98,
//' 19-30.
//'
//' @examples
//' library(curveDepth)
//' # Load digits and transform them to curves
//' data("mnistShort017")
//' n <- 10 # cardinality of each class
//' m <- 50 # number of points to sample
//' cst <- 1/10 # a threshold constant
//' alp <- 1/8 # a threshold constant
//' curves0 <- images2curves(mnistShort017$`0`[, , 1:n])
//' curves1 <- images2curves(mnistShort017$`1`[, , 1:n])
//' set.seed(1)
//' curves0Smpl <- sample.curves(curves0, 2 * m)
//' curves1Smpl <- sample.curves(curves1, 2 * m)
//' # Calculate depths
//' depthSpace = matrix(NA, nrow = n * 2, ncol = 2)
//' depthSpace[, 1] = depth.curve.Tukey(
//'   c(curves0Smpl, curves1Smpl), curves0Smpl,
//'   exactEst = TRUE, minMassObj = cst/m^alp)
//' depthSpace[, 2] = depth.curve.Tukey(
//'   c(curves0Smpl, curves1Smpl), curves1Smpl,
//'   exactEst = TRUE, minMassObj = cst/m^alp)
//' # Draw the DD-plot
//' plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
//'      xlab = paste("Depth w.r.t. '0'"),
//'      ylab = paste("Depth w.r.t. '1'"),
//'      main = paste("DD-plot for '0' vs '1'"))
//' grid()
//' # Draw the separating rule
//' dat1 <- data.frame(cbind(
//'   depthSpace, c(rep(0, n), rep(1, n))))
//' ddalpha1 <- ddalpha.train(X3 ~ X1 + X2, data = dat1,
//'                           depth = "ddplot",
//'                           separator = "alpha")
//' ddnormal <- ddalpha1$classifiers[[1]]$hyperplane[2:3]
//' pts <- matrix(c(0, 0, 1, ddnormal[1] / -ddnormal[2]),
//'               nrow = 2, byrow = TRUE)
//' lines(pts, lwd = 2)
//' # Draw the points
//' points(depthSpace[1:n, ],
//'        col = "red", lwd = 2, pch = 1)
//' points(depthSpace[(n + 1):(2 * n), ],
//'        col = "blue", lwd = 2, pch = 3)
// [[Rcpp::export(depth.curve.Tukey)]]
NumericVector depthCurveTukey(List objects, List data, int nDirs = 100,
                              bool subs = true,
                              double fracInt = 0.5, double fracEst = 0.5,
                              List subsamples = R_NilValue,
                              bool exactEst = true,
                              double minMassObj = 0, double minMassDat = 0){
  // 1. Read the input via Rcpp:
  // 1.a) Data sample ("data"):
  int n = data.size();
  Curve* curvesData = new Curve[n];
  for (int i = 0; i < n; i++){
    curvesData[i].type = VOXCOORDS; // voxels are specified directly
    NumericMatrix rcppCurveDataVals = as<List>(data(i))["coords"];
    // Convert to stl to obtain a pointer to the data
    std::vector<double> stlCurveDataVals = Rcpp::as<std::vector<double> >
      (transpose(rcppCurveDataVals));
    // Set voxel dimensions
    curvesData[i].n = rcppCurveDataVals.nrow();
    curvesData[i].d = rcppCurveDataVals.ncol();
    // Copy and double-index data
    curvesData[i].rawVals = new double[curvesData[i].n * curvesData[i].d];
    memcpy(curvesData[i].rawVals, stlCurveDataVals.data(),
           curvesData[i].n * curvesData[i].d * sizeof(double));
    curvesData[i].vals = asMatrix(curvesData[i].rawVals, curvesData[i].n,
                                  curvesData[i].d); // new d.*
    // DEBUG
    // Rcout << curvesData[i].n << " " << curvesData[i].d << "\n";
    // for (int j = 0; j < curvesData[i].n; j++){
    //   for (int k = 0; k < curvesData[i].d; k++){
    //     Rcout << curvesData[i].vals[j][k] << " ";
    //   }
    //   Rcout << "\n";
    // }
    // DEBUG (end)
  }
  // Convert the bunch of curves to an empirical distribution
  EmpDist genEmpDist(curvesData, n, 0.000001);
  // DEBUG
  // Rcout << "genEmpDist:\n";
  // Rcout << genEmpDist.n << " " << genEmpDist.d << "\n";
  // DEBUG (end)
  genEmpDist.updateWeights(true);
  // DEBUG
  // Rcout << "genEmpDist:\n";
  // Rcout << genEmpDist.n << " " << genEmpDist.d << "\n";
  // DEBUG (end)
  //  Rcout << genEmpDist.n << std::endl;
  /*
  Rcout << "genEmpDist" << std::endl;
  Rcout << "vals (weights)" << std::endl;
  for (int j = 0; j < genEmpDist.n; j++){
    for (int k = 0; k < genEmpDist.d; k++){
      Rcout << genEmpDist.vals[j][k] << " ";
    }
    Rcout << "(" << genEmpDist.weights[j] << ")";
    Rcout << std::endl;
  }
   */

  // 1.b) Curves for which depth should be computed ("objects"):
  int m = objects.size();
  Curve* curvesObj = new Curve[m];
  for (int i = 0; i < m; i++){
    curvesObj[i].type = VOXCOORDS;
    NumericMatrix rcppCurveVals = as<List>(objects(i))["coords"];
    // Convert to stl to obtain a pointer to the data
    std::vector<double> stlCurveVals = Rcpp::as<std::vector<double> >
      (transpose(rcppCurveVals));
    // Set voxel dimensions
    curvesObj[i].n = rcppCurveVals.nrow();
    curvesObj[i].d = rcppCurveVals.ncol();
    // Copy and double-index data
    curvesObj[i].rawVals = new double[curvesObj[i].n * curvesObj[i].d];
    memcpy(curvesObj[i].rawVals, stlCurveVals.data(),
           curvesObj[i].n * curvesObj[i].d * sizeof(double));
    curvesObj[i].vals = asMatrix(curvesObj[i].rawVals, curvesObj[i].n,
                                 curvesObj[i].d); // new d.*
  }
  // 2. Generate random directions if necessary
  int d = curvesObj[0].d;
  double* dirsRaw = 0;
  double** dirs = 0;
  if (!exactEst){
    dirsRaw = new double[nDirs * d];
    dirs = asMatrix(dirsRaw, nDirs, d); // new double*
    generateDirections(1, nDirs, d, dirs);
  }
  // 3. The loop for calculating depths (of each "object"):
  NumericVector depths(m);
  for (int i = 0; i < m; i++){
    // 3.a) Convert current object to two empirical distributions:
    EmpDist refEmpDist(curvesObj + i, 1, 0.000001);
    EmpDist curEmpDist(curvesObj + i, 1, 0.000001);
    refEmpDist.updateWeights(true);
    curEmpDist.updateWeights(true);
    /*
    Rcout << "refEmpDist[" << i << "]" << std::endl;
    Rcout << "vals (weights)" << std::endl;
    for (int j = 0; j < refEmpDist.n; j++){
      for (int k = 0; k < refEmpDist.d; k++){
        Rcout << refEmpDist.vals[j][k] << " ";
      }
      Rcout << "(" << refEmpDist.weights[j] << ")";
      Rcout << std::endl;
    }

    Rcout << "curEmpDist[" << i << "]" << std::endl;
    Rcout << "vals (weights)" << std::endl;
    for (int j = 0; j < curEmpDist.n; j++){
      for (int k = 0; k < curEmpDist.d; k++){
        Rcout << curEmpDist.vals[j][k] << " ";
      }
      Rcout << "(" << curEmpDist.weights[j] << ")";
      Rcout << std::endl;
    }
     */

    // 3.b) Subsample points of the object is necessary:
    if (subs){
      // Here, all points available in the object "EmpDist" are considered
      // 3.b.i) Subsample "refEmpDist":
      // Calculate number of points to calculate depth for
      int nFracInt = floor(refEmpDist.n * fracInt);
      // Draw subsample of this size sorted increasingly
      Environment base_env("package:base");
      Function base_sampleint = base_env["sample.int"];
      Function base_sort = base_env["sort"];
      IntegerVector inds = base_sort(base_sampleint(
        refEmpDist.n, nFracInt, false));
      // Set not chosen points' weights to zero
      int iInds = 0;
      // Rcout << "inds.size()=" << inds.size() << std::endl;
      for (int j = 0; j < refEmpDist.n; j++){
        // Rcout << "iInds=" << iInds << "; ";
        if(iInds < nFracInt && j == inds[iInds] - 1){
          iInds++;
        }else{
          refEmpDist.weights[j] = 0;
        }
      }
      // Rcout << std::endl;

      /*
      Rcout << "refEmpDist[" << i << "]" << std::endl;
      Rcout << "vals (weights)" << std::endl;
      for (int j = 0; j < refEmpDist.n; j++){
        for (int k = 0; k < refEmpDist.d; k++){
          Rcout << refEmpDist.vals[j][k] << " ";
        }
        Rcout << "(" << refEmpDist.weights[j] << ")";
        Rcout << std::endl;
      }
       */

      // 3.b.ii) Subsample "curEmpDist" taking "refEmpDist" into account:
      // Calculate number of points used for depth calculation
      int nFracEst = ceil(refEmpDist.n * fracEst);
      if (nFracEst > curEmpDist.n - nFracInt){ // if too much, then cut it down
        nFracEst = curEmpDist.n - nFracInt;
      }
      //Rcout << refEmpDist.n << " (" << nFracInt << ", " << nFracEst <<
      //  ")" << std::endl;
      // If all resting points should be used for depth estimation
      if (nFracEst == curEmpDist.n - nFracInt){
        // Set their weights to zero in "curEmpDist"
        for (int j = 0; j < refEmpDist.n; j++){
          // It is (clearly) assumed that curEmpDist.n = refEmpDist.n
          if (refEmpDist.weights[j] > 0){
            curEmpDist.weights[j] = 0;
          }
        }
      }else{
        // Determine the size of the subsample
        nFracEst = nFracEst / (curEmpDist.n - nFracInt);
        // Draw subsample of this size sorted increasingly
        Environment base_env("package:base");
        Function base_sampleint = base_env["sample.int"];
        Function base_sort = base_env["sort"];
        IntegerVector indsIgn = base_sort(base_sampleint(
          refEmpDist.n, nFracEst, false));
        // Set weights of the not chosen points' among the rest to zero
        int iZerosInDist = 0;
        int iIndsIgn = 0;
        for (int j = 0; j < curEmpDist.n; j++){
          // It is (clearly) assumed that curEmpDist.n = refEmpDist.n
          if (refEmpDist.weights[j] > 0){
            curEmpDist.weights[j] = 0;
          }else{
            if (iZerosInDist == indsIgn[iIndsIgn] - 1){
              iIndsIgn++;
            }else{
              curEmpDist.weights[j] = 0;
            }
            iZerosInDist++;
          }
        }
      }
      // Update weights of "refEmpDist" and of "curEmpDist" to sum up to one
      refEmpDist.updateWeights(true);
      curEmpDist.updateWeights(true);

      /*
      Rcout << "refEmpDist[" << i << "]" << std::endl;
      Rcout << "vals (weights)" << std::endl;
      for (int j = 0; j < refEmpDist.n; j++){
        for (int k = 0; k < refEmpDist.d; k++){
          Rcout << refEmpDist.vals[j][k] << " ";
        }
        Rcout << "(" << refEmpDist.weights[j] << ")";
        Rcout << std::endl;
      }

      Rcout << "curEmpDist[" << i << "]" << std::endl;
      Rcout << "vals (weights)" << std::endl;
      for (int j = 0; j < curEmpDist.n; j++){
        for (int k = 0; k < curEmpDist.d; k++){
          Rcout << curEmpDist.vals[j][k] << " ";
        }
        Rcout << "(" << curEmpDist.weights[j] << ")";
        Rcout << std::endl;
      }
       */

    }else{
      if (!Rf_isNull(subsamples)){ // whether subsamples indicated in input
        NumericVector refSubsample = subsamples(i);
        // If subsample length coincides with number of points in the object
        if (refSubsample.length() == refEmpDist.n){
//          Rcout << "Subsamples accepted." << std::endl;
          // Do manual subsampling
          for (int j = 0; j < refEmpDist.n; j++){
            switch ((int)refSubsample(j)){
            case 0:
              refEmpDist.weights[j] = 0;
              curEmpDist.weights[j] = 0;
              break;
            case 1:
              curEmpDist.weights[j] = 0;
              break;
            case 2:
              refEmpDist.weights[j] = 0;
              break;
            }
          }
          // Update weights of "refEmpDist" and of "curEmpDist" to sum up to one
          refEmpDist.updateWeights(true);
          curEmpDist.updateWeights(true);
        }
      }
    }
//    return depths;
    // 3.c) Calculate the depth
    if (exactEst){
      // Calculate exact depth of the curve's reference points
      depths(i) = calcOneDepth(refEmpDist, curEmpDist, genEmpDist, d,
             minMassObj, minMassDat);
    }else{
      // DEBUG
      // Rcout << "refEmpDist:\n";
      // Rcout << refEmpDist.n << " " << refEmpDist.d << "\n";
      // Rcout << "curEmpDist:\n";
      // Rcout << curEmpDist.n << " " << curEmpDist.d << "\n";
      // Rcout << "genEmpDist:\n";
      // Rcout << genEmpDist.n << " " << genEmpDist.d << "\n";
      // DEBUG (end)
      // Approximate depth of the curve's reference points
      depths(i) = approxOneDepth(refEmpDist, curEmpDist, genEmpDist, dirs,
             nDirs, d, minMassObj, minMassDat);
    }
  }
  // 4. Release memory and return the result
  delete[] curvesData;
  delete[] curvesObj;
  if (!exactEst){
    delete[] dirsRaw;
    delete[] dirs;
  }
  return depths;
}

//' Sample points on curves
//'
//' Samples points uniformly on curves interpolated as linear consequent
//' segments.
//'
//' @param curves A list where each element is a function being a list
//' containing a matrix \code{coords} (curves' values, d columns).
//'
//' @param ptsPerCurve A vector of numbers of points to be sampled on each
//' curve. If \code{length(ptsPerCurve) < length(curves)} then the first entry
//' of \code{ptsPerCurve} is considered only, and corresponds to the number of
//' points on a curve.
//'
//' @return A list of curves with each entry being a list constiting of [[1]]
//' the drawn curve being a matrix named \code{coords}, [[2]] length of the
//' curve as in \code{curves} named \code{length.init}, and [[3]] length of the
//' drawn curve named \code{length}.
//'
//' @references Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//' Depth for curve data and applications.
//'
//' @examples
//' library(curveDepth)
//' # Load digits and transform them to curves
//' data("mnistShort017")
//' n <- 10 # cardinality of each class
//' m <- 50 # number of points to sample
//' cst <- 1/10 # a threshold constant
//' alp <- 1/8 # a threshold constant
//' curves0 <- images2curves(mnistShort017$`0`[, , 1:n])
//' curves1 <- images2curves(mnistShort017$`1`[, , 1:n])
//' set.seed(1)
//' curves0Smpl <- sample.curves(curves0, 2 * m)
//' curves1Smpl <- sample.curves(curves1, 2 * m)
//' # Calculate depths
//' depthSpace = matrix(NA, nrow = n * 2, ncol = 2)
//' depthSpace[, 1] = depth.curve.Tukey(
//'   c(curves0Smpl, curves1Smpl), curves0Smpl,
//'   exactEst = TRUE, minMassObj = cst/m^alp)
//' depthSpace[, 2] = depth.curve.Tukey(
//'   c(curves0Smpl, curves1Smpl), curves1Smpl,
//'   exactEst = TRUE, minMassObj = cst/m^alp)
//' # Draw the DD-plot
//' plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
//'      xlab = paste("Depth w.r.t. '0'"),
//'      ylab = paste("Depth w.r.t. '1'"),
//'      main = paste("DD-plot for '0' vs '1'"))
//' grid()
//' # Draw the separating rule
//' dat1 <- data.frame(cbind(
//'   depthSpace, c(rep(0, n), rep(1, n))))
//' ddalpha1 <- ddalpha.train(X3 ~ X1 + X2, data = dat1,
//'                           depth = "ddplot",
//'                           separator = "alpha")
//' ddnormal <- ddalpha1$classifiers[[1]]$hyperplane[2:3]
//' pts <- matrix(c(0, 0, 1, ddnormal[1] / -ddnormal[2]),
//'               nrow = 2, byrow = TRUE)
//' lines(pts, lwd = 2)
//' # Draw the points
//' points(depthSpace[1:n, ],
//'        col = "red", lwd = 2, pch = 1)
//' points(depthSpace[(n + 1):(2 * n), ],
//'        col = "blue", lwd = 2, pch = 3)
// [[Rcpp::export(sample.curves)]]
List curvesSubsample(List curves, IntegerVector ptsPerCurve =
  IntegerVector::create(500)){
  // 0. Check input data
  // TODO: -||-
  // Rcout << ptsPerCurve << "\n";
  // 1. Draw samples on curves
  List res = List(); // the structure to collect all results
  // 1.1. Calculate length of each curve and the total length
  double lenTotal = 0;
  int maxNCurvePoints = 0;
  double* lengths = new double[curves.length()];
  NumericMatrix rcppCurveVals = as<List>(curves(0))["coords"];
  int dim = rcppCurveVals.ncol(); // space dimension
  double* tmpPoint = new double[dim];
  for (int i = 0; i < curves.length(); i++){
    NumericMatrix rcppCurveVals = as<List>(curves(i))["coords"];
    lengths[i] = 0;
    for (int j = 0; j < rcppCurveVals.nrow() - 1; j++){
      for (int k = 0; k < dim; k++){
        tmpPoint[k] = rcppCurveVals(j + 1, k) - rcppCurveVals(j, k);
      }
      lengths[i] += norm2(tmpPoint, dim);
    }
    lenTotal += lengths[i];
    // Update maximal length
    if (maxNCurvePoints < rcppCurveVals.nrow()){
      maxNCurvePoints = rcppCurveVals.nrow();
    }
  }
  // 1.2. Calculate number of points for each curve
  int* nsPoints = new int[curves.length()];
  if (ptsPerCurve.length() == curves.length()){ // if given, ...
    for (int i = 0; i < curves.length(); i++){
      nsPoints[i] = ptsPerCurve(i); // just copy
    }
  }else{ // if not given, then ...
    // first entry of 'ptsPerCurve' is the average number of points per curve
    for (int i = 0; i < curves.length(); i++){
      //nsPoints[i] = ptsPerCurve(0) * lengths[i] /
      //  (lenTotal / curves.length());
      nsPoints[i] = ptsPerCurve(0);
    }
  }
  // 1.3. Draw points on curves
  double* cumCurLengths = new double[maxNCurvePoints];
  for (int i = 0; i < curves.length(); i++){ // for each input
    // 1.3.1. Draw and sort the points
    Environment stats_env("package:stats");
    Environment base_env("package:base");
    Function stats_runif = stats_env["runif"];
    Function base_sort = base_env["sort"];
    NumericVector rcppPoints =
      base_sort(stats_runif(nsPoints[i], 0, lengths[i]));
    // 1.3.2. Extract the curve's points
    NumericMatrix rcppCurveVals = as<List>(curves(i))["coords"];
    int n = rcppCurveVals.nrow();
    // 1.3.3. Create structure for the returned points (and other info)
    List curCurve = List();
    NumericMatrix curPoints = NumericMatrix(nsPoints[i], dim);
    // 1.3.4. "Project" points onto the curve
    // 1.3.4.1. Calculate cumulative piecewise lengths
    cumCurLengths[0] = 0;
    for (int j = 1; j < n; j++){
      for (int k = 0; k < dim; k++){
        tmpPoint[k] = rcppCurveVals(j, k) - rcppCurveVals(j - 1, k);
      }
      cumCurLengths[j] = cumCurLengths[j - 1] + norm2(tmpPoint, dim);
    }
    // 1.3.4.2. Interpolate
    double curLength = 0;
    int k = 1; // index for the given curve points
    for (int j = 0; j < nsPoints[i]; j++){ // loop through drawn points
      while (k < n && cumCurLengths[k] < rcppPoints(j)){
        k++; // go to the (right)closest given point
      }
      // The interpolation step
      for (int l = 0; l < dim; l++){
        curPoints(j, l) = rcppCurveVals(k - 1, l) +
          (rcppPoints(j) - cumCurLengths[k - 1]) /
          (cumCurLengths[k] - cumCurLengths[k - 1]) *
          (rcppCurveVals(k, l) - rcppCurveVals(k - 1, l));
        if (j > 0){
          tmpPoint[l] = curPoints(j, l) - curPoints(j - 1, l);
        }
      }
      // Add to real length
      curLength += norm2(tmpPoint, dim);
    }
    // 1.3.5. Add the created structure to the return value
    curCurve.push_back(curPoints, "coords");
    curCurve.push_back(lengths[i], "length.init");
    curCurve.push_back(curLength, "length");
    res.push_back(curCurve);
  }
  // Release memory and return
  delete[] nsPoints;
  delete[] lengths;
  delete[] tmpPoint;
  delete[] cumCurLengths;
  return res;
}

//' Convert images to curves
//'
//' Converts images to curves with points sorted in traversing order.
//'
//' @param images A 3-dimensional array with each slice (matrix in first two
//' dimensions) corresponding to an image. Each (eps-strictly) positive
//' entry is regarded as an occupied pixel (one), otherwise it is regarded
//' as an empty pixel, of an image.
//'
//' @return A list of curves where each element is a function being a list
//' containing a matrix \code{coords} (curve's values, d columns).
//'
//' @references Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//' Depth for curve data and applications.
//'
//' @examples
//' library(curveDepth)
//' # Pixel-grid filling function for an image
//' plotGridImage <- function(dims){
//'   redDims1 <- dims[1] - 1
//'   redDims2 <- dims[2] - 1
//'   for (i in 1:(dims[1] - 1)){
//'     lines(c(i / redDims1 - 0.5 / redDims1,
//'             i / redDims1 - 0.5 / redDims1),
//'             c(0 - 0.5 / redDims2, 1 + 0.5 / redDims2),
//'             lwd = 1, lty = 3, col = "lightgray")
//'     lines(c(0 - 0.5 / redDims1, 1 + 0.5 / redDims1),
//'           c(i / redDims2 - 0.5 / redDims2,
//'             i / redDims2 - 0.5 / redDims2),
//'             lwd = 1, lty = 3, col = "lightgray")
//'   }
//'   rect(0 - 0.5 / redDims1, 0 - 0.5 / redDims2,
//'        1 + 0.5 / redDims1, 1 + 0.5 / redDims2)
//' }
//' # Pixel-grid filling function for a curve
//' plotGridCurve <- function(dims){
//'   for (i in 1:(dims[1] - 1)){
//'     lines(c(i / dims[1], i / dims[1]), c(0, 1),
//'           lwd = 1, lty = 3, col = "lightgray")
//'     lines(c(0, 1), c(i / dims[2], i / dims[2]),
//'           lwd = 1, lty = 3, col = "lightgray")
//'   }
//'   rect(0, 0, 1, 1)
//' }
//' # Load a digit and plot it
//' data("mnistShort017")
//'   aSevenDigit <- mnistShort017$`7`[, , 5]
//' image(as.matrix(rev(as.data.frame(aSevenDigit))),
//'       col = gray((255:0) / 256), asp = 1,
//'       xlim = c(0 - 1 / 27, 1 + 1 / 27),
//'       ylim = c(0 - 1 / 27, 1 + 1 / 27))
//' plotGridImage(dim(aSevenDigit)[1:2])
//' # Convert the digit to a curve and plot it
//' aSevenCurve <- images2curves(array(
//'   aSevenDigit, dim = c(28, 28, 1)))[[1]]
//' plot(cbind(aSevenCurve$coords[, 1],
//'            1 - aSevenCurve$coords[, 2]),
//'            type = "l", lwd = 3, asp = 1,
//'            xlim = c(0, 1), ylim = c(0, 1),
//'            xlab = "x", ylab = "y")
//'   plotGridCurve(dim(aSevenDigit)[1:2])
// [[Rcpp::export]]
List images2curves(arma::cube images){
  int d = 2; // images' dimension
  int* dims = new int[d];
  List res = List(); // the output structure
  // Go through all images, convert them, and output
  for (int i = 0; i < images.n_slices; i++){
    // Read input image
    dims[0] = images.slice(i).n_rows;
    dims[1] = images.slice(i).n_cols;
    ImageDensity imageDens(d, dims, images.slice(i).memptr());
    // Convert it to an ordered curve
    Curve* curve = imageDens.toOrderedCurve();
    // Extract point coordinates
    NumericMatrix coords = NumericMatrix(curve->n, d);
    for (int j = 0; j < curve->n; j++){
      for (int k = 0; k < d; k++){
        coords(j, k) = curve->vals[j][k];
      }
    }
    // Add to the output
    List curRcppCurve = List();
    curRcppCurve.push_back(coords, "coords");
    res.push_back(curRcppCurve);
    // Release memory
    delete curve;
  }
  return res;
}

//' Calculate Tukey curve depth for curves
//'
//' Calculates Tukey curve depth of each curve in \code{objects} w.r.t. the
//' sample of curves in \code{data}. First, \code{m} points are sampled from a
//' uniform distribution on a piecewise linear approximation of each of the
//' curves in \code{data} and \code{m / fracEst * (fracInt + fracEst)} points
//' on each of the curves in \code{objects}. Second, these samples are used to
//' calculate the Tukey curve depth.
//'
//' Calculation of partial depth of each single point can be either exact
//' or approximate. If exact, an extension of the method of
//' Dyckerhoff and Mozharovskyi (2016) is used; if approximate, approximation
//' is performed by projections on directions - points uniformly distributed on
//' the unit hypersphere.
//'
//' @param objects A list where each element is a multivariate curve being a
//' list containing a matrix \code{coords} (values, d columns).
//'
//' @param data A list where each element is a multivariate curve being a list
//' containing a matrix \code{coords} (values, d columns). The depths are
//' computed w.r.t. this data set.
//'
//' @param nDirs Number of directions used to inspect the space, drawn from the
//' uniform distribution on the sphere.
//'
//' @param subs Whether to split each object into two disjunctive subsets (one
//' for integrating and one for estimation) when computing the depth.
//'
//' @param m Number of points used for estimation.
//'
//' @param fracInt Portion of an object used for integrating.
//'
//' @param fracEst Portion of an object used for estimation,
//' maximum: \code{1 - fracInt}.
//'
//' @param exactEst Is calculation of depth for each reference point of the
//' curve exact (\code{TRUE}, by default) or approximate (\code{FALSE}).
//'
//' @param minMassObj Minimal portion of the \code{objects} distribution in the
//' halfspace to be considered when calculating depth.
//'
//' @param minMassDat minimal portion of the \code{data} distribution in the
//' halfspace to be considered when calculating depth.
//'
//' @return A vector of doubles having the same length as \code{objects}, whose
//' each entry is the depth of each element of \code{objects} w.r.t.
//' \code{data}.
//'
//' @references Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//' Depth for curve data and applications.
//'
//' Dyckerhoff, R. and Mozharovskyi P. (2016). Exact computation of the
//' halfspace depth. \emph{Computational Statistics and Data Analysis}, 98,
//' 19-30.
//'
//' @examples
//' library(curveDepth)
//' # Load digits and transform them to curves
//' data("mnistShort017")
//' n <- 10 # cardinality of each class
//' m <- 50 # number of points to sample
//' cst <- 1/10 # a threshold constant
//' alp <- 1/8 # a threshold constant
//' curves0 <- images2curves(mnistShort017$`0`[, , 1:n])
//' curves1 <- images2curves(mnistShort017$`1`[, , 1:n])
//' # Calculate depths
//' depthSpace = matrix(NA, nrow = n * 2, ncol = 2)
//' set.seed(1)
//' depthSpace[, 1] = depthc.Tukey(
//'   c(curves0, curves1), curves0, m = m,
//'   exactEst = TRUE, minMassObj = cst/m^alp)
//' depthSpace[, 2] = depthc.Tukey(
//'   c(curves0, curves1), curves1, m = m,
//'   exactEst = TRUE, minMassObj = cst/m^alp)
//' # Draw the DD-plot
//' plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
//'      xlab = paste("Depth w.r.t. '0'"),
//'      ylab = paste("Depth w.r.t. '1'"),
//'      main = paste("DD-plot for '0' vs '1'"))
//' grid()
//' # Draw the separating rule
//' dat1 <- data.frame(cbind(
//'   depthSpace, c(rep(0, n), rep(1, n))))
//' ddalpha1 <- ddalpha.train(X3 ~ X1 + X2, data = dat1,
//'                           depth = "ddplot",
//'                           separator = "alpha")
//' ddnormal <- ddalpha1$classifiers[[1]]$hyperplane[2:3]
//' pts <- matrix(c(0, 0, 1, ddnormal[1] / -ddnormal[2]),
//'               nrow = 2, byrow = TRUE)
//' lines(pts, lwd = 2)
//' # Draw the points
//' points(depthSpace[1:n, ],
//'        col = "red", lwd = 2, pch = 1)
//' points(depthSpace[(n + 1):(2 * n), ],
//'        col = "blue", lwd = 2, pch = 3)
// [[Rcpp::export(depthc.Tukey)]]
NumericVector depthCTukey(List objects, List data, int nDirs = 100,
                          bool subs = true, int m = 500,
                          double fracInt = 0.5, double fracEst = 0.5,
                          bool exactEst = true,
                          double minMassObj = 0, double minMassDat = 0){
  List tmpObj = curvesSubsample(objects, IntegerVector(1, m / fracEst *
    (fracInt + fracEst)));
  List tmpDat = curvesSubsample(data, IntegerVector(1, m));
  return depthCurveTukey(tmpObj, tmpDat, nDirs, subs, fracInt, fracEst,
                         exactEst, minMassObj, minMassDat);
}
