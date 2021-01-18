//*--------------------------------------------------------------------------*//
//  File:               voxelize.cpp
//  Created by:         Pavlo Mozharovskyi
//  First released:     01.11.2018
//
//  Contains R-level Rcpp-interfaced function(s) for transforming parametrized
//  functions into a voxel format.
//
//  Subsequent changes are listed below:
//  01.11.2018 (Pavlo Mozharovskyi): First version released.
//*--------------------------------------------------------------------------*//

#include "stdafx.h"

//' Voxelization of functions
//'
//' Convertes a pice-wise linear parametrized funtion into a discretized
//' voxel representation.
//'
//' @param f A parametrized function as a list containing a vector
//' "args" (arguments), and a matrix "vals" (values, d columns).
//'
//' @param from A vector of d numbers, each giving a starting discretization
//' point for one dimension.
//'
//' @param to A vector of d numbers, each giving a finishing discretization
//' point for one dimension.
//'
//' @param by A vector of d numbers, each giving discretization step for one
//' dimension.
//'
//' @return A list containing two matrices: "voxels" with rows being voxel
//' numbers, and "coords" with rows being coordinates of voxel centers.
//'
//' @references Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//' Depth for curve data and applications.
//'
//' @examples
//' library(curveDepth)
//' # Create some data based on growth curves
//' g1d <- dataf.growth()
//' g3d <- list("")
//' set.seed(1)
//' for (i in 1:length(g1d$dataf)){
//'   g3d[[i]] <- list(
//'     args = g1d$dataf[[1]]$args,
//'     vals = cbind(g1d$dataf[[i]]$vals,
//'                  g1d$dataf[[i]]$vals[length(g1d$dataf[[i]]$vals):1],
//'                  rnorm(length(g1d$dataf[[i]]$vals), sd = 1) +
//'                    rnorm(1, mean = 0, sd = 10)))
//' }
//' # Define voxels' bounds and resolution
//' from <- c(65, 65, -25)
//' to <- c(196, 196, 25)
//' steps <- 100
//' by <- (to - from) / steps
//' # Voxelize all curves
//' fs <- list("")
//' for (i in 1:length(g3d)){
//'   fs[[i]] <- voxelize(g3d[[i]], from, to, by)
//' }
//' \dontrun{
//' # Plot first 10 curves
//' library(rgl)
//' rgl.open()
//' rgl.bg(color = "white")
//' for (i in 1:10){
//'   spheres3d(fs[[i]]$voxels[, 1], fs[[i]]$voxels[, 2], fs[[i]]$voxels[, 3],
//'             col = "red", radius = 0.5)
//' }}
// [[Rcpp::export]]
List voxelize(List f, NumericVector from, NumericVector to, NumericVector by){
  // 0. Read the input via Rcpp
  NumericMatrix rcppFcnVals = f["vals"];
  NumericVector rcppFcnArgs = f["args"];
  std::vector<double> stlFcnVals = Rcpp::as<std::vector<double> >(rcppFcnVals);
  std::vector<double> stlFcnArgs = Rcpp::as<std::vector<double> >(rcppFcnArgs);
  int nPoints = rcppFcnVals.nrow();
  int nFcns = rcppFcnVals.ncol();
  double* rawFVals = stlFcnVals.data();
  double** fVals = new double*[nFcns];
  for (int i = 0; i < nFcns; i++){
    fVals[i] = rawFVals + i * nPoints;
  }
  // 1. Prepare inner temporary structures
  // Calculate number of grid discretes
  int* nDiscretes = new int[nFcns]; // number of discretes on each axis
  int totalNDiscretes = 0; // number of discretes on all axes
  for (int i = 0; i < nFcns; i++){
    nDiscretes[i] = ceil((to[i] - from[i]) / by[i]) + .5;
    totalNDiscretes += nDiscretes[i];
  }
  // Prepare the voxel structure: a matrix with a voxel in each row
  int nAllVoxels = 0; // curren amount of voxels
  int** allVoxels = new int*[totalNDiscretes * timesGrid];
  int* rawVoxels = new int[nFcns * totalNDiscretes * timesGrid];
  for (int i = 0; i < totalNDiscretes * timesGrid; i++){
    allVoxels[i] = rawVoxels + i * nFcns;
  }
  // Prepare binary hypermatrix to check for including repeated voxels
  binaryHypermatrix<unsigned long long> checkVoxels(nFcns, nDiscretes);
  // 2. Namely voxelize
  double* curPoint = new double[nFcns];
  int* curVoxel = new int[nFcns];
  for (int i = 0; i < nFcns; i++){ // for each functional coordinate
    for (int j = 0; j < nPoints - 1; j++){ // go through all line segments
      int jSmaller = fVals[i][j + 1] > fVals[i][j] ? j : j + 1;
      int jLarger = fVals[i][j + 1] > fVals[i][j] ? j + 1 : j;
      int iVoxelStart = (fVals[i][jSmaller] - from[i]) / by[i];
      int iVoxelEnd = (fVals[i][jLarger] - from[i]) / by[i];
      // Go through all joints of two neighboring voxels in the line segment
      for (int k = iVoxelStart + 1; k < iVoxelEnd; k++){
        // Evaluate function on the current hyperplane
        double curPortion = (from[i] + by[i] * k - fVals[i][jSmaller]) /
        (fVals[i][jLarger] - fVals[i][jSmaller]);
        for (int l = 0; l < nFcns; l++){
          curPoint[l] = fVals[l][jSmaller] +
            (fVals[l][jLarger] - fVals[l][jSmaller]) * curPortion;
        }
        // Determine the voxel where function interssects the hyperplane
        for (int l = 0; l < nFcns; l++){
          curVoxel[l] = (curPoint[l] - from[l]) / by[l];
        }
        // Save the two corresponding voxels: current and previous (a joint)
        if (checkVoxels.setIfNotSet(curVoxel)){
          memcpy(rawVoxels + nFcns * nAllVoxels++,
                 curVoxel, nFcns * sizeof(int));
        }
        if (curVoxel[i] > 0){
          curVoxel[i]--;
          if (checkVoxels.setIfNotSet(curVoxel)){
            memcpy(rawVoxels + nFcns * nAllVoxels++,
                   curVoxel, nFcns * sizeof(int));
          }
        }
      }
    }
  }
  // 3. Prepare R-resutls, release memory, and return the results
  // Copy voxel numbers
  IntegerMatrix voxels(nAllVoxels, nFcns);
  for (int i = 0; i < nAllVoxels; i++){
    for (int j = 0; j < nFcns; j++){
      voxels(i,j) = allVoxels[i][j] + 1;
    }
  }
  // Calculate voxel coordinates
  NumericMatrix coords(nAllVoxels, nFcns);
  for (int i = 0; i < nFcns; i++){
    NumericVector voxelCenters(nDiscretes[i]);
    for (int j = 0; j < nDiscretes[i]; j++){
      voxelCenters(j) = j * by[i] + by[i] / 2 + from[i];
    }
    for (int j = 0; j < nAllVoxels; j++){
      coords(j,i) = voxelCenters(voxels(j,i) - 1);
    }
  }
  // Release the memory
  delete[] nDiscretes;
  delete[] allVoxels;
  delete[] rawVoxels;
  delete[] curPoint;
  delete[] curVoxel;
  return List::create(_("voxels") = voxels,
                      _("coords") = coords);
}
