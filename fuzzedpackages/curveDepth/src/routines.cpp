//*--------------------------------------------------------------------------*//
//  File:               routines.cpp
//  Created by:         Pavlo Mozharovskyi
//  First released:     01.11.2018
//
//  Contains routines for computing Tukey curve depth and related functions.
//
//  For a description of the algorithms, see:
//  Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
//  Depth for curve data and applications.
//
//  Subsequent changes are listed below:
//  01.11.2018 (Pavlo Mozharovskyi): First version released.
//*--------------------------------------------------------------------------*//

#include "stdafx.h"

int Compare(RecSort &rs1, RecSort &rs2){
  return rs1.value < rs2.value;
}

void Swap(RecSort* rs1, RecSort* rs2){
  RecSort rsTmp = *rs1;
  *rs1 = *rs2;
  *rs2 = rsTmp;
}

int Compare(RecMatrix &rm1, RecMatrix &rm2) {
  return rm1.val < rm2.val;
}

void Swap(RecMatrix* rm1, RecMatrix* rm2) {
  RecMatrix rmTmp = *rm1;
  *rm1 = *rm2;
  *rm2 = rmTmp;
}

double norm2(double* x, int d){
  double sum = 0;
  for (int i = 0; i < d; i++){
    sum += pow(x[i], 2);
  }
  return sqrt(sum);
}

// Generates directions uniformly distributed on the unit sphere
// Args:
//   seed: The seed to initialize the random number generator; when 0, current
//     timestamp is used.
//   nDirs: Number of directions to be generated.
//   d: Dimension of the space.
//   dirs (return): Double-indexed array of directions.
// Returns: 0 in case of success.
// Remark: "dirs" shoud be allocated before the call!
int generateDirections(int seed, int nDirs, int d, double** dirs){
  // Prepare random genertor
  //boost::mt19937 rEngine;
  //if (seed != 0){
  //  rEngine.seed(seed);
  //}else{
  //  rEngine.seed(time(0));
  //}
  //boost::normal_distribution<> normDist(0., 1.);
  //boost::variate_generator<boost::mt19937&,
  //                         boost::normal_distribution<> >
  //  varNormDist(rEngine, normDist);
  GetRNGstate();
  for (unsigned int i = 0; i < nDirs; i++){
    double sqrSum = 0;
    for (unsigned int j = 0; j < d; j++){
      //dirs[i][j] = normDist(rEngine);
      dirs[i][j] = norm_rand();
      sqrSum += pow(dirs[i][j], 2);
    }
    sqrSum = sqrt(sqrSum);
    for (unsigned int j = 0; j < d; j++){
      dirs[i][j] = dirs[i][j] / sqrSum;
    }
  }
  PutRNGstate();
  return 0;
}

// Projects curve's voxels onto directions
// Args:
//   curve: A "Curve" with fields "vals", resp. "rawVals", filled.
//   nDirs: Number of directions.
//   dirs: Double-indexed directions.
//   prjs (return): Double-indexed projections.
// Returns: 0 in case of success.
// Remark: "prjs" shoud be allocated before the call!
int projectCurveVoxels(const Curve &curve, int nDirs, double** dirs,
                       double** prjs){
  // std::cout << curve.d << " " << curve.n << std::endl;
  //
  // for (int i = 0; i < curve.n; i++){
  //   for (int j = 0; j < curve.d; j++){
  //     std::cout << curve.vals[i][j] << " ";
  //   }
  // }
  // std::cout << std::endl;
  //
  // for (int i = 0; i < nDirs; i++){
  //   for (int j = 0; j < curve.d; j++){
  //     std::cout << dirs[i][j] << " ";
  //   }
  //   std::cout << std::endl;
  // }

  for (int i = 0; i < nDirs; i++){
    for(int j = 0; j < curve.n; j++){
      prjs[i][j] = 0;
      for (int k = 0; k < curve.d; k++){
        prjs[i][j] += curve.vals[j][k] * dirs[i][k];
      }
    }
  }
  return 0;
}

// Calculates depth of "curve" w.r.t. "curves" based on projections
// Args:
//   curve: "Curve" whose depth has to be calculated.
//   curvePrj: Double-indexed projection of "curve".
//   nDirs: Number of directions.
//   curves: An array of "Curve"s.
//   curvePrjs: An array of double-indexed projections of "curves".
//   n: number of curves in "curves".
// Returns: depth of "curve" w.r.t. "curves".
double calcOneDepth(Curve curve, double** curvePrj, int nDirs,
                    Curve* curves, double*** curvePrjs, int n){
  // Average over the entire curve
  double aveDepth = 0;
  for (int i = 0; i < curve.n; i++){
    // Minimize over directions
    double minDepth = 1000; // The depth value
    for (int j = 0; j < nDirs; j++){
      double fcnPortions = 0;
      // Calculate the cut off portion of each curve ...
      for (int k = 0; k < n; k++){
        fcnPortions += curvePortion(curvePrjs[k][j], curvePrj[j][i],
                                    curves[k].n);
      }
      fcnPortions /= (double)n; // ... and average
      // Calculate the cut off portion of the current function
      double curFcnPortion = curvePortion(curvePrj[j], curvePrj[j][i], curve.n);
      // Calculate and update the depth of the current point
      double tmpDepth = fcnPortions / curFcnPortion;
      minDepth = tmpDepth < minDepth ? tmpDepth : minDepth;
    }
    aveDepth += minDepth; // Augment average depth
  }
  // Return the depth of the single function
  return aveDepth / (double)curve.n;
}

// Calculates portion of the curve on projection
// Args:
//   curvePrj: Projection of the curve.
//   pointPrj: Projection of a single voxel.
//   nVoxels: Number of voxels.
// Returns: Portion of the curve in the positive halfspace cut off by
//   "pointPrj".
double curvePortion(double* curvePrj, double pointPrj, int nVoxels){
  // Sum positive projections
  int posProb = 0;
  for (int i = 0; i < nVoxels; i++){
    if (curvePrj[i] >= pointPrj){
      posProb++;
    }
  }
  return posProb / (double)nVoxels;
}

// Calculates portion of the image density on projection
// Args:
//   imagePrj: Projection of the image.
//   imageDns: Density of the image's voxels.
//   pointPrj: Projection of a single voxel.
//   nVoxels: Number of voxels.
// Returns: Portion of the curve's density in the positive halfspace cut off by
//   "pointPrj".
double imagePortion(double* imagePrj, double* imageDns, double pointPrj,
                    int nVoxels){
  // Sum positive and all projections
  double posProb = 0;
  double totProb = 0;
  for (int i = 0; i < nVoxels; i++){
    if (imageDns[i] > eps){
      totProb += imageDns[i];
      if (imagePrj[i] >= pointPrj){
        posProb += imageDns[i];
      }
    }
  }
  // std::cout << "posProb = " << posProb << std::endl;
  // std::cout << "totProb = " << totProb << std::endl;
  if (fabs(totProb) < eps){
    return 0;
  }
  return posProb / totProb;
}

// Calculates depth of "image" w.r.t. "image" based on projections
// Args:
//   object: "Image" whose depth has to be calculated.
//   objectPrj: Double-indexed projection of "object".
//   nDirs: Number of directions.
//   data: "Image" w.r.t. whom depth has to be calcualted.
//   dataPrj: An array of double-indexed projections of "data".
//   sbSmpl: A flag indicating whether to subsample each object's points when
//     estimating its depth.
// Returns: depth of "object" w.r.t. "data".
double calcOneDepth(ImageDensity &object, double** objectPrj, int nDirs,
                    ImageDensity &data, double** dataPrj, bool sbSmpl){
  // Average over the entire image
  double aveDepth = 0;
  double totDensity = 0;
  // std::cout << object.n << std::endl;
  // std::cout << "Starting to calculate the depth of a single curve." << std::endl;
  for (int i = 0; i < object.n; i++){
    if (sbSmpl){ // If we are subsampling - calculate depth for negative voxels
      if (object.body[i] > eps){
        // std::cout << "object.body[i] < eps == true" << std::endl;
        continue;
      }
    }else{
      // If voxel is empty - do nothing
      if (object.body[i] < eps){
        // std::cout << "object.body[i] < eps == true" << std::endl;
        continue;
      }
    }
    // Minimize over directions
    double minDepth = 1000; // the depth value
    for (int j = 0; j < nDirs; j++){
      // Calculate and update the depth of the current voxel
      double portion1 = imagePortion(dataPrj[j], data.body, objectPrj[j][i],
                                     data.n);
      // std::cout << "portion1 = " << portion1 << std::endl;
      double portion2 = imagePortion(objectPrj[j], object.body, objectPrj[j][i],
                                     object.n);
      // std::cout << "portion2 = " << portion2 << std::endl;
      double tmpDepth = 1002;
      if (fabs(portion2) < eps){
        tmpDepth = 1001;
      }else{
        tmpDepth = portion1 / portion2;
      }
      // std::cout << "tmpDepth = " << tmpDepth << std::endl;
      minDepth = tmpDepth < minDepth ? tmpDepth : minDepth;
    }
    // std::cout << "minDepth = " << minDepth << std::endl;
    // Augment average depth weighting by current voxel
    aveDepth += minDepth * object.body[i];
    totDensity += object.body[i];
  }
  // Return the depth of the single image
  // std::cout << aveDepth << "\t" << totDensity << std::endl;
  if (aveDepth == 0){
    return 0;
  }
  return aveDepth / totDensity;
}

// UNUSED FUNCTION!
// Calculates portion of the curve on projection
// Args:
//   curve: "Curve" with fields "vals", resp. "rawVals", filled.
//   point: Point by which the halfspace is cut off.
//   direction: Direction determining the positive halfspace.
// Returns: Portion of "curve" lying in the positive halfpace on "direction"
//   cut off by "point".
double curvePortion(Curve &curve, double* point, double* direction){
  // Project the point
  double pointPrj = 0;
  for (int i = 0; i < curve.d; i++){
    pointPrj += point[i] * direction[i];
  }
  // Iteratively project the function and sum positive iterations
  int posProb = 0;
  for (int i = 0; i < curve.n; i++){
    double curVoxelPrj = 0;
    for (int j = 0; j < curve.d; j++){
      curVoxelPrj += curve.vals[i][j] * direction[j];
    }
    // Add if the voxel center lies in the positive halfspace
    if (curVoxelPrj >= pointPrj){
      posProb++;
    }
  }
  // Return the portion of voxels in the positive halfspace
  return posProb / (double)curve.n;
}

// Converts an object of type "ImageDensity" into an object of type "Curve"
// Returns: an object of type "Curve" storing only voxels with nonzero values.
Curve* ImageDensity::toCurve(){
  // Create a buffer of exceeding size for voxels with nonzero values
  double* tmpRawVals = new double[n * d];
  int numVoxels = 0;
  // Copy nonzero voxels' coordinates into this buffer:
  // - the array of voxels' numbers
  int* counters = new int[d];
  for(int i = 0; i < d - 1; i++){
    counters[i] = 0;
  }
  counters[d - 1] = -1;
  // - go through all cells
  for (int i = 0; i < n; i++){
    // Choose current cell
    counters[d - 1]++;
    for (int j = d - 1; j > 0; j--){
      if (counters[j] == ns[j]){
        counters[j] = 0;
        counters[j - 1]++;
      }else{
        break;
      }
    }
    // for (int j = 0; j < d; j++){
    //   std::cout << counters[j] << " ";
    // }
    // std::cout << std::endl;
    // Copy the voxel if it is nozero
    if (body[i] > eps){
      for (int j = 0; j < d; j++){
        tmpRawVals[numVoxels * d + j] = vals[i][j];
      }
      numVoxels++;
      // std::cout << vals[i][j] << " ";
    }
  }
  Curve* curve = new Curve(d, numVoxels, tmpRawVals);
  // Release memory
  delete[] tmpRawVals;
  delete[] counters;
  // Return
  return curve;
}

// Calculates Hausdorff distance to another curve
// Args:
//   curve: an object of type "Curve".
// Returns: Hausdorff distance regarding a curve as a set of voxels.
double Curve::distHausdorff(const Curve &curve){
  // TODO: Check that dimensions coincide
  // Create the voxels' similarity matrix
  int nThis = n;
  int nThat = curve.n;
  double* rawSimilarity = new double[nThis * nThat];
  double** similarity = new double*[nThis]; // double-indexing
  for (int i = 0; i < nThis; i++){
    similarity[i] = rawSimilarity + i * nThat;
  }
  // Fill the similarity matrix, i.e. for each pair of pixels
  for (int i = 0; i < nThis; i++){
    for (int j = 0; j < nThat; j++){
      // Calculate the Euclidean distance between each pair of pixels
      similarity[i][j] = 0;
      for (int k = 0; k < d; k++){
        similarity[i][j] += pow(vals[i][k] - curve.vals[j][k], 2);
      }
      similarity[i][j] = sqrt(similarity[i][j]);
    }
  }
  // Calculate minimum for each row
  double* rowMins = new double[nThis];
  for (int i = 0; i < nThis; i++){
    rowMins[i] = DBL_MAX;
    for (int j = 0; j < nThat; j++){
      if (similarity[i][j] < rowMins[i]){
        rowMins[i] = similarity[i][j];
      }
    }
  }
  // Calculate minimum for each column
  double* colMins = new double[nThat];
  for (int j = 0; j < nThat; j++){
    colMins[j] = DBL_MAX;
    for (int i = 0; i < nThis; i++){
      if (similarity[i][j] < colMins[j]){
        colMins[j] = similarity[i][j];
      }
    }
  }
  // Calculate Hausdorff distance as maximum of minimums over rows and columns
  double dist = 0;
  for (int i = 0; i < nThis; i++){
    if (dist < rowMins[i]){
      dist = rowMins[i];
    }
  }
  for (int j = 0; j < nThat; j++){
    if (dist < colMins[j]){
      dist = colMins[j];
    }
  }
  // Release the memory
  delete[] similarity;
  delete[] rawSimilarity;
  delete[] rowMins;
  delete[] colMins;
  // Return
  return dist;
}

// Converts an object of type "ImageDensity" into an object of type "Curve"
// Returns: an object of type "Curve" storing only voxels with nonzero values
// ordered to be directly passed.
Curve* ImageDensity::toOrderedCurve(){
  //Rcout << "ns = " << ns[0] << " " << ns[1] << std::endl;
  // Find a non-zero cell:
  int rawICell = 0;
  // - go through all cells
  for (int i = 0; i < n; i++){
    // - stop if a non-zero cell has been found
    if (body[i] > eps){
      rawICell = i;
      break;
    }
  }
  // Double-indexing
  int* curCell = new int[2];
  curCell[1] = rawICell / ns[0];
  curCell[0] = rawICell % ns[0];
  //Rcout << "Starting point: " << curCell[0] << " " << curCell[1] << std::endl;
  // A matrix for storing seen entries
  binaryHypermatrix<unsigned long long> binImage(2, ns);
  // Prepare first- and second order cells' indices
  int iFirstOrderX[] = { 0, -1, 0, 1 }; // direct neighbors in X
  int iFirstOrderY[] = { -1, 0, 1, 0 }; // direct neighbors in Y
  int iSecondOrderX[] = { -1, -1, 1, 1 }; // diagonal neighbors in X
  int iSecondOrderY[] = { -1, 1, 1, -1 }; // diagonal neighbors in Y
  // Find the curve's fringe (beginning)
  bool beginFound = false;
  binImage.setIfNotSet(curCell);
  int* tmpCell = new int[2];
  while(!beginFound){
    bool found = false;
    // Check first-order neighbors
    for (int i = 0; i < 4; i++){
      // Position the 'tmpCell'
      tmpCell[0] = curCell[0] + iFirstOrderX[i];
      tmpCell[1] = curCell[1] + iFirstOrderY[i];
      //Rcout << "Tmp cell: " << curCell[0] << " " << curCell[1] << std::endl;
      // Check whether it belongs to image
      if (tmpCell[0] >= 0 && tmpCell[0] < ns[0] &&
          tmpCell[1] >= 0 && tmpCell[1] < ns[1]){
        //Rcout << "Cell's ontent: " << body[tmpCell[1] * ns[0] + tmpCell[0]] << std::endl;
        // Check whether the cell is filled
        if (body[tmpCell[1] * ns[0] + tmpCell[0]] > eps){
          // Check whether the cell has not been seen before
          bool isNew = binImage.setIfNotSet(tmpCell);
          if (isNew){
            //Rcout << "Cur cell: " << curCell[0] << " " << curCell[1] << std::endl;
            curCell[0] = tmpCell[0];
            curCell[1] = tmpCell[1];
            found = true;
            break;
          }
        }
      }
    }
    if (!found){
      // Check second-order neighbors
      for (int i = 0; i < 4; i++){
        // Position the 'tmpCell'
        tmpCell[0] = curCell[0] + iSecondOrderX[i];
        tmpCell[1] = curCell[1] + iSecondOrderY[i];
        // Check whether it belongs to image
        if (tmpCell[0] >= 0 && tmpCell[0] < ns[0] &&
            tmpCell[1] >= 0 && tmpCell[1] < ns[1]){
          // Check whether the cell is filled
          if (body[tmpCell[1] * ns[0] + tmpCell[0]] > eps){
            // Check whether the cell has not been seen before
            bool isNew = binImage.setIfNotSet(tmpCell);
            if (isNew){
              curCell[0] = tmpCell[0];
              curCell[1] = tmpCell[1];
              found = true;
              break;
            }
          }
        }
      }
    }
    if (!found){
      beginFound = true;
    }
  }
  // Construct the path of the curve accounting for points' order
  binImage.clear();
  int nPixels = getNNonzero();
  double* rawVals = new double[nPixels * this->d];
  // Traverse all the pixels starting in the beginning
  nPixels = 0;
  bool endFound = false;
  binImage.setIfNotSet(curCell);
  rawVals[nPixels * 2] = ((double)curCell[0] + 0.5) / (double)this->ns[0];
  rawVals[nPixels * 2 + 1] = ((double)curCell[1] + 0.5) / (double)this->ns[1];
  nPixels++;
  while(!endFound){
    bool found = false;
    // Check first-order neighbors
    for (int i = 0; i < 4; i++){
      // Position the 'tmpCell'
      tmpCell[0] = curCell[0] + iFirstOrderX[i];
      tmpCell[1] = curCell[1] + iFirstOrderY[i];
      // Check whether it belongs to image
      if (tmpCell[0] >= 0 && tmpCell[0] < ns[0] &&
          tmpCell[1] >= 0 && tmpCell[1] < ns[1]){
        // Check whether the cell is filled
        if (body[tmpCell[1] * ns[0] + tmpCell[0]] > eps){
          // Check whether the cell has not been seen before
          bool isNew = binImage.setIfNotSet(tmpCell);
          if (isNew){
            curCell[0] = tmpCell[0];
            curCell[1] = tmpCell[1];
            rawVals[nPixels * 2] = ((double)curCell[0] + 0.5) / (double)ns[0];
            rawVals[nPixels * 2 + 1] = ((double)curCell[1] + 0.5) / (double)ns[1];
            found = true;
            break;
          }
        }
      }
    }
    if (!found){
      // Check second-order neighbors
      for (int i = 0; i < 4; i++){
        // Position the 'tmpCell'
        tmpCell[0] = curCell[0] + iSecondOrderX[i];
        tmpCell[1] = curCell[1] + iSecondOrderY[i];
        // Check whether it belongs to image
        if (tmpCell[0] >= 0 && tmpCell[0] < ns[0] &&
            tmpCell[1] >= 0 && tmpCell[1] < ns[1]){
          // Check whether the cell is filled
          if (body[tmpCell[1] * ns[0] + tmpCell[0]] > eps){
            // Check whether the cell has not been seen before
            bool isNew = binImage.setIfNotSet(tmpCell);
            if (isNew){
              curCell[0] = tmpCell[0];
              curCell[1] = tmpCell[1];
              rawVals[nPixels * 2] = ((double)curCell[0] + 0.5) / (double)ns[0];
              rawVals[nPixels * 2 + 1] = ((double)curCell[1] + 0.5) / (double)ns[1];
              found = true;
              break;
            }
          }
        }
      }
    }
    if (!found){
      endFound = true;
    }else{
      nPixels++;
    }
  }
  //for (int i = 0; i < nPixels; i++){
  //  Rcout << round(rawVals[i * 2 + 0] * ns[0] - 0.5) << " " << round(rawVals[i * 2 + 1] * ns[1] - 0.5) << std::endl;
  //}
  Curve* curve = new Curve(2, nPixels, rawVals);
  if (fmax(fabs(round(rawVals[(0) * 2 + 0] * ns[0] - 0.5) -
      round(rawVals[(nPixels - 1) * 2 + 0] * ns[0] - 0.5)),
      fabs(round(rawVals[(0) * 2 + 1] * ns[1] - 0.5) -
        round(rawVals[(nPixels - 1) * 2 + 1] * ns[1] - 0.5))) < 1.5){
    curve->closed = true;
    //Rcout << "Closed curve detected." << std::endl;
  }
  // Release memory
  delete[] curCell;
  delete[] tmpCell;
  delete[] rawVals;
  return curve;
}

// Calculates distance from curve 1 to curve 2 in the curve space.
// 1 and 2, starting points, and order of processing points are kept constant.
// For this a path in the distance matrix is created with the smallest maximum
// value on it.
// Args:
//   n1: number of voxels in the first curve.
//   n2: number of voxels in the second curve.
//   dists: matrix of distances with n1 rows and n2 columns.
// Returns: (minimized) maximum value on the path through matrix.
double getMinMaxDist(int n1, int n2, double* dists){
  // Create structure describing the state of the matrix:
  // - borders of the filled space for rows and columns and thier number
  int* minIRows = new int[n1];
  int* maxIRows = new int[n1];
  int* nIRows = new int[n1];
  for (int i = 0; i < n1; i++){
    minIRows[i] = -1;
    maxIRows[i] = n2;
    nIRows[i] = 0;
  }
  int* minICols = new int[n2];
  int* maxICols = new int[n2];
  int* nICols = new int[n2];
  for (int i = 0; i < n2; i++){
    minICols[i] = -1;
    maxICols[i] = n1;
    nICols[i] = 0;
  }
  // - the indicator matrix for filled entries
  // g++ compiler bug
  int tmpArray1[] = {n1, n2};
  //binaryHypermatrix<unsigned long long> fulfillment(2, (int[]) { n1, n2 });
  binaryHypermatrix<unsigned long long> fulfillment(2, tmpArray1);
  fulfillment.clear();
  // Sort the elements (=distances) by their size
  RecMatrix* recDists = new RecMatrix[n1 * n2];
  for (int i = 0; i < n1; i++){
    for (int j = 0; j < n2; j++){
      recDists[i * n2 + j].val = dists[i * n2 + j];
      recDists[i * n2 + j].iRow = i;
      recDists[i * n2 + j].iCol = j;
      //Rcout << std::setprecision(3) << recDists[i * n2 + j].val << " ";
    }
    //Rcout << std::endl;
  }
  //Rcout << std::endl;
  quick_sort(recDists, 0, n1 * n2 - 1);
  //for (int i = n1 * n2 - 1; i >= 0; i--){
  //  Rcout << i << ": " << recDists[i].iRow << " " << recDists[i].iCol << " " <<
  //    recDists[i].val << " " << std::endl;
  //}
  // Main loop - determine the minimax
  double minMax = -1;
  // Print the matrix: begin
  //for (int i = 0; i < n1; i++){
  //  for (int j = 0; j < n2; j++){
  //    Rcout << (int)(fulfillment.isSet((int[]){i, j}));
  //  }
  //  Rcout << std::endl;
  //}
  // Print the matrix: end
  //Rcout << "Point 1 reached.\n" << std::endl;
  for (int i = n1 * n2 - 1; i >= 0; i--){
    // Remove the element from the matrix
    //Rcout << std::endl << recDists[i].iRow << " " << recDists[i].iCol << " " << recDists[i].val << std::endl;
    // g++ compiler bug (see most top comment)
    tmpArray1[0] = recDists[i].iRow;
    tmpArray1[1] = recDists[i].iCol;
    //if (fulfillment.setIfNotSet((int[]) { recDists[i].iRow,
    //                            recDists[i].iCol })){
    if (fulfillment.setIfNotSet(tmpArray1)){
      nIRows[recDists[i].iRow]++;
      nICols[recDists[i].iCol]++;
    }
    // Recalculate the row and column mins
    if (minIRows[recDists[i].iRow] + 1 == recDists[i].iCol){
      // Fill the gaps in the row from the left
      // g++ compiler bug (see most top comment)
      tmpArray1[0] = recDists[i].iRow;
      tmpArray1[1] = minIRows[recDists[i].iRow] + 1;
      // while(minIRows[recDists[i].iRow] + 1 < n2 &&
      //       fulfillment.isSet((int[]) { recDists[i].iRow,
      //                         minIRows[recDists[i].iRow] + 1 })){
      while(minIRows[recDists[i].iRow] + 1 < n2 &&
            fulfillment.isSet(tmpArray1)){
        minIRows[recDists[i].iRow]++;
        tmpArray1[0] = recDists[i].iRow;
        tmpArray1[1] = minIRows[recDists[i].iRow] + 1;
      }
    }
    if (minICols[recDists[i].iCol] + 1 == recDists[i].iRow){
      // Fill the gaps in the column from above
      // g++ compiler bug (see most top comment)
      tmpArray1[0] = minICols[recDists[i].iCol] + 1;
      tmpArray1[1] = recDists[i].iCol;
      // while(minICols[recDists[i].iCol] + 1 < n1 &&
      //       fulfillment.isSet((int[]) { minICols[recDists[i].iCol] + 1,
      //                         recDists[i].iCol })){
      while(minICols[recDists[i].iCol] + 1 < n1 &&
              fulfillment.isSet(tmpArray1)){
        minICols[recDists[i].iCol]++;
        tmpArray1[0] = minICols[recDists[i].iCol] + 1;
        tmpArray1[1] = recDists[i].iCol;
      }
    }
    // Refill the matrix due to route constraints
    if (minIRows[recDists[i].iRow] >= 0){
      // Refill downwards from a row
      for (int j = 0; j <= minIRows[recDists[i].iRow]; j++){
        for (int k = recDists[i].iRow; k < n1; k++){
          // Fill the column from above
          // g++ compiler bug (see most top comment)
          tmpArray1[0] = k;
          tmpArray1[1] = j;
          //if(fulfillment.setIfNotSet((int[]) { k, j })){
          if(fulfillment.setIfNotSet(tmpArray1)){
            nIRows[k]++;
            nICols[j]++;
          }
        }
      }
    }
    if (minICols[recDists[i].iCol] >= 0){
      // Refill rightwards from a column
      for (int j = 0; j <= minICols[recDists[i].iCol]; j++){
        for (int k = recDists[i].iCol; k < n2; k++){
          // Fill the column from above
          // g++ compiler bug (see most top comment)
          tmpArray1[0] = j;
          tmpArray1[1] = k;
          // if(fulfillment.setIfNotSet((int[]) { j, k })){
          if(fulfillment.setIfNotSet(tmpArray1)){
            nIRows[j]++;
            nICols[k]++;
          }
        }
      }
    }
    // Recalculate the row and column max
    if (maxIRows[recDists[i].iRow] - 1 == recDists[i].iCol){
      // Fill the row from the right
      // Fill the column from above
      // g++ compiler bug (see most top comment)
      tmpArray1[0] = recDists[i].iRow;
      tmpArray1[1] = maxIRows[recDists[i].iRow] - 1;
      // while(maxIRows[recDists[i].iRow] - 1 >= 0 &&
      //       fulfillment.isSet((int[]) { recDists[i].iRow,
      //                         maxIRows[recDists[i].iRow] - 1 })){
      while(maxIRows[recDists[i].iRow] - 1 >= 0 &&
            fulfillment.isSet(tmpArray1)){
        maxIRows[recDists[i].iRow]--;
        tmpArray1[0] = recDists[i].iRow;
        tmpArray1[1] = maxIRows[recDists[i].iRow] - 1;
      }
    }
    if (maxICols[recDists[i].iCol] - 1 == recDists[i].iRow){
      // Fill the column from bottom
      // Fill the row from the right
      // Fill the column from above
      // g++ compiler bug (see most top comment)
      tmpArray1[0] = maxICols[recDists[i].iCol] - 1;
      tmpArray1[1] = recDists[i].iCol;
      // while(maxICols[recDists[i].iCol] - 1 >= 0 &&
      //       fulfillment.isSet((int[]) { maxICols[recDists[i].iCol] - 1,
      //                         recDists[i].iCol })){
      while(maxICols[recDists[i].iCol] - 1 >= 0 &&
            fulfillment.isSet(tmpArray1)){
        maxICols[recDists[i].iCol]--;
        tmpArray1[0] = maxICols[recDists[i].iCol] - 1;
        tmpArray1[1] = recDists[i].iCol;
      }
    }
    // Refill the matrix due to route constraints
    if (maxIRows[recDists[i].iRow] < n2){
      // Refill upwards from a row
      for (int j = maxIRows[recDists[i].iRow]; j < n2; j++){
        for (int k = 0; k <= recDists[i].iRow; k++){
          // g++ compiler bug (see most top comment)
          tmpArray1[0] = k;
          tmpArray1[1] = j;
          // if(fulfillment.setIfNotSet((int[]) { k, j })){
          if(fulfillment.setIfNotSet(tmpArray1)){
            nIRows[k]++;
            nICols[j]++;
          }
        }
      }
    }
    if (maxICols[recDists[i].iCol] < n1){
      // Refill leftwards from a column
      for (int j = maxICols[recDists[i].iCol]; j < n1; j++){
        for (int k = 0; k <= recDists[i].iCol; k++){
          // g++ compiler bug (see most top comment)
          tmpArray1[0] = j;
          tmpArray1[1] = k;
          // if(fulfillment.setIfNotSet((int[]) { j, k })){
          if(fulfillment.setIfNotSet(tmpArray1)){
            nIRows[j]++;
            nICols[k]++;
          }
        }
      }
    }
    //Rcout << "Point 7 reached.\n" << std::endl;
    // Print the matrix: begin
    //for (int i = 0; i < n1; i++){
    //  for (int j = 0; j < n2; j++){
    //    Rcout << (int)(fulfillment.isSet((int[]){i, j}));
    //  }
    //  Rcout << std::endl;
    //}
    // Print the matrix: end
    // Check whether the matrix is blocked
    bool blocked = false;
    //if (fulfillment.isSet((int[]) { 0, 0 })){
      //Rcout << "Legal exit (0,0)" << std::endl;
      //blocked = true;
    //}
    //if (fulfillment.isSet((int[]) { n1 - 1, n2 - 1 })){
      //Rcout << "Legal exit (n1-1,n2-1)" << std::endl;
      //blocked = true;
    //}
    for (int j = 0; j < n1; j++){
      if (nIRows[j] > n2){
        Rcout << "!Alert: nIRows[" << j << "] > " << n2 << std::endl;
      }
      if (nIRows[j] == n2){
        //Rcout << "Legal exit (nIRows " << j << ")" << std::endl;
        blocked = true;
        break;
      }
    }
    for (int j = 0; j < n2; j++){
      if (nICols[j] > n1){
        Rcout << "!Alert: nICols[" << j << "] > " << n1 << std::endl;
      }
      if (nICols[j] == n1){
        //Rcout << "Legal exit (nICols " << j << ")" << std::endl;
        blocked = true;
        break;
      }
    }
    // Return the element if the matrix is blocked
    if (blocked){
      minMax = recDists[i].val;
      break;
    }

    //Rcout << "Point 10 reached.\n" << std::endl;
    // Print the matrix: begin
    //for (int i = 0; i < n1; i++){
    //  for (int j = 0; j < n2; j++){
    //    Rcout << (int)(fulfillment.isSet((int[]){i, j}));
    //  }
    //  Rcout << std::endl;
    //}
    // Print the matrix: end

    //if (i == (n1 * n2 - 10)){
      //Rcout << "Urgent exit" << std::endl;
    //  blocked = true;
    //  break;
    //}

  }
  //Rcout << "Point 11 reached.\n" << std::endl;
  // Print the matrix: begin
  //for (int i = 0; i < n1; i++){
  //  for (int j = 0; j < n2; j++){
  //    Rcout << (int)(fulfillment.isSet((int[]){i, j}));
  //  }
  //  Rcout << std::endl;
  //}
  // Print the matrix: end
  // Release memory
  delete[] minIRows;
  delete[] maxIRows;
  delete[] nIRows;
  delete[] minICols;
  delete[] maxICols;
  delete[] nICols;
  delete[] recDists;
  return minMax;
}

// Calculates distance in the space of curves to another curve
// Args:
//   curve: an object of type "Curve" with ordered voxels.
//   oneWay: whether curve should be traverced in one direction only
// Returns: metric-based distance regarding a curve as an ordered set of voxels.
double Curve::distCurve(const Curve &curve, bool oneWay){
  // If both curves are closed, employ the Hausdorff distance
  if (this->closed && curve.closed){
    return this->distHausdorff(curve);
  }
  // Initialize variables
  int n1 = n;
  int n2 = curve.n;
  double theDist = DBL_MAX;
  // Calculate distances:
  double* dists = new double[n1 * n2];
  double tmpDist = DBL_MAX;
  // - treat the closedness of one of the curves
  int nRotations = 1;
  int shift1 = 0;
  int shift2 = 0;
  nRotations = 1;
  if (this->closed){
    nRotations = n1;
  }
  if (curve.closed){
    nRotations = n2;
  }
  for (int sh = 0; sh < nRotations; sh++){
    // The closedness variables
    shift1 = 0;
    shift2 = 0;
    if (this->closed){
      shift1 = sh;
      shift2 = 0;
    }
    if (curve.closed){
      shift1 = 0;
      shift2 = sh;
    }
    // - this start at begin, the other start at begin
    //Rcout << std::endl;
    for (int i = 0; i < n1; i++){
      for (int j = 0; j < n2; j++){
        dists[i * n2 + j] = 0;
        for (int k = 0; k < d; k++){
          dists[i * n2 + j] += pow(vals[(i + shift1) % n1][k] -
            curve.vals[(j + shift2) % n2][k], 2);
        }
        dists[i * n2 + j] = sqrt(dists[i * n2 + j]);
        //Rcout << dists[i * n2 + j] << " ";
      }
      //Rcout << std::endl;
    }
    tmpDist = getMinMaxDist(n1, n2, dists);
    if (tmpDist < theDist){
      theDist = tmpDist;
    }
    //Rcout << tmpDist << std::endl;
    // If one-way only, quit here
    if (oneWay){
      break;
    }
    // - this start at end, the other start at begin
    for (int i = 0; i < n1; i++){
      for (int j = 0; j < n2; j++){
        dists[i * n2 + j] = 0;
        for (int k = 0; k < d; k++){
          dists[i * n2 + j] += pow(vals[(n1 - 1 - i + shift1) % n1][k] -
            curve.vals[(j + shift2) % n2][k], 2);
        }
        dists[i * n2 + j] = sqrt(dists[i * n2 + j]);
      }
    }
    tmpDist = getMinMaxDist(n1, n2, dists);
    if (tmpDist < theDist){
      theDist = tmpDist;
    }
    // - this start at begin, the other start at end
    for (int i = 0; i < n1; i++){
      for (int j = 0; j < n2; j++){
        dists[i * n2 + j] = 0;
        for (int k = 0; k < d; k++){
          dists[i * n2 + j] += pow(vals[(i + shift1) % n1][k] -
            curve.vals[(n2 - 1 - j + shift2) % n2][k], 2);
        }
        dists[i * n2 + j] = sqrt(dists[i * n2 + j]);
      }
    }
    tmpDist = getMinMaxDist(n1, n2, dists);
    if (tmpDist < theDist){
      theDist = tmpDist;
    }
    // - this start at end, the other start at end
    for (int i = 0; i < n1; i++){
      for (int j = 0; j < n2; j++){
        dists[i * n2 + j] = 0;
        for (int k = 0; k < d; k++){
          dists[i * n2 + j] += pow(vals[(n1 - 1 - i + shift1) % n1][k] -
            curve.vals[(n2 - 1 - j + shift2) % n2][k], 2);
        }
        dists[i * n2 + j] = sqrt(dists[i * n2 + j]);
      }
    }
    tmpDist = getMinMaxDist(n1, n2, dists);
    if (tmpDist < theDist){
      theDist = tmpDist;
    }
  }
  // Release memory
  delete[] dists;
  return theDist;
}

// Calculate depth of empirical measure 'curEmpDist' w.r.t. empirical measure
// 'genEmpDist' based on directions 'dirs'
double calcOneDepth(EmpDist &curEmpDist, EmpDist &genEmpDist, double** dirs,
             int nDirs, int d){
  // 1. Initialize strucrutes
  double* pDepths = new double[curEmpDist.n]; // point-wise depths
  RecSort* proj = new RecSort[curEmpDist.n + genEmpDist.n]; // projection
  double pProj = 0; // projection of the current point of 'curEmpDist'
  // 2. For each point of 'curEmpDist'
  for (int i = 0; i < curEmpDist.n; i++){
    //for (int l = 0; l < d; l++){
    //  Rcout << curEmpDist.vals[i][l] << " ";
    //}
    //Rcout << std::endl;
    //if (i > 4){
    //  break;
    //}
    // 2.a. Initialize point's depth
    pDepths[i] = 1; // ??? Check this
    // 2.b. For each direction
    for (int j = 0; j < nDirs; j++){
      // 2.b.i. Project 'curEmpDist'
      for (int k = 0; k < curEmpDist.n; k++){
        proj[k].value = 0;
        for (int l = 0; l < d; l++){
          proj[k].value += curEmpDist.vals[k][l] * dirs[j][l];
        }
        proj[k].index = 0; // 'curEmpDist'
        proj[k].weight = curEmpDist.weights[k];
      }
      // 2.b.ii. Project 'genEmpDist'
      for (int k = 0; k < genEmpDist.n; k++){
        proj[curEmpDist.n + k].value = 0;
        for (int l = 0; l < d; l++){
          proj[curEmpDist.n + k].value += genEmpDist.vals[k][l] * dirs[j][l];
        }
        proj[curEmpDist.n + k].index = 1; // 'genEmpDist'
        proj[curEmpDist.n + k].weight = genEmpDist.weights[k];
      }
      // 2.b.iii. Sort ptojection
      quick_sort(proj, 0, curEmpDist.n + genEmpDist.n - 1);
      for (int m = 0; m < curEmpDist.n + genEmpDist.n; m++){
        //Rcout << proj[m].value << " ";
        //Rcout << proj[m].weight << "(" << proj[m].index << ")" << " ";
      }
      //Rcout << std::endl;
      // 2.b.iv Project point
      pProj = 0;
      for (int l = 0; l < d; l++){
        pProj += curEmpDist.vals[i][l] * dirs[j][l];
        //Rcout << curEmpDist.vals[i][l] << " ";
      }
      //Rcout << std::endl;
      // 2.b.v. Treat dejenerate cases
      if ((pProj == proj[0].value && proj[0].value < proj[1].value) ||
          (pProj == proj[curEmpDist.n + genEmpDist.n - 1].value &&
          proj[curEmpDist.n + genEmpDist.n - 1].value >
          proj[curEmpDist.n + genEmpDist.n - 2].value)){
        pDepths[i] = 0;
        Rcout << "z";
        break;
      }
      // 2.b.vi. Check all possible halfspaces
      //int curPos = 0; // current position in the projection
      // Set up initial (pro)portions in the CLOSED halfspaces
      double wNomHfspLower = 0;
      double wNomHfspUpper = 1;
      double wDenHfspLower = 0;
      double wDenHfspUpper = 1;
      // Regard each point (accounting for ties) as the border of the halfspace
      for(int curPos = 0; curPos < curEmpDist.n + genEmpDist.n; curPos++){
        // Change the weights in the halfspaces:
        // - add weight to the lower halfspace
        if (proj[curPos].index == 0){
          wDenHfspLower += proj[curPos].weight; // 'curEmpDist'
        }else{
          wNomHfspLower += proj[curPos].weight; // 'genEmpDist'
        }
        // - subtract weight from the upper halfspace
        if (curPos > 0){
          if (proj[curPos - 1].index == 0){
            wDenHfspUpper -= proj[curPos - 1].weight; // 'curEmpDist'
          }else{
            wNomHfspUpper -= proj[curPos - 1].weight; // 'genEmpDist'
          }
        }
        // If coordinate is different w.r.t. the next one ...
        if (curPos < curEmpDist.n + genEmpDist.n - 1 &&
            proj[curPos + 1].value != proj[curPos].value &&
            pProj <= proj[curPos].value){
          // ... check lower halfspace
          if (wNomHfspLower == 0){
            pDepths[i] = 0; // minimum found
            break;
          }else{
            if (wDenHfspLower > 0){
              // non-degenerate case, update using lower halfspace
              double pDepthTmp = wNomHfspLower / wDenHfspLower;
              if (pDepthTmp < pDepths[i]){
                pDepths[i] = pDepthTmp;
                //Rcout << "(" <<
                  //pDepthTmp <<
                  //"," << curPos <<
                  //", " << wNomHfspLower <<
                  //"," << wDenHfspLower
                  //<< ") ";
              }
            }
          }
        }
        // If coordinate is different w.r.t. the previous one ...
        if (curPos > 0 &&
            proj[curPos - 1].value != proj[curPos].value &&
            pProj >= proj[curPos].value){
          // ... check upper halfspace
        }
        if (wNomHfspUpper == 0){
          pDepths[i] = 0; // minimum found
          Rcout << "z";
          break;
        }else{
          if (wDenHfspUpper > 0){
            // non-degenerate case, update using lower halfspace
            double pDepthTmp = wNomHfspUpper / wDenHfspUpper;
            if (pDepthTmp < pDepths[i]){
              //Rcout << "(" <<
                //pDepthTmp <<
                //"," << curPos <<
                  //"," << wNomHfspUpper <<
                    //"," << wDenHfspUpper
              //<< ") ";
              pDepths[i] = pDepthTmp;
            }
          }
        }
      }
      //// If zero depth already achieved on one of the projections, stop searching
      //if (pDepths[i] == 0){
      //  Rcout << "z";
      //  break;
      //}
    }
    // Augment the depth of 'cur'
  }
  // 3. Caclulate depth of 'curEmpDist' as the average over all its points
  double theDepth = 0; // the depth to be calculated by averaging
  double wis = 0;
  Rcout << std::endl;
  for (int i = 0; i < curEmpDist.n; i++){
    Rcout << pDepths[i] << " ";
    theDepth += curEmpDist.weights[i] * pDepths[i];
    wis += curEmpDist.weights[i];
  }
  Rcout << "WS: " << wis << " ";
  // Release memory
  delete[] pDepths;
  delete[] proj;
  return theDepth;
}

// Creates an empty (noise-filled) empirical distribution of a predefined size
EmpDist::EmpDist(int n, int d, bool zeroInit) : Curve() {
  // Create the structures
  this->d = d;
  this->n = n;
  this->rawVals = new double[n * d];
  this->vals = asMatrix(this->rawVals, n, d);
  this->weights = new double[n];
  // Fill with zeros if required
  if (zeroInit){
    for (int i = 0; i < n; i++){
      for (int j = 0; j < d; j++){
        this->vals[i][j] = 0;
      }
      this->weights[i] = 1 / (double)n;
    }
  }
}

// Creates an empirical distribution as a copy of another one
EmpDist::EmpDist(EmpDist &empDist) : Curve() {
  // Create the structures
  d = empDist.d;
  n = empDist.n;
  rawVals = new double[n * d];
  vals = asMatrix(rawVals, n, d);
  weights = new double[n];
  // Fill them by copying
  memcpy(rawVals, empDist.rawVals, n * d * sizeof(double));
  memcpy(weights, empDist.weights, n * sizeof(double));
}

// Creates an empirical distribution from a bunch of curves
EmpDist::EmpDist(Curve* curves, int nCurves, double tiesPrecision) : Curve() {
  // Estimate maximal (can be repetitions) number of points
  int nMax = 0; // this is also total number of points
  for (int i = 0; i < nCurves; i++){
    nMax += curves[i].n;
  }
  // DEBUG
  // Rcout << "nMax:" << nMax << "\n";
  // DEBUG (end)
  // Allocate temporary dimension, points and weights
  int d = curves[0].d;
  double* rawPoints = new double[nMax * d];
  double** points = asMatrix(rawPoints, nMax, d); // new double*
  double* tmpWeights = new double[nMax];
  int nPoints = 0; // current number of unique points
  // Collect unique points from all curves:
  for (int i = 0; i < nCurves; i++){
    for (int j = 0; j < curves[i].n; j++){
      // a) Check uniqueness of the point
      int counter = nPoints - 1 > -1 ? nPoints - 1 : 0;
      while (counter < nPoints){
        bool isRepeating = true;
        // DEBUG
        // Rcout << fabs(curves[i].vals[j][0] - points[counter][0]) << "\n";
        // DEBUG (end)
        for (int k = 0; k < d; k++){
          if (fabs(curves[i].vals[j][k] - points[counter][k]) > tiesPrecision){
            isRepeating = false;
            break;
          }
        }
        if (isRepeating){
          break;
        }else{
          counter++;
        }
      }
      // b) Add the point
      if (counter == nPoints){
        // If the point is unique - create a new point at the end
        memcpy(points[nPoints], curves[i].vals[j], d * sizeof(double));
        tmpWeights[nPoints] = 1 / (double)nMax;
        nPoints++;
      }else{
        // If a repeated point - just increase its weight
        //tmpWeights[counter] += 1 / (double)nMax;
      }
    }
  }
  // DEBUG
  // Rcout << "nPoints:" << nPoints << "\n";
  // DEBUG (end)
  // Copy the created data structures into this new curve
  this->d = d;
  this->n = nPoints;
  this->rawVals = new double[nPoints * d];
  memcpy(this->rawVals, rawPoints, nPoints * d * sizeof(double));
  this->vals = asMatrix(this->rawVals, this->n, this->d);
  this->weights = new double[nPoints];
  memcpy(this->weights, tmpWeights, nPoints * sizeof(double));
  // Release memory
  delete[] rawPoints;
  delete[] points;
  delete[] tmpWeights;
}

// Creates an empirical distribution from a bunch of images
// (Assumes that all images have the same resolution, this of the first one)
EmpDist::EmpDist(ImageDensity* imageDensities, int nImageDensities,
                 double tiesPrecision){
  // If more than one image provided, then fuse them
  int d, *ns; // dimension of the image
  imageDensities[0].getSize(&d, &ns); // new
  ImageDensity genPattern(d, ns); // generalized pattern
  for (int i = 0; i < nImageDensities; i++){
    genPattern.add(imageDensities[i]);
  }
  genPattern.standardize();
  Rcout << genPattern.getTotal() << " ";
  // Allocate temporary dimension, points and weights
  int nMax = 1; // maximumally possible number of pixels
  for (int i = 0; i < d; i++){
    nMax *= ns[i];
  }
  double* rawPoints = new double[nMax * d];
  double** points = asMatrix(rawPoints, nMax, d); // new
  double* weights = new double[nMax];
  int nPoints = 0; // current number of unique points
  // !Operations below repeat in a more complicated way some existing functions!
  // !This is to have comparison (with zero) in one place and to control it!
  double sum = 0;
  int n = genPattern.n;
  for(int i = 0; i < n; i++){
    // If pixel value different from zero, ...
    if (genPattern.body[i] > 0){
      // ... add it to the empirical distribution
      memcpy(points[nPoints], genPattern.vals[i], d * sizeof(double));
      weights[nPoints] = genPattern.body[i];
      nPoints++;
      // ... and sum it up
      sum += genPattern.body[i];
    }
  }
  Rcout << nPoints << " " << sum << " ";
  // Normalize to add up to 1
  for (int i = 0; i < nPoints; i++){
    weights[i] /= sum;
  }
  // Copy the created data structures into this new curve
  this->d = d;
  this->n = nPoints;
  this->rawVals = new double[nPoints * d];
  memcpy(this->rawVals, rawPoints, nPoints * d * sizeof(double));
  this->vals = asMatrix(this->rawVals, this->n, this->d);
  this->weights = new double[nPoints];
  memcpy(this->weights, weights, nPoints * sizeof(double));
  // Release memory
  delete[] rawPoints;
  delete[] points;
  delete[] weights;
}

void EmpDist::updateWeights(bool dropZeros){
  // Eliminate points with zero weights from the list if necessary
  if (dropZeros){
    // TODO: Write weight elimination to speed up the calculation
    // Count the number of non-zero weights
    int nNonZeros = 0;
    for (int i = 0; i < n; i++){
      if (weights[i] > 0){
        nNonZeros++;
      }
    }
    // Create new point and weihgt structures
    double* newRawVals = new double[nNonZeros * d];
    double** newVals = new double*[nNonZeros];
    double* newWeights = new double[nNonZeros];
    // Filter points and weights by removing entries with zero weight
    int counter = 0;
    for(int i = 0; i < n; i++){
      if (weights[i] > 0){
        // Copy point and its weight
        newVals[counter] = newRawVals + counter * d;
        memcpy(newVals[counter], vals[i], d * sizeof(double));
        newWeights[counter] = weights[i];
        // Go on to the next non-zero-weight point
        counter++;
      }
    }
    // Change the content
    n = nNonZeros;
    delete[] rawVals;
    delete[] vals;
    delete[] weights;
    rawVals = newRawVals;
    vals = newVals;
    weights = newWeights;
  }
  // Calculate the (new) sum of the weights
  double total = 0;
  for (int i = 0; i < n; i++){
    total += weights[i];
  }
  // Reweight if necessary
  if (total < 1){
    // Divide each weight through 'total'
    for (int i = 0; i < n; i++){
      weights[i] /= total;
    }
  }
}

// Calculates portion of an empirical distributoin on a projection
// Args:
//   empDistPrj: Projection of an empirical distribution.
//   empDistWeights: Weights of an empirical distribution.
//   pointPrj: Projection of a single voxel.
//   nVoxels: Number of voxels.
// Returns: Portion of the empirical distribution in the positive halfspace cut
//   away by "pointPrj".
double empDistPortion(double* empDistPrj, double* empDistWeights,
                      double pointPrj, int nVoxels){
  // Sum positive projections weighting them
  double posProb = 0;
  for (int i = 0; i < nVoxels; i++){ // for each voxel of the 'empDist':
    if (empDistPrj[i] - pointPrj >= - eps){ // if it is in a positive halfspace
      posProb += empDistWeights[i]; // then add its weight
    }
  }
  return posProb <= 1 ? posProb : 1;
}

// Calculates approximate depth of empirical measure 'curEmpDist' w.r.t.
// empirical measure 'genEmpDist' based on 'refEmpDist' projecting on
// directions 'dirs'
double approxOneDepth(EmpDist &refEmpDist, EmpDist &curEmpDist,
                      EmpDist &genEmpDist, double** dirs,
                      int nDirs, int d, double curMinMass, double genMinMass){
  // Rcout << "RefEmpDistN:" << refEmpDist.n << ", CurEmpDistN:" << curEmpDist.n << ", GenEmpDistN:" << genEmpDist.n << std::endl;
  double depth = 0;// the depth
  // 1. Initialize strucrutes
  double* pDepths = new double[refEmpDist.n]; // point-wise depths
  // 2. Project the empirical distributions on the given directions
  // 2.a) For 'genEmpDist'
  double* rawGenEmpPrjs = new double[nDirs * genEmpDist.n];
  double** genEmpPrjs = asMatrix(rawGenEmpPrjs, nDirs, genEmpDist.n); // new d.*
  projectCurveVoxels(genEmpDist, nDirs, dirs, genEmpPrjs); // do project
  // 2.b) For 'curEmpDist'
  double* rawCurEmpPrjs = new double[nDirs * curEmpDist.n];
  double** curEmpPrjs = asMatrix(rawCurEmpPrjs, nDirs, curEmpDist.n); // new d.*
  projectCurveVoxels(curEmpDist, nDirs, dirs, curEmpPrjs); // do project
  // 2.c) For 'refEmpDist'
  double* rawRefEmpPrjs = new double[nDirs * refEmpDist.n];
  double** refEmpPrjs = asMatrix(rawRefEmpPrjs, nDirs, refEmpDist.n); // new d.*
  projectCurveVoxels(refEmpDist, nDirs, dirs, refEmpPrjs); // do project
  // DEBUG
  // Rcout << "genEmpDist:\n";
  // for (int i = 0; i < nDirs; i++){
  //   for (int j = 0; j < genEmpDist.n; j++){
  //     Rcout << genEmpPrjs[i][j] << " ";
  //   }
  //   Rcout << "\n";
  // }
  // Rcout << "curEmpDist:\n";
  // for (int i = 0; i < nDirs; i++){
  //   for (int j = 0; j < curEmpDist.n; j++){
  //     Rcout << curEmpPrjs[i][j] << " ";
  //   }
  //   Rcout << "\n";
  // }
  // Rcout << "refEmpDist:\n";
  // for (int i = 0; i < nDirs; i++){
  //   for (int j = 0; j < refEmpDist.n; j++){
  //     Rcout << refEmpPrjs[i][j] << " ";
  //   }
  //   Rcout << "\n";
  // }
  // DEBUG (end)
  // 3. Calculate point-wise depths on 'refEmpDist'
  for (int i = 0; i < refEmpDist.n; i++){ // for each point in 'refEmpDist'
    // Ignore the point if its weight is equal to zero
    if (refEmpDist.weights[i] == 0){
      pDepths[i] = 0;
      //Rcout << pDepths[i] << ". ";
      continue;
    }
    pDepths[i] = 1002; // initialize with zero
    for (int j = 0; j < nDirs; j++){ // for each direction
      // Calculate the cut away portions of of distributions
      double genPortion = empDistPortion(genEmpPrjs[j], genEmpDist.weights,
                                         refEmpPrjs[j][i], genEmpDist.n);
      double curPortion = empDistPortion(curEmpPrjs[j], curEmpDist.weights,
                                         refEmpPrjs[j][i], curEmpDist.n);
      //if (i == 0 && j == 366){
      //  Rcout << "Portions: " << genPortion << " and " << curPortion << "\n";
      //  if (genPortion < 1.01){
      //    Rcout << "genPortion < 1.01\n";
      //  }
      //}
      // Check whether mass of 'empDist' in halfspace equals zero
//      if (genPortion == 0 || genPortion >= 1){
      if (genPortion == 0){
        //Rcout << "Zero depth found\n";
        pDepths[i] = 0;
        //Rcout << pDepths[i] << ".. ";
        break;
      }
      double tmpDepth = 1001; // depth of current point on current projection
      // Check halfspace occupacy condition
      if (genPortion > genMinMass && curPortion > curMinMass){
        if (curPortion > 0){
          tmpDepth = genPortion / curPortion;
//          if (tmpDepth < 0.642202){
//            Rcout << "D:";
//            for (int k = 0; k < d; k++){
//              Rcout << dirs[j][k] << " ";
//            }
//            Rcout << std::endl;
//            Rcout << "N:" << genPortion << ", " << "D:" << curPortion << std::endl;
//          }
        }
      }
      pDepths[i] = tmpDepth < pDepths[i] ? tmpDepth : pDepths[i]; // update
      if (pDepths[i] == 0){
        break;
      }
      // Check halfspace occupacy condition
//      if (1 - genPortion > genMinMass && 1 - curPortion > curMinMass){
//        if (1 - curPortion > 0){
//          tmpDepth = (1 - genPortion) / (1 - curPortion);
//        }
//      }
      //if (i == 0 && j == 366){
      //  Rcout << "tmpDepth = " << tmpDepth << "\n";
      //}
      //if (tmpDepth < 0 & tmpDepth < pDepths[i]){
      //  Rcout << " " << j << " ";
      //}
//      pDepths[i] = tmpDepth < pDepths[i] ? tmpDepth : pDepths[i]; // update
//      if (pDepths[i] == 0){
//        break;
//      }
    }
    //Rcout << pDepths[i] << " ";
    // Add to the depth of the entire distribution
    depth += pDepths[i] * refEmpDist.weights[i];
  }
  //Rcout << std::endl;
  // Release memory and return
  delete[] pDepths;
  delete[] rawGenEmpPrjs;
  delete[] genEmpPrjs;
  delete[] rawCurEmpPrjs;
  delete[] curEmpPrjs;
  delete[] rawRefEmpPrjs;
  delete[] refEmpPrjs;
  return depth;
}

double calcExPointDepth(double* point, EmpDist &curEmpDist,
                        EmpDist &genEmpDist,
                        double curMinMass, double genMinMass){
  // Make a copy of each measure
  EmpDist tmpCurEmpDist(curEmpDist);
  EmpDist tmpGenEmpDist(genEmpDist);
  // Subtract 'point'
  int d = tmpCurEmpDist.d;
  for (int i = 0; i < tmpCurEmpDist.n; i++){
    for (int j = 0; j < d; j++){
      tmpCurEmpDist.vals[i][j] -= point[j];
  }
    }
  for (int i = 0; i < tmpGenEmpDist.n; i++){
    for (int j = 0; j < d; j++){
      tmpGenEmpDist.vals[i][j] -= point[j];
    }
  }
  // Return the depth
  double* tmpPoint = new double[d];
  for (int j = 0; j < d; j++){
    tmpPoint[j] = 0;
  }
  // DEBUG
/*
  Rcout << "tmpCurEmpDist:" << std::endl;
  for (int i = 0; i < tmpCurEmpDist.n; i++){
    for (int j = 0; j < tmpCurEmpDist.d; j++){
      Rcout << tmpCurEmpDist.vals[i][j] << " ";
    }
    Rcout << "(" << tmpCurEmpDist.weights[i] << ")" << std::endl;
  }
  Rcout << std::endl;
  Rcout << "tmpGenEmpDist:" << std::endl;
  for (int i = 0; i < tmpGenEmpDist.n; i++){
    for (int j = 0; j < tmpGenEmpDist.d; j++){
      Rcout << tmpGenEmpDist.vals[i][j] << " ";
    }
    Rcout << "(" << tmpGenEmpDist.weights[i] << ")" << std::endl;
  }
  // (end) DEBUG
*/
  double depth = calcExPointDepthRec(tmpPoint, tmpCurEmpDist, tmpGenEmpDist,
                                     curMinMass, genMinMass);
  delete[] tmpPoint;
  return depth;
}

double calcExPointDepthRec(double* point, EmpDist &curEmpDist,
                           EmpDist &genEmpDist,
                           double curMinMass, double genMinMass){
  // If dimension 2 then call the two-dimensional procedure
  if (curEmpDist.d == 2){
    return calcExPointDepth2D(point, curEmpDist, genEmpDist,
                              curMinMass, genMinMass, 0, 0, 0, 0);
  }
  // Otherwise call this function recursively:
  // TODO: This function works currently for dimension 3 only
  // Create two (d-1)-dimensional measure templates used in all calls after
  EmpDist tmpCurEmpDist(curEmpDist.n, curEmpDist.d - 1, false);
  EmpDist tmpGenEmpDist(genEmpDist.n, genEmpDist.d - 1, false);
  // Temporary structures
  int d = curEmpDist.d;
//  Rcout << "d = " << d << ", curEmpDist.n = " << curEmpDist.n << ", genEmpDist.n = " << genEmpDist.n << std::endl;
  double* y = new double[d - 1]; // temporary point
  double* z = new double[d]; // normalization etalon
  double result = 1.1; // the variable (depth) to be minimized
  // Projecting loop through all points of 'curEmpDist'
  for (int i = 0; i < curEmpDist.n + genEmpDist.n; i++){ // for each point
//    if (i == 27 + 20){
//      for (int j = 0; j < d; j++){
//        Rcout << genEmpDist.vals[27][j] << " ";
//      }
//      Rcout << std::endl;
//    }
    // Find maximal coordinate
    int kMax = d;
    double xMax = 0;
    if (i < curEmpDist.n){
      // In 'curEmpDist'
      for (int k = 0; k < d; k++){
        if (fabs(curEmpDist.vals[i][k]) > xMax){
          xMax = fabs(curEmpDist.vals[i][k]);
          kMax = k;
        }
      }
    }else{
      // In 'genEmpDist'
      for (int k = 0; k < d; k++){
        if (fabs(genEmpDist.vals[i - curEmpDist.n][k]) > xMax){
          xMax = fabs(genEmpDist.vals[i - curEmpDist.n][k]);
          kMax = k;
        }
      }
    }
//    if (i == 27 + 20){
//      Rcout << xMax << std::endl;
//    }
    // If it is positive then project using this point
    if (xMax > eps){
      // Obtain orthogonal vector
      if (i < curEmpDist.n){
        // In 'curEmpDist'
        for (int k = 0; k < d; k++){
          z[k] = curEmpDist.vals[i][k] / curEmpDist.vals[i][kMax];
        }
      }else{
        // In 'genEmpDist'
        for (int k = 0; k < d; k++){
          z[k] = genEmpDist.vals[i - curEmpDist.n][k] /
            genEmpDist.vals[i - curEmpDist.n][kMax];
        }
      }
      // Initialize points' counts and masses
      //int nCurPos = 0;
      //int nCurNeg = 0;
      //int nGenPos = 0;
      //int nGenNeg = 0;
      double curPosMass = 0;
      double curNegMass = 0;
      double genPosMass = 0;
      double genNegMass = 0;
      int mCur = 0;
      int mGen = 0;
      // Project 'curEmpDist'
      for (int j = 0; j < curEmpDist.n; j++){
        // Obtain point's projection
        double alpha = -curEmpDist.vals[j][kMax];
        for (int k = 0; k < kMax; k++){
          y[k] = curEmpDist.vals[j][k] + alpha * z[k];
        }
        for (int k = kMax; k < d - 1; k++){
          y[k] = curEmpDist.vals[j][k + 1] + alpha * z[k + 1];
        }
        // Decide what to do with it based on its norm
        if (norm2(y, d - 1) > eps){
          // Copy 'y' into the new empirical distribution
          memcpy(tmpCurEmpDist.vals[mCur], y, (d - 1) * sizeof(double));
          tmpCurEmpDist.weights[mCur] = curEmpDist.weights[j];
          mCur++;
        }else{
          if (alpha > eps){
            // 'y' is on the positive side
            //nCurPos++;
            curPosMass += curEmpDist.weights[j];
          }else{
            if (alpha < -eps){
              // 'y' is on the negative side
              //nCurNeg++;
              curNegMass += curEmpDist.weights[j];
            }else{
              // 'y' is in the origin and in the (hyper)plane
              //memcpy(tmpCurEmpDist.vals[mCur], y, (d - 1) * sizeof(double));
              for (int k = 0; k < d - 1; k++){
                tmpCurEmpDist.vals[mCur][k] = 0;
              }
              tmpCurEmpDist.weights[mCur] = curEmpDist.weights[j];
              mCur++;
            }
          }
        }
      }
      // Project 'genEmpDist'
      for (int j = 0; j < genEmpDist.n; j++){
        // Obtain point's projection
        double alpha = -genEmpDist.vals[j][kMax];
        for (int k = 0; k < kMax; k++){
          y[k] = genEmpDist.vals[j][k] + alpha * z[k];
        }
        for (int k = kMax; k < d - 1; k++){
          y[k] = genEmpDist.vals[j][k + 1] + alpha * z[k + 1];
        }
        // Decide what to do with it based on its norm
        if (norm2(y, d - 1) > eps){
          // Copy 'y' into the new empirical distribution
          memcpy(tmpGenEmpDist.vals[mGen], y, (d - 1) * sizeof(double));
          tmpGenEmpDist.weights[mGen] = genEmpDist.weights[j];
          mGen++;
        }else{
          if (alpha > eps){
            // 'y' is on the positive side
            //nGenPos++;
            genPosMass += genEmpDist.weights[j];
          }else{
            if (alpha < -eps){
              // 'y' is on the negative side
              //nGenNeg++;
              genNegMass += genEmpDist.weights[j];
            }else{
              // 'y' is in the origin and in the (hyper)plane
              //memcpy(tmpGenEmpDist.vals[mGen], y, (d - 1) * sizeof(double));
              for (int k = 0; k < d - 1; k++){
                tmpGenEmpDist.vals[mGen][k] = 0;
              }
              tmpGenEmpDist.weights[mGen] = genEmpDist.weights[j];
              mGen++;
            }
          }
        }
      }
      // Update point numbers
      tmpCurEmpDist.n = mCur;
      //tmpCurEmpDist.updateWeights(false);
      tmpGenEmpDist.n = mGen;
      //tmpGenEmpDist.updateWeights(false);
      // Update the depth by the recursive call of this ptocedure
      double tmpResult = calcExPointDepth2D(point, tmpCurEmpDist, tmpGenEmpDist,
                                            curMinMass, genMinMass,
                                            curPosMass, genPosMass,
                                            curNegMass, genNegMass);

      // DEBUG
      /*
      double* dirsRaw = 0;
      double** dirs = 0;
      int nDirs = 10101;
      dirsRaw = new double[nDirs * (d - 1)];
      dirs = asMatrix(dirsRaw, nDirs, (d - 1)); // new double*
      generateDirections(1, nDirs, (d - 1), dirs);
      EmpDist tmpRefEmpDist(1, d - 1, true);
      double tmpAResult = approxOneDepth(tmpRefEmpDist,
                                         tmpCurEmpDist, tmpGenEmpDist, dirs,
                                         nDirs, d - 1, curMinMass, genMinMass);
      delete[] dirs;
      delete[] dirsRaw;
      if (fabs(tmpResult - tmpAResult) > eps){
        Rcout << "i:" << i << ", E:" << tmpResult << ", A:" << tmpAResult << std::endl;
      }
       */
      // (end) DEBUG

      //double tmpResult = -1;
      if (tmpResult < result){
        result = tmpResult;
      }
      if (result == 0){
        break;
      }
    }
  }
  // Clear the memory
  delete[] y;
  delete[] z;
  return result;
}

// Update the depth
int updateDepth(double* nhMassCur, double nzMassCur, double curPosMass,
                double curNegMass, double* nhMassGen, double nzMassGen,
                double genPosMass, double genNegMass,
                double curMinMass, double genMinMass, double* depth){
  // Correct if eps-smaller than zero
  if (*nhMassCur <= 0){
    *nhMassCur = 0;
  }
  if (*nhMassGen <= 0){
    *nhMassGen = 0;
  }
  // Check favorable degenerate case of zero nominator
  if (*nhMassGen + nzMassGen + genPosMass == 0){
    *depth = 0;
    return 0;
  }
  if (*nhMassGen + nzMassGen + genNegMass == 0){
    *depth = 0;
    return 0;
  }
  if (*nhMassGen + nzMassGen + genPosMass + genNegMass == 0){
    *depth = 0;
    return 0;
  }
  double tmpResult = 1.1; // initialize against the compiler warning
  // Checks are against division by zero only
  if (*nhMassGen + nzMassGen + genPosMass > genMinMass &&
      *nhMassCur + nzMassCur + curPosMass > curMinMass){
    // Include positive-side points into halfspace
    tmpResult = (*nhMassGen + nzMassGen + genPosMass) /
      (*nhMassCur + nzMassCur + curPosMass);
    if (tmpResult < *depth){
      *depth = tmpResult;
    }
  }
  if (*nhMassGen + nzMassGen + genNegMass > genMinMass &&
      *nhMassCur + nzMassCur + curNegMass > curMinMass){
    // Include negative-side points into halfspace
    tmpResult = (*nhMassGen + nzMassGen + genNegMass) /
      (*nhMassCur + nzMassCur + curNegMass);
    if (tmpResult < *depth){
      *depth = tmpResult;
    }
  }
  if (*nhMassGen + nzMassGen + genPosMass + genNegMass > genMinMass &&
      *nhMassCur + nzMassCur + curPosMass + curNegMass > curMinMass){
    // Include points on both sides into halfspace
    tmpResult = (*nhMassGen + nzMassGen + genPosMass + genNegMass) /
      (*nhMassCur + nzMassCur + curPosMass + curNegMass);
    if (tmpResult < *depth){
      *depth = tmpResult;
    }
  }
  return 0;
}

// Calculates exact depth of a point, work in dimension 2 only
// Closely follows:
// Rousseeuw, P.J.and Ruts, I. (1996).
// Algorithm AS 307: bivariate location depth.
// Journal of the Royal Statistical Society, C, 45, 516--526.
double calcExPointDepth2D(double* point, EmpDist &curEmpDist,
                          EmpDist &genEmpDist,
                          double curMinMass, double genMinMass,
                          double curPosMass, double genPosMass,
                          double curNegMass, double genNegMass){
  // Initialization
//  Rcout << point[0] << " " << point[1] << std::endl;
  int nzCur = 0; // number of zero points in curEmpDist (denominator)
  int nzGen = 0; // number of zero points in genEmpDist (nominator)
  int nhCur = 0; // number of points in the halfplane from curEmpDist
  int nhGen = 0; // number of points in the halfplane from curEmpDist
  double nzMassCur = 0; // probability mass in the denominator for zeros
  double nzMassGen = 0; // probability mass in the nominator for zeros
  double nhMassCur = 0; // probability mass in the denominator for halfplane
  double nhMassGen = 0; // probability mass in the nominator for halfplane
  // Prepare angles structure
  RecSort* angles = new RecSort[curEmpDist.n + genEmpDist.n];
  // Compute angles and count points in the halfplane (for curEmpDist)
  for (int i = 0; i < curEmpDist.n; i++){
    if (hypot(curEmpDist.vals[i][0] - point[0],
              curEmpDist.vals[i][1] - point[1]) <= eps){
      // If in zero, do not compute the angle
      nzCur++;
      nzMassCur += curEmpDist.weights[i];
    }else{
      angles[i - nzCur].index = 0; // for curEmpDist
      angles[i - nzCur].value = atan2(curEmpDist.vals[i][1] - point[1],
                                      curEmpDist.vals[i][0] - point[0]);
      angles[i - nzCur].weight = curEmpDist.weights[i];
      if (angles[i - nzCur].value > M_PI - eps ||
          angles[i - nzCur].value < -M_PI + eps){
        angles[i - nzCur].value = -M_PI;
      }
      if (angles[i - nzCur].value <= eps){ // if in lower halfplane ...
        nhCur++; // count it
        nhMassCur += angles[i - nzCur].weight;
      }
    }
  }
  // Compute angles and count points in the halfplane (for genEmpDist)
  for (int i = 0; i < genEmpDist.n; i++){
    if (hypot(genEmpDist.vals[i][0] - point[0],
              genEmpDist.vals[i][1] - point[1]) <= eps){
      // If in zero, do not compute the angle
      nzGen++;
      nzMassGen += genEmpDist.weights[i];
    }else{
      angles[curEmpDist.n + i - nzCur - nzGen].index = 1; // for genEmpDist
      angles[curEmpDist.n + i - nzCur - nzGen].value =
        atan2(genEmpDist.vals[i][1] - point[1],
              genEmpDist.vals[i][0] - point[0]);
      angles[curEmpDist.n + i - nzCur - nzGen].weight = genEmpDist.weights[i];
      if (angles[curEmpDist.n + i - nzCur - nzGen].value > M_PI - eps ||
          angles[curEmpDist.n + i - nzCur - nzGen].value < -M_PI + eps){
        angles[curEmpDist.n + i - nzCur - nzGen].value = -M_PI;
      }
      if (angles[curEmpDist.n + i - nzCur - nzGen].value <= eps){
        // If in lower halfplane ...
        nhGen++; // count it
        nhMassGen += angles[curEmpDist.n + i - nzCur - nzGen].weight;
      }
    }
  }
  // Generalize and sort
  int nn = (curEmpDist.n + genEmpDist.n) - (nzCur + nzGen);
  quick_sort(angles, 0, nn - 1);
  double result = 1.25; // the initial depth of 'point'
  updateDepth(&nhMassCur, nzMassCur, curPosMass, curNegMass, &nhMassGen,
              nzMassGen, genPosMass, genNegMass, curMinMass, genMinMass,
              &result);
//  Rcout << "[" << result << "]" << std::endl;
//  Rcout << result << std::endl;
//  Rcout << "Cur: " << nhMassCur << "(" << nzMassCur << "); " << "gen: " << nhMassGen << "(" << nzMassGen << ")" << std::endl;
  // Minimizing loop through all possible halfplane depths
  if (result > 0){ // otherwise minimum already achieved
    int j = nhCur + nhGen; // current position
//    Rcout << j << std::endl;
//    Rcout << nn << std::endl;
    // Loop for lower halfplane
    for (int i = 0; i <= nhCur + nhGen; i++){
      // Pick the point(s) from the upper halfplane into the lower one
//      Rcout << angles[i].value << std::endl;
//      break;
      // TODO: Condition is a dangerous code!!!
      while (j < nn && ((i == nhCur + nhGen && angles[j].value <= M_PI) ||
             angles[j].value - M_PI <= angles[i].value + eps)){
        // Update (increase) the masses in the halfplane
        if (angles[j].index == 0){
          nhMassCur += angles[j].weight;
        }else{
          nhMassGen += angles[j].weight;
        }
//        Rcout << angles[j].value << "(j) ";
        // If next point to include is eps-close (a tie) then add it as well
        if (j < nn - 1 && angles[j + 1].value - angles[j].value <= eps){
          j++;
          continue;
        }
        // Update the depth if not a tie
        // TODO: Condition is a dangerous code!!!
        if (i == nhCur + nhGen ||
            angles[j].value - M_PI < angles[i].value - eps){
//          Rcout << "Cur: " << nhMassCur << "(" << nzMassCur << "); " << "gen: " << nhMassGen << "(" << nzMassGen << ")" << std::endl;
          updateDepth(&nhMassCur, nzMassCur, curPosMass, curNegMass, &nhMassGen,
                      nzMassGen, genPosMass, genNegMass, curMinMass, genMinMass,
                      &result);
//          Rcout << "[" << result << "]" << std::endl;
        }
        j++; // rotate halfplane (counter-clock)
      }
      // If at the end - no points to throw away from the halfplane
      if (i == nhCur + nhGen){
        updateDepth(&nhMassCur, nzMassCur, curPosMass, curNegMass, &nhMassGen,
                    nzMassGen, genPosMass, genNegMass, curMinMass, genMinMass,
                    &result);
//        Rcout << "[" << result << "]" << std::endl;
        break;
      }
//      Rcout << "Cur: " << nhMassCur << "(" << nzMassCur << "); " << "gen: " << nhMassGen << "(" << nzMassGen << ")" << std::endl;
      // Throw away points from the lower halfplane (continue rotating):
      // Update (reduce) the masses in the halfplane
      if (angles[i].index == 0){
        nhMassCur -= angles[i].weight;
      }else{
        nhMassGen -= angles[i].weight;
      }
//      Rcout << angles[i].value << "(i) ";
      // If next point to exclude is eps-close (a tie) then exclude it as well
      if (i < nhCur + nhGen - 1 && angles[i + 1].value - angles[i].value <= eps){
        continue;
      }
//      Rcout << "Cur: " << nhMassCur << "(" << nzMassCur << "); " << "gen: " << nhMassGen << "(" << nzMassGen << ")" << std::endl;
      updateDepth(&nhMassCur, nzMassCur, curPosMass, curNegMass, &nhMassGen,
                  nzMassGen, genPosMass, genNegMass, curMinMass, genMinMass,
                  &result);
//      Rcout << "[" << result << "]" << std::endl;
    }
//    Rcout << "[" << result << "]" << std::endl;
    // Calculate again mass in the (upper) halfspace
    nhCur = 0;
    nhGen = 0;
    nhMassCur = 0;
    nhMassGen = 0;
    for (int i = 0; i < nn; i++){
      if (angles[i].value < -M_PI + eps){
        angles[i].value = M_PI;
      }
      if (angles[i].value >= -eps){
        if (angles[i].index == 0){
          nhCur++;
          nhMassCur += angles[i].weight;
        }else{
          nhGen++;
          nhMassGen += angles[i].weight;
        }
      }
    }
    quick_sort(angles, 0, nn - 1);
    updateDepth(&nhMassCur, nzMassCur, curPosMass, curNegMass, &nhMassGen,
                nzMassGen, genPosMass, genNegMass, curMinMass, genMinMass,
                &result);
//    Rcout << "[" << result << "]" << std::endl;
    j = 0; // current position
//    Rcout << nhCur + nhGen << std::endl;
//    Rcout << nn << std::endl;
//    Rcout << "Cur: " << nhMassCur << "(" << nzMassCur << "); " << "gen: " << nhMassGen << "(" << nzMassGen << ")" << std::endl;
    // Loop for upper halfplane
    for (int i = nn - (nhCur + nhGen); i <= nn; i++){
      // Pick the point(s) from the lower halfplane into the upper one
      // TODO: Condition is a dangerous code!!!
      while (j < nn - (nhCur + nhGen) && ((i == nn && angles[j].value <= eps) ||
             angles[j].value + M_PI <= angles[i].value + eps)){
        // Update (increase) the masses in the halfplane
        if (angles[j].index == 0){
          nhMassCur += angles[j].weight;
        }else{
          nhMassGen += angles[j].weight;
        }
//        Rcout << angles[j].value << "(j) ";
        // If next point to include is eps-close (a tie) then add it as well
        if (j < nn - (nhCur + nhGen) - 1 &&
            angles[j + 1].value - angles[j].value <= eps){
          j++;
          continue;
        }
        // Update the depth if not a tie
        // TODO: Condition is a dangerous code!!!
        if (i == nn || angles[j].value + M_PI < angles[i].value - eps){
//          Rcout << "Cur: " << nhMassCur << "(" << nzMassCur << "); " << "gen: " << nhMassGen << "(" << nzMassGen << ")" << std::endl;
          updateDepth(&nhMassCur, nzMassCur, curPosMass, curNegMass, &nhMassGen,
                      nzMassGen, genPosMass, genNegMass, curMinMass, genMinMass,
                      &result);
//          Rcout << "[" << result << "]" << std::endl;
        }
        j++; // rotate halfplane (counter-clock)
      }
      // If at the end - no points to throw away from the halfplane
      if (i == nn){
        updateDepth(&nhMassCur, nzMassCur, curPosMass, curNegMass, &nhMassGen,
                    nzMassGen, genPosMass, genNegMass, curMinMass, genMinMass,
                    &result);
//        Rcout << "[" << result << "]" << std::endl;
        break;
      }
//      Rcout << "Cur: " << nhMassCur << "(" << nzMassCur << "); " << "gen: " << nhMassGen << "(" << nzMassGen << ")" << std::endl;
//      Rcout << "[" << result << "]" << std::endl;
      // Throw away points from the upper halfplane (continue rotating):
      // Update (reduce) the masses in the halfplane
      if (angles[i].index == 0){
        nhMassCur -= angles[i].weight;
      }else{
        nhMassGen -= angles[i].weight;
      }
//      Rcout << angles[i].value << "(i) ";
      // If next point to exclude is eps-close (a tie) then exclude it as well
      if (i < nn - 1 && angles[i + 1].value - angles[i].value <= eps){
        continue;
      }
//      Rcout << "Cur: " << nhMassCur << "(" << nzMassCur << "); " << "gen: " << nhMassGen << "(" << nzMassGen << ")" << std::endl;
      updateDepth(&nhMassCur, nzMassCur, curPosMass, curNegMass, &nhMassGen,
                  nzMassGen, genPosMass, genNegMass, curMinMass, genMinMass,
                  &result);
//      Rcout << "[" << result << "]" << std::endl;
    }
//    Rcout << "[" << result << "]" << std::endl;
  }
//  Rcout << "[" << result << "]" << std::endl;
  // Release memory and return the results
  delete[] angles;
  return result;
}

// Calculates exact depth of empirical measure 'curEmpDist' w.r.t. empirical
// measure 'genEmpDist' based on 'refEmpDist' projecting on directions 'dirs'
double calcOneDepth(EmpDist &refEmpDist, EmpDist &curEmpDist,
                    EmpDist &genEmpDist, int d,
                    double minMassObj, double minMassDat){
//  Rcout << "RefEmpDistN:" << refEmpDist.n << ", CurEmpDistN:" << curEmpDist.n << ", GenEmpDistN:" << genEmpDist.n << std::endl;
//  Rcout << "refEmpDist.n = " << refEmpDist.n << std::endl;
  // DEBUG
  double sum = 0;
  for (int i = 0; i < refEmpDist.n; i++){
    sum += refEmpDist.weights[i];
  }
//  Rcout << "refEmpDist: mass = " << sum << std::endl;
  sum = 0;
  for (int i = 0; i < curEmpDist.n; i++){
    sum += curEmpDist.weights[i];
  }
//  Rcout << "curEmpDist: mass = " << sum << std::endl;
  sum = 0;
  for (int i = 0; i < genEmpDist.n; i++){
    sum += genEmpDist.weights[i];
  }
//  Rcout << "genEmpDist: mass = " << sum << std::endl;
  // DEBUG (end)
  double depth = 0;// the depth
  // 1. Initialize strucrutes
  double* pDepths = new double[refEmpDist.n]; // point-wise depths
  // 2. Calculate point-wise depths on 'refEmpDist'
  for (int i = 0; i < refEmpDist.n; i++){ // for each point in 'refEmpDist'
    // Ignore the point if its weight is equal to zero
    if (refEmpDist.weights[i] == 0){
      pDepths[i] = 0;
      //Rcout << pDepths[i] << ". ";
      continue;
    }
    // Calculate depth
    pDepths[i] = calcExPointDepth(refEmpDist.vals[i], curEmpDist, genEmpDist,
                                  minMassObj, minMassDat);
    //pDepths[i] = calcExPointDepth2D(refEmpDist.vals[i], curEmpDist, genEmpDist,
    //                              minMassObj, minMassDat, 0, 0, 0, 0);
    // Add to the depth of the entire distribution
    depth += pDepths[i] * refEmpDist.weights[i];
  }
  // Release memory and return
  delete[] pDepths;
  return depth;
}
