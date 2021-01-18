#include "stdafx.h"

const double eps_band = 1e-10;

void BandDepth(T3DMatrix x, T3DMatrix X, int m, int n, int t, int d, bool modif, 
               int J, double* depths){
  //Rcout << "m = " << m << endl;
  //Rcout << "n = " << n << endl;
  //Rcout << "t = " << t << endl;
  //Rcout << "d = " << d << endl;
  //Rcout << "modif = " << modif << endl;
  //Rcout << "J = " << J << endl;
  // Prepare data structures for the loop through all J combinations
  double* b = new double[d + 1]; b[d] = 1;
  double* z = new double[d + 1];
  int* counters = new int[d + 1];
  TDMatrix A = newM(d + 1, d + 1);
  unsigned long long div0 = choose(n, d + 1); // num simplices
  // Loop for all observations to compute depth for
  for (int iObs = 0; iObs < m; iObs++){
    unsigned long long theCounter = 0;
    //unsigned long long numSimplicesChecked = 0;
    // Loop to check all combinations of J functions out of n
    for (int i = 0; i < d; i++){ counters[i] = i; } counters[d] = d - 1;
    while (counters[0] != n - (d + 1)){
      int i = d;
      while (i > 0 && counters[i] == n - (d + 1) + i){ i--; }
      counters[i]++; int j = i + 1;
      while (j < d + 1){ counters[j] = counters[j - 1] + 1; j++; }
      // Execute logic for a single (d+1)-tuple of functoins:
      bool isInBand = true;
      // Loop for all time points
      for (int iTime = 0; iTime < t; iTime++){
        // Check whether current function is inside simplex for this time point
        for (int j = 0; j < d; j++){
          for (int k = 0; k < d + 1; k++){
            A[j][k] = X[counters[k]][iTime][j];
          }
        }
        for (int k = 0; k < d + 1; k++){
          A[d][k] = 1;
        }
        memcpy(b, x[iObs][iTime], d * sizeof(double)); b[d] = 1;
        if (solveUnique(A, b, z, d + 1)){
          bool isInside = true;
          for (int j = 0; j < d + 1; j++){
            if (z[j] < -eps_band){ isInside = false; break; }
          }
          if (isInside){ // if inside simplex
            if (modif){
              theCounter++;
            }
          }else{ // if outside simplex
            if (!modif){
              isInBand = false;
              break;
            }
          }
        }
        //Rcout << " " << isInBand;
      }
      if (!modif){ // if not modified version
        theCounter += isInBand; // add 1 once only (not each time point)
      }
      //Rcout << " = " << theCounter << endl;
      //if (modif){ // if modified version
      //  numSimplicesChecked += t;
      //}else{
      //  numSimplicesChecked++;
      //}
    }
    if (modif){ // if modified version
      depths[iObs] = (double)theCounter / (div0 * t);
    }else{ // if not modified version
      depths[iObs] = (double)theCounter / div0;
    }
  }
  // Release memory
  delete[] b;
  delete[] z;
  delete[] counters;
  deleteM(A);
}
