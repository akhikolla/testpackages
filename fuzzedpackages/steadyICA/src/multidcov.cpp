// includes from the plugin
 
#include <Rcpp.h>
#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif
#ifndef END_RCPP
#define END_RCPP
#endif
using namespace Rcpp;

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <R_ext/Arith.h>
#include <Rinternals.h>

// user includes
using namespace std;
// declarations
extern "C" {
  SEXP dcovustatC( SEXP x, SEXP y, SEXP alpha);
  SEXP gradmdcov( SEXP Z1, SEXP Z2) ;
  }

// Source cpp script for dCov.Ustat() function
SEXP dcovustatC( SEXP x, SEXP y, SEXP alpha ){
  BEGIN_RCPP

  Rcpp::NumericMatrix X(x);              // creates Rcpp matrix X from SEXP x
  Rcpp::NumericMatrix Y(y);              // creates Rcpp matrix Y from SEXP y
  double Alpha = as<double>(alpha);	 // convert from SEXP to double
  int n_X = X.nrow(), d_X = X.ncol();    // find number of rows n_X and cols d_X in matrix X
  int n_Y = Y.nrow(), d_Y = Y.ncol();    // find number of rows n_Y and cols d_Y in matrix Y
  double n = X.nrow();

  if( n_X != n_Y ){
    return Rcpp::wrap("Sample sizes must agree!");    // check that n_X == n_Y
  }

  double* pX = X.begin();    // creates pointer for X, matrix is indexed as vector, in column order 
  double* pY = Y.begin();    // creates pointer for Y, matrix is indexed as vector, in column order

  double T1  = 0;         // creates T1  scaler for dCov
  double T2x = 0;         // creates T2x scaler for dCov, T2 = T2x*T2y
  double T2y = 0;         // creates T2y scaler for dCov, T2 = T2x*T2y
  double T3  = 0;         // creates T3  scaler for dCov
  double* pT1  = &T1;     // creates pointer for T1
  double* pT2x = &T2x;    // creates pointer for T2x
  double* pT2y = &T2y;    // creates pointer for T2y
  double* pT3  = &T3;     // creates pointer for T3
  double dX = 0;          // creates scaler for |X_k - X_l| distance dX
  double dY = 0;          // creates scaler for |Y_k - Y_l| distance dY
  double* pdX = &dX;      // creates pointer for dX
  double* pdY = &dY;      // creates pointer for dY

  Rcpp::NumericVector Xsum(n_X);    // creates vector for single sum over dX
  Rcpp::NumericVector Ysum(n_Y);    // creates vector for single sum over dX
  double* pXsum = Xsum.begin();     // creates pointer for Xsum 
  double* pYsum = Ysum.begin();     // creates pointer for Ysum

  int i,j,k,l = 0;                // creates indices for looping  
  for (k = 1; k < n_X; k++){      // begin first loop over index k = 2:n
    for (l = 0; l < k; l++){      // begin second loop over index l < k

// calculate |X_k - X_l| distance
      *pdX = 0;                   // reset dX distance to 0 via pointer
      i = 0;                      // reset column counter for X to 0
      for(i = 0; i < d_X; i++){
        *pdX += std::pow( (pX[n_X*i+k] - pX[n_X*i+l]), 2);    // calculate sum of squared element-wise differences 
      }
      *pdX = std::pow( sqrt(*pdX), Alpha);    // calculate Euclidean distance |X_k - X_l|^Alpha

// calculate |Y_k - Y_l| distance
      *pdY = 0;                   // reset dY distance to 0 via pointer
      j = 0;                      // reset column counter for Y to 0
      for(j = 0; j < d_Y; j++){
        *pdY += std::pow( (pY[n_Y*j+k] - pY[n_Y*j+l]), 2);    // calculate sum of squared element-wise differences 
      }
      *pdY = std::pow( sqrt(*pdY), Alpha);    // calculate Euclidean distance |Y_k - Y_l|^Alpha

// update T1, and k_th element of Xsum and Ysum
      *pT1 += (*pdX) * (*pdY);    // update T1 via pointer
    //BRisk edits 25 September 2012: changed from Xsum[k] to pXsum[k];
    //appears to be slightly faster. 
      pXsum[k] += *pdX;    // update Xsum, for updating T2x and T3  
      pYsum[k] += *pdY;    // update Ysum, for updating T2y and T3
      pXsum[l] += *pdX;    // update Xsum, for updating T3
      pYsum[l] += *pdY;    // update Ysum, for updating T3

    }    // end second loop over index l

// update T2x, T2y
    *pT2x += pXsum[k] ;    // update T2x
    *pT2y += pYsum[k] ;    // update T2y

  }    // end first loop over index k

// calculate T3
  for (k = 0; k < n_X; k++){                // begin final loop over index k
    *pT3  += pXsum[k] * pYsum[k] / 3;    // update T3 via pointer
  }                                         // end final loop over index k
// remove extra T1 type terms
    *pT3 -= 2*(*pT1)/3;

// Rf_PrintValue(Rcpp::wrap(T3));    // example print statement for debugging 

// update T1, T2x, T2y, T3
  *pT1  = 2 * (*pT1) * ((n) / (n-1));           // update T1 via pointer
  *pT2x = 2 * (*pT2x) / (n-1);                  // update T2x via pointer
  *pT2y = 2 * (*pT2y) / (n-1);                  // update T2y via pointer
  *pT3  = 6 * (*pT3) * ((n) / ((n-1)*(n-2)));   // update T3 via pointer
//  *pT1  = *pT1  / (n*(n-1)/2);         // update T1 via pointer
//  *pT2x = *pT2x / (n*(n-1)/2);         // update T2x via pointer
//  *pT2y = *pT2y / (n*(n-1)/2);         // update T2y via pointer
//  *pT3  = *pT3  / (n*(n-1)*(n-2)/6);   // update T3 via pointer

// calculate T2
  double T2 = *pT2x * *pT2y ;    // creates T2 element of dCov, T2 = T2x*T2y

// calculate sample distance covariance V squared
  double V2 = (*pT1 + T2 - (*pT3));

//return Rcpp::List::create(
//Rcpp::Named("T1") = Rcpp::wrap(T1),
//Rcpp::Named("T2x") = Rcpp::wrap(T2x),
//Rcpp::Named("T2y") = Rcpp::wrap(T2y),
//Rcpp::Named("T3") = Rcpp::wrap(T3),
//Rcpp::Named("T2") = Rcpp::wrap(T2),
//Rcpp::Named("V2") = Rcpp::wrap(V2)
//) ;

return Rcpp::wrap(V2);    // coerced to a comparable R type
END_RCPP
}






SEXP gradmdcov( SEXP Z1, SEXP Z2 ){
BEGIN_RCPP
  
  NumericVector Sq(Z1);
  NumericMatrix Snq(Z2);
  int nrow = Sq.length();
  int* pnrow = &nrow;
  int nrowsq =  *pnrow * *pnrow;
  int* pnrowsq = &nrowsq;
  int ncol = Snq.ncol()+1;
  int* pncol = &ncol;
  int ncolm1 = *pncol - 1;
  int* pncolm1 = &ncolm1;
  double* pSq = Sq.begin(); //pointer to Sq indexed as a vector.
  double* pSnq = Snq.begin();


  //********************
  //Start of the program that will be used in R:  
  //Declarations:

  double coefT1 = 2.0 * *pnrow / (*pnrow-1); //modified for input S = S/n 
  double coefT2 = 2.0 / (*pnrow-1); //modified for input S = S/n
  double coefT3 = 2.0 * *pnrow / ((*pnrow-1) * (*pnrow-2)); //modified for input S = S/n
  double* pcoefT1 = &coefT1;
  double* pcoefT2 = &coefT2;
  double* pcoefT3 = &coefT3;

  int index;
  int* pindex = &index;
  int tindex;
  int* ptindex = &tindex;
  double diffSp = 0;
  double* pdiffSp = &diffSp;
  double dSnq = 0;
  double* pdSnq = &dSnq;

  NumericVector distSq(*pnrow * *pnrow);
  double* pdistSq = distSq.begin();

  IntegerVector gDistSq(*pnrowsq);
  int* pgDistSq = gDistSq.begin();

  NumericVector distSnq(*pnrowsq);
  double* pdistSnq = distSnq.begin();

  NumericVector colSumsDistSnq(*pnrow);
  double* pcolSumsDistSnq = colSumsDistSnq.begin();
  NumericVector colSumsDistSq(*pnrow);
  double* pcolSumsDistSq = colSumsDistSq.begin();

  double T2x = 0;
  double T2y = 0;
  double* pT2x = &T2x;
  double* pT2y = &T2y;

  double gDistSpNq;
  double* pgDistSpNq = &gDistSpNq;
  double gT1p;
  double* pgT1p = &gT1p;
  double gT3p;
  double* pgT3p = &gT3p;

  double sgSp;
  double* psgSp = &sgSp;
  NumericVector gSq(*pnrow);
  double* pgSq = gSq.begin();
  NumericVector gSnq(*pnrow * *pncolm1);
  double* pgSnq = gSnq.begin();
  
  //Create distance matrices, gDistSq, and column sums of distances
  for(int i = 1; i < *pnrow; i++) {
    for(int j = 0; j < i; j++) {
      *pindex = *pnrow * i + j;
      *ptindex = *pnrow * j + i;
      *pdiffSp = pSq[i] - pSq[j];
      pgDistSq[*pindex] = -1L * (*pdiffSp < 0) + 1L * (*pdiffSp > 0);
      pgDistSq[*ptindex] = -1L * pgDistSq[*pindex];
      *pdiffSp = std::abs(*pdiffSp);
      pdistSq[*pindex] = *pdiffSp;
      pdistSq[*ptindex] = *pdiffSp;

      *pdSnq = 0;    
      for(int p = 0; p < *pncolm1; p++) { 
        *pdSnq += std::pow(pSnq[*pnrow * p + i] - pSnq[*pnrow * p + j],2);
        } 
      *pdSnq = sqrt(*pdSnq);
      pdistSnq[*pindex] = *pdSnq;
      pdistSnq[*ptindex] = *pdSnq;

      pcolSumsDistSnq[i] += *pdSnq;
      pcolSumsDistSnq[j] += *pdSnq;

      pcolSumsDistSq[i] += *pdiffSp;
      pcolSumsDistSq[j] += *pdiffSp;

      *pT2y += *pdSnq;
      *pT2x += *pdiffSp;
  }
}

  *pT2x *= *pcoefT2; 
  *pT2y *= *pcoefT2;

  //Gradient of Sq:
  for(int i = 0; i < *pnrow; i++) {
    for(int j = *pnrow-1; j > i; j--) {
      *pindex = i * *pnrow + j; //corresponds to matrix indices [j,i]
      //pgDistSq[*pindex] TO DO: make a reference
      *pgT1p =  pgDistSq[*pindex] * pdistSnq[*pindex];
      *pgT3p = pgDistSq[*pindex] * (pcolSumsDistSnq[i] + pcolSumsDistSnq[j]) - 2.0 * *pgT1p;
      *psgSp = ((*pgT1p * *pcoefT1) + (pgDistSq[*pindex] * *pT2y * *pcoefT2) - (*pgT3p * *pcoefT3)) / *pnrow;  
//NOTE: need to divide by *pnrow because scaling S by 1/n does not affect pgDistSq
//TO DO: Bring multiplication outside of j loop; would take more memory. 
      pgSq[i] += *psgSp; // add to [i,q]
      pgSq[j] -= *psgSp; //subtract from [j,q]
    }
}

  //Gradient of Sp for p ne q:
  for (int p = 0; p < *pncolm1; p++) {
    for (int i = 0; i < *pnrow; i++) {  
     for (int j = *pnrow-1; j > i; j--) {
        *pindex = i * *pnrow + j;
	      *pdiffSp = pSnq[p * *pnrow + i] - pSnq[p * *pnrow + j];
        /*note: we have already calculated this diffSp, but we need the
        denominator before we can use it; thus, we would have to store
        all these p * nrow * nrow terms*/

	      if (pdistSnq[*pindex] == 0) 
          *pgDistSpNq = 0;
	      else 
          *pgDistSpNq = *pdiffSp / pdistSnq[*pindex];
          
          *pgT1p = pdistSq[*pindex] * *pgDistSpNq;
        //TO DO: make reference pgT2p = *pgDistSpnq;

	      *pgT3p = *pgDistSpNq * (pcolSumsDistSq[i] + pcolSumsDistSq[j]) - 2.0 * *pgT1p;
       
        *pindex = p * *pnrow + i; //corresponds to matrix indices [i,p]
        *ptindex = p * *pnrow + j; //corresponds to matrix indices [j,p]
        *psgSp = ((*pgT1p * *pcoefT1) + (*pgDistSpNq  * *pT2x * *pcoefT2) - (*pgT3p * *pcoefT3)) / *pnrow;
//NOTE: Need to divide by 1/n because scaling S by (1/n) does not affect pgDistSpNq
        pgSnq[*pindex] += *psgSp; // add to [i,p]
        pgSnq[*ptindex] -= *psgSp; //subtract from [j,p]    
      }
     }
    }
return List::create(Named("gSnq") = gSnq, Named("gSq") = gSq);

END_RCPP
}
