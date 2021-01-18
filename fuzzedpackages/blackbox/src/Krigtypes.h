#ifndef KRIGTYPES_H
#define KRIGTYPES_H

/*** the sensitive part of the computations are the sums involving 1/d_i where eigenvalue d_i can be close to 0
 therefore
* the variables where such summations are made should have maximum precision;
the Kriging predictor itself is so
* the computation of the d_i themselves sould have high precisions
i.e. the eigensystem computation ***/
typedef long double internalTypeEigen; // gives type to 'Real' dans jama_eig.h which has been 'scanned for numerical precision'
typedef double CQRtype; ///DOUBLE: gives type to 'QRtype' dans qr.h.

/*** other computations may be less important. In particular, the result of the above summations could itself be stored in lower precision;
in the end only a double is considered by R. **/
typedef double ioType; // for input/output communication with R. no need for template here, but it will be easy to add template if necess.

/*** minimization/maximization operations  ***/
/*** unclear status operations ***/
/// float covTypedef compiles but fails...
typedef  double covTypedef; // the value of the covType typename : for elements of covariance matricex, Bessel fns in particular
// float generates an error (nicely handled by error catching code...)
// although long double \sim double in win32, long double seems faster... !
// floats below do not seem to make things faster...
typedef double brentTypedef; // internal computations of brent one-dimensional minimization

extern int fittedparamnbr;
extern bool verbose,batchDebug;

#endif
