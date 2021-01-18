
#include "BinTree.h"
#include "Jump.h"
#include "TriArray.h"
#include "Bounds.h"

#include "Debug.h"

#include <Rversion.h>
#if defined(R_VERSION) && R_VERSION >= R_Version(3, 3, 0)
  /* nothing */
#else
  #define NO_C_HEADERS true // disables the including of stdlib.h, stdio.h, limits.h, math.h by R.h
  #include <cstdlib>        // manually loading of cpp versions of disabled headers
  #include <cstdio>
  #include <climits>
  #include <cmath>
  #include <cstddef>
  using std::size_t;
#endif

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/***************
* class Step
* virtual class allowing to fit step functions to data
* Thomas Hotz, 2007-2011
***************/

class Step {
  public:
    unsigned int N; // length of input data
    double* lb; // lower bound
    double* ub; // upper bound

    // initialization
    Step(unsigned int n); // constructor for n data points
    Step(unsigned int n, double* xlb, double* xub); // constructor for n data points

    // needed for all distributions
    // needs to be declared as pure virtual (= 0) here (or else vtable and typeinfo are not generated in object file, as visible with "nm -C Step.o")
    virtual double cost(unsigned int startIndex, unsigned int endIndex) const = 0; // calculate cost of a block
//     virtual double cost(unsigned int startIndex, unsigned int endIndex) const; // calculate cost of a block
    virtual double costBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const; // calculate cost of a block given bounds
    virtual double estBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const; // corresponding estimate

    // computations
    Jump findCandidate(const Jump& prev, const Jump& next) const; // find a candidate jump
    SEXP forward(unsigned int maxNum) const; // find candidates through forward selection
    void flattenTree(BinTree<Jump>* bt, int* number, int* depth, int* rightIndex, double* improve, double &totalcost) const;

    SEXP path(unsigned int maxNum) const; // compute solution path
    SEXP bounded(Bounds& B) const; // compute optimal solution with minimal jumps fulfilling bounds
};
