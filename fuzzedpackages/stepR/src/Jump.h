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

/***************
* class Jump
* describes a jump in a block signal
* Thomas Hotz, 2007-2008
***************/

class Jump {
  public:
    Jump(int n, int ri, double im); // constructor
    Jump(); // default constructor initialises to "before" the data

    int number; // the (ordinal) number of the jump
    int rightIndex; // the index of the right end of the block described
    double improve; // the improvement that the inclusion of this jump brought
};
