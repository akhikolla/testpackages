/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#ifndef nro_INCLUDED
#define nro_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <ctime>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <Rcpp.h>
#include "medusa.h"
#include "abacus.h"
#include "punos.h"
#include "koho.h"
#include "scriptum.h"

using namespace std;
using namespace Rcpp;
using namespace medusa;
using namespace abacus;
using namespace punos;
using namespace koho;
using namespace scriptum;

/*
 *
 */
namespace nro {

  /*
   *
   */
  extern vector<vector<mdreal> > matrix2reals(const SEXP&, const mdreal);
  extern NumericMatrix reals2matrix(const vector<vector<mdreal> >&);
  extern NumericVector reals2vector(const vector<mdreal>&);
  extern vector<mdreal> vector2reals(const SEXP&);
  extern vector<mdsize> vector2sizes(const SEXP&);
  extern punos::Topology reals2topology(const vector<vector<mdreal> >&,
					const mdreal);
}

using namespace nro;

#endif /* nro_INCLUDED */

