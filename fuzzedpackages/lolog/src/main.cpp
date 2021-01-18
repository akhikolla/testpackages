//#define INSIDE
#ifdef INSIDE

#include <vector>
#include <Rcpp.h>
#include "BinaryNet.h"
#include <RInside.h> // for the embedded R via RInside
#include "Stat.h"
#include "Stats.h"
#include "StatController.h"
#include "Offset.h"
#include "Offsets.h"
#include "Constraint.h"
#include "Constraints.h"
#include "Model.h"
#include "VarAttrib.h"
#include <memory>
#include <boost/shared_ptr.hpp>
#include <assert.h>
#include "tests.h"
#include <boost/container/flat_set.hpp>
#include <ctime>
#undef NDEBUG
using namespace Rcpp;
using namespace std;
using namespace lolog;


RcppExport SEXP _rcpp_module_boot_lolog(){return NULL;}
/*!
 * Main entry point when called from outside R (with R embedded via RInside).
 */
int main(int argc, char *argv[]) {
    RInside Rins(argc, argv);
    initStats();
    tests::runLologTests();
    exit(0);
    return 0;
}


#endif
