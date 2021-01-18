#ifndef __sqp_types_included__        // if include guard for 'sqp_types.h' is undefined
#define __sqp_types_included__        // define include guard for 'sqp_types.h'



// Config
#include "sqp_config.h"

#if SQP_USE_SUPERLU == 1
#define ARMA_USE_SUPERLU 1
#endif

// External headers
#include "RcppArmadillo.h"
#include "stdlib.h"
#include "limits.h"


// Conditional external headers
#if SQP_USE_EIGEN == 1
#include "RcppEigen.h"
#include "Eigen/SparseQR"
#include "Eigen/SparseLU" 
#include "Eigen/SparseCholesky"
#endif // SQP_USE_EIGEN == 1


#if SQP_USE_VIENNACL == 1
#define VIENNACL_WITH_ARMADILLO 1
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
#include "viennacl/linalg/row_scaling.hpp"
#include "viennacl/linalg/detail/ilu/chow_patel_ilu.hpp"
#include "viennacl/linalg/detail/ilu/ilut.hpp"
#include "viennacl/linalg/detail/ilu/ilu0.hpp"
#include "viennacl/linalg/direct_solve.hpp"
#include "viennacl/linalg/lu.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/linalg/mixed_precision_cg.hpp"
#include "viennacl/io/matrix_market.hpp"
#endif // SQP_USE_VIENNACL == 1


// Package headers
#include "sqp_bfgs.h"
#include "sqp_misc.h"
#include "sqp_solvers.h"



#endif                                 // end of include guard for 'sqp_types.h'
