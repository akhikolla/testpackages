#ifndef __sqp_config_included__        // if include guard for 'sqp_config.h' is undefined
#define __sqp_config_included__        // define include guard for 'sqp_config.h'


/*---------------------------------------------------------
  Optional: Availability of solvers from external libraries
  --------------------------------------------------------*/

#define SQP_USE_EIGEN 1                //  'Eigen'     library: http://eigen.tuxfamily.org/
#define SQP_USE_SUPERLU 0              //  'SuperLU'   library: https://portal.nersc.gov/project/sparse/superlu/
#define SQP_USE_VIENNACL 0             //  'ViennaCL'  library: http://viennacl.sourceforge.net/

// #define ARMA_DONT_USE_WRAPPER 0       //  Armadillo:  Disable wrappers. You will need to directly link with BLAS, LAPACK
#define ARMA_USE_LAPACK 1             //  Armadillo:  Use Lapack 
#define ARMA_USE_BLAS 1               //  Armadillo:  Use BLAS 
#define ARMA_BLAS_LONG 1              //  Armadillo:  Use long number formats for BLAS 



#endif                                 // end of include guard for 'sqp_config.h'
