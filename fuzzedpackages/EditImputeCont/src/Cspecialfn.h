// #ifndef SpecialDistributionsH
// #define SpecialDistirbutionsH

#ifndef CSpecialFnH   // HANG changed this on May 16, 2015
#define CSpecialFnH

#ifdef use_namespace
using namespace RBD_LIBRARIES;
using namespace NEWRAN;
#endif
using namespace std;

double logdet(LowerTriangularMatrix& lchol);

ReturnMatrix exp_ColumnVector(ColumnVector &X);
ReturnMatrix log_ColumnVector( ColumnVector &X); 
ReturnMatrix Submatrix_elem( SymmetricMatrix &Sigma, ColumnVector &id_rows);
ReturnMatrix Submatrix_elem_2( SymmetricMatrix &Sigma, ColumnVector &id_rows, ColumnVector &id_cols ); // ADDED 2015/01/27
ReturnMatrix whichOnes( ColumnVector &ind_vector );
ReturnMatrix whichZeros( ColumnVector &ind_vector );
ColumnVector subvector(ColumnVector &Y, ColumnVector &S);
ColumnVector subvector_by_index(ColumnVector &Y, ColumnVector &Index);

ReturnMatrix rMVN_fn(ColumnVector &mu, LowerTriangularMatrix &LSigma );
double log_MVN_fn(ColumnVector &x, ColumnVector &mu, LowerTriangularMatrix &LSigma_i, double logdetandmore);
double log_MVN_fn(ColumnVector &x, ColumnVector &mu, LowerTriangularMatrix &LSigma_i);

int rdiscrete_fn(ColumnVector &Prob);
int runifdiscrete_fn(int n_var);

ReturnMatrix rIW_w_pd_check_fn( int nu, LowerTriangularMatrix &LPhi );

double rgamma_fn( double a, double b );
double rbeta_fn( double a, double b );

#endif
