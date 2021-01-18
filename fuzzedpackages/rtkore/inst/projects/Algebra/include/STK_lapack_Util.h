/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Algebra
 * created on: 10 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_lapack_Util.h
 *  @brief In this file we define and implement utilities classes and methods
 *  for the interface with lapack.
 **/


#ifndef STK_LAPACK_UTIL_H
#define STK_LAPACK_UTIL_H

#ifdef STKUSELAPACK

extern "C"
{
/** lapack routine computing the optimal block size in least square problems */
int ilaenv_(int *, char *, char *, int *, int *,int *, int *);

#ifdef STKREALAREFLOAT

/** LAPACK routine in float computing SVD decomposition */
void sgesdd_( char *jobz, int *M, int *N, float *A, int *lda, float *S, float *U, int* ldu
            , float *vt, int *ldvt, float *work, int *lWork, int *iwork, int *info);
/** LAPACK routine in float computing QR decomposition */
extern void sgeqrf_(int* M, int* N, float* A, int* lda, float* TAU, float* work, int* lWork, int* info );
/** LAPACK routine in float computing eigenvalues decomposition */
extern void ssyevr_( char *, char *, char *, int *, float *, int *, float *,
                     float *, int *, int *, float *, int *, float *, float *, int *,
                    int *, float *, int *, int *, int *, int *);
/** LAPACK routine in float computing least square solution */
extern int sgelsd_( int *m, int *n, int *nrhs
                  , float *a, int *lda, float *b, int *ldb, float *s, float *rcond
                  , int *rank, float *work, int *lWork
                  , int *iwork, int *info);

#else /* double */

/** LAPACK routine in double computing SVD decomposition */
void dgesdd_( char *jobz, int *M, int *N, double *A, int *lda, double *S, double *U, int* ldu
            , double *vt, int *ldvt, double *work, int *lWork, int *iwork, int *info);
/** LAPACK routine in double to computing QR decomposition */
extern void dgeqrf_(int* M, int* N, double* A, int* lda, double* TAU, double* work, int* lWork, int* info );
/** LAPACK routine in double to compute the eigenvalues */
extern void dsyevr_( char *, char *, char *, int *, double *, int *, double *,
                     double *, int *, int *, double *, int *, double *, double *, int *,
                    int *, double *, int *, int *, int *, int *);
/** LAPACK routine in double computing least square solution */
extern int dgelsd_( int *m, int *n, int *nrhs
                  , double *a, int *lda, double *b, int *ldb, double *s, double *rcond
                  , int *rank, double *work, int *lWork
                  , int *iwork, int *info);

#endif

} // extern "C"

#endif // STKUSELAPACK

namespace STK
{
namespace lapack
{
/** @ingroup Algebra
 *  @brief wrapper of the LAPACK GELSD routine: computes the minimum-norm
 *  solution to a real linear least squares problem.
 *
 *  GELSD computes the minimum-norm solution to a real linear least squares
 *  problem: minimize 2-norm(| b - A*x |) using the singular value decomposition
 *  (SVD) of A. A is an M-by-N matrix which may be rank-deficient.
 *
 *  @param[in] m  The number of rows of A. <em> m >= 0 </em>.
 *  @param[in] n  The number of columns of A. <em> n>= 0 </em>.
 *  @param[in] nrhs The number of right hand sides, i.e., the number of
 *  columns of the matrices B and X. <em>nrhs >= 0</em>.
 *  @param[in] a On entry, the M-by-N matrix A. On exit, A has been destroyed.
 *  @param[in] lda The leading dimension of the array A.  <em>lda >= max(1,m)</em>.
 *  @param[in,out] b On entry, the M-by-NRHS right hand side matrix B.
 *  @verbatim
 *    On exit, B is overwritten by the N-by-NRHS solution matrix X.
 *    If m >= n and RANK = n, the residual sum-of-squares for the solution in
 *    the i-th column is given by the sum of squares of elements n+1:m in that
 *    column.
 *  @endverbatim
 *  @param[in] ldb The leading dimension of the array B.
 *  <em> ldb >= max(1,max(m,n)) </em>.
 *
 *  @param[out] s The singular values of A in decreasing order.
 *  @verbatim
 *    The condition number of A in the 2-norm = s(1)/s(min(m,n)).
 *  @endverbatim
 *
 *  @param[out] rcond used to determine the effective rank of A.
 *  @verbatim
 *    Singular values s(i) <= RCOND*s(1) are treated as zero.
 *    If RCOND < 0, machine precision is used instead.
 *  @endverbatim
 *
 *  @param[out] rank The effective rank of A, i.e., the number of singular values
 *  which are greater than RCOND*S(1).
 *
 *  @param[out] work On exit, if info = 0, work(1) returns the optimal lWork.
 *
 *  @param[in] lWork The dimension of the array @c work. @c lWork must be at least 1.
 *  @verbatim
 *    The exact minimum amount of workspace needed depends on M, N and NRHS.
 *    As long as lWork is at least 12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
 *    if M is greater than or equal to N or 12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
 *    if M is less than N, the code will execute correctly.
 *    SMLSIZ is returned by ILAENV and is equal to the maximum size of
 *    the subproblems at the bottom of the computation tree (usually about 25),
 *    and NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
 *    For good performance, lWork should generally be larger.
 *    If lWork = -1, then a workspace query is assumed; the routine only
 *    calculates the optimal size of the work array, returns this value as the
 *    first entry of the work array, and no error message related to lWork is
 *    issued by XERBLA.
 *  @endverbatim
 *
 *  @param iwork  array of integer of dimension (MAX(1,LiWork))
 *  @verbatim
 *    LiWork >= 3 * MINMN * NLVL + 11 * MINMN, where MINMN = MIN( M,N ).
 *    On exit, if info = 0, iWork(1) returns the minimum LiWork
 *  @endverbatim
 *
 *  @return info
 *  @verbatim
 *    = 0:  successful exit
 *    < 0:  if info = -i, the i-th argument had an illegal value.
 *    > 0:  the algorithm for computing the SVD failed to converge;
 *    if info = i, i off-diagonal elements of an intermediate bidiagonal
 *    form did not converge to zero.
 *  @endverbatim
 *
 *  @verbatim
 *    Further Details
 *    ===============
 *      Based on contributions by
 *      Ming Gu and Ren-Cang Li, Computer Science Division, University of
 *      California at Berkeley, USA
 *      Osni Marques, LBNL/NERSC, USA
 * @endverbatim
 **/
inline int gelsd( int m, int n, int nrhs
                , Real * a, int lda
                , Real * b, int ldb
                , Real * s
                , Real *rcond, int *rank
                , Real *work, int lWork, int* iwork)
{
  int info = 1;
#ifdef STKUSELAPACK
#ifdef STKREALAREFLOAT
  sgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, rcond, rank, work, &lWork, iwork, &info);
#else
  dgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, rcond, rank, work, &lWork, iwork, &info);
#endif
#endif
  return info;
}

/** @ingroup Algebra
 *  @brief wrapper of the LAPACK DGEQRF routine: computes the Qr decomposition
 *  of a matrix.
 *
 * @param[in] m The number of rows of the matrix A.  M >= 0.
 * @param[in] n The number of columns of the matrix A.  N >= 0.
 *
 * @param[in,out] a Real array, dimension (lda, N)
 * \verbatim
 *     On entry, the M-by-N matrix A.
 *     On exit, the elements on and above the diagonal of the array
 *     contain the min(M,N)-by-N upper trapezoidal matrix R (R is
 *     upper triangular if m >= n); the elements below the diagonal,
 *     with the array TAU, represent the orthogonal matrix Q as a
 *     product of min(m,n) elementary reflectors (see Further Details).
 * \endverbatim
 *
 * @param[in] lda The leading dimension of the array A.  lda >= max(1,M).
 *
 * @param[out] tau Real array, dimension min(M,N)
 * The scalar factors of the elementary reflectors (see Further Details).
 *
 * @param[in,out] work Real array, dimension (MAX(1,lWork))
 * \verbatim
 *   On exit, if info = 0, work(1) returns the optimal lWork.
 * \endverbatim
 *
 * @param[in] lWork The  dimension  of  the array work
 * \verbatim
 *  lWork >= max(1,N).
 *  For optimum performance lWork >= N*NB, where NB is the optimal blocksize.
 *
 *  If lWork = -1, then a workspace query is assumed; the routine
 *  only calculates the optimal size of the work array, returns
 *  this value as the first entry of the work array, and no error
 *  message related to lWork is issued by XERBLA.
 * \endverbatim
 *
 * @return info
 * \verbatim
 *  = 0:  successful exit
 *  < 0:  if info = -i, the i-th argument had an illegal value
 * \endverbatim
 *
 * @verbatim
 *  Further Details
 *  ===============
 *
 *  The matrix Q is represented as a product of elementary reflectors
 *
 *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
 *
 *  Each H(i) has the form
 *
 *     H(i) = I - tau * v * v'
 *
 *  where tau is a real scalar, and v is a real vector with
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
 *  and tau in TAU(i).
 * @endverbatim
 */
inline int geqrf(int m, int n, Real* a, int lda, Real* tau, Real *work, int lWork)
{
  int info = 1;

#ifdef STKUSELAPACK
#ifdef STKREALAREFLOAT
  sgeqrf_(&m, &n, a, &lda, tau, work, &lWork, &info);
#else
  dgeqrf_(&m, &n, a, &lda, tau, work, &lWork, &info);
#endif
#endif

  return info;
}

/** @ingroup Algebra
 *  @brief wrapper of the LAPACK DGESDD routine.
 *  DGESDD computes the Svd decomposition of a matrix.
 *
 *  @verbatim
 *   DGESDD computes the singular value decomposition (SVD) of a real
 *   m-by-n matrix a, optionally computing the left and right singular
 *   vectors.  If singular vectors are desired, it uses a
 *   divide-and-conquer algorithm.
 *
 *   The SVD is written
 *
 *      a = u * SIGMA * transpose(v)
 *
 *   where SIGMA is an m-by-n matrix which is zero except for its
 *   min(m,n) diagonal elements, u is an m-by-m orthogonal matrix, and
 *   v is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
 *   are the singular values of a; they are real and non-negative, and
 *   are returned in descending order. The first min(m,n) columns of
 *   u and v are the left and right singular vectors of a.
 *  @endverbatim
 *
 *  @note routine returns vt = v**T, not v.
 *
 *
 *  Arguments:
 *  ==========
 *
 * @param[in] jobz
 * @verbatim
 *          jobz is Char*1
 *          Specifies options for computing all or part of the matrix u:
 *          = 'A':  all m columns of u and all n rows of v**T are
 *                  returned in the arrays u and vt;
 *          = 'S':  the first min(m,n) columns of u and the first
 *                  min(m,n) rows of v**T are returned in the arrays U
 *                  and vt;
 *          = 'O':  If m >= n, the first n columns of u are overwritten
 *                  on the array a and all rows of v**T are returned in
 *                  the array vt;
 *                  otherwise, all columns of u are returned in the
 *                  array u and the first m rows of v**T are overwritten
 *                  in the array a;
 *          = 'N':  no columns of u or rows of v**T are computed.
 * @endverbatim
 *
 * @param[in] m
 * @verbatim
 *          m is Integer
 *          The number of rows of the input matrix a.  m >= 0.
 * @endverbatim
 *
 * @param[in] n
 * @verbatim
 *          n is Integer
 *          The number of columns of the input matrix a.  n >= 0.
 * @endverbatim
 *
 * @param[in,out] a
 * @verbatim
 *          a is STK::Real array, dimension (lda,n)
 *          On entry, the m-by-n matrix a.
 *          On exit,
 *          if jobz = 'O',  a is overwritten with the first n columns
 *                          of u (the left singular vectors, stored
 *                          columnwise) if m >= n;
 *                          a is overwritten with the first m rows
 *                          of v**T (the right singular vectors, stored
 *                          rowwise) otherwise.
 *          if jobz .ne. 'O', the contents of a are destroyed.
 * @endverbatim
 *
 * @param[in] lda
 * @verbatim
 *          lda is Integer
 *          The leading dimension of the array a.  lda >= max(1,m).
 * @endverbatim
 *
 * @param[out] s
 * @verbatim
 *          s is STK::Real array, dimension (min(m,n))
 *          The singular values of a, sorted so that s[i] >= s[i+1].
 * @endverbatim
 *
 * @param[out] u
 * @verbatim
 *          u is STK::Real array, dimension (ldu,ucol)
 *          ucol = m if jobz = 'A' or jobz = 'O' and m < n;
 *          ucol = min(m,n) if jobz = 'S'.
 *          If jobz = 'A' or jobz = 'O' and m < n, u contains the m-by-m
 *          orthogonal matrix u;
 *          if jobz = 'S', u contains the first min(m,n) columns of u
 *          (the left singular vectors, stored columnwise);
 *          if jobz = 'O' and m >= n, or jobz = 'N', u is not referenced.
 * @endverbatim
 *
 * @param[in] ldu
 * @verbatim
 *          ldu is Integer
 *          The leading dimension of the array U.  ldu >= 1; if
 *          jobz = 'S' or 'A' or jobz = 'O' and m < n, ldu >= m.
 * @endverbatim
 *
 * @param[out] vt
 * @verbatim
 *          vt is STK::Real array, dimension (ldvt,n)
 *          If jobz = 'A' or jobz = 'O' and m >= n, vt contains the
 *          N-by-N orthogonal matrix v**T;
 *          if jobz = 'S', vt contains the first min(m,n) rows of
 *          v**T (the right singular vectors, stored rowwise);
 *          if jobz = 'O' and m < n, or jobz = 'N', vt is not referenced.
 * @endverbatim
 *
 * @param[in] ldvt
 * @verbatim
 *          ldvt is Integer
 *          The leading dimension of the array vt.  ldvt >= 1; if
 *          jobz = 'A' or jobz = 'O' and m >= n, ldvt >= n;
 *          if jobz = 'S', ldvt >= min(m,n).
 * @endverbatim
 *
 * @param[out] work
 * @verbatim
 *          work is STK::Real array, dimension (MAX(1,lWork))
 *          On exit, if info = 0, work(1) returns the optimal lWork;
 * @endverbatim
 *
 * @param[in] lWork
 * @verbatim
 *          lWork is Integer
 *          The dimension of the array work. lWork >= 1.
 *          If jobz = 'N',
 *            lWork >= 3*min(m,n) + max(max(m,n),7*min(m,n)).
 *          If jobz = 'O',
 *            lWork >= 3*min(m,n) +
 *                     max(max(m,n),5*min(m,n)*min(m,n)+4*min(m,n)).
 *          If jobz = 'S' or 'A'
 *            lWork >= min(m,n)*(6+4*min(m,n))+max(m,n)
 *          For good performance, lWork should generally be larger.
 *          If lWork = -1 but other input arguments are legal, work(1)
 *          returns the optimal lWork.
 * @endverbatim
 *
 * @param[out] iWork
 * @verbatim
 *          iWork is Integer array, dimension (8*min(m,n))
 * @endverbatim
 *
 * @return info
 * @verbatim
 *          info is Integer
 *          = 0:  successful exit.
 *          < 0:  if info = -i, the i-th argument had an illegal value.
 *          > 0:  DBDSDC did not converge, updating process failed.
 * @endverbatim
 *
 *  Authors:
 *  ========
 *
 * @author Univ. of Tennessee
 * @author Univ. of California Berkeley
 * @author Univ. of Colorado Denver
 * @author NAG Ltd.
 *
 * Contributors:
 * ============
 *
 *     Ming Gu and Huan Ren, Computer Science Division, University of
 *     California at Berkeley, USA
 *
 **/
inline int gesdd( char jobz, int m, int n, Real *a, int lda
                , Real *s, Real *u, int ldu, Real *vt, int ldvt
                , Real *work, int lWork, int *iWork
                )
{
  int info = 1;
#ifdef STKUSELAPACK
#ifdef STKREALAREFLOAT
  sgesdd_( &jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lWork, iWork, &info);
#else
  dgesdd_( &jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lWork, iWork, &info);
#endif
#endif
  return info;
}

 /** @ingroup Algebra
  * @brief wrapper of the LAPACK SYEVR routine.
  * SYEVR computes the eigenvalues of a symmetric square matrix.
  *
  *  @param[in] jobz
  * @verbatim
  *  CHARACTER*1
  *  = 'N':  Compute eigenvalues only;
  *  = 'V':  Compute eigenvalues and eigenvectors.
  * @endverbatim
  *
  *  @param[in] range
  *  @verbatim
  *  CHARACTER*1
  *  = 'A': all eigenvalues will be found.
  *  = 'V': all eigenvalues in the half-open interval  (VL_,VU_]  will be found.
  *  = 'I': the IL_-th through IU_-th eigenvalues will be found.
  * @endverbatim
  *
  *  @param[in] uplo
  * @verbatim
  * CHARACTER*1
  * = 'U':  Upper triangle of A is stored;
  * = 'L':  Lower triangle of A is stored.
  * @endverbatim
  *
  * @param[in] n The order of the matrix A.  N >= 0.
  *
  * @param[in,out] a Real array, dimension (lda, N)
  * @verbatim
  * On entry, the symmetric matrix A.  If UPLO_ = 'U', the leading
  * N-by-N upper triangular part of A contains the upper  triangular part
  * of the  matrix  A.   If  UPLO_  = 'L', the leading N-by-N lower triangular
  * part of A contains the lower triangular part of the matrix A.  On
  * exit, the lower triangle (if UPLO_='L') or the upper triangle
  * (if UPLO_='U') of A, including the diagonal, is destroyed.
  * @endverbatim
  *
  * @param[in] lda The leading dimension of the array A.  lda >= max(1,N).
  *
  * @param[in] vl,vu
  * @verbatim
  *  Real If RANGE_='V', the lower and  upper  bounds
  *  of  the  interval to be searched for eigenvalues. VL_ < VU_.  Not
  *  referenced if RANGE_ = 'A' or 'I'.
  * @endverbatim
  *
  * @param[in] il, iu
  * @verbatim
  *  Integer If RANGE_='I', the indices (in ascending order)
  *  of the smallest and largest eigenvalues to be returned.
  *  1 <= IL_ <= IU_ <= NL, if NL > 0; IL_ = 1 and IU_ = 0 if NL = 0. Not
  *  referenced if RANGE_ = 'A' or 'V'.
  * @endverbatim
  *
  * @param[in] abstol
  * @verbatim
  *  The  absolute error tolerance for the eigenvalues.  An approximate
  *  eigenvalue is accepted as converged when it is  determined
  *  to lie in an interval [a,b] of width less than or equal to
  *  ABSTOL + EPS *   max( |a|,|b| ) ,
  *  where  EPS is the machine precision.  If ABSTOL is less than or
  *  equal to zero, then  EPS*|T|  will be used in its place,  where
  *  |T|  is the 1-norm of the tridiagonal matrix obtained by reducing A
  *  to tridiagonal form.
  *  If high relative accuracy is important, set ABSTOL  to  SLAMCH(
  *  'Safe minimum' ).  Doing so will guarantee that eigenvalues are
  *  computed to high relative  accuracy  when  possible  in  future
  *  releases.   The current code does not make any guarantees about
  *  high relative accuracy, but future releases will. See J. Barlow
  *  and J. Demmel, "Computing Accurate Eigensystems of Scaled Diagonally
  *  Dominant Matrices", LAPACK Working Note #7, for  a  discussion of
  *  which matrices define their eigenvalues to high relative accuracy.
  * @endverbatim
  *  @see "Computing Small Singular  Values  of  Bidiagonal  Matrices
  *  with  Guaranteed  High Relative Accuracy," by Demmel and Kahan,
  *  LAPACK Working Note #3.
  *
  * @param[out] m
  * @verbatim
  *   The  total number of eigenvalues found.  0 <= M <= NL.  If RANGE_
  *   = 'A', M = NL, and if RANGE_ = 'I', M = IU_-IL_+1.
  * @endverbatim
  *
  * @param[out] w
  * @verbatim
  *  array, dimension (NL)
  *  The first  M  elements  contain  the  selected  eigenvalues  in
  *  ascending order.
  * @endverbatim
  *
  * @param[out] z
  * @verbatim
  *  array, dimension (LDZ, max(1,M))
  *  If  jobz_ = 'V', then if info = 0, the first M columns of Z contain
  *  the orthonormal eigenvectors of the matrix A corresponding
  *  to  the selected eigenvalues, with the i-th column of Z holding
  *  the eigenvector associated with W(i).  If jobz_ = 'N', then Z is
  *  not  referenced.   Note:  the  user  must  ensure that at least
  *  max(1,M) columns are supplied in the array Z; if RANGE_  =  'V',
  *  the exact value of M is not known in advance and an upper bound
  *  must be used.  Supplying N columns is always safe.
  * @endverbatim
  *
  * @param[in] ldz
  * @verbatim
  *  The leading dimension of the array Z.  LDZ >= 1, and if jobz_  =
  *  'V', LDZ >= max(1,N).
  * @endverbatim
  *
  * @param[out] isuppz array, dimension ( 2*max(1,M) )
  * @verbatim
  *  The  support  of the eigenvectors in Z, i.e., the indices indicating
  *  the nonzero elements  in  Z.  The  i-th  eigenvector  is
  *  nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ( 2*i ).
  * @endverbatim
  *
  * @param[in,out] work Real array, dimension (MAX(1,lWork))
  * @verbatim
  *   On exit, if info = 0, work(1) returns the optimal lWork.
  * @endverbatim
  *
  * @param[in] lWork The  dimension  of  the array work
  * @verbatim
  *  lWork >= max(1,26*N). For optimal efficiency, lWork >= (NB+6)*N, where
  *  NB is the max of the  blocksize for SSYTRD and SORMTR returned by ILAENV.
  *  If lWork = -1, then a workspace query is assumed; the routine only
  *  calculates  the  optimal  sizes  of  the work and iWork arrays,
  *  returns these values as the first entries of the work and iWork
  *  arrays,  and  no  error  message  related to lWork or LiWork is
  *  issued by XERBLA.
  * @endverbatim
  *
  * @param[in,out] iwork array, dimension (MAX(1,LiWork))
  * @verbatim
  *  On exit, if info = 0, iWork(1) returns the optimal lWork.
  * @endverbatim
  *
  * @param[in] liwork The dimension of the array iWork.
  * @verbatim
  *  LiWork >=  max(1,10*N). If LiWork  =  -1, then a workspace query is
  *  assumed; the routine only calculates the optimal sizes of the work and
  *  iWork arrays, returns these values as the first entries of the work and
  *  iWork arrays, and no error message related  to  lWork  or  LiWork  is
  *  issued by XERBLA.
  * @endverbatim
  *
  * @return info
  * @verbatim
  *  = 0:  successful exit
  *  < 0:  if info = -i, the i-th argument had an illegal value
  *  > 0:  Internal error
  * @endverbatim
  */
inline int syevr( char jobz, char range, char uplo
                , int n, Real* a, int lda
                , Real vl, Real vu, int il, int iu
                , Real abstol, int *m, Real *w
                , Real *z, int ldz, int *isuppz
                , Real *work, int lWork, int *iwork, int liwork
                )
{
  int info = 0;
#ifdef STKUSELAPACK
#ifdef STKREALAREFLOAT
  ssyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il,
          &iu, &abstol, m, w, z, &ldz, isuppz, work,
          &lWork, iwork, &liwork, &info);
#else
  dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il,
          &iu, &abstol, m, w, z, &ldz, isuppz, work,
          &lWork, iwork, &liwork, &info);
#endif
#endif
  return info;
}

} // namespace lapack

} // namespace STK


#endif /* STK_LAPACK_UTIL_H */
