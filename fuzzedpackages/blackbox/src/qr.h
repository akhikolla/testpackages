#ifndef H_QR
#define H_QR

#include <R_ext/Linpack.h>
#include "RcppEigen.h"
#include "R.h"


/** we use QR decomp to obtain (a part of) the Q matrix Q_2 that is combined with
covTypedef arguments to yield the matrix whose eigensystem will be used. That means that
QR output should be covTypedef; and the QR method (Householder transformations) is numericallly stable, so it should not
require a higher precision than covTypedef

Thus QRtype = covTypedef should be OK. But it is constrained to be double.

**/


template <typename QRtype>
class CQR {


   /** Array for internal storage of decomposition.
   @serial internal array storage.
   */

   QRtype* a; //linpack qr in vector form
   int nrow,ncol;
   bool allocated;

   int lda,k; //linpack variables... typically nrow and ncol
   QRtype* qraux; //linpack qraux
   QRtype* y; //linpack y
   QRtype* qy; //will contain output Q.y
   QRtype* qty; //will contain output Q'.y
   QRtype* b; //will contain solution of least square problem y=Q.b+ error
   /////////////////////std::vector<covTypedef> result; // outer world (covType) copy of qy, qty or b
   /////////////////////std::vector<std::vector<covTypedef> >resultm; // outer world copy of several qy, qty or b computations
   QRtype* rsd; //would contain solution of other problem, dim not known
   QRtype* ab; //would contain solution of other problem, dim not known
   int* jpvt; //dqrdc arg
   QRtype* work; //dqrdc arg
public:
	CQR() {allocated=false;};
//dqrdc
/**
	Create a QR factorization object for A.

	@param A rectangular (m>=n) matrix.
*/
	void reinit(std::vector<std::vector<covTypedef> > &A) {
	  if (allocated && ( ! (nrow==int(A.size()) && ncol==int(A[0].size())))) deallocate();
	  if ( ! allocated) {
         nrow=lda=int(A.size());//nb lignes
         ncol=k=int(A[0].size());//nb cols

      ///////////////////////      result.resize(nrow,0);
      /////////////////////////////resultm.resize(nrow);

         a=new QRtype[nrow*ncol];// linpack's vector format for the qr variable
         qraux=new QRtype[ncol];
         jpvt=new int[ncol];
         work=new QRtype[ncol];

         y=new QRtype[nrow];
         qy=new QRtype[nrow];
         qty=new QRtype[nrow];
         b=new QRtype[nrow];
         allocated=true;
	} /** else {
         if ( ! (nrow==A.size() && ncol==A[0].size())) {
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From CQR(<matrix>): trying to reuse a CQR instantiation for a matrix of different size.";
           if (batchDebug) std::cin.get();
           exit(-1);
#else
        Rf_error("(!) From CQR(<matrix>): trying to reuse a CQR instantiation for a matrix of different size.\n");
#endif
         }
	} **/

      for (int jj=0;jj<ncol;jj++) jpvt[jj]=jj;  // inoperant with only one column.
      for (int jj=0;jj<nrow;jj++) {
          b[jj]=qty[jj]=qy[jj]=0; //this only so that valgrind does not complain about uninitialized values
      }
      for (int ii=0;ii<nrow;ii++) {
         for (int jj=0;jj<ncol;jj++) {
               a[jj*nrow+ii]=A[ii][jj]; // successive elements are different row for a given column
         }
      }
      dqrdc ( 0 ); // computes linpack's qr and qraux as in R
   }
   void deallocate() {
      if (allocated) {
          delete[] a;// linpack's vector format for the qr variable
          delete[] qraux;
          delete[] jpvt;
          delete[] work;
          delete[] y;
          delete[] qy;
          delete[] qty;
          delete[] b;
          allocated=false;
      }
   }

	CQR(std::vector<std::vector<covTypedef> > &A) {		/* constructor */
	    allocated=false;
	    this->reinit(A);
	};
   ~CQR() {
      deallocate();
   };

   int dqrdc (int job=0) {// cFR->FR temporary wrapper to F77_CALL(dqrdc); to be replaced by Rcpp::QR
    F77_CALL(dqrdc)(a, &lda, &nrow, &ncol, qraux, jpvt, work, &job);
    //  Parameters:
    //
    //    Input/output, QRtype A(LDA,P).  On input, the N by P matrix
    //    whose decomposition is to be computed.  On output, A contains in
    //    its upper triangle the upper triangular matrix R of the QR
    //    factorization.  Below its diagonal A contains information from
    //    which the orthogonal part of the decomposition can be recovered.
    //    Note that if pivoting has been requested, the decomposition is not that
    //    of the original matrix A but that of A with its columns permuted
    //    as described by JPVT.
/// A is stored in columns in a !
    //
    //    Input, int LDA, the leading dimension of the array A.  LDA must
    //    be at least N.
    //
    //    Input, int N, the number of rows of the matrix A.
    //
    //    Input, int P, the number of columns of the matrix A.
    //
    //    Output, QRtype QRAUX[P], contains further information required
    //    to recover the orthogonal part of the decomposition.
    //
    //    Input/output, integer JPVT[P].  On input, JPVT contains integers that
    //    control the selection of the pivot columns.  The K-th column A(*,K) of A
    //    is placed in one of three classes according to the value of JPVT(K).
    //      > 0, then A(K) is an initial column.
    //      = 0, then A(K) is a free column.
    //      < 0, then A(K) is a final column.
    //    Before the decomposition is computed, initial columns are moved to
    //    the beginning of the array A and final columns to the end.  Both
    //    initial and final columns are frozen in place during the computation
    //    and only free columns are moved.  At the K-th stage of the
    //    reduction, if A(*,K) is occupied by a free column it is interchanged
    //    with the free column of largest reduced norm.  JPVT is not referenced
    //    if JOB == 0.  On output, JPVT(K) contains the index of the column of the
    //    original matrix that has been interchanged into the K-th column, if
    //    pivoting was requested.
    //
    //    Workspace, QRtype WORK[P].  WORK is not referenced if JOB == 0.
    //
    //    Input, int JOB, initiates column pivoting.
    //    0, no pivoting is done.
    //    nonzero, pivoting is done.
    //
      return 0;
    } // end dqrdc (in template)

    double logabsdet() {
        double lad=0;
//std::cout<<std::endl<<">";
        for (unsigned int ii=0;ii<std::min(ncol,nrow);ii++) {
            lad+=log(std::abs(a[ii*(nrow+1)])); /// A stored by columns in a !
//std::cout<<a[ii*(nrow+1)]<<" ";
        }
//std::cout<<"<"<<std::endl;
        return(lad);
    }

    int dqrsl (int job )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    DQRSL computes transformations, projections, and least squares solutions.
    //
    //  Discussion:
    //
    //    DQRSL requires the output of DQRDC.
    //
    //    For K <= min(N,P), let AK be the matrix
    //
    //      AK = ( A(JPVT[0]), A(JPVT(2)), ..., A(JPVT(K)) )
    //
    //    formed from columns JPVT[0], ..., JPVT(K) of the original
    //    N by P matrix A that was input to DQRDC.  If no pivoting was
    //    done, AK consists of the first K columns of A in their
    //    original order.  DQRDC produces a factored orthogonal matrix Q
    //    and an upper triangular matrix R such that
    //
    //      AK = Q * (R)
    //               (0)
    //
    //    This information is contained in coded form in the arrays
    //    A and QRAUX.
    //
    //    The parameters QY, QTY, B, RSD, and AB are not referenced
    //    if their computation is not requested and in this case
    //    can be replaced by dummy variables in the calling program.
    //    To save storage, the user may in some cases use the same
    //    array for different parameters in the calling sequence.  A
    //    frequently occuring example is when one wishes to compute
    //    any of B, RSD, or AB and does not need Y or QTY.  In this
    //    case one may identify Y, QTY, and one of B, RSD, or AB, while
    //    providing separate arrays for anything else that is to be
    //    computed.
    //
    //    Thus the calling sequence
    //
    //      dqrsl ( a, lda, n, k, qraux, y, dum, y, b, y, dum, 110, info )
    //
    //    will result in the computation of B and RSD, with RSD
    //    overwriting Y.  More generally, each item in the following
    //    list contains groups of permissible identifications for
    //    a single calling sequence.
    //
    //      1. (Y,QTY,B) (RSD) (AB) (QY)
    //
    //      2. (Y,QTY,RSD) (B) (AB) (QY)
    //
    //      3. (Y,QTY,AB) (B) (RSD) (QY)
    //
    //      4. (Y,QY) (QTY,B) (RSD) (AB)
    //
    //      5. (Y,QY) (QTY,RSD) (B) (AB)
    //
    //      6. (Y,QY) (QTY,AB) (B) (RSD)
    //
    //    In any group the value returned in the array allocated to
    //    the group corresponds to the last member of the group.
    //
    //  Modified:
    //
    //    07 June 2005
    //
    //  Author:
    //
    //    C++ translation by John Burkardt.
    //
    //  Reference:
    //
    //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
    //    LINPACK User's Guide,
    //    SIAM, (Society for Industrial and Applied Mathematics),
    //    3600 University City Science Center,
    //    Philadelphia, PA, 19104-2688.
    //    ISBN 0-89871-172-X
    //
    //  Parameters:
    //
    //    Input, QRtype A[LDA*P], contains the output of DQRDC.
    //
    //    Input, int LDA, the leading dimension of the array A.
    //
    //    Input, int N, the number of rows of the matrix AK.  It must
    //    have the same value as N in DQRDC.
    //
    //    Input, int K, the number of columns of the matrix AK.  K
    //    must not be greater than min(N,P), where P is the same as in the
    //    calling sequence to DQRDC.
    //
    //    Input, QRtype QRAUX[P], the auxiliary output from DQRDC.
    //
    //    Input, QRtype Y[N], a vector to be manipulated by DQRSL.
    //
    //    Output, QRtype QY[N], contains Q * Y, if requested.
    //
    //    Output, QRtype QTY[N], contains Q' * Y, if requested.
    //
    //    Output, QRtype B[K], the solution of the least squares problem
    //      minimize norm2 ( Y - AK * B),
//qr.coef
    //    if its computation has been requested.  Note that if pivoting was
    //    requested in DQRDC, the J-th component of B will be associated with
    //    column JPVT(J) of the original matrix A that was input into DQRDC.
    //
    //    Output, QRtype RSD[N], the least squares residual Y - AK * B,
    //    if its computation has been requested.  RSD is also the orthogonal
    //    projection of Y onto the orthogonal complement of the column space
    //    of AK.
    //
    //    Output, QRtype AB[N], the least squares approximation Ak * B,
//(FR): le produit A.B et non B. That is, qr.fitted et non qr.coef
    //    if its computation has been requested.  AB is also the orthogonal
    //    projection of Y onto the column space of A.
    //
    //    Input, integer JOB, specifies what is to be computed.  JOB has
    //    the decimal expansion ABCDE, with the following meaning:
    //
    //      if A != 0, compute QY.
    //      if B != 0, compute QTY.
    //      if C != 0, compute QTY and B.
    //      if D != 0, compute QTY and RSD.
    //      if E != 0, compute QTY and AB.
    //
    //    Note that a request to compute B, RSD, or AB automatically triggers
    //    the computation of QTY, for which an array must be provided in the
    //    calling sequence.
    //
    //    Output, int DQRSL, is zero unless the computation of B has
    //    been requested and R is exactly singular.  In this case, INFO is the
    //    index of the first zero diagonal element of R, and B is left unaltered.
    //
    {
      int info=0;
      //      dqrsl ( a, lda, n, k, qraux, y, dum, y, b, y, dum, 110, info )
      F77_CALL(dqrsl)(a, &lda, &nrow, &k, qraux, y, //inputs
                      qy, qty , b, rsd, ab,  //outputs
                      &job, &info);
      return info;
    } // end dqrsl
    /// VECTOR version of Qy...
    template <typename inputType,typename outputType>
    int Qy(std::vector<inputType>& yv,std::vector<outputType>& locresult) { //stands for column vector // '&' yv is useful
        if (int(yv.size())!=nrow) {
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From qy(): y vector of wrong size";
           if (batchDebug) std::cin.get();
           exit(-1);
#else
           Rf_error("(!) From qy(): y vector of wrong size\n");
#endif
        }
        for (int ii=0;ii<nrow;ii++) y[ii]=yv[ii];
        dqrsl(10000); //qy
        locresult.resize(0);
        for (int ii=0;ii<nrow;ii++) locresult.push_back(qy[ii]); //qy is a C array...
        return(0);
    }
    /// MATRIX version of Qy...
    template <typename inputType,typename outputType>
    int QY(std::vector<std::vector<inputType> >& ym,std::vector<std::vector<outputType> >& locresult) { //stands for column vector // '&' yv is useful
        if (int(ym.size())!=nrow) {
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From QY(): y matrix of wrong size";
           if (batchDebug) std::cin.get();
           exit(-1);
#else
           Rf_error("(!) From QY(): y matrix of wrong size\n");
#endif
        }
        int ncoly=ym[0].size();
////////////////        for (int ii=0;ii<nrow;ii++) resultm[ii].resize(ncoly);
        //filling by column makes it uneasy to use same container as input and output
        locresult.resize(nrow);
        for ( int ii=0;ii<nrow;ii++) locresult[ii].resize(0);
if (true) {  /// version based on calling dqrsl on each column
        for ( int col=0;col<ncoly;col++) {
            for ( int ii=0;ii<nrow;ii++) y[ii]=ym[ii][col]; //input column
            dqrsl(10000); //qy
            for ( int ii=0;ii<nrow;ii++) locresult[ii].push_back(qy[ii]); //output as column
        }
} else { /// trying a more direct code... slower...
    QRtype* it;
    QRtype* jt;
    int j;
    int jj;
    QRtype t;
    QRtype temp;
    int ju = std::min ( k, nrow-1 );
    for ( int col=0;col<ncoly;col++) { // loop over columns of y matrix
        for ( int ii=0;ii<nrow;ii++) qy[ii]=ym[ii][col]; //input column
        for ( jj = 1; jj <= ju; jj++ ) {
          j = ju - jj + 1;

          if ( qraux[j-1] != 0.0 ) {
            temp = a[j-1+(j-1)*lda];
            a[j-1+(j-1)*lda] = qraux[j-1];
                it=a+j-1+(j-1)*lda;
                jt=qy+j-1;
                t=0;
                /// ddot:
                for(;it!=a+j-1+(j-1)*lda+nrow-j+1;it++,jt++) {
                    t -=(*it)*(*jt);
                };
                t /= a[j-1+(j-1)*lda];
                it=a+j-1+(j-1)*lda;
                jt=qy+j-1;
                /// daxpy:
                for(;it!=a+j-1+(j-1)*lda+nrow-j+1;it++,jt++) {
                    (*jt)+=t*(*it);
                } /// so that the code for the first value of jj changes a series of values in qy...
            a[j-1+(j-1)*lda] = temp; // restore original a value
          }
        }
        for ( int ii=0;ii<nrow;ii++) locresult[ii].push_back(qy[ii]); //output as column
    }
}
        return(0);
    }
    template <typename outputType>
    int getR(std::vector<std::vector<outputType> >& locresult) {
        /// may be valid only in nrow>=ncol...
        locresult.resize(ncol);
        for ( int ii=0;ii<ncol;ii++) {
            locresult[ii].resize(ncol,0);
            for ( int jj=ii;jj<ncol;jj++) {
                locresult[ii][jj]=a[ii+jj*lda]; // stockage en colonne...
            }
        }
        return(0);
    }
    template <typename outputType>
    int getQ(std::vector<std::vector<outputType> >& locresult) {
        /// may be valid only in nrow>=ncol...
        locresult.resize(nrow);
        for ( int ii=0;ii<nrow;ii++) locresult[ii].resize(0);
        int j;
        int jj;
        QRtype t;
        QRtype temp;
        int ju = std::min ( ncol, nrow-1 );
        for ( int col=0;col<ncol;col++) {
            for ( int ii=0;ii<nrow;ii++) qy[ii]=0;
            qy[col]=1;
            for ( jj = 1; jj <= ju; jj++ ) {
              j = ju - jj+1;

              if ( qraux[j-1] != 0.0 ) {
                temp = a[j-1+(j-1)*lda];
                a[j-1+(j-1)*lda] = qraux[j-1];
                if (j-2<int(col)) { // otherwise t=0...
                   t = -ddot ( nrow-j+1, a+j-1+(j-1)*lda, 1, qy+j-1, 1 ) / a[j-1+(j-1)*lda];
                   daxpy ( nrow-j+1, t, a+j-1+(j-1)*lda, 1, qy+j-1, 1 );  // changes qy from qy[j-1] onwards
                }
                a[j-1+(j-1)*lda] = temp; // restore original a value
              }
            }
            for ( int ii=0;ii<nrow;ii++) locresult[ii].push_back(qy[ii]); //output as column
        }
        return(0);
    }
    /// VECTOR version of Qty...
    template <typename inputType,typename outputType>
    void Qty(std::vector<outputType> yv) { //stands for column vector
        if (int(yv.size())!=nrow) {
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From Qty(): y vector of wrong size";
           if (batchDebug) std::cin.get();
           exit(-1);
#else
           Rf_error("(!) From Qty(): y vector of wrong size\n");
#endif
        }
        for ( int ii=0;ii<nrow;ii++) y[ii]=yv[ii];
        dqrsl(1000); //qty
        for ( int ii=0;ii<nrow;ii++) {
            yv[ii]=qty[ii];
        }
     }
    /// MATRIX version of Qty...
    template <typename Type>
    void QtY(std::vector<std::vector<Type> >& ym) {
        if (int(ym.size())!=nrow) {
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From Qty(): y matrix of wrong size";
           if (batchDebug) std::cin.get();
           exit(-1);
#else
           Rf_error("(!) From Qty(): y matrix of wrong size\n");
#endif
        }
         int ncoly=ym[0].size();
////////////////        for (int ii=0;ii<nrow;ii++) resultm[ii].resize(ncoly);
        //filling by column makes it uneasy to use same container as input and output
        for ( int col=0;col<ncoly;col++) {
            for ( int ii=0;ii<nrow;ii++) y[ii]=ym[ii][col]; //input column
            dqrsl(1000); //qty
            for ( int ii=0;ii<nrow;ii++) ym[ii][col]=qty[ii]; //output column
        }
//for (int ii=0;ii<nrow;ii++) {ostream_vector(resultm[ii],std::cout);getchar();}
////////////////        return resultm;
    }
    template <typename inputType,typename outputType>
    Eigen::MatrixXd Qtyt(std::vector<std::vector<inputType> >& ym) {
        // input x rows -> output x cols
        /** note transposition input / output;  hence if input and output containers are identical,
        one cannot write successive columns in output; hence anotther container must be defined and copied **/
         int ncoly=ym.size();
        if (int(ym[0].size())!=nrow) {
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From Qtyt(): y matrix of wrong size";
           if (batchDebug) std::cin.get();
           exit(-1);
#else
           Rf_error("(!) From Qtyt(): y matrix of wrong size\n");
#endif
        }
        Eigen::MatrixXd resultm(nrow,ncoly);
        for ( int col=0;col<ncoly;col++) {
            for ( int ii=0;ii<nrow;ii++) y[ii]=ym[col][ii]; //input row
//for (int ii=0;ii<nrow;ii++) {std::cout<<y[ii]<<" ";}
            dqrsl(1000); //qty
            for ( int ii=0;ii<nrow;ii++) resultm(ii,col)=qty[ii]; //output column
        }
//for (int ii=0;ii<nrow;ii++) {ostream_vector(resultm[ii],std::cout);getchar();}
        return resultm;
    }
    template <typename inputType,typename outputType>
    std::vector<std::vector<inputType> > Qyt(std::vector<std::vector<outputType> >& ym) {
        // input x rows -> output x cols
        int ncoly=ym.size();
        if (int(ym[0].size())!=nrow) {
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From Qyt(): y matrix of wrong size";
           if (batchDebug) std::cin.get();
           exit(-1);
#else
           Rf_error("(!) From Qyt(): y matrix of wrong size\n");
#endif
        }
        std::vector<std::vector<outputType> > resultm(0);
        resultm.resize(nrow);
        for ( int ii=0;ii<nrow;ii++) resultm[ii].resize(ncoly);
        for ( int col=0;col<ncoly;col++) {
            for ( int ii=0;ii<nrow;ii++) y[ii]=ym[col][ii]; //input row
//for (int ii=0;ii<nrow;ii++) {std::cout<<y[ii]<<" ";}
            dqrsl(10000); //computes Qyt using qty return after calling qy computation...
            for ( int ii=0;ii<nrow;ii++) resultm[ii][col]=qty[ii]; //output column
        }
//for (int ii=0;ii<nrow;ii++) {ostream_vector(resultm[ii],std::cout);getchar();}
        return resultm;
    }
    template <typename typeforES> /// serves to compute the coef_fixed !
    void coef(std::vector<typeforES>& yv) { //stands for column vector  /// here pass by reference is necessary as this returns void!
        if (int(yv.size())!=nrow) {
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From coef(): y vector of wrong size";
           if (batchDebug) std::cin.get();
           exit(-1);
#else
           Rf_error("(!) From coef(): y vector of wrong size\n");
#endif
        }
        for ( int ii=0;ii<nrow;ii++) y[ii]=yv[ii];
        dqrsl(100); //b
        for ( int ii=0;ii<nrow;ii++) {
            yv[ii]=b[ii];
        }
    }

}; //end TEMPLATE class CQR

#endif
// QR__H

