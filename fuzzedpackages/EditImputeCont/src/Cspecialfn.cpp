#include "CHeader.h"  // this code uses c++ library "newmat11" and "Newran03" by R B Davies
#define LOG_2_PI 1.83787706640935

double rbeta_fn( double a, double b ){
  	//Tracer tr("rbeta_fn");
    Gamma randGamma_a(a);
    double temp_a = randGamma_a.Next();

    Gamma randGamma_b(b);
    double temp_b = randGamma_b.Next();

    double out = ( temp_a ) / ( temp_a + temp_b );
    return out;
}

double log_beta_fn( double x, double a, double b ){
		//Tracer tr("log_beta_fn");
    double out1 = ln_gamma( a + b ) - ln_gamma( a ) - ln_gamma( b ) ;
    double out2 = (a-1.0) * log(x) + (b-1.0) * log(1-x) ;
    double out = out1 + out2 ;
    return out;
}


ReturnMatrix Submatrix_elem(SymmetricMatrix &Sigma, ColumnVector &id_rows){
  //Tracer tr("Submatrix_elem");
    // ex. id_rows (1,0,1) , id_columns (1,1,1,0)
    int i;
    int n_row = Sigma.nrows();
    int n_sub_row = id_rows.sum();
    int *index = new int[n_sub_row];
    int count = 0;
    for (i = 1; i <= n_row; i++) {
      if (id_rows(i) > 0) {
        index[count++] = i;
      }
    }
    SymmetricMatrix Y_sub(n_sub_row);
    for (i =1; i <= n_sub_row; i++) {
      for (int j = i; j <= n_sub_row; j++) {
        Y_sub(i,j) = Sigma(index[i-1],index[j-1]);
      }
    }
    delete [] index;
    Y_sub.Release();return Y_sub;
}

// ADDED 01/27/2015
ReturnMatrix Submatrix_elem_2( SymmetricMatrix &Sigma, ColumnVector &id_rows, ColumnVector &id_cols ){
  Tracer tr("Submatrix_elem_2");
  // ex. id_rows (1,0,1) , id_cols (1,1,1,0)

  int n_row = Sigma.nrows() ; int n_col = n_row ;
  int i;
  int n_sub_row = id_rows.sum() ; int n_sub_col = id_cols.sum() ;
  int *index_row = new int[n_sub_row]; int *index_col = new int[n_sub_col];
  int count_row = 0 ;
  for (i = 1; i <= n_row; i++) {
    if (id_rows(i) > 0) {
      index_row[count_row++] = i;
    }
  }
  int count_col = 0 ;
  for (i = 1; i <= n_col; i++) {
    if (id_cols(i) > 0) {
      index_col[count_col++] = i;
    }
  }
  Matrix Sigma_temp(n_row,n_col) ; Sigma_temp << Sigma ;
  for (i =1; i <= n_row; i++) {
      for (int j = i; j <= n_col; j++) {
        Sigma_temp(i,j) = Sigma_temp(j,i);
      }
  }
  Matrix Y_sub(n_sub_row,n_sub_col) ;
  for (i =1; i <= n_sub_row; i++) {
      for (int j = 1; j <= n_sub_col; j++) {
        Y_sub(i,j) = Sigma_temp(index_row[i-1],index_col[j-1]);
      }
  }

  delete [] index_row;
  delete [] index_col;
  Y_sub.Release(); Sigma_temp.Release();
  return Y_sub;
}
// ADDED 01/27/2015

ReturnMatrix log_ColumnVector( ColumnVector &X ){
  //Tracer tr("log_ColumnVector");
	ColumnVector logX = X ;
	for (int i=1; i<=X.nrows(); i++){
		logX(i) = log(X(i)) ;
	}
	logX.Release();return logX ;
}

ReturnMatrix whichOnes( ColumnVector &ind_vector ){
  // Elements of ind_vector have values 0 or 1
	ColumnVector index = ind_vector;
  int count = 0;
  for (int i=1; i<=ind_vector.nrows(); i++){
      if (ind_vector(i)==1){
			  index(++count) = i;
      }
  }
  index = index.rows(1,count);
  index.release(); return index;
}

ReturnMatrix whichZeros( ColumnVector &ind_vector ){
  ColumnVector index = ind_vector;
  int count = 0;
  for (int i=1; i<=ind_vector.nrows(); i++){
      if (ind_vector(i)==0){
  		  index(++count) = i;
      }
  }
  index = index.rows(1,count);
  index.release(); return index;
}

ColumnVector subvector_by_index(ColumnVector &Y, ColumnVector &Index) {
    int n_var_sub = Index.nrows();
    ColumnVector Y_sub(n_var_sub) ;
    for (int i = 1; i <= n_var_sub; i++) {
      Y_sub(i) = Y(Index(i));
    }
    return Y_sub;
}
ColumnVector subvector(ColumnVector &Y, ColumnVector &S ){
    // ex. S = (1,0,0,1)
    int n_var_sub = S.sum();
    ColumnVector Y_sub(n_var_sub) ;
    for (int i_var=1, count = 1; count <=n_var_sub; i_var++) {
        if (S(i_var)==1){
            Y_sub(count) = Y(i_var);
            count++;
        }
    }
    return Y_sub;
}

ReturnMatrix exp_ColumnVector( ColumnVector &X){
	ColumnVector expX = X;
	for (int i=1; i<=X.nrows(); i++){
		expX(i) = exp(X(i));
	}
	expX.release(); return expX;
}

ReturnMatrix rMVN_fn( ColumnVector &mu, LowerTriangularMatrix &LSigma ){
  // LowerTriangularMatrix LSigma = Cholesky(Sigma);
  //Tracer tr("rMVN_fn");
  int n_var = LSigma.nrows();
  Normal randNorm;
  ColumnVector z(n_var);
  for (int j=1; j<=n_var; j++){
      z(j)=randNorm.Next();
  }
  ColumnVector out = mu + (LSigma * z);
  out.release(); return out;
}

int rdiscrete_fn(ColumnVector &Prob){
  //Tracer tr("rdiscrete_fn");
  int n_comp = Prob.nrows();
  double* Prob_in = new double[n_comp];
  for (int k=1; k<=n_comp; k++) Prob_in[k-1]=Prob(k);
  DiscreteGen D(n_comp,Prob_in);
  int out = D.Next(); out = out + 1;
	delete [] Prob_in;
  return out;
}

int runifdiscrete_fn(int n_var ){
  //Tracer tr("runifdiscrete_fn");
  ColumnVector temp_prob(n_var) ; temp_prob = 1.0 / n_var ;
	int out = rdiscrete_fn(temp_prob) ;
  return out;
}

//////
// Note: use eigenvalue decoposition is less precise than Cholesky
//       adjust n.p.d. covariacne matrix
//////
ReturnMatrix rIW_w_pd_check_fn( int nu, LowerTriangularMatrix &LPhi ){
    Tracer tr("rIW_w_pd_check_fn");
    // LowerTriangularMatrix LPhi = Cholesky(Phi);  // error if argument is not SymmetricMatrix
    int n_var = LPhi.nrows();
    LowerTriangularMatrix invLPhi = LPhi.i();

    Normal randNorm;
    Matrix X(nu,n_var);
    ColumnVector z_i(n_var), x_i(n_var);

    for (int i=1; i<=nu; i++){
        for (int j=1; j<=n_var; j++){
            z_i(j)=randNorm.Next();
        }
        x_i = invLPhi.t() * z_i;
        X.row(i) = x_i.t() ;
    }
    Matrix A = X.t()*X ;
		SymmetricMatrix A_inv(n_var) ; A_inv = 0.0 ;
		A_inv << A.i() ;

		// added to avoid non-positive definite matrix
		DiagonalMatrix D(n_var) ; Matrix V(n_var,n_var) ;
		Jacobi(A_inv,D,V) ;

		int is_zero_exist = 0 ;
		for (int i_var=1; i_var<=n_var; i_var++){
			if ( D(i_var) < 1e-7 ){
				D(i_var) = 1e-7 ;
				is_zero_exist = 1 ;
			}
		} // for (int i_var=1; i_var<=n_var; i_var++)

		if ( is_zero_exist == 1 ){
			A_inv << V * D * V.t() ;
		}
		// added to avoid non-positive definite matrix

    LowerTriangularMatrix out = Cholesky(A_inv) ;
    out.release(); return out;
}

double rgamma_fn( double a, double b ){
  	//Tracer tr("rgamma_fn");
    Gamma randGamma(a);
    double out = randGamma.Next();
    out = out / b ;
    return out;
}

double logdet(LowerTriangularMatrix& lchol) {
  double log_det_sigma = 0;
  int dim = lchol.Ncols();
	for (int i = 1; i <= dim; i++) {
		log_det_sigma += log(lchol(i,i));
	}
	return log_det_sigma;
}

double log_MVN_fn(ColumnVector &x, ColumnVector &mu, LowerTriangularMatrix &LSigma_i) {
  double logdetandmore = -0.5*mu.nrows()* LOG_2_PI + logdet(LSigma_i);
  return log_MVN_fn(x, mu, LSigma_i, logdetandmore);
}
double log_MVN_fn(ColumnVector &X, ColumnVector &Mu, LowerTriangularMatrix &LSigma_i, double logdetandmore){
    int i,j;
    int dim = Mu.nrows();
    double *x = X.data();
    double *mu = Mu.data();
    double *xx = new double[dim];
  	for (i = 0; i < dim; i++) {
  		xx[i] = x[i] - mu[i];
  	}
  	double discrim = 0;
  	double *s = LSigma_i.data();
  	for ( i = 0; i < dim; i++) {
  		double sum = 0;
  		for (j = 0; j <= i; j++) {
  			sum += *s++ * xx[j];
  		}
  		discrim += sum * sum;
  	}
  	delete [] xx;
    return  logdetandmore - 0.5*discrim;
    /*
    ColumnVector xi = LSigma_i * ( x - mu );
    Matrix temp2 = xi.t() * xi ;
    return  logdetandmore - 0.5*temp2.as_scalar();
    */
}
