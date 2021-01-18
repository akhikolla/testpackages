//
// constructor and destructor for class Clmbr
//


#include "lmbr.h"




Clmbr::Clmbr(  NumericVector  yR,  NumericMatrix  xR,  NumericMatrix  wR,  int model_num,
						int  inv,  int  var_k )
// constructor
{
	w_in = x_in = y_in = xs = NULL;
	is = NULL;
	rS = irS = Q = tau = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
    f01 = f0x = NULL;
    B = C = NULL;
	px = psig1 = psigx = pv1h = pxh = NULL;
	nan_m1 = pnse1 = pnuse1 = pusen = NULL;
	nan_m = puqe1 = puqen = puqx = NULL;
	ps1 = psx = NULL;
	pq1 = pqx = NULL;
	pmq1 = pm1h = NULL; 
	py = psy = NULL;
	pqy = NULL;


	int i, j;

	n = yR.size();

	xrank = xR.ncol();

	model_in = model_num;

	variance_unknown= !( static_cast<bool>(var_k) );

	inverse= static_cast<bool>(inv);

	bool  cov_matrix_I = true;
	cov_matrix_diagonal = true;

	if( wR(0,0) > -0.5 )  {		// flag for null 'weights'
		if( wR.ncol()==1 )  {
			for (i=0;i<n;i++)  if( fabs( wR(i,0) - 1. ) > zero_eq )  cov_matrix_I = false;
		}  else  {
			for (i=0;i<n;i++) for (j=0;j<n;j++) {
				if (i==j &&  fabs( wR(i,j) - 1. )>zero_eq) cov_matrix_I = false;
				if (i!=j && fabs( wR(i,j) )>zero_eq) {cov_matrix_I = false; cov_matrix_diagonal = false;}
			}
		}
	}

	vectorS = false;
	matrixS = false;
	if( !cov_matrix_I )  {
		if( cov_matrix_diagonal )  vectorS = true;  else  matrixS = true;
	}


	y_in = Calloc( n, double );  
	x_in = Calloc( n*xrank, double ); 
	if( vectorS )  w_in = Calloc( n, double );  
	if( matrixS )  w_in = Calloc( n*n, double );  
 

// store input values

	for (i=0;i<n;i++) {
		y_in[i] =  yR[i];
		for(j=0; j< xrank; j++)  *(x_in+j*n+i) = xR(i,j);
		if( vectorS )  { if( wR.ncol()==1 )  *(w_in+i) = wR(i,0);  else  *(w_in+i) = wR(i,i); }
		if( matrixS )  for(j=0; j< n; j++)  *(w_in+j*n+i) = wR(i,j);
	}



	initialize();
}





Clmbr::Clmbr( const Clmbr  &initM )
//copy constructor
{
	w_in = x_in = y_in = xs = NULL;
	is = NULL;
	rS = irS = Q = tau = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
    f01 = f0x = NULL;
    B = C = NULL;
	px = psig1 = psigx = pv1h = pxh = NULL;
	nan_m1 = pnse1 = pnuse1 = pusen = NULL;
	nan_m = puqe1 = puqen = puqx = NULL;
	ps1 = psx = NULL;
	pq1 = pqx = NULL;
	pmq1 = pm1h = NULL; 
	py = psy = NULL;
	pqy = NULL;

	
	model_in = initM.model_in;

	variance_unknown= initM.variance_unknown;

	inverse= initM.inverse;

	n = initM.n;

	xrank = initM.xrank;

	cov_matrix_diagonal = initM.cov_matrix_diagonal;
	vectorS = initM.vectorS;
	matrixS = initM.matrixS;

	y_in = Calloc( n, double );  
	x_in = Calloc( n*xrank, double ); 
	if( vectorS )  w_in = Calloc( n, double );  
	if( matrixS )  w_in = Calloc( n*n, double );  


// store input values

	for (int i=0;i<n;i++) {
		y_in[i] = (*initM.py)[i];
		for(int j=0; j< xrank; j++)  *(x_in+j*n+i) = *(initM.x_in+j*n+i);
		if( vectorS )  w_in[i] = initM.w_in[i];
		if( matrixS )  for(int j=0; j< n; j++)  *(w_in+j*n+i) = *(initM.w_in+j*n+i);
	}


	initialize();
}




Clmbr::~Clmbr()
// destructor
{
	const Vector<double>  free(0);
	*px = *pv1h = *pxh = *psig1 = *psigx = *nan_m1 = *nan_m = free;
	*pnse1 = *pnuse1 = *pusen = *puqe1 = *puqen = *puqx = free;
	*py = *psy = *pqy = free;
	for(int i=0; i<ns+1; i++) {
		ps1[i]= free;  psx[i]= free;  pq1[i]= free;  pqx[i]= free;  pmq1[i]= free;
	}
	if(Model==M3)  *pm1h = free;

	Free( w_in );  Free( x_in );  Free( y_in );  Free( xs );
	Free( px );
	Free( rS );  Free( irS );  Free( Q ); Free( tau );
	Free( is );
	Free( q11 );  Free( qx1 );  Free( qxx );  Free( ck );  Free( qff );
	Free( q10 );  Free( qx0 );  Free( a0 );  Free( b0 );
	Free( f01 );  Free( f0x );
	Free( B );  Free( C );
	Free( psig1 );  Free( psigx );  Free( pv1h );  Free( pxh );
	Free( nan_m1 );  Free( pnse1 );  Free( pnuse1 );  Free( pusen );
	Free( nan_m );  Free( puqe1 );  Free( puqen );  Free( puqx );
	Free( ps1 );  Free( psx );
	Free( pq1 );  Free( pqx );
	Free( pmq1 ); Free( pm1h );
	Free( py );  Free( psy );
	Free( pqy );

}


