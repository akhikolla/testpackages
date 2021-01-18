//
//

#include "lmbr.h"



void Clmbr::initialize( void )
//  initialize variables 
//  precalculate working variables
{

	if(model_in==1)  Model= M1;  else 
		if(model_in==2 || model_in== -2)  Model= M2;			  // treat Model -2 as M2 internally
			else  if(model_in== 3 || model_in== -3)  Model= M3;	  // treat Model -3 as M3 internally

	if(Model==M1)  { m= n-2-(xrank-2);  m1= m+2;  k1=  1; } 
	if(Model==M2)  { m= n-1-(xrank-2);  m1= m+1;  k1=  0; }
	if(Model==M3)  { m= n-(xrank-1);    m1= m;    k1= -1; }



	px = Calloc( 1, Vector<double> );
	set_x();


//  'S' = 'Sigma'  is the covariate or covariance matrix,  such that
//   errors  ~ N( 0, var * Sigma ) ,   Sigma = inverse( 'weights' )

	if( vectorS || matrixS )  {

		if(vectorS)  {
			rS = Calloc( n, double );
			irS = Calloc( n, double );
		}  else  {
			rS = Calloc( n*n, double );
			irS = Calloc( n*n, double );
		}

		set_Sigma();
	}


//  The computations for a multiple-regression model are the same as for a simple, 
//  except that the bottom 'm1' rows of Q*inv-root-Sigma replaces inv-root-Sigma.
//  So, if multiple, set the 'cov_matrix' flag to non-diagonal to invoke 
//  more general routines, even if 'weights' is the default identity matrix.
 
	if( m1 < n )  cov_matrix_diagonal = false;


	Q = Calloc( n*xrank, double );
	tau = Calloc( xrank, double );

	set_Q();


// allocate memory for pre-calculated variables
	pv1h = Calloc( 1, Vector<double> );
	pxh = Calloc( 1, Vector<double> );
	psig1 = Calloc( 1, Vector<double> );
	psigx = Calloc( 1, Vector<double> );
	nan_m1 = Calloc( 1, Vector<double> );
	nan_m = Calloc( 1, Vector<double> );
	pnse1 = Calloc( 1, Vector<double> );
	pnuse1 = Calloc( 1, Vector<double> );
	pusen = Calloc( 1, Vector<double> );
	puqe1 = Calloc( 1, Vector<double> );
	puqen = Calloc( 1, Vector<double> );
	puqx = Calloc( 1, Vector<double> );

// array sizes that depend on 'ns'
	is = Calloc( ns, int );
	xs = Calloc( ns,  double );
	ps1 = Calloc( ns+1, Vector<double> );
	psx = Calloc( ns+1, Vector<double> );
	pq1 = Calloc( ns+1, Vector<double> );
	pqx = Calloc( ns+1, Vector<double> );
	pmq1 =  Calloc( ns+1, Vector<double> );
	if(Model==M3)  pm1h =  Calloc( 1, Vector<double> );
	q11 = Calloc( ns+1, double );
	qx1 = Calloc( ns+1, double );
	qxx = Calloc( ns+1, double );
	ck = Calloc( ns+1, double );
	qff = Calloc( ns+1, double );

	pre_calc();


	py = Calloc( 1, Vector<double> );
	psy = Calloc( 1, Vector<double> );
	pqy = Calloc( 1, Vector<double> );

	set_y();	// calls 'set_sy'

	sety_called = false;


	q10 = Calloc( ns+1, double );
	qx0 = Calloc( ns+1, double );
	a0 = Calloc( ns+1, double );
	b0 = Calloc( ns+1, double );
	f01 = Calloc( ns+1, double );
	f0x = Calloc( ns+1, double );
	B = Calloc( ns+1, double );

	const double th_0 = xs[1];
	th0 = th_0 + 1;
	th0MC = xs[ns-1] + 1;
	set_theta0(th_0, INIT);		// initializes 'z', 'w'

	const double a0 = (*py)[ is[1] ];
	alpha0 = a0 + 1;
	set_alpha0(a0, INIT);

	prev_SL= -1; 
	set_SL();

	set_tol();

	C = Calloc( 3, double );
	C[0]= get_C(m-2); C[1]= get_C(m-1); C[2]= get_C(m); 


	old_th = prev_th =  Inf; 
	ah = a_low = a_high =  0; 


	return;
}

