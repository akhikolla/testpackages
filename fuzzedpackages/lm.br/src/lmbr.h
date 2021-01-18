//
//  define class  Clmbr  for  data-and-model objects 
//  which are applications of the broken-line model to observed data sets
//


#if !defined  CLMBR_H__			//prevents compiler from repeating this code in other files
#define  CLMBR_H__


#include  "globals.h"


class Clmbr
{

private:

	MODEL Model;
	int model_in, n, m, m1, ns, xrank;
	bool variance_unknown, inverse, vectorS, matrixS;
	bool cov_matrix_diagonal, trivial, sety_called;

	bool th0ex;
	int k1, k0, subints;
	double s11, sx1, sxx, n1, se1sq, y1, yx, sysq, qysq, omega, c, th0, alpha0, z, w;
	double prime_z, g0u1, g0u2, c1, c2, Lgamma, lambdasq, lambda, th0MC, th0a0;

	double ah, old_th, prev_th, a_low, a_high, rel_print_eps;
	double SL, prev_SL, cFex, cCHIex, cF, cCHI, x_vu_ex, x_vk_ex, x_vu, x_vk;
	double tol_rho, tol_sl_abs, tol_sl_rel, tol_xb, tol_yb, inc_x, inc_y, xinc;

	double *x_in, *y_in, *w_in;
	double *rS, *irS, *Q, *tau;
	int *is;
	double *xs;
	double *q11, *qx1, *qxx, *ck, *qff;
	double *q10, *qx0, *a0, *b0;
	double *f01, *f0x;
	double *B, *C;

	Vector<double> *px, *psig1, *psigx, *pv1h, *pxh;
	Vector<double> *nan_m1, *pnse1, *pnuse1, *pusen;
	Vector<double> *nan_m, *puqe1, *puqen, *puqx;
	Vector<double> *ps1, *psx;
	Vector<double> *pq1, *pqx;
	Vector<double> *pmq1; 
	Vector<double> *pm1h; 

	Vector<double> *py, *psy;
	Vector<double> *pqy;



// function prototypes:

//initializing functions
	void initialize( void );
	void set_x( void );
	void set_Sigma( void );
	void set_Q( void );
	void pre_calc( void );
	void set_y( void );
	void set_theta0( double th_0, METHOD met =INIT );
	void set_alpha0( double a_0, METHOD met =INIT );
	void set_SL( double cSL =0.05);
	void set_tol( double tol =0.001);

// gamma and f functions
	Vector<double>  gam(double th, int data_interval) const;
	Vector<double>  gfr(double th, int data_interval) const;
	Vector<double>  gsm(double th, int data_interval) const;
	Vector<double>  gbar(double th, int data_interval) const;
	Vector<double>  gbar_prime(double th, int data_interval) const;
	Vector<double>  q_f(double th, int data_interval) const;
	Vector<double>  sf(double th, int data_interval) const;
	Vector<double>  sfc(double th, int data_interval) const;
	double  ff(double th, int data_interval) const;

//in file 'rho_etc':
	double rho(double th) const;
	double rho(double theta, int data_interval) const;
	double rhosq(double th, int data_interval) const;
	double drho(double th, int data_interval) const;
	double drhosq(double th, int data_interval) const;
	double dgsq(double th, int data_interval) const;
	double rho_inv(double s, int data_interval, int hi_lo =1) const;

//in file 'fm_Fm'
	double get_C(int k) const;
	double fk(int k, double arg) const;
	double F(int k, double arg) const;
	double sF(int k, double arg) const;

//in file 'bisect'
	double bisect(double x1, double x2, double (Clmbr::*fn)(double,int), int k, double value, double crit);
	double bisect(double x1, double x2, double (Clmbr::*fn)(double,int) const, int k, double value, double crit) const;
	double bisect_sl(double x1, double x2, METHOD met, double crit);

//in files  'geo'  and  'geo_ex'
	double geo(double th2, double *err) const;
	double geo_vu_D(double th2, double *err) const;
	double geo_vu_ND(double th2, double *err) const;
	double geo_vu_NDab(int data_interval, double th_a, double th_b, int hi_lo, double *err) const;
	double geo_vk_D(double th2, double *err) const;
	double geo_vk_ND(double th2, double *err) const;
	double geo_vk_NDab(int data_interval, double th_a, double th_b, int hi_lo, double *err) const;
	double geo_ex(void) const;
	double geo_vu_ex(void) const;
	double geo_vk_ex(void) const;

//in file 'geo_i'
	double amu_by_Omega(double rho_, int data_interval) const;
	double Emupr(double theta, int data_interval) const;
	double Emupr_vk(double theta, int data_interval) const;

// sl computing functions
	double sl_af(int mode =0) const;
	double sl_af2(void) const;
	double sl_geo(double *err=0);
	double sl_geo2(double *err=0);
	double prden(double xi, double *err);
	double sl_mc(void) const;
	double sl_mc2(void) const;
	bool m_ge_w(double wsq, const Vector<double> &s) const;

//'ci' and 'cr' sub-subroutines
	int ci_geo( METHOD met, double increments, double *bounds =0);
	int ci_af( METHOD met, double *bounds =0);
	double a_sl(METHOD met, double th, int high_low);
	double a_af(double th, int high_low);
	double ahigh(METHOD met, double th);
	double sl_a(double alpha, int k);

//private engines for interface functions
	double sl(double theta0, METHOD met = GEO, bool verbose =true);
	double sl(double theta0, double alpha0, METHOD met = GEO, bool verbose =true);
	int ci(METHOD met =GEO, double increments = -1, bool verbose =true, double *bounds =0);
	int cr(METHOD met =GEO, double increments = -1, bool verbose =true, double *bounds =0);
	double mle( bool verbose =true, double *max_gamma_dot_Qy_sq =NULL, double *param = NULL ) const;
	void set_sy(double *irsy, METHOD met =INIT);


public:
// constructors
	Clmbr( NumericVector  yR,  NumericMatrix  xR,  NumericMatrix  wR,  int  model_num,
			 int  inverse,  int  var_known );
	Clmbr( const Clmbr  &initM );	// copy constructor
	~Clmbr();		// destructor


// interface functions
	void sl3R( int met, double tol, double theta0 );
	void sl4R( int met, double tol, double theta0, double alpha0 );
	double sl5R( int met, int verboseR, int valueR, double tol, 
				double theta0 );
	double sl6R( int met, int verboseR, int valueR, double tol, 
				double theta0, double alpha0 ); 
	void ciR( double CL, int met );
	void cr3R( double CL, int met, double incr );
	NumericMatrix cr4R( double CL, int met, double incr, int verboseR );
	void MLE(void) const;
	NumericVector PARAM(void) const;	// internal
	void SET_rWy( NumericVector rWy );


// friend functions, integrands to 'Rdqag' routines in files 'geo' and 'sl_geo'
	friend void igeo(double *x, const int n, void *const ex);
	friend void igeo2(double *x, const int n, void *const ex);

};



#endif


