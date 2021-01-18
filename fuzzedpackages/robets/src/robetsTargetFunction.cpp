#include <cmath>
#include <R.h>
#include "robetsTargetFunction.h"

#include <R_ext/Print.h>

const int NONE = 0;
const int ADD = 1;
const int MULT = 2;
const int DAMPED = 1;
const double TOL = 1.0e-10;
const double HUGEN = 1.0e10;
const double NA  = -99999.0;

const double LAMBDA_SIGMA = 0.1;

void RobetsTargetFunction::init(std::vector<double> & p_y, int p_errortype,
		int p_trendtype, int p_seasontype, bool p_damped,
		std::vector<double> & p_lower, std::vector<double> & p_upper, std::string p_opt_crit,
		int p_nmse, std::string p_bounds, int p_m,
		bool p_optAlpha, bool p_optBeta, bool p_optGamma, bool p_optPhi, bool p_optSigma0 ,bool p_optInit, bool p_optK,
		bool p_givenAlpha, bool p_givenBeta, bool p_givenGamma, bool p_givenPhi, bool p_givenSigma0, bool p_givenInit, bool p_givenK,
		double alpha, double beta, double gamma, double phi, double p_sigma0, std::vector<double> & p_initstate, double k) {

	this->y = p_y;
	this->n = this->y.size();
  this->nstate = 2 + ((p_trendtype != NONE) ? 1 : 0) + ((p_seasontype != NONE) ? 1 : 0)*p_m;
	this->errortype = p_errortype;

	this->trendtype = p_trendtype;
	this->seasontype = p_seasontype;
	this->damped = p_damped;

	this->lower = p_lower;
	this->upper = p_upper;

	this->opt_crit = p_opt_crit;
	this->nmse = p_nmse;
	this->bounds = p_bounds;

	this->m = p_m;

	this->optAlpha = p_optAlpha;
	this->optBeta = p_optBeta;
	this->optGamma = p_optGamma;
	this->optPhi = p_optPhi;

	this->givenAlpha = p_givenAlpha;
	this->givenBeta = p_givenBeta;
	this->givenGamma = p_givenGamma;
	this->givenPhi = p_givenPhi;

	this->alpha = alpha;
	this->beta = beta;
	this->gamma = gamma;
	this->phi = phi;
  
  this->sigma0 = p_sigma0;
  this->initstate = p_initstate;
  if(initstate.size() < (unsigned) nstate ){
    if(seasontype != NONE) {
  		double sum=0;
  		for(int i=0;i<(m-1);i++) {
  			sum += initstate[(1+((trendtype != NONE) ? 1 : 0)) + i];
  		}
  
  		double new_state = m*((seasontype == MULT) ? 1 : 0) - sum;
  
  		initstate.push_back(new_state);
  	}
  }
  
  this->k = k;
  
  this->optSigma0 = p_optSigma0;
  this->optInit = p_optInit;
	this->optK = p_optK;
  
  this->givenSigma0 = p_givenSigma0;
  this->givenInit = p_givenInit;
  this->givenK = p_givenK;

	this->lik = 0;
	this->objval = 0;

	//	for(int i=0; i < 10; i++) this->amse.push_back(0);
	//	for(int i=0; i < n; i++) this->e.push_back(0);
	this->amse.resize(30, 0);
	this->e.resize(n, 0);

}

void RobetsTargetFunction::eval(const double* p_par, int p_par_length) {

	bool equal=true;
  
	// Check if the parameter configuration has changed, if not, just return.
	if((unsigned)p_par_length != this->par.size()) {
		equal=false;
	} else {
		for(int j=0;j < p_par_length;j++) {
			if(p_par[j] != this->par[j]) {
				equal=false;
				break;
			}
		}
	}

	if(equal) return;

	this->par.clear();

	for(int j=0;j < p_par_length;j++) {
		this->par.push_back(p_par[j]);
	}

	int j=0;
	if(optAlpha) this->alpha = par[j++];
	if(optBeta) this->beta = par[j++];
	if(optGamma) this->gamma = par[j++];
	if(optPhi) this->phi = par[j++];

	if(!this->check_params()) {
		this->objval = R_PosInf;
		return;
	}

	this->state.clear();
  if(optSigma0) this->state.push_back(par[j++]);
  if(givenSigma0) this->state.push_back(sigma0);
  
  if(optInit){ 
    for(int i=(optAlpha+optBeta+optGamma+optPhi+optSigma0); i < (p_par_length-optK); i++) {
      this->state.push_back(par[i]);
    }
  } 
  
  if(givenInit){
    for(unsigned i=0; i < initstate.size(); i++) {
      this->state.push_back(initstate[i]);
    }
  }
  
  if(state.size() < (unsigned) nstate ){
    if(seasontype != NONE) {
      double sum=0;
    
    	for(unsigned i=(2+((trendtype != NONE) ? 1 : 0));i<state.size();i++) {
    		sum += state[i];
    	}
    
    	double new_state = m*((seasontype == MULT) ? 1 : 0) - sum;
    
  		state.push_back(new_state);
    }
  }
  
  if(optK){
    this->k = par[p_par_length-1];
  }	

	for(int i=0; i < nstate*n; i++) state.push_back(0);

	robetscalc();


	// Avoid perfect fits
	if (this->lik < -1e10) this->lik = -1e10;

	// isnan() is a C99 function
	//if (isnan(this->lik)) this->lik = 1e8;
	if (ISNAN(this->lik)) this->lik = R_PosInf;

	if(fabs(this->lik+99999) < 1e-7) this->lik = R_PosInf;

	if(this->opt_crit=="lik") this->objval = this->lik;
	else if(this->opt_crit=="mse") this->objval = this->amse[0];
	else if(this->opt_crit=="amse") {

		//return(mean(e$amse[1:nmse]))
		double mean=0;
		for(int i=0;i < this->nmse;i++) {
			mean+=amse[i]/this->nmse;
		}
		this->objval=mean;

	}
	else if(this->opt_crit=="sigma") {
		//return(mean(e$e^2))
		double mean=0;
		int ne=e.size();
		for(int i=0;i<ne;i++) {
			mean+=e[i]*e[i]/ne;
		}
		this->objval=mean;

	}
	else if (this->opt_crit=="mae") {
		//return(mean(abs(e$e)))

		double mean=0;
		int ne=e.size();
		for(int i=0;i<ne;i++) {
			mean+=fabs(e[i])/ne;
		}
		this->objval=mean;

	}
  else if (this->opt_crit=="roblik") {
		this->objval=roblik;
	}
  else if (this->opt_crit=="tau2") {
  	this->objval=tau2;
	}
  
    //	---------show params----------
		//Rprintf("par: ");
		//for(int j=0;j < p_par_length;j++) {
		//	Rprintf("%f ", p_par[j]);
		//}
		//Rprintf(" objval: %f | ", this->objval);
	  //	---------show params----------
  
}

bool RobetsTargetFunction::check_params() {

	if(bounds != "admissible")
	{
		if(optAlpha)
		{
			if(alpha < lower[0] || alpha > upper[0])
				return(false);
		}
		if(optBeta)
		{
			if(beta < lower[1] || beta > alpha || beta > upper[1])
				return(false);
		}
		if(optPhi)
		{
			if(phi < lower[3] || phi > upper[3])
				return(false);
		}
		if(optGamma)
		{
			if(gamma < lower[2] || gamma > 1-alpha || gamma > upper[2])
				return(false);
		}
	}
	if(bounds != "usual")
	{
		if(!admissible()) return(false);
	}
	return(TRUE);

}

bool RobetsTargetFunction::admissible() {

	if(phi < 0 || phi > 1+1e-8) return(false);

	//If gamma was set by the user or it is optimized, the bounds need to be enforced
	if(!optGamma && !givenGamma) {
		if(alpha < 1-1/phi || alpha > 1+1/phi) return(false);

		if(optBeta || givenBeta)
		{
			if(beta < alpha * (phi-1) || beta > (1+phi)*(2-alpha)) return(false);
		}
	}

	else if(m > 1) //Seasonal model
	{

		if(!optBeta && !givenBeta) beta = 0;


		//max(1-1/phi-alpha,0)
		double d = 1-1/phi-alpha;
		if(gamma < ((d > 0) ? d : 0) || gamma > 1+1/phi-alpha) return(false);

		if(alpha < 1-1/phi-gamma*(1-m+phi+phi*m)/(2*phi*m)) return(false);

		if(beta < -(1-phi)*(gamma/m+alpha)) return(false);

		// End of easy tests. Now use characteristic equation

		std::vector<double> opr;
		opr.push_back(1);
		opr.push_back(alpha+beta-phi);

		for(int i=0;i<m-2;i++) opr.push_back(alpha+beta-alpha*phi);

		opr.push_back(alpha+beta-alpha*phi+gamma-1);
		opr.push_back(phi*(1-alpha-gamma));

		int degree = opr.size()-1;

		std::vector<double> opi;
		opi.resize(opr.size(),0);

		std::vector<double> zeror(degree);
		std::vector<double> zeroi(degree);

		Rboolean fail;

		cpolyroot(&opr[0], &opi[0], &degree, &zeror[0], &zeroi[0], &fail);

		double max = 0;
		for(unsigned i=0;i<zeror.size();i++) {
		  double abs_val = std::sqrt(zeror[i]*zeror[i] + zeroi[i]*zeroi[i]);
		  if(abs_val>max) max = abs_val;
		}

		//Rprintf("maxpolyroot: %f\n", max);

		if(max > 1+1e-10) return(false);

		// P <- c(phi*(1-alpha-gamma),alpha+beta-alpha*phi+gamma-1,rep(alpha+beta-alpha*phi,m-2),(alpha+beta-phi),1)
		// roots <- polyroot(P)
		// if(max(abs(roots)) > 1+1e-10) return(false);
	}

	//Passed all tests
	return(true);
}

void RobetsTargetFunction::robetscalc(){
  int i, j, maxvalue;
  double oldsigma, sigma, oldl, l, oldb, b, olds[24], s[24], f[30], lik2, tmp, ydown;
  std::vector<double> absYpred;
    
  if((m > 24) & (seasontype > NONE)){
    return;
  }else if(m < 1){
    m = 1;
  }

  if(nmse > 30){
    nmse = 30;
  }
  
  // Copy initial state components
  sigma = state[0];
  l = state[1];
  if(trendtype > NONE){
    b = state[2];
  }
        
  if(seasontype > NONE){
    for(j=0; j<m; j++){
      s[j] = state[(trendtype>NONE)+j+2];
    }
  }
  lik = 0.0;
  lik2 = 0.0;
  
  
  for(j=0; j<nmse; j++){
    amse[j] = 0.0;
  }

  for (i=0; i<n; i++){
        // COPY PREVIOUS STATE
        oldsigma = sigma;
        oldl = l;
        if(trendtype > NONE)
            oldb = b;
        if(seasontype > NONE)
        {
            for(j=0; j<m; j++)
                olds[j] = s[j];
        }

        maxvalue = m; //std::max(nmse,m);
        // ONE STEP FORECAST
        forecast(oldl, oldb, olds, f, maxvalue);
        
        absYpred.push_back(std::abs(f[0]));
        
        if(fabs(f[0]-NA) < TOL){
            lik = NA;
            return;
        }
        
        if(errortype == ADD)
            e[i] = y[i] - f[0];
        else
            e[i] = (y[i]- f[0])/f[0];
        for(j=0; j<nmse; j++)
        {
            if(i+j<n)
            {
                tmp = y[i+j]-f[j];
                amse[j] += (tmp*tmp)/n;
            }
        } 
        
          // ROBUSTNESS STEP
        sigma = std::sqrt( LAMBDA_SIGMA*rhobiweight(e[i]/oldsigma)*oldsigma*oldsigma + (1.0-LAMBDA_SIGMA)*oldsigma*oldsigma); 
        if(errortype==ADD){
          ydown = sigma*psihuber(e[i]/sigma) + f[0];
        }else{
          ydown = f[0] * sigma*psihuber(e[i]/sigma) + f[0];
        }
        
        //if(f[0] < TOL){
        //    lik = HUGE;
        //    roblik = HUGE;
        //    return;
        //}
        
        // UPDATE STATE
        update(&oldl, &l, &oldb, &b, &olds[0], &s[0], ydown);

        // STORE NEW STATE
        state[nstate*(i+1)] = sigma;
        state[nstate*(i+1)+1] = l;
        if(trendtype > NONE)
            state[nstate*(i+1)+2] = b;
        if(seasontype > NONE)
        {
           for(j=0; j<m; j++){
             state[(trendtype>NONE)+nstate*(i+1)+j+2] = s[j];
           }
          
        }
        lik = lik + e[i]*e[i];
        lik2 += log(fabs(f[0]));
    }
    
    
    //char c[10];
    //sprintf(c,"%lf",lik);
    //Rprintf("%s ",c);
    lik = n * log(lik);
    //sprintf(c,"%lf",lik);
    //Rprintf("%s ",c);
    if(errortype == MULT)
      lik += 2*lik2;
    //sprintf(c,"%lf",lik2);
    //Rprintf("%s | ",c);
    
    if(errortype == ADD){
      tau2 = computeTau2(e);
      roblik = n*log(n*tau2); 
    }else{ // errortype = MULT
      tau2 = computeTau2(e); 
      roblik = n*log(n*tau2)+2*n*log(median(absYpred));  //  n*log(n*tau2)+2*lik2;  // 
    }

}

void RobetsTargetFunction::forecast(double& l,  double& b, double *s, double *f, int& h)
{
    int i,j;
    double phistar;

    phistar = phi;

    // FORECASTS
    for(i=0; i<h; i++)
    {
        if(trendtype == NONE)
            f[i] = l;
        else if(trendtype == ADD)
            f[i] = l + phistar*b;
        else if(b<0)
            f[i] = NA;
        else
            f[i] = l * pow(b,phistar);
        j = m-1-i;
        while(j < 0)
             j += m;
        if(seasontype == ADD)
            f[i] = f[i] + s[j];
        else if(seasontype == MULT)
            f[i] = f[i] * s[j];
        if(i < (h-1))
        {
            if(fabs(phi-1.0) < TOL)
                phistar = phistar + 1.0;
            else
                phistar = phistar + pow(phi, (double) (i+1));
        }
    }
}

void RobetsTargetFunction::update(double* oldl, double* l, double* oldb, double* b, double *olds, double *s, double& y)
{
    int j;
    double q, p, r, t, phib;

    // NEW LEVEL
    if(trendtype==NONE)
  {
        q = *oldl;                 // l(t-1)
  	phib = 0;
	}
    else if(trendtype==ADD)
    {
        phib = phi*(*oldb);
        q = *oldl + phib;          // l(t-1) + phi*b(t-1)
    }
    else if(fabs(phi-1.0) < TOL)
    {
        phib = *oldb;
        q = *oldl * (*oldb);       // l(t-1)*b(t-1)
    }
    else
    {
        phib = pow(*oldb,phi);
        q = (*oldl) * phib;          // l(t-1)*b(t-1)^phi
    }
    if(seasontype==NONE)
        p = y;
    else if(seasontype==ADD)
        p = y - olds[m-1];         // y[t] - s[t-m]
    else
    {
        if(fabs(olds[m-1]) < TOL)
            p = HUGEN;
        else
            p = y / olds[m-1];     // y[t]/s[t-m]
    }
    *l = q + alpha*(p-q);

    // NEW GROWTH
    if(trendtype > NONE)
    {
        if(trendtype==ADD)
           r = (*l) - (*oldl);       // l[t]-l[t-1]
        else //if(trendtype==MULT)
        {
            if(fabs(*oldl) < TOL)
                r = HUGEN;
            else
                r = (*l)/(*oldl);    // l[t]/l[t-1]
        }
        *b = phib + (beta/alpha)*(r - phib);   // b[t] = phi*b[t-1] + beta*(r - phi*b[t-1])
                                               // b[t] = b[t-1]^phi + beta*(r - b[t-1]^phi)
    }

    // NEW SEASON
    if(seasontype > NONE)
    {
        if(seasontype==ADD)
            t = y - q;
        else //if(seasontype==MULT)
        {
            if(fabs(q) < TOL)
                t = HUGEN;
            else
                t = y / q;
        }
        s[0] = olds[m-1] + gamma*(t - olds[m-1]); // s[t] = s[t-m] + gamma*(t - s[t-m])
        for(j=1; j<m; j++)
            s[j] = olds[j-1];                     // s[t] = s[t]
    }
}

void RobetsTargetFunction::oneEval(std::vector<double> & p_y, int p_errortype,
  	int p_trendtype, int p_seasontype, bool p_damped,
		int p_nmse, int p_m,
		double alpha, double beta, double gamma, double phi, std::vector<double> & initstate, double k){
    this->y = p_y;
    this->n = this->y.size();
    this->nstate = 2 + ((p_trendtype != NONE) ? 1 : 0) + ((p_seasontype != NONE) ? 1 : 0)*p_m ;
    this->errortype = p_errortype;
      
    this->trendtype = p_trendtype;
    this->seasontype = p_seasontype;
    this->damped = p_damped;
      
    this->nmse = p_nmse;
      
    this->m = p_m;
      
    this->alpha = alpha;
    this->beta = beta;
    this->gamma = gamma;
    this->phi = phi;
    
    this->k = k;
      
    this->lik = 0;
        
    this->amse.resize(30, 0);
    this->e.resize(n, 0);

	this->state.clear();

	for(unsigned i=0; i < initstate.size(); i++) {
    this->state.push_back(initstate[i]);
	}
  
  if(state.size() < (unsigned)nstate ){
    if(seasontype != NONE) {
      double sum=0;
  
    	for(unsigned i=(2+((trendtype != NONE) ? 1 : 0));i<state.size();i++) {
    		sum += state[i];
    	}
    
    	double new_state = m*((seasontype == MULT) ? 1 : 0) - sum;
    
    	state.push_back(new_state);
    }
  }

  for(unsigned i=0; i < nstate*this->y.size(); i++) state.push_back(0);    
  robetscalc();
}

double RobetsTargetFunction::rhobiweight(double x){
  double rho = 1.0;
  if(std::abs(x)<k){
    rho = 1.0-pow((1.0-pow(x/k,2)),3);
  }
  return rho/Erho() ;
}

double RobetsTargetFunction::Erho(){
  double dnormk = R::dnorm(k, 0.0, 1.0, 0); //exp(-k*k/2.0)/sqrt(2.0*M_PI);
  double pnormk = R::pnorm(k, 0.0, 1.0, 1, 0); //0.5*erfc(-k/sqrt(2));
  double d1 = pnormk-0.5-k*dnormk;
  double d2 = 3.0*d1-pow(k,3)*dnormk;
  double d3 = 5.0*d2-pow(k,5)*dnormk;
  double Erho = (6.0/pow(k,2))*d1 - (6.0/pow(k,4))*d2 + (2.0/pow(k,6))*d3 + 2.0*(1.0-pnormk);
  return Erho;
}

double RobetsTargetFunction::psihuber(double x){
  if(std::abs(x)<k){
    return x;
  }else{
    return copysign(k,x);
  }
}

double RobetsTargetFunction::median(std::vector<double> x) {
  std::size_t n = x.size() / 2;
  std::nth_element(x.begin(), x.begin() + n, x.end());

  if (x.size() % 2) return x[n]; 
  return (x[n] + *std::max_element(x.begin(), x.begin() + n)) / 2.;
}

double RobetsTargetFunction::computeTau2(std::vector<double>& x){
  // compute Erho
  const double kt = 3.0;//2.0; 
  const double dnormk = R::dnorm(kt, 0.0, 1.0, 0); //exp(-kt*kt/2.0)/sqrt(2.0*M_PI);
  const double pnormk = R::pnorm(kt, 0.0, 1.0, 1, 0); //0.5*erfc(-kt/sqrt(2.0));
  const double d1 = pnormk-0.5-kt*dnormk;
  const double d2 = 3.0*d1-pow(kt,3)*dnormk;
  const double d3 = 5.0*d2-pow(kt,5)*dnormk;
  const double Erho = (6.0/pow(kt,2))*d1 - (6.0/pow(kt,4))*d2 + (2.0/pow(kt,6))*d3 + 2.0*(1.0-pnormk);
  
  std::vector<double> x2;
  for(unsigned i = 0 ; i<x.size() ; ++i){
    x2.push_back(x[i]*x[i]);
  }
  double sc = 1.482602*sqrt(median(x2));
  double tauscale = 0.0;
  double rho;
  for(unsigned i = 0 ; i < x.size() ; ++i){
    rho = 1.0;
    if(std::abs(x[i]/sc)<kt){
      rho = 1.0-pow((1.0-pow((x[i]/sc)/kt,2)),3);
    }
    tauscale = tauscale + rho/Erho;
  }
  return (1.0/x.size())*pow(sc,2)*tauscale;
}
