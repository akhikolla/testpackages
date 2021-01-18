#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/Linpack.h>
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )

using namespace Rcpp;

extern "C" void setulb_(int *n, int *m, double *x, double *l, double *u,
			int *nbd, double *f, double *g, double *factr, double *pgtol,
			double *wa, int *iwa, int *itask, int *iprint,
			int *icsave, int *lsave, int *isave, double *dsave);

typedef double optimfn(int n, double *par, void *ex);

typedef void optimgr(int n, double *par, double *gr, void *ex);

List lbfgsb3Cinfo;

extern "C" void lbfgsb3C_(int n, int lmm, double *x, double *lower,
			  double *upper, int *nbd, double *Fmin, optimfn fn,
			  optimgr gr, int *fail, void *ex, double factr,
			  double pgtol, int *fncount, int *grcount,
			  int maxit, char *msg, int trace, int iprint,
			  double atol, double rtol, double *g){
  // Optim compatible interface
  int itask= 2;
  // *Fmin=;
  double *lastx = new double[n];
  std::copy(&x[0],&x[0]+n,&lastx[0]);
  int nwa = 2*lmm*n + 11*lmm*lmm + 5*n + 8*lmm;
  double *wa= new double[nwa];
  int niwa = 3*n;
  int *iwa= new int[niwa];
  int icsave = 0;
  int lsave[4] = {0};
  int isave[44] = {0};
  int i=0;
  double dsave[29]= {0};
  // Initial setup
  int doExit=0;
  fncount[0]=0;
  grcount[0]=0;
  int itask2=0;
  CharacterVector taskList(28);
  taskList[0]="NEW_X";
  taskList[1]="START";
  taskList[2]="STOP";
  taskList[3]="FG";//,  // 1-4
  taskList[4]="ABNORMAL_TERMINATION_IN_LNSRCH";
  taskList[5]="CONVERGENCE"; //5-6
  taskList[6]="CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL";//7
  taskList[7]="CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH";//8
  taskList[8]="ERROR: FTOL .LT. ZERO"; //9
  taskList[9]="ERROR: GTOL .LT. ZERO";//10
  taskList[10]="ERROR: INITIAL G .GE. ZERO"; //11
  taskList[11]="ERROR: INVALID NBD"; // 12
  taskList[12]="ERROR: N .LE. 0"; // 13
  taskList[13]="ERROR: NO FEASIBLE SOLUTION"; // 14
  taskList[14]="ERROR: STP .GT. STPMAX"; // 15
  taskList[15]="ERROR: STP .LT. STPMIN"; // 16
  taskList[16]="ERROR: STPMAX .LT. STPMIN"; // 17
  taskList[17]="ERROR: STPMIN .LT. ZERO"; // 18
  taskList[18]="ERROR: XTOL .LT. ZERO"; // 19
  taskList[19]="FG_LNSRCH"; // 20
  taskList[20]="FG_START"; // 21
  taskList[21]="RESTART_FROM_LNSRCH"; // 22
  taskList[22]="WARNING: ROUNDING ERRORS PREVENT PROGRESS"; // 23
  taskList[23]="WARNING: STP .eq. STPMAX"; // 24
  taskList[24]="WARNING: STP .eq. STPMIN"; // 25
  taskList[25]="WARNING: XTOL TEST SATISFIED"; //
  taskList[26] = "CONVERGENCE: Parameters differences below xtol";
  taskList[27] = "Maximum number of iterations reached";
  while (true){
    if (trace >= 2){
      Rprintf("\n================================================================================\nBefore call f=%f task number %d, or \"%s\"\n", *Fmin, itask, (as<std::string>(taskList[itask-1])).c_str());
    }
    if (itask==3) doExit=1;
    setulb_(&n, &lmm, x, lower, upper, nbd, Fmin, g, &factr, &pgtol,
	  wa, iwa, &itask, &iprint, &icsave, lsave, isave, dsave);
    if (trace > 2) {
      Rprintf("returned from lbfgsb3 \n");
      Rprintf("returned itask is %d or \"%s\"\n",itask,(as<std::string>(taskList[itask-1])).c_str());
    }
    switch (itask){
    case 4:
    case 20:
    case 21:
      if (trace >= 2) {
	Rprintf("computing f and g at prm=\n");
	NumericVector xv(n);
	std::copy(&x[0],&x[0]+n,&xv[0]);
	print(xv);
      }
      // Calculate f and g
      Fmin[0] = fn(n, x, ex);
      fncount[0]++;
      gr(n, x, g, ex);
      grcount[0]++;
      if (trace > 0) {
	Rprintf("At iteration %d f=%f ", isave[33], *Fmin);
	if (trace > 1) {
	  double tmp = fabs(g[n-1]);
	  for (unsigned int j=n-1; j--;){
	    if (tmp > fabs(g[j])){
	      tmp = fabs(g[j]);
	    }
	  }
	  Rprintf("max(abs(g))=%f",tmp);
	}
	Rprintf("\n");
      }
      break;
    case 1:
      // New x;
      if (maxit <= fncount[0]){
      	itask2=28;
	doExit=1;
      	itask=3; // Stop -- gives the right results and restores gradients
	if (trace > 2){
	  Rprintf("Exit becuase maximum number of function calls %d met.\n", maxit);
	}
      } else {
      	bool converge=fabs(lastx[n-1]-x[n-1]) < fabs(x[n-1])*rtol+atol;
      	if (converge){
      	  for (i=n-1;i--;){
      	    converge=fabs(lastx[i]-x[i]) < fabs(x[i])*rtol+atol;
      	    if  (!converge){
      	      break;
      	    }
      	  }
      	}
      	if (converge){
      	  itask2=27;
      	  itask=3; // Stop -- gives the right results and restores gradients
	  if (trace > 2){
	    Rprintf("CONVERGENCE: Parameters differences below xtol.\n", maxit);
	  }
	  doExit=1;
      	}
      }
      std::copy(&x[0],&x[0]+n,&lastx[0]);
      break;
    default:
      doExit=1;
    }
    if (doExit) break;
  }
  if (itask2){
    itask=itask2;
  }
  LogicalVector lsaveR(4);
  NumericVector dsaveR(29);
  IntegerVector isaveR(44);
  std::copy(&lsave[0],&lsave[0]+4, &lsaveR[0]);
  std::copy(&dsave[0],&dsave[0]+29, dsaveR.begin());
  std::copy(&isave[0],&isave[0]+44, isaveR.begin());;
  CharacterVector taskR(1);
  taskR[0] = taskList[itask-1];
  lbfgsb3Cinfo = List::create(_["task"] = taskR,
			      _["itask"]= IntegerVector::create(itask),
			      _["lsave"]= lsaveR,
			      _["icsave"]= IntegerVector::create(icsave),
			      _["dsave"]= dsaveR,
			      _["isave"] = isaveR);
  // info <- list(task = task, itask = itask, lsave = lsave,
  //      icsave = icsave, dsave = dsave, isave = isave)
  fail[0]= itask;
  delete[] wa;
  delete[] iwa;
  delete[] lastx;
}

Environment grho;

CharacterVector gnames;

List ev;

double gfn(int n, double *x, void *ex){
  Rcpp::NumericVector par(n);
  std::copy(&x[0], &x[0]+n, &par[0]);
  Function fn = as<Function>(ev["fn"]);
  par.attr("names") = ev["pn"];
  double ret = as<double>(fn(par, grho));
  return ret;
}

void ggr(int n, double *x, double *gr, void *ex){
  Rcpp::NumericVector par(n), ret(n);
  std::copy(&x[0], &x[0]+n, &par[0]);
  Function grad = as<Function>(ev["gr"]);
  par.attr("names") = ev["pn"];
  ret = grad(par, grho);
  std::copy(&ret[0], &ret[0]+n, &gr[0]);
}

//[[Rcpp::export]]
Rcpp::List lbfgsb3cpp(NumericVector par, Function fn, Function gr, NumericVector lower, NumericVector upper, List ctrl, Environment rho){
  Rcpp::List ret;
  ev["fn"] = fn;
  ev["gr"] = gr;
  ev["pn"] = par.attr("names");
  Rcpp::NumericVector g(par.size());
    // CONV in 6, 7, 8; ERROR in 9-19; WARN in 23-26
  IntegerVector traceI = as<IntegerVector>(ctrl["trace"]);
  if (traceI.size() != 1) stop("trace has to have one element in it.");
  int trace = traceI[0];
  NumericVector factrN = as<NumericVector>(ctrl["factr"]);
  if (factrN.size() != 1) stop("factr has to have one element in it.");
  double factr = factrN[0];
  NumericVector pgtolN = as<NumericVector>(ctrl["pgtol"]);
  if (pgtolN.size() != 1) stop("pgtol has to have one element in it.");
  double pgtol = pgtolN[0];
  NumericVector atolN = as<NumericVector>(ctrl["abstol"]);
  if (atolN.size() != 1) stop("abstol has to have one element in it.");
  double atol = atolN[0];
  NumericVector rtolN = as<NumericVector>(ctrl["reltol"]);
  if (atolN.size() != 1) stop("reltol has to have one element in it.");
  double rtol = rtolN[0];
  LogicalVector infoN = as<LogicalVector>(ctrl["info"]);
  if (infoN.size() != 1) stop("info has to have one element in it.");
  bool addInfo = infoN[0];
  IntegerVector lmmN = as<IntegerVector>(ctrl["lmm"]);
  if (lmmN.size() != 1) stop("lmm has to have one element in it.");
  int lmm = lmmN.size();
  int n = par.size();
  IntegerVector maxitN = as<IntegerVector>(ctrl["maxit"]);
  if (maxitN.size() != 1) stop("maxit has to have one element in it.");
  int maxit = maxitN[0];
  IntegerVector iprintN = as<IntegerVector>(ctrl["iprint"]);
  if (iprintN.size() != 1) stop("iprint has to have one element in it.");
  int iprint = iprintN[0];
  // double *g = new double[par.size()];
  double *low = new double[par.size()];
  if (lower.size() == 1){
    std::fill_n(&low[0],par.size(),lower[0]);
  } else if (lower.size() == par.size()){
    std::copy(lower.begin(),lower.end(),&low[0]);
  } else {
    delete [] low;
    stop("Lower bound must match the size of par or only have one element.");
  }
  double *up = new double[par.size()];
  if (upper.size() == 1){
    std::fill_n(&up[0],par.size(),upper[0]);
  } else if (upper.size() == par.size()){
    std::copy(upper.begin(),upper.end(),&up[0]);
  } else {
    delete [] low;
    delete [] up;
    stop("Upper bound must match the size of par or only have one element.");
  }
  double *x = new double[par.size()];
  std::copy(par.begin(),par.end(),&x[0]);
  int *nbd = new int[par.size()];
  int i;
  for (i = par.size();i--;){
    /*
	   nbd(i)=0 if x(i) is unbounded,
		  1 if x(i) has only a lower bound,
		  2 if x(i) has both lower and upper bounds,
		  3 if x(i) has only an upper bound.
    */    
    nbd[i] = 0;
    if (R_FINITE(low[i])) nbd[i] = 1;
    if (R_FINITE(up[i]))  nbd[i] = 3 - nbd[i];
  }
  double fmin=std::numeric_limits<double>::max();
  int fail = 0, fncount=0, grcount=0;
  grho=rho;
  //void *ex = (void*)rho; //Should work but use global instead.
  void *ex =NULL;
  char msg[120];
  lbfgsb3C_(n, lmm, x, low, up, nbd, &fmin, gfn, ggr,
	    &fail, ex, factr, pgtol, &fncount,
	    &grcount, maxit, msg, trace, iprint , atol, rtol, &g[0]);
  NumericVector parf(par.size());
  std::copy(&x[0],&x[0]+par.size(),parf.begin());
  parf.attr("names")=ev["pn"];
  g.attr("names")=ev["pn"];
  ret["par"]=parf;
  ret["grad"]=g;
  ret["value"] = fmin;
  IntegerVector cnt = IntegerVector::create(fncount,grcount);
  ret["counts"] = cnt;
  switch (fail){
  case 6:
  case 7:
  case 8:
  case 27:
    ret["convergence"]=0;
    break;
  case 28:
    ret["convergence"]=1;
    break;
  case 23:
  case 24:
  case 25:
  case 26:
    ret["convergence"] = 51;
    break;
  case 9:
  case 10:
  case 11:
  case 12:
  case 13:
  case 14:
  case 15:
  case 16:
  case 17:
  case 18:
  case 19:
    ret["convergence"] = 52;
    break;
  }
  CharacterVector taskList(28);
  taskList[0]="NEW_X";
  taskList[1]="START";
  taskList[2]="STOP";
  taskList[3]="FG";//,  // 1-4
  taskList[4]="ABNORMAL_TERMINATION_IN_LNSRCH";
  taskList[5]="CONVERGENCE"; //5-6
  taskList[6]="CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL";//7
  taskList[7]="CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH";//8
  taskList[8]="ERROR: FTOL .LT. ZERO"; //9
  taskList[9]="ERROR: GTOL .LT. ZERO";//10
  taskList[10]="ERROR: INITIAL G .GE. ZERO"; //11
  taskList[11]="ERROR: INVALID NBD"; // 12
  taskList[12]="ERROR: N .LE. 0"; // 13
  taskList[13]="ERROR: NO FEASIBLE SOLUTION"; // 14
  taskList[14]="ERROR: STP .GT. STPMAX"; // 15
  taskList[15]="ERROR: STP .LT. STPMIN"; // 16
  taskList[16]="ERROR: STPMAX .LT. STPMIN"; // 17
  taskList[17]="ERROR: STPMIN .LT. ZERO"; // 18
  taskList[18]="ERROR: XTOL .LT. ZERO"; // 19
  taskList[19]="FG_LNSRCH"; // 20
  taskList[20]="FG_START"; // 21
  taskList[21]="RESTART_FROM_LNSRCH"; // 22
  taskList[22]="WARNING: ROUNDING ERRORS PREVENT PROGRESS"; // 23
  taskList[23]="WARNING: STP .eq. STPMAX"; // 24
  taskList[24]="WARNING: STP .eq. STPMIN"; // 25
  taskList[25]="WARNING: XTOL TEST SATISFIED"; //
  taskList[26] = "CONVERGENCE: Parameters differences below xtol";
  taskList[27] = "Maximum number of iterations reached";      

  ret["message"]= CharacterVector::create(taskList[fail-1]);
  if (addInfo) ret["info"] = lbfgsb3Cinfo;
  delete [] x;
  delete [] low;
  delete [] up;
  delete [] nbd;
  return ret;
}
