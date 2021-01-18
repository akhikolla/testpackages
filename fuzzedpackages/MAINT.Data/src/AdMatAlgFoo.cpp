#include "AdMatAlgFoo.h"
#include <math.h>

bool pdsolve(const mat& M, mat& MInv, double* logDet)
{
  mat MSr,MSrInv;

  if (!chol(MSr,M))  return false;  // Maybe replace this by a better rank-reaveling criteron
  if(!inv(MSrInv,trimatu(MSr))) return false;
  MInv = MSrInv.t() * MSrInv;   
  if (logDet) {
    *logDet = -log(MInv(0,0));
    for (int i=1;i<M.n_rows;i++) {
      *logDet -= log(MInv(i,i));
    }
  }	
  return true;
}

bool pdsolve(const mat& M, mat& rhs, mat& res, double* logDet)  // Try replacing pointer to logDet by a constant passed by reference
{
  double DetSign;

  if (!solve(res,M,rhs,solve_opts::likely_sympd)) return false;
  if (logDet) log_det(*logDet,DetSign,M);
  return true; 
}

void SetZero(vec& v,const int n,bool cheksize)
{ 
  if (cheksize && v.n_elem!=n) v.zeros(n); 
  else v.zeros();
}


void SetZero(mat& M,const int m,const int n,bool cheksize)
{ 
  if ( cheksize && (M.n_rows!=m || M.n_cols!=n) ) M.zeros(m,n); 
  else M.zeros();
}

bool chcksing(const mat& M, double& logDet, double& viol, double& minegv, double& maxegv,
  const double minlndet, const double maxlnk2, const bool scale)
{
  int p(M.n_rows);

  static mat M1;
  static vec diagsrqt;
  const mat* Mpnt;

  for (int i=0;i<p;i++) for (int j=0;j<p;j++) 
    if (!std::isfinite(M(i,j))) return false;  

  if (scale) {  
    if (M1.n_rows!=p|| M1.n_cols!=p) M1.set_size(p,p);
    if (diagsrqt.n_elem!=p) diagsrqt.set_size(p);
    for (int i=0;i<p;++i) {
      if (M(i,i) <= 0.) return false;
      diagsrqt(i) = sqrt(M(i,i));
    }
    for (int i=0;i<p;++i) {
      M1(i,i) = 1.;
      for (int j=0;j<i;++j) M1(i,j) = M1(j,i) = M(i,j)/(diagsrqt(i)*diagsrqt(j));
    }
    Mpnt = &M1;
  } else {
    Mpnt = &M;
  }   

  double DetSign;
  log_det(logDet,DetSign,*Mpnt);
  if (logDet < minlndet) {
    viol = minlndet-logDet;
    return false;
  }  
  double tr = trace(*Mpnt);
  if (tr < std::numeric_limits<double>::min()) {
    viol = std::numeric_limits<double>::infinity();
    return false;
  }  

//  double logtrace = log(tr);
//  if (p*logtrace-logDet > maxlnk2) {    // Things to do: reactive this test when eigenvalues do not have to be returned (use eigenvalues pointers with and test if they are NULL)
  vec egval = eig_sym(*Mpnt); 
  minegv = egval(0);    
  maxegv = egval(p-1);  
  if (minegv < MINEGVAL) {
    viol = std::numeric_limits<double>::infinity();
    return false;
  }  
  double logk2 = log(maxegv/minegv); 
  if (logk2 > maxlnk2) {
    viol = logk2 - maxlnk2;
    return false;
  }  

  if (scale) for (int i=0;i<p;++i) logDet += 2*log(diagsrqt(i));

  return true;
}

bool safepdsolve(const mat& M, mat& MInv, double& logDet, double& viol, double& minegv, double& maxegv, 
  const double minlndet, const double maxlnk2, const bool scale) 
{

  int p(M.n_rows);
  static mat M1,Mscl;
  static vec diagsrqt;
  const mat* Mpnt;
  if (scale) {  
    if (M1.n_rows!=p || M1.n_cols!=p) M1.set_size(p,p);
    if (Mscl.n_rows!=p || Mscl.n_cols!=p) Mscl.set_size(p,p);
    if (diagsrqt.n_elem!=p) diagsrqt.set_size(p);
    for (int i=0;i<p;++i) diagsrqt(i) = sqrt(M(i,i));
    for (int i=0;i<p;++i) for (int j=0;j<=i;++j) Mscl(i,j) = Mscl(j,i) = diagsrqt(i)*diagsrqt(j);
    M1 = M / Mscl;
    Mpnt = &M1;
  } else {
    Mpnt = &M;
  }   

  if ( !chcksing(*Mpnt, logDet, viol, minegv, maxegv, minlndet, maxlnk2) )  return false;
  if (!inv_sympd(MInv,*Mpnt)) {
    viol =  std::numeric_limits<double>::infinity();
    return false;
  }
  if (scale) {
    MInv = MInv / Mscl;
    for (int i=0;i<p;++i) logDet += 2*log(diagsrqt(i));
  }  

  return true;
}

bool safepdsolve(const mat& M, const mat& rhs, mat& res, double& logDet, double& viol, double& minegv, double& maxegv,  
  const double minlndet, const double maxlnk2, const bool scale)
{
  int p(M.n_rows);
  int q(rhs.n_cols);
  static mat M1,rhs1;
  static vec diagsrqt;
  const mat* Mpnt;
  const mat* rhspnt;
  if (scale) {  
    if (M1.n_rows!=p|| M1.n_cols!=p) M1.set_size(p,p);
    if (diagsrqt.n_elem!=p) diagsrqt.set_size(p);
    if (rhs1.n_rows!=p || rhs1.n_cols!=q) rhs1.set_size(p,q);
    for (int i=0;i<p;++i) diagsrqt(i) = sqrt(M(i,i));
    for (int i=0;i<p;++i) {
      M1(i,i) = 1.;
      for (int j=0;j<i;++j) M1(i,j) = M1(j,i) = M(i,j)/(diagsrqt(i)*diagsrqt(j));
      for (int j=0;j<q;++j) rhs1(i,j) = rhs(i,j)/diagsrqt(i);  
    }
    Mpnt = &M1;
    rhspnt = &rhs1;
  } else {
    Mpnt = &M;
    rhspnt = &rhs;
  }   

  if ( !chcksing(*Mpnt, logDet, viol, minegv, maxegv, minlndet, maxlnk2) ) return false;
  if (!solve(res,*Mpnt,*rhspnt,solve_opts::likely_sympd)) {
    viol =  std::numeric_limits<double>::infinity();
    return false;
  }
  if (scale) for (int i=0;i<p;++i) {
    for (int j=0;j<q;++j) res(i,j) /= diagsrqt(i); 
    logDet += 2*log(diagsrqt(i));
  }

  return true; 
}




