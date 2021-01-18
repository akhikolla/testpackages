// MIT License
// Copyright (c) 2019 Vincent Runge

#include<iostream>
#include "Costs.h"
#include "math.h"

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

Costs::Costs(){}

//####### slopeCost #######////####### slopeCost #######////####### slopeCost #######//
//####### slopeCost #######////####### slopeCost #######////####### slopeCost #######//

double Costs::slopeCost(double& u, double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& S2t, double& S2T, double& SPt, double& SPT)
{
  return(S2T-S2t - (2.0/(T-t))*((T*u-t*v)*(S1T-S1t) + (v-u)*(SPT-SPt)) + (v*v-u*u)/2.0 + (T-t)*(u*u + u*v + v*v)/3.0 + (v-u)*(v-u)/(6.0*(T-t)));
}

//####### vhat #######////####### vhat #######////####### vhat #######//
//####### vhat #######////####### vhat #######////####### vhat #######//

double Costs::vhat(double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& SPt, double& SPT)
{
  return((6.0/((T-t-1)*(2.0*(T-t)-1)))*(T*(S1T-S1t) - (SPT-SPt)) - v*(T-t+1)/(2.0*(T-t)-1));
}

//####### closestStateIndex #######////####### closestStateIndex #######////####### closestStateIndex #######//
//####### closestStateIndex #######////####### closestStateIndex #######////####### closestStateIndex #######//

unsigned int Costs::closestStateIndex(double& v, double* states, unsigned int p)
{
  if(v <= states[0]){return(0);}
  if(v >= states[p - 1]){return(p - 1);}

  // binary search
  unsigned int i = 0;
  unsigned int j = p;
  unsigned int mid = 0;
  while(i < j)
  {
    mid = (i + j)/2;
    if(states[mid] == v){return(mid);}

    if(v < states[mid])
    {
      if(mid > 0 && v > states[mid - 1])
        {if(states[mid - 1] + states[mid] > 2*v){return(mid - 1);}else{return(mid);}}
      j = mid;
    }
    else
    {
      if (mid < (p - 1) && v < states[mid + 1])
        {if(states[mid] + states[mid + 1] > 2*v){return(mid);}else{return(mid + 1);}}
      i = mid + 1;
    }
  }
  return(mid);
}




//####### fillCoeffsAG #######////####### fillCoeffsAG #######////####### fillCoeffsAG #######//
//####### fillCoeffsAG #######////####### fillCoeffsAG #######////####### fillCoeffsAG #######//

void Costs::fillCoeffsAG(double** coeffs, double* sumY, unsigned int n)
{
  ///fill coeffs[t] = CoeffsAG[t] for t = 1 to n-1 (0 and n un-necessary)
  double* gtT = new double[n + 1];
  for(unsigned int t = 1; t < (n - 1); t++) /// Fill coeffs[t] using gtT[t + 1]... gtT[n]
  {
    gtT[t + 1] = 0;
    for(unsigned int T = t + 1; T < n; T++) //update gtT from gtT[t + 2] to gtT[n-1]
    {
      gtT[T + 1] = ((T-t)*gtT[T] + (sumY[T]-sumY[t]))/(T-t+1);
    }
    linReg(coeffs[t], gtT, t, n); //consider gtT only from t + 1 to n in this function
    //std::cout << t << " : " << coeffs[t][0] << " "<< coeffs[t][1] << " "<< coeffs[t][2] << " "<< coeffs[t][3] << std::endl;
  }
  coeffs[n - 1][0] = 0;
  coeffs[n - 1][1] = 0;
  coeffs[n - 1][2] = 0;
  coeffs[n - 1][3] = 0;

  coeffs[n][0] = 0;
  coeffs[n][1] = 0;
  coeffs[n][2] = 0;
  coeffs[n][3] = 0;

  delete [] gtT;
  gtT = NULL;
}


//####### linReg #######////####### linReg #######////####### linReg #######//
//####### linReg #######////####### linReg #######////####### linReg #######//

void Costs::linReg(double* coeff, double *gtT, unsigned int t, unsigned int n)
{
  double sum_gtT = 0;
  double sum_igtT = 0;
  for(unsigned int i = t + 2; i < n + 1; i++) // 0 in t + 1
  {
    sum_gtT = sum_gtT + gtT[i];
    sum_igtT = sum_igtT + i * gtT[i];
  }
  double prod3 = 1.0*(n-t)*(n-t-1)*(n-t+1);
  double sum3 = 1.0*(n+t+1);

  double slope = (-6*sum3/prod3)*sum_gtT + (12/prod3)*sum_igtT;
  double intercept =((1/(1.0*(n-t))) + 3*sum3/prod3)*sum_gtT - (6*sum3/prod3)*sum_igtT;

  double residuMax = 0;
  double residuMin = 0;

  for(unsigned int i = t + 1; i < n + 1; i++) // 0 in t + 1
  {
    if(gtT[i] - (slope*i + intercept) > residuMax){residuMax = gtT[i] - (slope*i + intercept);}
    if(gtT[i] - (slope*i + intercept) < residuMin){residuMin = gtT[i] - (slope*i + intercept);}
  }

  ///matrix of coefficients: Ap, Am, Gp, Gm
  coeff[0] = slope;
  coeff[1] = slope;
  coeff[2] = intercept + residuMin;
  coeff[3] = intercept + residuMax;
}


//####### pruningTest #######////####### pruningTest #######////####### pruningTest #######//
//####### pruningTest #######////####### pruningTest #######////####### pruningTest #######//

bool Costs::pruningTest(double r, double s, unsigned int tau, unsigned int t, unsigned int n, double S, double Ap, double Am, double Gp, double Gm)
{
  bool response = false;

  if(r == s || t == n){response = true;}
  else if(s > r)
  {
    double state1 = (r+2*s)/6.0;
    double tempGp = tau*state1 - (s-r)/(12.0*(t-tau)) + S/(1.0*(t-tau)) + Gp;
    double ApC = Ap - state1;
    if((ApC*(t+1) + tempGp > 0) && (ApC*n + tempGp > 0)){response = true;}
  }
  else /// s < r
  {
    double state2 = (r+2*s)/6.0;
    double tempGm = tau*state2 - (s-r)/(12.0*(t-tau)) + S/(1.0*(t-tau)) + Gm;
    double AmC = Am - state2;
    if((AmC*(t+1) + tempGm < 0) && (AmC*n + tempGm < 0)){response = true;}
  }

  //response = false; /// TO REMOVE
  return(response);
}


//####### angleTest #######////####### angleTest #######////####### angleTest #######//
//####### angleTest #######////####### angleTest #######////####### angleTest #######//

bool Costs::angleTest(unsigned int& t1, unsigned int& t2, unsigned int& t3, double& v1, double& v2, double& v3, double& minAngle)
{
  bool response = false;
  double cosAngleRad = ((1.0*t1-1.0*t2)*(1.0*t3-1.0*t2) + (v1-v2)*(v3-v2))/(sqrt(((1.0*t1-1.0*t2)*(1.0*t1-1.0*t2) + (v1-v2)*(v1-v2))*((1.0*t3-1.0*t2)*(1.0*t3-1.0*t2) + (v3-v2)*(v3-v2))));
  double theta = acos(cosAngleRad) *180.0 / M_PI; // in degree

  if(theta >= minAngle){response = true;}
  if((t1 == t2) && (v1 == v2)){response = true;}

  return(response);
}

