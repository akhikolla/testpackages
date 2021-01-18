// Implementation of chopthin in C++
// Copyright (C)2015 Axel Gandy
// Version 0.2.1

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <limits>

inline double myfloor(double x){
  double res=floor(x);
  if (res==x&&res>=2.) res-=1.;
  return(res);
}

template<bool checks> void chopthin_internal(std::vector<double>& w, unsigned int N, std::vector<double>& wres, std::vector<int>& ires, double eta=5.828427,bool normalise=true){
  if (checks){
    if (eta==std::numeric_limits<double>::infinity()){
      if (N>w.size())
	chopthin_error("Number of particles cannot be increased just by thinning.");
    }
    if (eta<4) chopthin_error("eta less than 4");
    if (N<=0) chopthin_error("N must be >0");
    for (std::vector<double>::iterator wi = w.begin(); wi!=w.end(); wi++){
      if (*wi<(double)0.) chopthin_error("Negative weights not allowed.");
    }
  }
  if (ires.size()!=N) ires.resize(N);
  if (wres.size()!=N) wres.resize(N);
  int n=w.size();
  std::vector<double> vl= w;
  std::vector<double>::iterator vli, vui, j;

  std::vector<double> vu=vl;
  double sl=0,cm=0,su=0,cu=0;
  double a, b;
  double afinal=-1.;
  double h,sltmp,sutmp;
  int cutmp, cmtmp;
  while(vl.size()>0||vu.size()>0){
    if (vl.size()>=vu.size()){
      a=vl[intrand(vl.size())];
      b=a*eta/2.;
    }else{
      b=vu[intrand(vu.size())];
      a=2.*b/eta;
    }
    cmtmp=0;
    sltmp=0.;
    for (vli=vl.begin(); vli!=vl.end(); ++vli){
      if (*vli<=a) sltmp+=*vli;
      else cmtmp++;
    }
    sutmp=0.;
    cutmp=0;
    for (vui=vu.begin(); vui!=vu.end(); ++vui){
      if (*vui>=b) {sutmp+=(*vui);cutmp++;};
    }
    if (a<=0){ //have I picked particle with non-positive weight?
      if (cm+cmtmp==0){
	chopthin_error("No positive weights");
      }else{
	h=N+1;
      }
    }else{
      h=((cm+cmtmp)-(cu+cutmp))+(sl+sltmp)/a+(su+sutmp)/b;
    }
    if (h==N) {
      sl+=sltmp;
      su+=sutmp;
      cu+=cutmp;
      cm+=cmtmp;
      afinal=a;
      break;
    }
    if (h>N){
      sl+=sltmp;
      j=vl.begin();
      for (vli=vl.begin(); vli!=vl.end(); ++vli){
	if (*vli>a) {
	  *j=*vli;
	  ++j;
	}
      }
      vl.resize(j-vl.begin());
      j=vu.begin();
      for (vui=vu.begin(); vui!=vu.end(); ++vui){
	if (*vui>b){
	  *j=*vui;
	  ++j;
	}
      }
      vu.resize(j-vu.begin());
    }else{
      j=vl.begin();
      for (vli=vl.begin(); vli!=vl.end(); ++vli){
	if (*vli<a) {
	  *j=*vli;
	  ++j;
	}else cm++;
      }
      vl.resize(j-vl.begin());
      su+=sutmp;
      cu+=cutmp;
      j=vu.begin();
      for (vui=vu.begin(); vui!=vu.end(); ++vui){
	if (*vui<b){
	  *j=*vui;
	  ++j;
	}
      }
      vu.resize(j-vu.begin());
    }
  }
  std::vector<int>::iterator irespos=ires.begin();
  std::vector<double>::iterator wrespos=wres.begin();
  double wtot=0;

  if (afinal<0){
    if (sl==0. &&su==0. &&(N-cm+cu==0.)){      // no action necessary, directly copy results, omitting 0s
      for (int posv=0; posv<n; posv++){	
	if (w[posv]>0){
	  if (checks) if (irespos==ires.end()) chopthin_error("Internal Error: Trying to produce too many particles in special case.");
	  *(irespos++)=posv+1;
	  *(wrespos++)=w[posv];	  
	  wtot+=w[posv];
	}
      }
      if (checks) if (irespos!=ires.end()) chopthin_error("Internal error: Not enough particles produced in special case.");
      if (normalise){
	for (unsigned int count=0; count<N; count++){
	  wres[count]*=((double)N)/wtot;
	}
      }
      return;
    }
    a=(sl+2.*su/eta)/(N-cm+cu);
    if (eta==std::numeric_limits<double>::infinity()){
      b=std::numeric_limits<double>::infinity();
    }else{
      b=a*eta/2;
    }
  }

  //ensuring that total weight proportions are maintained on average
  // Start by thinning
  double u = myrunif();
  for (int posv=0; posv<n; posv++){
    if (w[posv]<a){
      u-=w[posv]/a;
    }else {
      if (w[posv]==a){
	if (a!=0){
	  u-=1.;
	}else{
	  if (N>cm)
	    u-=((double)(N-cm))/((double)n-cm);
	}
      }
    }
    if (u<0){
      u+=(double)1.;
      *(irespos++)=posv+1;
      *(wrespos++)=a;
      wtot+=a;
    }
  }
  // Correct for thinning by changing the overall weight
  double facw=1.;
  double faccount=1.;
  double difflow=sl-(irespos-ires.begin())*a;
  if (difflow!=0.&&su>0){ // Only needed if thinning happens and if weights above b exist.
    //need to compute sum of residuals
    double sures=0;
    for (std::vector<double>::iterator wi = w.begin(); wi!=w.end(); wi++){
      if (*wi>=b){
	double tmp=*wi/b;
	sures+=tmp-myfloor(tmp);
      }
    }    
    sures*=b;
    facw=(difflow+sures)/sures;
    faccount=(difflow/a+sures/b)/(sures/b);
    if (faccount <0){
      if (checks){
	if (faccount < -1e-5){
	  chopthin_error("Internal Error: faccount<-10^-5");
	}
      }
      faccount=0.;
    }
    if (facw<0) {
      if (checks){
	chopthin_error("Internal Error: facw<0");
      }else{
	facw=0.;
      }
    }	
  }
  // Now chop
  u = myrunif(); //new random variable
  for (int posv=0; posv<n; posv++){
    if (w[posv]>a){
      if (w[posv]>=b){
	double wint=myfloor(w[posv]/b)*b;
	double wres=w[posv]-wint;	
	double countakt=(wres*faccount+wint)/b;
	double wakt=wres*facw+wint;
	u-=countakt;
	if (u<0){
	  int ndes=ceil(-u);
	  u+=(double) ndes;
	  double wdes=wakt/ndes;
	  for (int j=0; j<ndes; j++){
	    if (checks) if (irespos==ires.end()) chopthin_error("Internal Error: Trying to produce too many particles.");
	    *(irespos++)=posv+1;
	    *(wrespos++)=wdes;
	    wtot+=wdes;
	  }
	}
      }else{ //a<w[posv]<b
	if (checks) if (irespos==ires.end()) chopthin_error("Internal Error: Trying to produce too many particles.");
	*(irespos++)=posv+1;
	*(wrespos++)=w[posv];
	wtot+=w[posv];
      }
    }
  }
  if (checks) if (irespos!=ires.end()) chopthin_error("Internal error: Not enough particles produced.");
  if (normalise){
    for (unsigned int count=0; count<N; count++){
      wres[count]*=((double)N)/wtot;
    }
  }
}
