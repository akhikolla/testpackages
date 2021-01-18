
#include <Rcpp.h>
#include <vector>
    
#include<stdio.h>
#include<stdlib.h>
#include<memory.h>

#include <math.h>
#include <assert.h>
#include <stdlib.h>

using namespace Rcpp;

#define close(x1, y1, x2, y2, ir) (( ((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)) )<=(ir)*(ir))

RcppExport SEXP rcpp_calculate_areas(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5) {
    std::vector<double> x=Rcpp::as<std::vector<double> >(arg1);
    std::vector<double> y=Rcpp::as<std::vector<double> >(arg2);
    std::vector<double> ir=Rcpp::as<std::vector<double> >(arg3);
    std::vector<double> rand1=Rcpp::as<std::vector<double> >(arg4);
    std::vector<double> rand2=Rcpp::as<std::vector<double> >(arg5);
    
    double accuracy=rand1.size();
    
    // For point in the input file, compute its cache vector
    
    std::vector<std::vector<double> > out(x.size(), std::vector<double>(x.size()));
    
    // Alternative implementation without bitwise ops. 
    for(unsigned int a=0; a<x.size(); a++) {
      for(unsigned int b=0; b<x.size(); b++) {
        long ptsina=0;
        long ptsinab=0;
        
        for(unsigned int i=0; i<accuracy; i++) {
          // Simulate points in unit square. 
          double simx=rand1[i];  //simu[i].first;   //gen();
          double simy=rand2[i];  //simu[i].second;  //gen();
          
          // Scale to point a. 
          simx=simx*2*ir[a]+x[a]-ir[a];
          simy=simy*2*ir[a]+y[a]-ir[a];
          if(close(x[a], y[a], simx, simy, ir[a])) {
            ptsina++;
            if(close(x[b], y[b], simx, simy, ir[b])) {
              ptsinab++;
            }
          }
        }
        out[a][b]=((double)ptsinab)/((double)ptsina);
      }
    }
    
    return Rcpp::wrap(out);
}
