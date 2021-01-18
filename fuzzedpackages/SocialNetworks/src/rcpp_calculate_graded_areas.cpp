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

RcppExport SEXP rcpp_calculate_graded_areas(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5) {
    std::vector<double> x=Rcpp::as<std::vector<double> >(arg1);
    std::vector<double> y=Rcpp::as<std::vector<double> >(arg2);
    std::vector<double> ir=Rcpp::as<std::vector<double> >(arg3);
    std::vector<double> rand1=Rcpp::as<std::vector<double> >(arg4);
    std::vector<double> rand2=Rcpp::as<std::vector<double> >(arg5);
    
    double accuracy=rand1.size();
    
    std::vector<std::vector<double> > out(x.size(), std::vector<double>(x.size()));
    
    for(unsigned int a=0; a<x.size(); a++) {
        for(unsigned int b=0; b<x.size(); b++) {
            long ptsinab=0;
            
            for(unsigned int i=0; i<accuracy; i++) {
                // Simulate points in unit square. 
                double simx=rand1[i];
                double simy=rand2[i];
                
                // Scale to point a. 
                simx=simx*0.5*ir[a]+x[a];
                simy=simy*0.5*ir[a]+y[a];
                
                if(close(x[b], y[b], simx, simy, ir[b])) {
                    ptsinab++;
                }
            }
            out[a][b]=((double)ptsinab)/((double)rand1.size());
        }
    }
    
    
    
    
    return Rcpp::wrap(out);
}
