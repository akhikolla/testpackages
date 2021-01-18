
#include <Rmath.h>
#include <R.h>
#include <math.h>
#include <Rcpp.h>

#include <cmath>


using namespace Rcpp;
using namespace std;
 //x1 < u1, x1+x2 >= u2, x1 from 0 to n1, x2 from 0 to n2-n1,x3 from 0 to n3-n2
//x1<u1,x1+x2<u2,x1+x2+x3>=u3


double ztest(double r1, double r2, double n1,double n2 ){

	double z;
	if ((r1+r2)/(n1+n2)*(1-(r1+r2)/(n1+n2))*(1/n1+1/n2)==0)
		z=0;

	z=(r1/n1-r2/n2)/sqrt((r1+r2)/(n1+n2)*(1-(r1+r2)/(n1+n2))*(1/n1+1/n2));


	return(z);

}


// [[Rcpp::export]]
double crejprob(double pe, double ps, double x1, double n1, double n2,double zalpha)
{
	if(n2==0)return 1;
	double value=0.0;
    int x2,y2;
	double b,c;
	//Rcpp::NumericVector zalpha(1);
	//zalpha[0]=1-alpha;
	
        for(x2 = 0; x2 <= n2; x2++) {
            
            b=Rf_dbinom(x2, n2,pe,0);
            
            for(y2 = 0; y2 <= n2; y2++) {
                c=Rf_dbinom(y2, n2,ps,0);
                
                value = value+ b*c*(ztest(x1+x2,y2,n1+n2,n2)>zalpha);
                
            }
            
        }
	return value;
	
}



// [[Rcpp::export]]
double rejprob(double pe, double ps, double r1, double n1, double n2,double zalpha)
{
		double value=0.0;
	    for(int x1 = r1+1; x1 <= n1; x1++) {
			value = value + Rf_dbinom(x1, n1,pe,0)*crejprob( pe,  ps,  x1,  n1,  n2, zalpha);
			
		}
		return value;
}
