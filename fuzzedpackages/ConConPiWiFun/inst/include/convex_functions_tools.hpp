/*
 * convex_functions_tools.hpp
 *
 *  Created on: 16 avr. 2013
 *      Author: robin
 */

#ifndef CONVEX_FUNCTIONS_TOOLS_HPP_
#define CONVEX_FUNCTIONS_TOOLS_HPP_

bool isincreasing(Rcpp::NumericVector arg){
	int length=arg.size();
	  bool res=true;
	   for (int n=1; n<(length); n++)
		  if (arg[n]<=arg[n-1]){
			  res=false;
			  break;
		  }
	  return res;
}



// here polynom is 1/2 ax^2+bx+c
double getSlope(std::pair<double,double> Coefficients,double val){
// returns the slope at val given Coefficients a and b.f
	// Coefficients are (a,b)
   if (val==-std::numeric_limits<double>::infinity()&&Coefficients.first!=0){
  	 if (Coefficients.first<0){
			 return(std::numeric_limits<double>::infinity());
		 }else{
			 return(-std::numeric_limits<double>::infinity());
		 }
	 }else{
		 if (val==std::numeric_limits<double>::infinity()&&Coefficients.first!=0){
			 if (Coefficients.first<0){
				 return(-std::numeric_limits<double>::infinity());
			 }else{
				 return(std::numeric_limits<double>::infinity());
			 }
		 }else
		 {
			if (Coefficients.first==0){
				return(Coefficients.second);
			}
			else if (Coefficients.first==-std::numeric_limits<double>::infinity()){
				if (val<0)
				{
					return(std::numeric_limits<double>::infinity());
				}else
				{
					return(-std::numeric_limits<double>::infinity());
				}

			}else if (Coefficients.first==std::numeric_limits<double>::infinity())
			{
				if (val<0)
				{
					return(-std::numeric_limits<double>::infinity());
				}else
				{
					return(std::numeric_limits<double>::infinity());
				}
			}else
			{
				return(Coefficients.first*val+Coefficients.second);
			}
		 }
	 }
}

double getVal(std::pair<double,double> Coefficients,double val){
// returns the val at val given Coefficients a and b
	 if (val==-std::numeric_limits<double>::infinity()&&Coefficients.first!=0){
		 if (Coefficients.first<0){
			 return(std::numeric_limits<double>::infinity());
		 }else{
			 return(-std::numeric_limits<double>::infinity());
		 }
	 }else{
		 if (val==std::numeric_limits<double>::infinity()&&Coefficients.first!=0){
			 if (Coefficients.first<0){
				 return(-std::numeric_limits<double>::infinity());
			 }else{
				 return(std::numeric_limits<double>::infinity());
			 }
		 }else{
			 return(Coefficients.first/2*val*val+Coefficients.second*val);
		 }
	 }
}

double getXetoile(std::pair<double,double> Coefficients){
	 if (Coefficients.first==0){
		 if (Coefficients.second==0)
		 {
			 return(0);
		 }
		 else{
			 if (Coefficients.second<0){
				 return(std::numeric_limits<double>::infinity());
			 }else{
				 return(-std::numeric_limits<double>::infinity());
			 }
		 }
	 }else{
		 return(-Coefficients.second/Coefficients.first);
	 }
}

std::pair<double,double> Slopes2Coeffs(double Slope0,double Slope1){
  // returns the a and b coefficient of 1/2 ax^2+bx+c given the slopes in zero and the slopes in 1
  // a= S1-S0
  std::pair<double,double> res;
  res.first=Slope1-Slope0;
  res.second=Slope0;
  return(res);
}



#endif /* CONVEX_FUNCTIONS_TOOLS_HPP_ */
