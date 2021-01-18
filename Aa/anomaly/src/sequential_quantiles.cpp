#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

/*
Online quantile estimation - methodology from
EFFICIENT RECURSIVE PROCEDURE FOR ESTIMATING A QUANTILE OF AN UNKNOWN DISTRIBUTION - LUKE TIERNEY - 1983
*/

/*
ONLINE QUANTILE ESTIMATION FUNCTIONS
helpers
update_quantile
*/

// helper functions
int Z(double x , double y)
{

  if (x <= y){
    return(1);
  }
  return(0);
    
}

int I(double x , double y , double z)
{

  if ( std::abs(x - y) <= z){
    return(1);
  }
  return(0);
  
}


// update \alpha^th quantile estimator with new observed data point x
std::tuple<double,double,double,double,int> update_quantile( std::tuple<double,double,double,double,int> state , const double& x,const double& alpha)
{
 
  double a = 0.25;
  double zeta = std::get<0>(state);
  double fn = std::get<1>(state);
  double d = std::get<2>(state);
  double d0 = std::get<3>(state);
  // counter
  int i = std::get<4>(state);

  zeta = zeta - ( d/(i+1) ) * ( Z( x , zeta ) - alpha );
  fn = ( 1.0/(i+1) )*( i*fn + I( zeta , x , 1.0/sqrt((double)(i+1)) )/(2.0/sqrt((double)(i+1))) );
  d = std::min( 1.0/fn , d0*pow( i+1 , a ) );
     
  int counter = i+1;
  state = std::make_tuple( zeta , fn , d , d0 , counter );
  return(state);
  
}


/*
SEQUENTIAL QUANTILE
input 
-data 
-n 
-burnin
output
-mu vecotr
-sigma vector
*/
std::tuple<std::vector<double>,std::vector<double>> sequential_ests(const std::vector<double>& data, int n, int burnin, std::tuple<double,double> lqinit, std::tuple<double,double> medinit, std::tuple<double,double> uqinit)
{

  // initial estimates from burnin
  double muburnin = std::get<0>(medinit);
  double sigburnin = (std::get<0>(uqinit) - std::get<0>(lqinit))/1.349;
     
  std::vector<double> mu(n, muburnin);
  std::vector<double> sigma(n, sigburnin);

  // initialise states
  double d0 = 1/( (std::get<0>(uqinit) - std::get<0>(lqinit)) );
  double d = d0;
  int counter = 0;

  double zeta = std::get<0>(lqinit); 
  double fn = std::get<1>(lqinit);
  std::tuple<double,double,double,double,int> lqstate = std::make_tuple(zeta, fn, d, d0, counter);

  zeta = std::get<0>(medinit); 
  fn = std::get<1>(medinit);
  std::tuple<double,double,double,double,int> medstate = std::make_tuple(zeta, fn, d, d0, counter);

  zeta = std::get<0>(uqinit);
  fn = std::get<1>(uqinit); 
  std::tuple<double,double,double,double,int> uqstate = std::make_tuple(zeta, fn, d, d0, counter);
  
  for(int i = burnin; i < n; i++)
    {
      lqstate = update_quantile( lqstate , data[i] , 0.25 );  
      medstate = update_quantile( medstate , data[i]  , 0.5 );
      uqstate = update_quantile( uqstate , data[i]  , 0.75 );
      mu[i] = std::get<0>(medstate);
      sigma[i] = (std::get<0>(uqstate) - std::get<0>(lqstate))/1.349;
    }

  std::tuple<std::vector<double>,std::vector<double>> params = std::make_tuple(mu, sigma);
  return(params);
  
}
