
#include <cmath>
#include <list>
#include <vector>
#include <tuple>
#include <boost/math/distributions/normal.hpp>
#include <algorithm>
#include <thread>
#include <chrono>
#include <future>
#include <numeric>

bool disjoint(const std::tuple<int,int>& x,
	      const std::tuple<int,int>& y)
{
    if(std::get<0>(x) > std::get<0>(y))
      {
	return disjoint(y,x);
      }
    if(std::get<1>(x) < std::get<0>(y))
      {
	return true;
      }
    return false;
}

std::tuple<std::list<std::tuple<int,int> >, std::list<double> > pass(const std::list<std::vector<double> >& Xi,
								     const int& Lmax,
								     const int& Lmin,
								     const int& alpha,
								     const double& lambda,
								     std::shared_future<void> futureobj)
{


  try
    {

  
  // size of data
  auto T = (*Xi.begin()).size();
  // create cumsums
  std::vector<std::vector<double> > cumsums(Xi.size());
  transform(Xi.begin(),Xi.end(),cumsums.begin(),[&T](auto& X)
	    {
	      std::vector<double> cumsum(T+1,0.0);
	      std::partial_sum(X.begin(),X.end(),cumsum.begin()+1);
	      return cumsum;
	    });
  // build J
  std::list<std::tuple<int,int> > J;
  for(int k = Lmin; k < Lmax+1; k++)
    {
      for(unsigned int j = 0; j < T-k; j++)
	{
	  auto p = std::make_tuple(j,j+k);
	  J.push_back(p);
	}
    }
    
  // calculate VNJ
  boost::math::normal norm;
  std::list<std::vector<double> >  pJi;
  auto N = Xi.size(); 
  auto sqrtN = std::sqrt(N);
  double dN = (double)N;
  std::vector<int> index(N);
  for(unsigned int i = 0; i < N; i++)
    {
      index[i]=i+1;
    }


  
  std::vector<double> VNJ(J.size());
  int i = 0;
  int ncalcs = 1;
  auto candidate = [&ncalcs,&i,futureobj,&cumsums,&N,&index,&norm,&sqrtN,&dN,&alpha](auto& p)
    { 
      if(i++ % ncalcs == 0 && futureobj.valid() && !(futureobj.wait_for(std::chrono::milliseconds(0)) == std::future_status::timeout))
	{
	  throw std::exception();
	}
      auto a = std::get<0>(p);
      auto b = std::get<1>(p);
      std::vector<double> pJ(N);
      std::transform(cumsums.begin(),cumsums.end(),pJ.begin(),[&a,&b,&norm](auto& C)
    {
      return 2*cdf(complement(norm, fabs((C[b+1]-C[a])/sqrt((double)(b-a+1)))));
    });
      sort(pJ.begin(),pJ.end());
      std::transform(pJ.begin(),pJ.end(),index.begin(),pJ.begin(),[&sqrtN,&dN](auto& p, auto& i)
    {
      return sqrtN*(i/dN-p)/std::sqrt(p*(1.0-p));
    });
      return *std::max_element(pJ.begin()+alpha-1,pJ.begin()+(N/2));
    }; 


  auto start = std::chrono::high_resolution_clock::now();
  if(J.size() > 1)
    {
      transform(J.begin(),std::next(J.begin(),1),VNJ.begin(),candidate);
    }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  ncalcs = int(1.0/elapsed.count()) + 1;
  transform(J.begin(),J.end(),VNJ.begin(),candidate);
  
  // prune candidates according to threshold value
  std::list<std::tuple<double,std::tuple<int,int> > > I;
  auto itVNJ = VNJ.begin();
  auto itJ = J.begin();
  while(itVNJ != VNJ.end())
    {
      if(*itVNJ > lambda)
	{
	  I.push_back(std::make_tuple(*itVNJ,*itJ));
	}
      itVNJ++;
      itJ++;
    }

  // create Ihat
  std::list<std::tuple<int,int> > Ihat;
  std::list<double> xstar;  
  if(I.size() > 0)
    {
      while(true)
	{
	  if(futureobj.valid() && !(futureobj.wait_for(std::chrono::milliseconds(0)) == std::future_status::timeout))
	    {
	      throw std::exception();
	    }
	  std::vector<double> V(I.size());
	  transform(I.begin(),I.end(),V.begin(),[](auto& p) { return std::get<0>(p); });
	  auto pos = std::distance(V.begin(), std::max_element(V.begin(), V.end()));
	  auto itI = I.begin();
	  std::advance(itI,pos);
	  auto valitI = *itI;
	  Ihat.push_back(std::get<1>(*itI));
	  xstar.push_back(*std::max_element(V.begin(), V.end()));
	  I.remove_if([&valitI](auto& p){return !disjoint(std::get<1>(p),std::get<1>(valitI));});
	  if(I.size()==0)
	    {
	      break;
	    }
	}
    }

  
  return std::make_tuple(Ihat,xstar);

    }
  catch(const std::exception& e)
    {
      std::tuple<std::list<std::tuple<int,int> >, std::list<double> > empty;
      return(empty); // interrupted - send back empty container      
    }

  
}
