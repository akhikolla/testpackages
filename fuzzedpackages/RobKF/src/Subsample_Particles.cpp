#include "Particle.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>

using namespace std;

std::list < struct Particle > Subsample_Particles(std::list < struct Particle > & candidates, const int &N)
{

	
	std::list<Particle>::iterator it = candidates.begin();

	double max_log_prob = it->log_prob;

	for (it = candidates.begin(); it != candidates.end(); ++it)
	{

		if (it->log_prob > max_log_prob)
		{ 

    			max_log_prob = it->log_prob;

		}	

	}
	
	list < Particle >  New_sample;

	int Size = candidates.size();

	Rcpp::NumericVector Probs(Size);
	Rcpp::IntegerVector Index(Size);
	
	int ii = 0;

	double dummy;

	it = candidates.begin();
	
  	while(it != candidates.end())
    	{
      		
		Index[ii] = ii;
		Probs[ii] = std::exp(it->log_prob - max_log_prob);
		ii++;
		it++;

    	}

	Rcpp::IntegerVector Indexsample = Rcpp::sample(Index,N,true,Probs);

	std::list<Particle>::iterator pos = candidates.begin();
	
	for (int jj = 0; jj < N; jj++)
	{

		advance(pos,Indexsample[jj]);		
		New_sample.push_back(*pos);
		pos = candidates.begin();

	}

	return New_sample;

};
