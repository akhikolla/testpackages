#ifndef TEST_H_
#define TEST_H_

#include <vector>
#include <string>
#include <set>
#include <iostream>

#include "functions.h"

struct Rank
{
	std::vector<int> rank;
	bool isPartial;
	std::vector<int> missingIndex;
	std::set<int> missingNumber;
};

/**
 * Simulate a sample of mixture of ISR
 * @param simul sample will be modified
 * @param n size of sample
 * @param m size of rank
 * @param mu reference rank
 * @param p dispersion parameter: probability of a godd comparaison
 * @param prop proportion of the mixture
 */
void simulMixtureISR(std::vector<std::vector<int> > &simul,int const& n,int const& m,std::vector<std::vector<int> > const& mu,std::vector<double> const& p,std::vector<double> const& prop);

/**
 * khi2 adequation test
 * @param data univariate rank
 * @param p dispersion parameter: probability of a godd comparaison
 * @param prop proportion of the mixture
 * @param mu reference rank
 * @param nBoot number of iteration for estimation
 * @return estimated pvalue of khi2 adequation test
 */
double khi2(std::vector<std::vector<int> > const& data,std::vector<double> const& p,std::vector<double> const& prop,std::vector<std::vector<int> > const& mu,int const& nBoot);

/**
 * khi2 adequation test for partial rank
 * @param data univariate partial rank
 * @param p dispersion parameter: probability of a godd comparaison
 * @param prop proportion of the mixture
 * @param mu reference rank
 * @param nBoot number of iteration for estimation
 * @return estimated pvalue of khi2 adequation test
 */
double khi2partial(std::vector<Rank > &data,std::vector<double> const& p,std::vector<double> const& prop,std::vector<std::vector<int> > const& mu,int const& nBoot);


/**
 * Updating the kullback leibler divergence
 * @param divKL actual kullback leibler divergence
 * @param p1 probabilities of each rank for the first set of parameters
 * @param p2 probabilities of each rank for the second set of parameters
 * @param d dimension
 * @param g number of cluster
 * @param proportion1 proportion of the mixture of the first set of parameters
 * @param proportion2 proportion of the mixture of the second set of parameters
 */
void updateD(double &divKL,std::vector<int> &index, std::vector<std::vector<std::vector<double> > > const& p1,std::vector<std::vector<std::vector<double> > >  const& p2,int const& d,int const& g,
		std::vector<double> const& proportion1,std::vector<double> const& proportion2);
/**
 * Recursive function for updating the index for compute probabilities in kullback 
 * @param index current index
 * @param i actual dimension
 * @param factm factorial of the size of rank for each dimension
 * @param stop if true stop the recursivity
 */
void updateIndex(std::vector<int> &index,int i,std::vector<int> const& factm,bool &stop);

/**
 * Compute probabilities of each rank for each set of parameters
 * @param divKL actual kullback leibler divergence
 * @param p probabilities of each rank for the first set of parameters
 * @param q probabilities of each rank for the second set of parameters
 * @param mu1 reference rank for the first set of parameters
 * @param mu2 reference rank for the second set of parameters
 * @param p1 dispersion parameter of the first set of parameters
 * @param p2 dispersion parameter of the second set of parameters
 * @param m size of rank for each dimension
 * @param d dimension
 * @param g number of cluster
 */
void computePQ(std::vector<std::vector<std::vector<double> > > &p, std::vector<std::vector<std::vector<double> > > &q,std::vector<std::vector<std::vector<int> > > const& mu1,
		std::vector<std::vector<std::vector<int> > > const& mu2,std::vector<std::vector<double> > const& p1,std::vector<std::vector<double> > const& p2,std::vector<int> const& m,int d, int g);

/**
 * Compute the kullback leibler divergence between the first and the second set of parameters
 * @param m size of rank for each dimension
 * @param mu1 reference rank for the first set of parameters
 * @param mu2 reference rank for the second set of parameters
 * @param p1 dispersion parameter of the first set of parameters
 * @param p2 dispersion parameter of the second set of parameters
 * @param proportion1 proportion of the mixture of the first set of parameters
 * @param proportion2 proportion of the mixture of the second set of parameters
 * @return kullback leibler divergence
 */
double divKL(std::vector<int> const& m,std::vector<std::vector<std::vector<int> > > const& mu1,std::vector<std::vector<std::vector<int> > > const& mu2,
		std::vector<std::vector<double> > const& p1, std::vector<std::vector<double> >  const& p2,std::vector<double> const& proportion1,std::vector<double> const& proportion2);


#endif /* TEST_H_ */

