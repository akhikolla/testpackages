/*
 * functions.h
 *
 *  Created on: 1 mars 2013
 *      Author: grimonprez
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include<vector>
#include <Rmath.h>

/**r random number generator, code from : http://gallery.rcpp.org/articles/r-function-from-c++/*/
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

/**
 * search the position of i in x
 * @param x rank
 * @param i integer which we want the position in x
 * @return the position of i in x
 */
int positionRank(std::vector<int> const& x,int const& i);

/**
 * compute  A(x,y) and G(x,y,mu)
 * @param x rank
 * @param y order of presentation of x
 * @param mu order of reference
 * @return a vector of 2 elements (A(x,y),G(x,y,mu))
 *
 * A(x,y)=total number of comparaison in the insertion sorting algorithm
 * G(x,y,mu)= total number of good comparaison according to mu in the insertion sorting algorithm
 *
 */
std::vector<int> comparaison(std::vector<int> const& x,std::vector<int> const& y,std::vector<int> const& mu);

/**
 * compute the conditional probability  p(x|y;mu,p)
 * @param x rank
 * @param y order of presentation of x
 * @param mu order of reference
 * @param p probability of making a good comparaison
 * @return p(x|y;mu,p)
 */
double probaCond(std::vector<int> const& x,std::vector<int> const& y,std::vector<int> const& mu,double const& p);

/**
 * factorial function (recursive)
 * @param m an integer >0
 * @return m!
 */
int factorial(int const& nombre);

/**
 * compute all the factorial from 1! to m!
 * @param m a positive integer
 * @return a vector of size m containing all the factorial from 1! to m!
 */
std::vector<int> tab_factorial(int const& m);

/**
 *-----------rang->index------------
 * convert a rank into an integer
 * @param rang rank
 * @param  tabFactorial \see tab_factorial
 * @return index
 */
int rank2index(std::vector<int> const& rang,std::vector<int> const& tabFactorial);

/**
 *-----------rang->index------------
 * @see rank2index for a vector of rank
 * @param listeRang vector of rank
 * @param  tabFactorial @see tab_factorial
 * @return vector of index
 */
std::vector<int> rank2index(std::vector<std::vector<int> > const& listeRang,std::vector<int> const& tabFact);

/**
 * conversion index->rank
 * @param index index of the rank
 * @param m size of the rank
 * @param tabFactorial vector containing 1! to m! (@see tab_factorial)(optional)
 * @return le rang correspondant à l'index
 * index=rank
 * 0=1 2 ... m
 *
 * m!= m m-1 .. 1
 */
std::vector<int> index2rank(int index,int const& m,std::vector<int> const& tabFactorial);
std::vector<int> index2rankb(int index,int const& m,std::vector<int> const& tabFactorial);
std::vector<int> index2rank(int index,int const& m);


/**
 * Return a vector containing the index of order of presentation y for which we have to compute probability
 * @param m  size of the rank
 * @param tabFactorial vector containing 1! to m! (@see tab_factorial)
 * @return a vector containing the index of order of presentation y for which we have to compute probability
 *
 * p(x|y;mu,p) doesn't change if the 2 first element of y are inverted
 * So we compute probabilities for the half of y
 */
std::vector<int> listeSigma(int const& m,std::vector<int> const& tabFactorial);

/**
 * invert a rank (1 2 3 become 3 2 1)
 * @param rang rank to invert
 * @param m size of the rank
 */
void inverseRang(std::vector<int> &rang);
void inverseRang(std::vector<int> &rang,int const& m);

/**
 * Compute the BIC
 * @param loglikelihood the loglikelihhod of the data
 * @param nbDonnees total number of data
 * @return BIC
 */
double BIC(double loglikelihood,int nbDonnees,int nbParam);

/**
 * compute the Rand index of 2 partitions
 * @param z1 partition
 * @param z2 partition
 * @return a double, the Rand index
 */
double indiceRand(std::vector<int> const& z1,std::vector<int>  const& z2);

/**
 * conversion from ordering representation to ranking representation
 * @param x a rank
 * @param m size of the rank
 * @return a rank : the ranking representation of x
 */
std::vector<int> order2rank(std::vector<int> const& x,int const& m);

/**
 *  Compute the Kendall's distance between 2 ranks (in ordering representation)
 * @param x a rank in ordering representation
 * @param y a rank in ordering representation of the same size than x
 * @return an integer, the kendall distance between x and y
 */
int distanceKendall(std::vector<int> const& x,std::vector<int> const& y);

/**
 * sort the parameters in order that the first cluster is the cluster with the more little ndex of mu
 * @param mu index of the rank of the first dimension of listeMu
 * @param p parameter of the ISR
 * @param prop proportion of the mixture model
 * @param listeMu reference rank
 * @param z partition
 * @param g number of cluster
 * @param d number of dimension
 * @param n number of individual
 * listeMu, p, prop and z are modify if necessary
 *
 */
void tri_insertionMulti(std::vector<int> &mu,std::vector<double> &prop,std::vector<std::vector<double> > &p,
		std::vector<std::vector<std::vector<int> > > &listeMu,std::vector<int> &z,int const& g,int const& d,int const& n);


/**
 * simulation of a n-sample of ISR(mu,p)
 * @param n size of the sample
 * @param m size of the rank
 * @param mu rank
 * @param p probability of a good comapraison
 * @return a n-sample of ISR(mu,p)
 */
std::vector<std::vector<int> > simulISR(int const& n,int const& m,std::vector<int> const& mu, double const& p);

/**
 * compute frequence of a data set
 * @param listeRang data
 * @return a pair with the unique rank and the frequence of each unique rank
 */
std::pair<std::vector<std::vector<std::vector<int> > >,std::vector<int> >
freqMulti(std::vector<std::vector<std::vector<int> > > const& listeRang);

/**
 * compute probability of x according to multivariate ISR
 * @param x rank for each dimension for compute probability
 * @param mu reference rank for each dimension
 * @param pi dispersion parameter for each dimension
 * @return p(x;mu,pi)
 */
double proba(std::vector<std::vector<int> > const& x, std::vector<std::vector<int> > const& mu, std::vector<double> const& pi);

#endif /* FUNCTIONS_H_ */
