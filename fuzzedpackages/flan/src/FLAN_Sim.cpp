/**
  * FLAN Software
  *
  * @author 2015-2020 Adrien Mazoyer  <adrien.mazoyer@imag.fr>
  * @see The GNU Public License (GPL)
  */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "FLAN_Sim.h"

using namespace Rcpp;

    // --------------------------------------------------------------------------------------------------------------
    //////////////////////////////////////////// Samples computation ////////////////////////////////////////////////
    // --------------------------------------------------------------------------------------------------------------

    /* Generates a sample of size n of couples (mutant count ; final count)
     * The final counts are sampled with the log-normal
     * distribution with coefficient of variation cvfn and mean mfn.
     *
     * n: number of experiments
     * mfn: mean of final counts
     * Cfn: variation number of final counts
     * alphapi: mean number of mutations probability of mutation (depends of cvfn)
     * rho: fitness parameter
     * death: probability of death
     */


List FLAN_Sim::computeSamplesMutantsFinalsNumber(int n)  {


  RNGScope rngScope;

  NumericVector mutantCount(n);

  if (mCvfn>0) {
    NumericVector finalCount(n);
      if(mDistfn.compare("lnorm") == 0){
        double sdLog2=log(1.+mCvfn*mCvfn);
        double sdLog=sqrt(sdLog2);
        double meanLog=log(mMfn)-sdLog2/2;

        finalCount=rlnorm(n,meanLog,sdLog);
      } else if(mDistfn.compare("gamma") == 0){
        double shape=1./(mCvfn*mCvfn);
        double scale=mMfn/shape;

        finalCount=rgamma(n,shape,scale);
      }
      mutantCount=computeSampleMutantsNumber(n,finalCount);
      return List::create(_["mc"]=mutantCount,
			  _["fn"]=finalCount);


  } else {
      mutantCount=computeSampleMutantsNumber(n);
      return List::create(_["mc"]=mutantCount,
			  _["fn"]=R_NilValue);
  }

}



    /* Generates a sample of size n of mutant counts
     * n: sample size
     * alpha: mean number of mutations
     * rho: fitness parameter
     * death: probability of death
     */



NumericVector FLAN_Sim::computeSampleMutantsNumber(int n)  {

    double s,sj;
    // set the size of the mutants
    NumericVector mutantCount=rpois(n,mMut);
    int mc;
    int i=0,j;
    bool testneg;
    for (NumericVector::iterator it = mutantCount.begin(); it != mutantCount.end(); ++it,i++) {
        // simulate a poisson number

	  mc=(int)(*it);

        if (mc>0) {
            //Poisson sum of Yule variables
            NumericVector sample=mClone->computeSample(mc);

            s=0;
	    j=0;
	    testneg=false;
	    while(j<mc && !testneg){
	      sj=sample[j];
	      if(sj < 0){
		testneg=true;
		s=sj;
	      } else {
		s+=sj;
		j++;
	      }
	    }
	    *it=s;
        } else
            *it=0;
    }

    return mutantCount;
}


    /* Generates a sample of size n of mutant counts, given a sample of final counts.
     *
     * n: number of experiments
     * pi: mutation probability
     * rho: fitness parameter
     * death: probability of death
     * finalCount:
     */

NumericVector FLAN_Sim::computeSampleMutantsNumber(int n,
					  NumericVector& finalCount)  {

    std::vector <double> mutantCount(n);

    // poisson number of mutation

//     NumericVector sample;
    double s,sj;
    // set the size of the mutants

    double lambda;
    int mc,j;
    bool testneg;
//     int i=0;
    NumericVector::iterator itfn=finalCount.begin();
    for (std::vector <double>::iterator itmc = mutantCount.begin();
	  itmc != mutantCount.end(); ++itmc, ++itfn) {
        // simulate a poisson number
      lambda=mMut*(*itfn);

      mc=(int)(rpois(1,lambda)[0]);

      if (mc>0) {
            //Poisson sum of Yule variables
            NumericVector sample=mClone->computeSample(mc);
            //s = sum_j sample[j]
//             int m=sample.getSize();
            s=0;
	    j=0;
	    testneg=false;
	    while(j<mc && !testneg){
	      sj=sample[j];
	      if(sj < 0){
		testneg=true;
		s=sj;
	      } else {
		s+=sj;
		j++;
	      }
	    }
	    *itmc=s;
        } else
            *itmc=0;
    }

    return NumericVector(mutantCount.begin(),mutantCount.end());

}


/* Class for simulation of inhomogeneous models */

List FLAN_SimInhomogeneous::computeSamplesMutantsFinalsNumber(int n)  {


  RNGScope rngScope;

  NumericVector mutantCount(n);

  if (mCvfn>0) {
      double sdLog2=log(1.+mCvfn*mCvfn);
      double sdLog=sqrt(sdLog2);
      double meanLog=log(mMfn)-sdLog2/2;

      NumericVector finalCount=rlnorm(n,meanLog,sdLog);

      mutantCount=computeSampleMutantsNumber(n,finalCount);
      return List::create(_["mc"]=mutantCount,
			  _["fn"]=finalCount);


  } else {
      mutantCount=computeSampleMutantsNumber(n);
      return List::create(_["mc"]=mutantCount,
			  _["fn"]=R_NilValue);
  }

}



    /* Generates a sample of size n of mutant counts
     * n: sample size
     * alpha: mean number of mutations
     * rho: fitness parameter
     * death: probability of death
     */



NumericVector FLAN_SimInhomogeneous::computeSampleMutantsNumber(int n)  {

    NumericVector mutantCount=rpois(n,mMut);
    int mc;
    NumericVector sample;
    double sumS;
    double tp;
    double s;
    double muInf=as<double>((*mMU)(0.,R_PosInf));
//     double f0=as<double>((*mF)(0.));
//     double fn=as<double>((*mF)(R_PosInf));
    double power=mFitness*(1-2*mDeath);

    for (NumericVector::iterator itmc = mutantCount.begin(); itmc != mutantCount.end(); ++itmc) {
        // simulate a poisson number of mutations

      mc=(int)(*itmc);

      if (mc>0) {
	// simulate a NHPP with expectation the derivative of mF
// 	sample.resize(mc);
	sample=runif(mc,0.,1.);
	sumS=0;
	for (NumericVector::iterator itS = sample.begin(); itS != sample.end(); ++itS) {
	  tp=(*itS);
// 	  std::cout<<"U ="<<tp<<std::endl;
// 	  tp=pow(tp*(pow(fn/f0,power)-1)+1,1./power)*f0;
	  tp*=(exp(power*muInf)-1);
	  tp++;
	  tp=log(tp);
	  tp/=power;

// 	  std::cout<<"Après transfo 1 ="<<tp<<std::endl;

// 	  else s=as<double>((*mFinv)(tp));
// 	  if(tp > muInf) sumS++;
// 	    s=-1;

// 	  else {
	    s=as<double>((*mMUinv0)(tp));

	    tp=mClone->computeSample(1,s)[0];

	    sumS+=tp;
// 	  }
	}
	*itmc=sumS;
      } else *itmc=0;
    }

    return mutantCount;
}


NumericVector FLAN_SimInhomogeneous::computeSampleMutantsNumber(int n,
					  NumericVector& finalCount)  {


    std::vector <double> mutantCount(n);
    int mc;
    double mut;
    NumericVector sample;
    double sumS;
    // NumericVector TP;
    double tp;
    double s;
    double muInf=as<double>((*mMU)(0.,R_PosInf));
//     double f0=as<double>((*mF)(0.));
//     double fn=as<double>((*mF)(R_PosInf));
//     std::cout<<"f0 ="<<f0<<std::endl;
//     std::cout<<"fn ="<<fn<<std::endl;
    double power=mFitness*(1-2*mDeath);

    NumericVector::iterator itfn=finalCount.begin();
    for (std::vector <double>::iterator itmc = mutantCount.begin();
	  itmc != mutantCount.end(); ++itmc, ++itfn) {
        // simulate a poisson number

      mut=mMut*(*itfn);
      mc=(int)(rpois(1,mut)[0]);

      if (mc>0) {
	// simulate a NHPP with exptation the derivative of mF
	sample=runif(mc,0.,1.);
	sumS=0;
	for (NumericVector::iterator itS = sample.begin(); itS != sample.end(); ++itS) {
	  tp=(*itS);
// 	  tp=pow(tp*(pow(fn/f0,power)-1)+1,1./power)*f0;
	  tp*=(exp(power*muInf)-1);
	  tp++;
	  tp=log(tp);
	  tp/=power;
// 	  if(tp < 0) s=0;
// 	  else s=as<double>((*mFinv)(tp));
// 	  else
	  s=as<double>((*mMUinv0)(tp));
// 	  else {
// 	    tpf=(*mFinv)(tp);
// 	    s=as<double>(tpf);
// 	  }

	  tp=mClone->computeSample(1,s)[0];
// 	  std::cout<<"Clone size ="<<tp<<std::endl;
	  sumS+=tp;
	}
	*itmc=sumS;
      } else *itmc=0;

    }

    return NumericVector(mutantCount.begin(),mutantCount.end());
}
