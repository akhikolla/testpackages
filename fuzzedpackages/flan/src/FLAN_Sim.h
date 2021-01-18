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


#ifndef FLAN_SIM_H
#define FLAN_SIM_H


#include "FLAN_Clone.h"
using namespace Rcpp ;



/* Simulation class */

class FLAN_Sim {
public:

    /*! \brief  create an object */

    FLAN_Sim(){
      if(!mDist) delete mDist;
      mDist=NULL;

      if(!mClone) delete mClone;
      mClone=NULL;
    };

    FLAN_Sim(List args){

      mMut=as<double>(args["mutations"]);
      mFitness=as<double>(args["fitness"]);
      mDeath=as<double>(args["death"]);

      List dist=args["dist"];

      if(!mDist) delete mDist;
      mDist=new FLAN_Dist(dist);          // Lifetime Distribution

      std::string distfn=args["distfn"]; // Final Count Distribution
      mDistfn=distfn;
      // std::cout<<"mDistfn :"<<mDistfn<<std::endl;

      mMfn=as<double>(args["mfn"]);
      mCvfn=as<double>(args["cvfn"]);

      if(!mClone) delete mClone;
      mClone=new FLAN_SimClone(mFitness,mDeath,mDist);

    };

    // DESTRUCTORS

    /*! \brief  destroy an object.
     */
    ~FLAN_Sim(){
      if(!mClone) delete mClone;
      if(!mDist) delete mDist;
    };



private:

  double mMut;
  double mFitness;
  double mDeath;
  FLAN_SimClone *mClone=NULL;
  FLAN_Dist *mDist=NULL;
  std::string mDistfn;
  double mMfn;
  double mCvfn;

    // --------------------------------------------------------------------------------------------------------------
    //////////////////////////////////////////// Samples computation ////////////////////////////////////////////////
    // --------------------------------------------------------------------------------------------------------------
    /* Generates a sample of size n of mutant counts
     * n: sample size
     * alpha: mean number of mutations
     * rho: fitness parameter
     * death: probability of death
     */

    NumericVector computeSampleMutantsNumber(int n) ;

    /* Generates a sample of size n of mutant counts, given a sample of final counts.
     *
     * n: number of experiments
     * pi: mutation probability
     * rho: fitness parameter
     * death: probability of death
     * finalCount:
     */
    NumericVector computeSampleMutantsNumber(int n,
				    NumericVector& finalCount) ;
public:
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
    List computeSamplesMutantsFinalsNumber(int n)  ;

};


/* Simulation class for inhomoegenous models */

class FLAN_SimInhomogeneous {

private:

  double mMut;
  double mFitness;
  double mDeath;
  FLAN_SimInhomogeneousClone *mClone=NULL;
  Function* mMU=NULL;
  Function* mMUinv0=NULL;
  double mMfn;
  double mCvfn;

    // --------------------------------------------------------------------------------------------------------------
    //////////////////////////////////////////// Samples computation ////////////////////////////////////////////////
    // --------------------------------------------------------------------------------------------------------------
    /* Generates a sample of size n of mutant counts
     * n: sample size
     * alpha: mean number of mutations
     * rho: fitness parameter
     * death: probability of death
     */

    NumericVector computeSampleMutantsNumber(int n) ;

    /* Generates a sample of size n of mutant counts, given a sample of final counts.
     *
     * n: number of experiments
     * pi: mutation probability
     * rho: fitness parameter
     * death: probability of death
     * finalCount:
     */
    NumericVector computeSampleMutantsNumber(int n,
				    NumericVector& finalCount) ;

public:

    /*! \brief  create an object */

    FLAN_SimInhomogeneous(){
      if(!mMU) delete mMU;
      mMU=NULL;
      if(!mMUinv0) delete mMUinv0;
      mMUinv0=NULL;
      if(!mClone) delete mClone;
      mClone=NULL;
    };

    FLAN_SimInhomogeneous(List args){

      mMut=as<double>(args["mutations"]);
      mFitness=as<double>(args["fitness"]);
      mDeath=as<double>(args["death"]);

      List muih=args["muih"];
      if(!mMU) delete mMU;
      if(!mMUinv0) delete mMUinv0;

      mMU=new Function("identity");
      mMUinv0=new Function("identity");

//       std::cout<<"Read fih"<<std::endl;
      *mMU=muih["mu"];
      *mMUinv0=muih["muinv0"];
                // Lifetime Distribution


      mMfn=as<double>(args["mfn"]);
      mCvfn=as<double>(args["cvfn"]);

      if(!mClone) delete mClone;
//       std::cout<<"Call constructor of IhClone"<<std::endl;
      mClone=new FLAN_SimInhomogeneousClone(mFitness,mDeath,mMU);

    };

    // DESTRUCTORS

    /*! \brief  destroy an object.
     */
    ~FLAN_SimInhomogeneous(){
      if(!mMU) delete mMU;
      if(!mMUinv0) delete mMUinv0;
      if(!mClone) delete mClone;
    };


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
    List computeSamplesMutantsFinalsNumber(int n)  ;

};

#endif
