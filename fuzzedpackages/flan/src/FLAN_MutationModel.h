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

#ifndef FLAN_MUTATION_MODEL_H
#define FLAN_MUTATION_MODEL_H

#include "FLAN_Clone.h"
using namespace Rcpp ;


class FLAN_MutationModel {

private:

    FLAN_Clone* mClone=NULL;    // Clone object

protected:

    double mMutNumber;     // Mean number of mutations
    double mFitness;       // Relative fitness
    double mDeath;         // Death probability
    double mPlateff;       // Plating efficiency

    double mMuinf;       // Plating efficiency

//     bool mLT;              // logical: if TRUE probabilities are P[X <= x]
			   //   otherwise, P[X > x]

    NumericVector mSample;
    double mMfn,mCvfn;
    double mScale;


public:

    //  Create an object
    FLAN_MutationModel(){};

    // Create object for GF method
    FLAN_MutationModel(double death,std::string model){

      mDeath=death;
      FLAN_Clone* clone=NULL;

      if(model.compare("H") == 0) clone=new FLAN_DiracClone(death);
      else clone=new FLAN_ExponentialClone(death);


      mClone=clone;
    //   mZ4=0.55;
//       mLT=true;
    };

    // Create object for GF method
    FLAN_MutationModel(double rho,double death,std::string model){

      mFitness=rho;
      mDeath=death;

      FLAN_Clone* clone=NULL;
      if(model.compare("H") == 0) clone=new FLAN_DiracClone(rho,death);
      else clone=new FLAN_ExponentialClone(rho,death);


      mClone=clone;

//       mLT=true;
    };

    // Create object for distribution
    FLAN_MutationModel(List args) {

      if(!Rf_isNull(args["mutations"]))  mMutNumber=as<double>(args["mutations"]);
      if(!Rf_isNull(args["fitness"])) mFitness=as<double>(args["fitness"]);
      if(!Rf_isNull(args["death"])) mDeath=as<double>(args["death"]);
      if(!Rf_isNull(args["muinf"])) mMuinf=as<double>(args["muinf"]);
      if(!Rf_isNull(args["plateff"])) mPlateff=as<double>(args["plateff"]);

      std::string model=args["model"];

      // if(args.size()==3){
      // 	mMutNumber=as<double>(args["mutations"]);
      // 	mFitness=as<double>(args["fitness"]);
      // 	mDeath=as<double>(args["death"]);
      //
      // } else if(args.size()==5){
	    //  mMutNumber=as<double>(args["mutations"]);
	    //  mFitness=as<double>(args["fitness"]);
	    //  mDeath=as<double>(args["death"]);
      //
      //  mMuinf=as<double>(args["muinf"]);
      //
    	// std::string model=args["model"];

      if(model.compare("N") != 0) {
      	FLAN_Clone* clone=NULL;
      	// if(model.compare("H") == 0) clone=new FLAN_DiracClone(mFitness,mDeath);
        // // else clone=new FLAN_ExponentialClone(mFitness,mDeath,args["integrands"]);
      	// else if(model.compare("LD") == 0) clone=new FLAN_ExponentialClone(mFitness,mDeath);
        // else clone=new FLAN_InhomogeneousClone(mFitness,mDeath,mMuinf);

        if(model.compare("H") == 0) clone=new FLAN_DiracClone(args);
        // else clone=new FLAN_ExponentialClone(mFitness,mDeath,args["integrands"]);
      	else if(model.compare("LD") == 0) clone=new FLAN_ExponentialClone(args);
        else clone=new FLAN_InhomogeneousClone(args);

      	mClone=clone;
      }
// 	mLT=as<bool>(args["lt"]);

      // } else if(args.size()==6){
      // 	mMutNumber=as<double>(args["mutations"]);
      // 	mFitness=as<double>(args["fitness"]);
      // 	mDeath=as<double>(args["death"]);
      // 	mPlateff=as<double>(args["plateff"]);
      //
      //   mMuinf=as<double>(args["muinf"]);
      //
      // 	std::string model=args["model"];
      //
      // 	FLAN_Clone* clone=NULL;
      // 	if(model.compare("H") == 0) clone=new FLAN_DiracClone(mFitness,mDeath);
      //   // else clone=new FLAN_ExponentialClone(mFitness,mDeath,args["integrands"]);
      // 	// else if(model.compare("LDpef") == 0) clone=new FLAN_ExponentialClone(mFitness,mDeath,mPlateff);
      // 	else if(model.compare("LD") == 0) clone=new FLAN_ExponentialClone(mFitness,mDeath,mPlateff);
      //   else if(model.compare("I") == 0) clone=new FLAN_InhomogeneousClone(mFitness,mDeath,mPlateff,mMuinf);
      //
      // 	mClone=clone;
      //
      // } else if (args.size()==9){

        if(!Rf_isNull(args["mc"])) mSample=args["mc"];
      	if(!Rf_isNull(args["mfn"])){
      	  mMfn=as<double>(args["mfn"]);
      	  mCvfn=as<double>(args["cvfn"]);
      	} else {
      	  mMfn=-1;
      	  mCvfn=-1;
      	}
      //
      // 	mFitness=as<double>(args["fitness"]);
      // 	mDeath=as<double>(args["death"]);
      // 	mPlateff=as<double>(args["plateff"]);
      //
      //   mMuinf=as<double>(args["muinf"]);
      //
      // 	std::string model=args["model"];
      //
    	// FLAN_Clone* clone=NULL;
    	// if(model.compare("H") == 0)clone=new FLAN_DiracClone(mFitness,mDeath);
      // else if(model.compare("LD") == 0) clone=new FLAN_ExponentialClone(mFitness,mDeath);
      // else clone=new FLAN_InhomogeneousClone(mFitness,mDeath,mMuinf);
      //
    	// mClone=clone;

    	if(!Rf_isNull(args["scale"])) mScale=as<double>(args["scale"]);

      // }

    };

    // Destructor
    ~FLAN_MutationModel(){};


    // Set attributes
    void setMutNumber(double alpha) {
      mMutNumber=alpha;
    };
    void setFitness(double rho){
      mFitness=rho;
    };
    void setDeath(double death){
      mDeath=death;
    };
    void setClone(FLAN_Clone* clone){
      mClone=clone;
    };


    // Get attributes

    double getMutNumber(){
      return mMutNumber;
    };
    double getFitness(){
      return mFitness;
    };
    double getDeath(){
      return mDeath;
    };
    FLAN_Clone* getClone(){
      return mClone;
    };



    // --------------------------
    // Probability methods
    //---------------------------


    NumericVector computeProbability(int m) ;

    NumericVector deduceProbability(int m,NumericVector& pClone) ;


    List computeProbability1DerivativeAlpha(int m)  ;
    List deduceProbability1DerivativeAlpha(int m,
					    NumericVector& pClone)  ;


    List computeProbability1DerivativeRho(int m)  ;
    List deduceProbability1DerivativeRho(int m,
					  NumericVector& pClone,
					  NumericVector& dpClone_r) ;


    List computeProbability1DerivativesAlphaRho(int m)  ;
    List deduceProbability1DerivativesAlphaRho(int m,
						NumericVector& pClone,
						NumericVector& dpClone_r)  ;


    NumericVector computeCumulativeFunction(int m,bool lower_tail) ;

//     // --------------------
//     // GF ESTIMATION covariance methods
//     // -------------------

    List MutationGFEstimation(bool WithSd);

    double computeGeneratingFunction(double z);

//
    double covariance2(double z1, double z2) ;
//
    NumericVector covariance(double z1,double z2,double z3) ;
//
//
//     // Unbiased estimation of pi and its standart deviation if fluctuation of final counts
    List unbiasPiEstimation(double sd, double z,
			    double mfn,double cvfn);

};
#endif
