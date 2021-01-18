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
#ifndef FLAN_CLONE_H
#define FLAN_CLONE_H

#include "MATH_Function.h"

using namespace Rcpp ;

// Lifetime Distributions class
class FLAN_Dist {

private:

  std::string mName;                // Name of the distribution
  std::vector<double> mParams;      // Parameter(s) of the distribution

public:

  FLAN_Dist(){};

  FLAN_Dist(List dist){
    std::string name=dist["name"];
    mName=name;
    if(mName.compare("lnorm") == 0) {
      mParams.resize(2);
      mParams[0]=dist["meanlog"];mParams[1]=dist["sdlog"];
    } else if(mName.compare("gamma") == 0) {
      mParams.resize(2);
      mParams[0]=dist["shape"];mParams[1]=dist["scale"];
    // } else if(mName.compare("exp") == 0) {
    //   mParams.resize(1);
    //   mParams[0]=dist["rate"];
    // } else if(mName.compare("dirac") == 0) {
    //   mParams.resize(1);
    //   mParams[0]=dist["location"];
    }
  };

  FLAN_Dist(std::string name, std::vector<double> params) {
    mName=name;
    mParams=params;
  }

  ~FLAN_Dist(){};

  std::string getDistName() {
    return mName;
  }

  std::vector<double> getDistParams() {
    return mParams;
  }

};


  /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_SimClone class to compute sample of clone size
  */

class FLAN_SimClone {

private:

    double mFitness;
    double mDeath;
    FLAN_Dist *mDist=NULL;   // Lifetime distribution
    FLAN_Dist *mDistN=NULL;  // Lifetime distribution for normal cells (drawing function)

protected:
  static const double DEATH_EPS_SIM;     // Threshold for death

public:
    FLAN_SimClone(){
      if(!mDist) delete mDist;
      mDist=NULL;
      if(!mDistN) delete mDistN;
      mDistN=NULL;
    };

    ~FLAN_SimClone(){
      if(!mDist) delete mDist;
      if(!mDistN) delete mDistN;
    };

    FLAN_SimClone(double rho,double death, List dist){

      mFitness=rho;
      mDeath=death;

      if(!mDist) delete mDist;
      mDist= new FLAN_Dist(dist);
      if(!mDistN) delete mDistN;
      mDistN=NULL;

    };

    FLAN_SimClone(double rho,double death, List distn, List distm){

      mFitness=rho;
      mDeath=death;

      if(!mDist) delete mDist;
      if(!mDistN) delete mDistN;

      mDist= new FLAN_Dist(distm);
      mDistN= new FLAN_Dist(distn);

    };
      // create object to generate samples
    FLAN_SimClone(double rho,double death, FLAN_Dist *dist){

      mFitness=rho;
      mDeath=death;

      if(!mDist) delete mDist;
      if(!mDistN) delete mDistN;

      mDist=dist;
      mDistN=NULL;

    };

    // -----------------------
    // Sample computation
    // ------------------------

    NumericVector computeSample(int n);

    int splitTimes(double t);


};

/*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_SimInhomogeneousClone class for samples of inhomogeneous lifetimes
  */

class FLAN_SimInhomogeneousClone {

private:

    double mFitness;
    double mDeath;

    Function* mMU=NULL;

protected:
  static const double DEATH_EPS_SIM;     // Threshold for death

public:
    FLAN_SimInhomogeneousClone(){

      if(!mMU) delete mMU;

      mMU=NULL;

    };

    ~FLAN_SimInhomogeneousClone(){

      if(!mMU) delete mMU;

      mMU=NULL;

    };


    FLAN_SimInhomogeneousClone(double rho,double death, List muih){

      mFitness=rho;
      mDeath=death;

      if(!mMU) delete mMU;

      // allocations
      mMU=new Function("identity");

      *mMU=muih["mu"];


    };

    FLAN_SimInhomogeneousClone(double rho,double death, Function* muih){

      mFitness=rho;
      mDeath=death;

      if(!mMU) delete mMU;

      // allocations
      mMU=new Function("identity");

      mMU=muih;


    };

    // -----------------------
    // Sample computation
    // ------------------------

    NumericVector computeSample(int n,double s);

};



  /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class for the distribution of a clone size
  */


class FLAN_Clone {

private:


protected:
    static const double DEATH_EPS_DIST;     // Threshold for death
    double mFitness;    // Realtive fitness
    double mDeath;      // Death probability
//     double mPlateff;

    FLAN_Clone(){};

    // create object for GF method

    FLAN_Clone(List params){

      // mDeath=as<double>(params["death"]);

      // if(params.size() >= 2) {
      if(!Rf_isNull(params["death"]))  mDeath=as<double>(params["death"]);
	    if(!Rf_isNull(params["fitness"]))  mFitness=as<double>(params["fitness"]);
      // }
      // std::cout<<"Fitness ="<<mFitness<<std::endl;
      // std::cout<<"Death ="<<mDeath<<std::endl;
    };

    FLAN_Clone(double death){
      mDeath=death;
    };

    // create object to compute distribution
    FLAN_Clone(double rho,double death){
      mFitness=rho;
      mDeath=death;
    };


    ~FLAN_Clone(){};


    // Set attributes
    void setFitness(double rho){
      mFitness=rho;
    }
    void setDeath(double death){
      mDeath=death;
    }


    // Get attributes

    double getFitness(){
      return mFitness;
    }
    double getDeath(){
      return mDeath;
    }



public:

    virtual NumericVector computeProbability(int m) = 0;


    virtual List computeProbability1DerivativeRho(int m) = 0;


    /*/////////////////////////////////////////////////////////////////////////////////////////////////////
     * Function for GF estimates
     */////////////////////////////////////////////////////////////////////////////////////////////////////


    virtual double computeGeneratingFunction(double z) {
      std::vector<double> Z (1);
      Z[0]=z;
      return computeGeneratingFunction2(mFitness,Z)[0];
    };

    virtual std::vector<double> computeGeneratingFunction2(double rho,std::vector<double> Z) = 0;

    virtual double computeGeneratingFunction1DerivativeRho(double z) = 0;

};



 /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the lifetime model is suppsoed to be exponential
  */


class FLAN_ExponentialClone : public FLAN_Clone {

  protected:

    double mPlateff;

  private:

    MATH_Integration* mIntegrator=NULL;

    void init() {
      List info=Environment::base_namespace().get(".Machine");
      double flantol=info["double.eps"];
      flantol=sqrt(flantol);
      int flansubd=1000;

      if(!mIntegrator) delete mIntegrator;
      mIntegrator=new MATH_Integration(flantol,flansubd);
    }

  public:

    FLAN_ExponentialClone():FLAN_Clone(){
      init();
    };

    FLAN_ExponentialClone(List params):FLAN_Clone(params) {
      // if(params.size() == 3) mPlateff=as<double>(params["plateff"]);
      // else mPlateff=1;
      if(!Rf_isNull(params["plateff"]))  mPlateff=as<double>(params["plateff"]);
      init();
    };
    FLAN_ExponentialClone(double death):FLAN_Clone(death) {
      mPlateff=1;
      init();
    };
    FLAN_ExponentialClone(double rho,double death):FLAN_Clone(rho,death) {
      mPlateff=1;
      init();

    };

    FLAN_ExponentialClone(double rho,double death, double plateff):FLAN_Clone(rho,death) {
      mPlateff=plateff;
      init();

    };
    ~FLAN_ExponentialClone(){
      if(!mIntegrator) delete mIntegrator;
    };


    double testIntegral(double a, double b);

    NumericVector computeProbability(int m);

    List computeProbability1DerivativeRho(int m) ;

    std::vector<double> computeGeneratingFunction2(double rho,std::vector<double> Z);

    double computeGeneratingFunction1DerivativeRho(double z)  ;

};

//
//  /*//////////////////////////////////////////////////////////////////////////////////////////////
//   * FLAN_Clone class when the lifetime model is suppsoed to be constant
//   */


class FLAN_DiracClone : public FLAN_Clone {

private:

  void init_death(){
      mPol=MATH_Polynom();
    };

protected:

    MATH_Polynom mPol;

public:

    FLAN_DiracClone():FLAN_Clone() {};
    FLAN_DiracClone(double death):FLAN_Clone(death) {
      if(death > 0) init_death();
    };
    FLAN_DiracClone(double rho,double death):FLAN_Clone(rho,death) {
      if(death > 0) init_death();
    };
    FLAN_DiracClone(List params):FLAN_Clone(params) {
      if(mDeath > 0) init_death();
    };
    ~FLAN_DiracClone(){};

    NumericVector computeProbability(int m);

    List computeProbability1DerivativeRho(int m);

    std::vector<double> computeGeneratingFunction2(double rho,std::vector<double> Z);

    double computeGeneratingFunction1DerivativeRho(double z)  ;

};


/*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the model is supposed to be inhomogeneous
  */

class FLAN_InhomogeneousClone : public FLAN_Clone {

protected:

  double mPlateff;
private:

  MATH_Integration* mIntegrator=NULL;
  double mMuinf;

  void init() {
    List info=Environment::base_namespace().get(".Machine");
    double flantol=info["double.eps"];
    flantol=sqrt(flantol);
    int flansubd=1000;

    if(!mIntegrator) delete mIntegrator;
    mIntegrator=new MATH_Integration(flantol,flansubd);
  }

public:

  FLAN_InhomogeneousClone():FLAN_Clone(){
    if(!mIntegrator) delete mIntegrator;
    mIntegrator=NULL;
  };

  FLAN_InhomogeneousClone(List params):FLAN_Clone(params) {
    // if(params.size() == 4) mPlateff=as<double>(params["plateff"]);
    if(!Rf_isNull(params["plateff"]))  mPlateff=as<double>(params["plateff"]);
    if(!Rf_isNull(params["muinf"]))  mMuinf=as<double>(params["muinf"]);
    init();
  };

  FLAN_InhomogeneousClone(double death,double muinf):FLAN_Clone(death) {
    init();
    mMuinf=muinf;
  };

  FLAN_InhomogeneousClone(double rho,double death,double muinf):FLAN_Clone(rho,death) {
    init();
    mMuinf=muinf;
  };

  FLAN_InhomogeneousClone(double rho,double death,double plateff,double muinf):FLAN_Clone(rho,death) {
    init();
    mPlateff=plateff;
    mMuinf=muinf;
  };


  ~FLAN_InhomogeneousClone(){
    if(!mIntegrator) delete mIntegrator;
  };

  NumericVector computeProbability(int m);

  List computeProbability1DerivativeRho(int m) ;

  std::vector<double> computeGeneratingFunction2(double rho,std::vector<double> Z);

  double computeGeneratingFunction1DerivativeRho(double z)  ;

};
#endif
