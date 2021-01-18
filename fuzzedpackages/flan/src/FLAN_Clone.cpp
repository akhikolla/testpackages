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


#include "FLAN_Clone.h"

const double FLAN_SimClone::DEATH_EPS_SIM=1.e-4;

NumericVector FLAN_SimClone::computeSample(int n) {

  std::string name=mDist->getDistName();

  NumericVector sample=rexp(n,mFitness);


  if(name.compare("dirac") == 0){

    // specialized method
    if (mDeath<DEATH_EPS_SIM) {
        // case without mDeath
      for(NumericVector::iterator it = sample.begin(); it != sample.end(); ++it) {
	  *it = pow(2.,floor(*it/log(2.)));
        }
    } else {
        // case with death
        double t,a=log(2.*(1-mDeath));
	for(NumericVector::iterator it = sample.begin(); it != sample.end(); ++it) {
	  double cg=0;
	  int m=1;
	  int nAlives=0,nDeaths=0;
	  t=*it;
            while ((cg<t-a) && (m>0)) {
		nDeaths=(int)(rbinom(1,m,mDeath)[0]);
                // the number of alives among m
                nAlives=m-nDeaths;
                // the new alive cells
                m+=nAlives-nDeaths;
                // current generation
                cg+=a;
            }
	      *it=m;
        }
    }
   // end if dirac distribution
  } else if(name.compare("exp") == 0) {
      double p,up;
      double geoD,beD;

      if (mDeath<DEATH_EPS_SIM) {
  //         for (int k=0;k<n;k++) {
	for(NumericVector::iterator it = sample.begin();it != sample.end(); ++it) {
	      p=exp(-(*it));

	      geoD=rgeom(1,p)[0];

	      if(geoD>=0) *it=geoD+1;
	      else if (geoD<0 || geoD!=geoD)*it=-1.e-5;
	  }
      } else {
	for(NumericVector::iterator it = sample.begin();it != sample.end(); ++it) {

	    p=exp(-(*it));

	    up=(1-2*mDeath)/(1-mDeath*(1+p));

	    beD=rbinom(1,1,up)[0];

	    if (beD==1){

	      geoD=rgeom(1,p*up)[0];

	      if(geoD>=0) *it=geoD+1;
	      else if (geoD<0 || geoD!=geoD) *it=-1.e-5;

	    } else *it=0;
	}
      }
    // end if exponential distribution
  } else {

    // transform the sample

//     std::vector<double> splitTimesList;
    int nc=0;
    for(NumericVector::iterator it = sample.begin();
	  it != sample.end(); ++it) {
	     nc=splitTimes(*it);
        *it=nc;
    }

  }

  return sample;
}


int FLAN_SimClone::splitTimes(double t) {

    std::string name=mDist->getDistName();
    std::vector<double> params=mDist->getDistParams();

    // return the clone size
    int cloneSize=0;

    //current generation
    int g=0;

    // make a sample of size 1
    //division times of current generation
    std::vector<double> st (1);
    if(name.compare("lnorm") == 0) st=as< std::vector<double> >(rlnorm(1,params[0],params[1]));

    if(name.compare("gamma") == 0) st=as< std::vector<double> >(rgamma(1,params[0],params[1]));

    // take only the times if before t
    // ng is the number of simulated value less than t
    // number of alive cells
    int ng=0, nActiveAlives;
    std::vector<double> gtf;

    // If no cells death
    if(mDeath<DEATH_EPS_SIM){
      if(st[0]<t){
        ng=1;
      } else {
	// no division before t
	ng=0;
	cloneSize++;
      }
      nActiveAlives=ng;
      int i=0;
//       std::cout<<"Main loop"<<std::endl;
      while ((ng>0) && (ng<1.e+6)) {
	g++;
	i++;
	gtf.resize(2*ng);
	if(name.compare("lnorm") == 0) gtf=as< std::vector<double> >(rlnorm(2*ng,params[0],params[1]));
	if(name.compare("gamma") == 0) gtf=as< std::vector<double> >(rgamma(2*ng,params[0],params[1]));

	//  split times of daughters with duplicating split times
        for (int i=0;i<ng;i++) {
            gtf[2*i]+=st[i];
            gtf[2*i+1]+=st[i];
        }

        ng*=2;
        // keep only values less than t
        // daughters that still divide
        // keep their split times

	nActiveAlives=0;
	st.resize(ng);
        for (int i=0;i<ng;i++) {
            if (gtf[i]<t) {
                // division at  t of an alive cell
                st[nActiveAlives]=gtf[i];
                nActiveAlives++;
            } else {
                // no division before t
                cloneSize++;
            }
        }
//         splitTimesList.insert(splitTimesList.end(), st.begin(), st.end());
        // ng is the size of st
        //dividing cells in next generation
        ng=nActiveAlives;
      }
    // If cells death with prob mDeath
    } else {

      NumericVector u=runif(1,0.,1.);
      if((st[0]<t) && (u[0]>mDeath)){
        ng=1;
      } else if (st[0]>=t){
	// no division before t
	ng=0;
	cloneSize++;
      } else if(u[0] <= mDeath) {
	// the cell dies out
	ng=0;
      }
      nActiveAlives=ng;

      while ((ng>0) && (ng<1.e+6)) {
	g++;
// 	NumericVector gtf;
	gtf.resize(2*ng);
	if(name.compare("lnorm") == 0) gtf=as< std::vector<double> >(rlnorm(2*ng,params[0],params[1]));
	if(name.compare("gamma") == 0) gtf=as< std::vector<double> >(rgamma(2*ng,params[0],params[1]));

	//  split times of daughters with duplicating split times
        for (int i=0;i<ng;i++) {
            gtf[2*i]+=st[i];
            gtf[2*i+1]+=st[i];
        }

        ng*=2;

	u=runif(ng,0.,1.);

	nActiveAlives=0;
	st.resize(ng);
	for (int i=0;i<ng;i++) {

            if ((gtf[i]<t) && (u[i]>mDeath)) {
                // division at  t of an alive cell
                st[nActiveAlives]=gtf[i];
                nActiveAlives++;
            } else if (gtf[i]>=t) {
                // no division before t
                cloneSize++;
            }
        }
        // ng is the size of st
        //dividing cells in next generation
        ng=nActiveAlives;
      }
  }

    return cloneSize;
}


/*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_SimInhomogeneousClone class for samples of inhomogeneous lifetimes
  */

const double FLAN_SimInhomogeneousClone::DEATH_EPS_SIM=1.e-4;

NumericVector FLAN_SimInhomogeneousClone::computeSample(int n,double s) {


  NumericVector sample(n);
  double p;

  //double fs=as<double>((*mF)(s));
  //double fn=as<double>((*mF)(R_PosInf));
  double mu=as<double>((*mMU)(s,R_PosInf));
  double geoD,beD;
  double up;
//   p=pow(fs/fn,(1-2*mDeath));
  p=exp(-mu*(1-2*mDeath));

    if (mDeath<DEATH_EPS_SIM) {
//         for (int k=0;k<n;k++) {
      for(NumericVector::iterator it = sample.begin();it != sample.end(); ++it) {
	  geoD=rgeom(1,p)[0];
	  if(geoD>=0) *it=geoD+1;
	  else if (geoD<0 || geoD!=geoD) *it=-1.e-5;
	}
    } else {
      for(NumericVector::iterator it = sample.begin();it != sample.end(); ++it) {

	up=(1-2*mDeath)/(1-mDeath*(1+p));

	beD=rbinom(1,1,up)[0];

	if (beD==1){

	  geoD=rgeom(1,p*up)[0];

	  if(geoD>=0) *it=geoD+1;
	  else if (geoD<0 || geoD!=geoD) *it=-1.e-5;
	} else *it=0;
      }
    }

  return sample;
}



/*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class for the distribution of a clone size
  */
const double FLAN_Clone::DEATH_EPS_DIST=1.e-4;

 /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the lifetime model is supposed to be exponential
  */


// create object for GF method
//
// double FLAN_ExponentialClone::testIntegral(double a, double b){
//   MATH_Params params;
//   params.rho=mFitness;
//   double d1=mDeath/(1-mDeath);
//   params.delta=d1;
//
//   mIntegrator->setFunction("TEST",&params);
//   //integrate the function in [0,1]
//   double I;
//   // m=0 probability
//   // std::cout<<"Call of computeIntegral"<<std::endl;
//   I=mIntegrator->computeIntegral(a,b);
//
//   return I;
// }

NumericVector FLAN_ExponentialClone::computeProbability(int m){

//   std::cout<<"Compute CloneExp probabilities"<<std::endl;

  std::vector<double> P(m+1);

  if(mPlateff < 1){

    MATH_Params params;
    params.rho=mFitness;
    params.delta=mDeath;
    params.zeta=mPlateff;
    params.k=0.;
    mIntegrator->setFunction("CLONE_P0_WD_WPEF",&params);

    //integrate the function in [0,1]
    double I;
    I=mIntegrator->computeIntegral(0.,1.);

    P[0]=I*mFitness;

    if(m > 0){

//       mIntegrator->setFunctionName("CLONE_P1_WD_WPEF");
//
//       I=mIntegrator->computeIntegral(0.,1.);

//       P[1]=mFitness*I;

      std::vector<double>::iterator it=P.begin()+1;
//       double pkm1=P[1];
//       int m1=m;
//       if (m1>=m_max) m1=m_max;

      for (int k=1;k<=m;k++,++it) {
	params.k=k;
	mIntegrator->setFunction("CLONE_PK_WD_WPEF",&params);
	I=mIntegrator->computeIntegral(0.,1.);
	*it=mFitness*I;
      }
    }
  } else {
    if(mDeath<DEATH_EPS_DIST){
//     std::cout<<"Compute CloneExp probabilities (death=0)"<<std::endl;

      P[0]=0;
      if (m > 0){
	int k=1;
	for(std::vector<double>::iterator it=P.begin()+1 ; it!= P.end() ; ++it,k++) *it=mFitness*R::beta(mFitness+1,k);
      }
    } else {
  //         std::cout<<"Compute CloneExp probabilities (death>0)"<<std::endl;

      double d1=mDeath/(1-mDeath);
      int m_max=1000;
      // std::cout<<"Set the integrand and the parameters"<<std::endl;
      MATH_Params params;
      params.rho=mFitness;
      params.delta=d1;
      params.zeta=1;
      params.k=0.;
      mIntegrator->setFunction("CLONE_P0_WD",&params);
      //integrate the function in [0,1]
      double I;
      // m=0 probability
      // std::cout<<"Call of computeIntegral"<<std::endl;
      I=mIntegrator->computeIntegral(0.,1.);
      // std::cout<<"DONE : I ="<<I<<std::endl;
      P[0]=I*d1*mFitness;

      if(m > 0){

	double d2=(1.-2.*mDeath)/(1-mDeath);
	d2*=d2;

	//m>0 probability
	int m1=m;
	if (m1>=m_max) m1=m_max;
	std::vector<double>::iterator it=P.begin()+1;
	// MATH_Params params;
	// params.rho=mFitness;
	// params.delta=d1;
	for (int k=1;k<=m1;k++,++it) {
	params.k=k;
	mIntegrator->setFunction("CLONE_PK_WD",&params);
	      I=mIntegrator->computeIntegral(0.,1);
	      *it=I*d2*mFitness;
	}

	// equivalent computation
	double a=pow(d2,(1.-mFitness)/2.)*mFitness*R::gammafn(mFitness+1);
	for (int k=m1+1;k<=m;k++,++it) {
	  *it=a*pow((double)(k),-mFitness-1);
	}
      }
  }
 }
  return NumericVector(P.begin(),P.end());

}

List FLAN_ExponentialClone::computeProbability1DerivativeRho(int m){

  std::vector<double> P(m+1);
  std::vector<double> dP_dr(m+1);
//   double pk;

  if(mPlateff < 1){

        MATH_Params params;
    params.rho=mFitness;
    params.delta=mDeath;
    params.zeta=mPlateff;
    params.k=0.;
    mIntegrator->setFunction("CLONE_P0_WD_WPEF",&params);

    //integrate the function in [0,1]
    double I,Idr;
    I=mIntegrator->computeIntegral(0.,1.);
    mIntegrator->setFunctionName("CLONE_dP0_dr_WD_WPEF");
    Idr=mIntegrator->computeIntegral(0.,1.);


    dP_dr[0]=I+mFitness*Idr;
    P[0]=I*mFitness;


    if(m > 0){

      std::vector<double>::iterator it=P.begin()+1;
      std::vector<double>::iterator itdP=dP_dr.begin()+1;
//       double pkm1=P[1];
//       int m1=m;
//       if (m1>=m_max) m1=m_max;

      for (int k=1;k<=m;k++,++it,++itdP) {
	params.k=k;
	mIntegrator->setFunction("CLONE_PK_WD_WPEF",&params);
	I=mIntegrator->computeIntegral(0.,1.);

	mIntegrator->setFunctionName("CLONE_dPK_dr_WD_WPEF");
	Idr=mIntegrator->computeIntegral(0.,1.);

	*itdP=I+mFitness*Idr;
	*it=mFitness*I;
      }
    }
  } else {
    if(mDeath<DEATH_EPS_DIST){

      P[0]=0;
      dP_dr[0]=0;
      if(m == 0) return List::create(_["P"]=P[0],_["dP_dr"]=dP_dr[0]);

  //     int k=1;
      int k=1;
      std::vector<double>::iterator itdP=dP_dr.begin()+1 ;
      double dg=R::digamma(mFitness+1);
      for(std::vector<double>::iterator itP=P.begin()+1 ; itP!= P.end(); ++itP,++itdP,k++){
	*itP=mFitness*R::beta(mFitness+1,k);
	*itdP=(*itP)*(1/mFitness+dg-R::digamma(mFitness+k+1));
      }
    } else {
      double d1=mDeath/(1-mDeath);

      double I;

      int m_max=1000;
      // set the function type
      MATH_Params params;
      params.rho=mFitness;
      params.delta=d1;
      params.zeta=1;
      params.k=0.;
      mIntegrator->setFunction("CLONE_P0_WD",&params);
      I=mIntegrator->computeIntegral(0.,1.);
      P[0]=I*d1*mFitness;

      mIntegrator->setFunctionName("CLONE_dP0_dr_WD");
      //integrate the function in [0,1]
      // m=0 probability
      I=mIntegrator->computeIntegral(0.,1.);
      dP_dr[0]=I*d1*mFitness+P[0]/mFitness;

      if (m==0) return List::create(_["P"]=P[0],_["dP_dr"]=dP_dr[0]);


      double d2=(1.-2.*mDeath)/(1-mDeath);
      d2*=d2;
      //m>0 probability
      int m1=m;
      if (m1>=m_max) m1=m_max;
      std::vector<double>::iterator itP=P.begin()+1 ;
      std::vector<double>::iterator itdP=dP_dr.begin()+1 ;
      for (int k=1;k<=m1;k++,++itP,++itdP) {
	params.k=k;
	  mIntegrator->setFunction("CLONE_PK_WD",&params);
	  I=mIntegrator->computeIntegral(0.,1);
	  *itP=I*d2*mFitness;

	  mIntegrator->setFunctionName("CLONE_dPK_dr_WD");
	  I=mIntegrator->computeIntegral(0.,1.);
	  *itdP=(*itP)/mFitness+I*d2*mFitness;
      }

      // equivalent computation
      double gfn=R::gammafn(mFitness+1);
      double rhodg=mFitness*R::digamma(mFitness+1);
      double a=pow(d2,(1-mFitness)/2.);
      double b=a*mFitness*R::gammafn(mFitness+1);
      double c=-0.5*log(d2)*mFitness+1.;
      double kpr;
      double kd;
      for(int k=m1+1; k<=m ;k++ , ++itP,++itdP){
	kd=(double)(k);
	kpr=pow(kd,-mFitness-1);
	*itP=a*kpr;

	*itdP=b*kpr*(
	      gfn*(c-mFitness*log(kd))
	      +rhodg);
      }
    }
  }
  return List::create(_["P"]=NumericVector(P.begin(),P.end()),
		      _["dP_dr"]=NumericVector(dP_dr.begin(),dP_dr.end()));

}


std::vector<double> FLAN_ExponentialClone::computeGeneratingFunction2(double rho,std::vector<double> Z){

  double eps =1e-8;

  std::vector<double> H(Z.size());
  std::vector<double>::iterator itH=H.begin();

  for (std::vector<double>::iterator itZ = Z.begin() ; itZ != Z.end() ; ++itZ, ++itH){
    if(fabs(*itZ) <= eps) *itH=0.;
    else if(fabs(1-(*itZ)) <= eps) *itH=1.;
    else {
      double dstar,zstar,I;

      dstar=mDeath/(1-mDeath);
      zstar=((*itZ)-dstar)/(1-(*itZ));
      MATH_Params params;
      params.rho=rho;
      params.delta=zstar;
      params.zeta=1.;
      params.k=0.;
      mIntegrator->setFunction("CLONE_PGF",&params);
      I=mIntegrator->computeIntegral(0.,1.);

      *itH=dstar+I*zstar*(1-dstar)*rho;


    }
  }

  return H;

}


double FLAN_ExponentialClone::computeGeneratingFunction1DerivativeRho(double z){

  double eps =1e-8;

  if(fabs(z) <= eps) return 0.;
  else if(fabs(1-z) <= eps) return 0.;
  else {
    double dstar,zstar,I1,I2,h;

    dstar=mDeath/(1-mDeath);
    zstar=(z-dstar)/(1-z);
    MATH_Params params;
    params.rho=mFitness;
    params.delta=zstar;
    params.zeta=1.;
    params.k=0.;

    mIntegrator->setFunction("CLONE_PGF",&params);
    I1=mIntegrator->computeIntegral(0.,1.);

    mIntegrator->setFunctionName("CLONE_dPGF_dr");
    I2=mIntegrator->computeIntegral(0.,1.);

    h=zstar*(1-dstar)*(I1+I2*mFitness);

    return h;
  }

}





/*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the lifetime model is suppsoed to be constant
  */

NumericVector FLAN_DiracClone::computeProbability(int m){

  std::vector<double> P(m+1);

  if(mDeath<DEATH_EPS_DIST){

    P[0]=0;

    if(m > 0){
      int index=1;
      int n=floor(log((double)(m))/log(2.));
      double tp = (1-pow(2.,-mFitness));
        //take only the indices index=2^ind for ind in [0,n] where index < m
      for (int k=0;k<=n;k++) {
	P[index]=tp*pow(index,-mFitness);
	index*=2;
      }
    }

  } else {

    // case with death >0
    // double eps=1.e-6;
    double deps=1.e-8;

    int itmax=19;  // Arbitrary choice

    // initiate the polynome X
    mPol.setDegree(1);

    // minimum number of iterations
    int i=0;
    double umd=1-mDeath;
    double a=log(2.*umd);
    double t=exp(-mFitness*a);
    double ti=1;//t^i
    int Pol_degree=1;

    double sumP;
    // double sumP_old;


    // initialize P
    P[0]=0;P[1]=1;
    sumP=1;
    // double err = 1;
    int dmax;
    int k;
    // loop
    std::vector<double>::iterator itP;
    // while (((Pol_degree<=m) || (err>eps)) && (i<itmax)) {
    while(i < itmax){
	// next iteration
	i++;
	ti*=t;
	// sumP_old=sumP;

	// update P

	mPol.square_fft();

	mPol*=umd;
	mPol+=mDeath;

	mPol.reduce(deps);
	Pol_degree=mPol.getDegree();
	// update pk=sum_{i=1}^{i=number of polynoms computer} Pi[k] t^i
	dmax=(m<Pol_degree)?m:Pol_degree;
	sumP=0;
	itP=P.begin();
	for (k=0;k<=dmax;k++,++itP) {
	    *itP+=mPol[k]*ti;
	    sumP+=*itP;
	}

	// next iteration
	// err = fabs(1-sumP_old/sumP);
  // std::cout<<"écart relatif (it n°"<<i<<") ="<<err<<std::endl;
    } //end loop

    // multiply by (1-t)
    for(itP=P.begin() ; itP!=P.end();++itP) (*itP)*=(1-t);

  }

  return NumericVector(P.begin(),P.end());

}


List FLAN_DiracClone::computeProbability1DerivativeRho(int m){

  std::vector<double> P(m+1);
  std::vector<double> dP_dr(m+1);

  if(mDeath<DEATH_EPS_DIST){

    P[0]=0;
    dP_dr[0]=0;
    if(m == 0) return List::create(_["P"]=P[0],_["dP_dr"]=dP_dr[0]);


    int n=floor(log((double)(m))/log(2.));
    int index=1;
    double tp=(1-pow(2.,-mFitness));
    double logindex=0;
    double log2=log(2.);
    double tp2=log2*(1/tp-1);

    for(int k=0;k<=n;k++) {
//       index=(int)(pow(2,k));
      P[index]=tp*pow(index,-mFitness);
//       dP_dr[index]=log(2.)*(pow(2.,-mFitness*(k+1))-k*P[index]);
      dP_dr[index]=P[index]*(tp2-logindex);
      index*=2;
      logindex+=log2;
    }

  } else {
     // case with death >0
    // double eps=1.e-6;
    double deps=1.e-8;

    int itmax=20;

    // initiate the polynome X
    mPol.setDegree(1);


    // minimum number of iterations
    int i=0;
    double umd=1-mDeath;
    double a=log(2*umd);
    double t=exp(-mFitness*a);
    double umt=1-t;
    double ti=1;//t^i
    double iti=0; // i*t^i
    int Pol_degree=1;

    double sumP;
    // double sumP_old;


    // initialize P
    P[0]=0;P[1]=1;
    dP_dr[0]=0;dP_dr[1]=t-umt;
    sumP=1;
    // double err = 1;
    int dmax;
    int k;
    // loop
    std::vector<double>::iterator itP;
    std::vector<double>::iterator itdP;
    // while (((Pol_degree<=m) || (err>eps)) && (i<itmax)) {
    while(i < itmax){
	// next iteration
	i++;
	ti*=t;
	iti=i*ti;
	// sumP_old=sumP;

	// update P
	mPol.square_fft();

	// update P
	mPol*=(1-mDeath);
	mPol+=mDeath;
	mPol.reduce(deps);

	Pol_degree=mPol.getDegree();

	dmax=(m<Pol_degree)?m:Pol_degree;
	sumP=0;

	itP=P.begin();
	itdP=dP_dr.begin();
	for (k=0;k<=dmax;k++,++itP,++itdP) {

	  *itP+=mPol[k]*ti;
	  sumP+=*itP;

	  *itdP-=mPol[k]*iti;
	}

      // err = fabs(1-sumP_old/sumP);
    }

    itdP=dP_dr.begin();
    for(itP=P.begin() ; itP!=P.end() ; ++itP,++itdP){
      (*itdP)*=umt;
      (*itdP)+=t*(*itP);
      (*itP)*=umt;
    }

  }

  return List::create(_["P"]=NumericVector(P.begin(),P.end()),
		      _["dP_dr"]=NumericVector(dP_dr.begin(),dP_dr.end()));

}


std::vector<double> FLAN_DiracClone::computeGeneratingFunction2(double rho,std::vector<double> Z){


    std::vector<double> H (Z.size());
    std::vector<double>::iterator itH=H.begin();
    double eps=1.e-8;

    for(std::vector<double>::iterator itZ=Z.begin() ; itZ != Z.end() ; ++itZ, ++itH){
      // z=0 return 0
      if (fabs((*itZ))<eps) *itH=0.;

      // z=1 return 1
      if (fabs(1-(*itZ))<eps) *itH=1.;

      // otherwize
      double s=0;
      if (mDeath<DEATH_EPS_DIST) {
	  double a=pow(2.,-rho);
	  int n=floor(4.-log(fabs(log((*itZ))))/log(2.))+1;

	  for (int k=0;k<=n;k++) {
	      s+=pow((*itZ),pow(2.,k))*pow(a,k);
	  }
	  s*=(1-a);
      } else {
	  double a=log(2.*(1.-mDeath));
	  double dstar=mDeath/(1.-mDeath);
	  int n=floor(-log(eps)/(rho*a))+1;
	  double tp=exp(-rho*a);
	  double tpi=1;
	  double bi=(*itZ);
	  s=(*itZ);
	  for (int i=1;i<=n;i++) {
	      bi=mDeath+(1-mDeath)*bi*bi;
	      tpi*=tp;
	      s+=tpi*bi;
	  }
	  s*=(1-tp);
	  s+=dstar*tpi*tp;
      }
      *itH=s;
    }

  return H;

}
//
//
double FLAN_DiracClone::computeGeneratingFunction1DerivativeRho(double z) {

    double eps=1.e-8;
    // z=0 return 0
    if (fabs(z)<eps) return 0;
    // z=1 return 1
    if (fabs(1-z)<eps) return 0;
    // otherwise

    if (mDeath<DEATH_EPS_DIST) {
        double a=pow(2,-mFitness);
        int n=floor(4.-log(fabs(log(z)))/log(2.))+1;
// 	n++;
        double s1=0,s2=0,t=0;
        for (int k=0;k<=n;k++) {
            t=pow(z,pow(2.,k))*pow(a,k);
            s1+=t;
            s2+=k*t;

        }
        return log(2.)*(a*s1-(1-a)*s2);
    } else {
        double a=log(2.*(1.-mDeath));
        double dstar=mDeath/(1-mDeath);
        int n=floor(-log(eps)/(mFitness*a))+1;

	double tp=exp(-mFitness*a);
        double tpi=1;
        double bi=z;
        double s=z;
        double ds=0;
        for (int i=1;i<=n;i++) {
            bi=mDeath+(1-mDeath)*bi*bi;
            tpi*=tp;
            s+=tpi*bi;
            ds-=i*tpi*bi;
        }
        s=a*(ds*(1-tp)+s*tp);
        // no usefull when n is big
	s+=dstar*tpi*tp;
        return s;
    }
}
//
//
//



/*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the model is supposed to be inhomogeneous
  */
  NumericVector FLAN_InhomogeneousClone::computeProbability(int m){

    // std::cout<<"HERE"<<std::endl;
    std::vector<double> P(m+1);
    std::vector<double>::iterator it;
    int k=1;

    double emuinf=exp(-mMuinf);

    if(mPlateff < 1){
      // std::cout<<"PEF"<<std::endl;
      MATH_Params params;
      params.rho=mFitness;
      params.delta=mDeath;
      params.zeta=mPlateff;
      params.k=0.;
      mIntegrator->setFunction("CLONE_P0_WD_WPEF",&params);

      double I;
      double cste=mFitness/(1-pow(emuinf,mFitness));
      I=mIntegrator->computeIntegral(emuinf,1.);

      P[0]=I*cste;

      if(m > 0){
        std::vector<double>::iterator it=P.begin()+1;
        for (int k=1;k<=m;k++,++it) {
        	params.k=k;
        	mIntegrator->setFunction("CLONE_PK_WD_WPEF",&params);
        	I=mIntegrator->computeIntegral(emuinf,1.);
        	*it=cste*I;
        }
      }
    } else {
      // std::cout<<"NOPEF"<<std::endl;
      if(mDeath<DEATH_EPS_DIST){
        // std::cout<<"NODEATH"<<std::endl;
        P[0]=0;
        if (m > 0){
          double cste=mFitness/(1-pow(emuinf,mFitness));
          for(it=P.begin()+1 ; it!= P.end() ; ++it, k++) {
            // std::cout<<"k = "<<k<<std::endl;
            *it=cste*R::beta(mFitness+1,k)*(1-R::pbeta(emuinf,mFitness+1,k,1,0));
          }
        }
      } else {
        // std::cout<<"DEATH"<<std::endl;
        emuinf=pow(emuinf,1-2*mDeath);
        double d1=mDeath/(1-mDeath);
        int m_max=1000;
        double cste=mFitness/(1-pow(emuinf,mFitness));

        MATH_Params params;
        params.rho=mFitness;
        params.delta=d1;
        params.k=0.;
        mIntegrator->setFunction("CLONE_P0_WD",&params);
        //integrate the function in [mEmuinf,1]
        double I;
        // m=0 probability
        I=mIntegrator->computeIntegral(emuinf,1.);
        P[0]=I*d1*cste;

        if(m > 0){

          double d2=(1.-2.*mDeath)/(1-mDeath);
          d2*=d2;
          //m>0 probability
    //       int m1=m;
          if (m>=m_max) {
    	for(it=P.begin()+1 ; k<=m_max ; ++it,k++)  {
    	  params.k=k;
    	  mIntegrator->setFunction("CLONE_PK_WD",&params);
    	  I=mIntegrator->computeIntegral(emuinf,1.);
    	  *it=I*d2*cste;
    	}

    	// equivalent computation
    	double a=pow(d2,(1.-mFitness)/2.)*R::gammafn(mFitness+1);
            for (k=m_max+1 ; k<=m ; k++,++it) {
    	  *it=a*pow((double)(k),-mFitness-1);
    	  params.k=k;
    	  mIntegrator->setFunction("CLONE_PK_WD",&params);
    	  I=mIntegrator->computeIntegral(0.,emuinf);
    	  (*it)-=d2*I;
    	  (*it)*=cste;
    	}
          } else {
    	for(it=P.begin()+1 ; it!=P.end() ; ++it,k++)  {
    	  params.k=k;
    	  mIntegrator->setFunction("CLONE_PK_WD",&params);
    	  I=mIntegrator->computeIntegral(emuinf,1.);
    	  *it=I*d2*cste;
    	}

          }
        }
     }
   }
  return NumericVector(P.begin(),P.end());

}



  List FLAN_InhomogeneousClone::computeProbability1DerivativeRho(int m){

    std::vector<double> P(m+1);
    std::vector<double> dP_dr(m+1);
    std::vector<double>::iterator itP;
    std::vector<double>::iterator itdP=dP_dr.begin()+1 ;
  //   double pk;
    double emuinf=exp(-mMuinf);
    double enuinf=pow(emuinf,mFitness);
    double cste1,cste2,tp;
  //   double cste2=1./mFitness+log(mEmuinf)*enuinf/(1-enuinf);
    double I,Idr;
    int k=1;

    if(mPlateff < 1){

      MATH_Params params;
      params.rho=mFitness;
      params.delta=mDeath;
      params.zeta=mPlateff;
      params.k=0.;
      mIntegrator->setFunction("CLONE_P0_WD_WPEF",&params);

      //integrate the function in [0,1]
      cste1=mFitness/(1-enuinf);
      cste2=-mMuinf*enuinf/(1-enuinf);
      I=mIntegrator->computeIntegral(emuinf,1.);
      mIntegrator->setFunctionName("CLONE_dP0_dr_WD_WPEF");
      Idr=mIntegrator->computeIntegral(emuinf,1.);

      P[0]=I*cste1;
      dP_dr[0]=cste2*P[0]+cste1*(I/mFitness+Idr);

      if(m > 0){
        std::vector<double>::iterator it=P.begin()+1;
        std::vector<double>::iterator itdP=dP_dr.begin()+1;

        for (int k=1;k<=m;k++,++it,++itdP) {
        	params.k=k;
        	mIntegrator->setFunction("CLONE_PK_WD_WPEF",&params);
        	I=mIntegrator->computeIntegral(emuinf,1.);

        	mIntegrator->setFunctionName("CLONE_dPK_dr_WD_WPEF");
        	Idr=mIntegrator->computeIntegral(emuinf,1.);

          tp=cste1*I;
          *it=tp;
        	*itdP=cste2*tp+cste1*(I/mFitness+Idr);
        }
      }
    } else {
      if(mDeath<DEATH_EPS_DIST){

        P[0]=0;
        dP_dr[0]=0;
        if(m == 0) return List::create(_["P"]=P[0],_["dP_dr"]=dP_dr[0]);
        cste1=mFitness/(1-enuinf);
        cste2=1./mFitness-mMuinf*enuinf/(1-enuinf);
    //     int k=1;

    //     double I;
        double dg=R::digamma(mFitness+1);
        double bk;

        MATH_Params params;

        params.rho=mFitness;
        params.delta=0;

        for(itP=P.begin()+1 ; itP!= P.end(); ++itP,++itdP,k++){
          params.k=k;
    //       *itP=mFitness*R::beta(mFitness+1,k);
    //       *itdP=(*itP)*(1/mFitness+dg-R::digamma(mFitness+k+1));
    //       mIntegrator->setFunctionName("CLONE_PK");
    //       I=mIntegrator->integralFunction(mEmuinf,1.,mFitness,0,k);
    //       *itP=cste1*I;
          bk=R::beta(mFitness+1,k);
            *itP=cste1*bk*(1-R::pbeta(emuinf,mFitness+1,k,1,0));

          mIntegrator->setFunction("CLONE_PK_dr",&params);
          I=mIntegrator->computeIntegral(0.,emuinf);
          *itdP=(*itP)*cste2+cste1*(bk*(dg-R::digamma(mFitness+k+1))-I);
          //       *itdP=(*itP)*cste2+cste1*I;
        }
      } else {

        emuinf=pow(emuinf,1-2*mDeath);
        enuinf=pow(enuinf,1-2*mDeath);
        cste1=mFitness/(1-enuinf);
        cste2=1./mFitness-(1-2*mDeath)*mMuinf*enuinf/(1-enuinf);
        double d1=mDeath/(1-mDeath);
    //     double enuinf=pow(mEmuinf,mFitness);
    //     double cste1=mFitness/(1-enuinf);
    //     double cste2=1./mFitness+log(mEmuinf)*enuinf/(1-enuinf);

    //     double I;
        MATH_Params params;

        params.rho=mFitness;
        params.delta=d1;
        params.k=0;
        int m_max=1000;
        // set the function type
        mIntegrator->setFunction("CLONE_P0_WD",&params);
        I=mIntegrator->computeIntegral(emuinf,1.);
        P[0]=I*d1*cste1;

        mIntegrator->setFunctionName("CLONE_dP0_dr_WD");
        //integrate the function in [0,1]
        // m=0 probability
        I=mIntegrator->computeIntegral(emuinf,1.);
        dP_dr[0]=I*d1*cste1+P[0]*cste2;

        if (m==0) return List::create(_["P"]=P[0],_["dP_dr"]=dP_dr[0]);


        double d2=(1.-2.*mDeath)/(1-mDeath);
        d2*=d2;
        //m>0 probability
    //     int m1=m;
    //     if (m1>=m_max) m1=m_max;
    //     std::vector<double>::iterator itP=P.begin()+1 ;
    //     int k=1;
    //     for (int k=1;k<=m1;k++,++itP,++itdP) {
        if(m>=m_max){
          for(itP=P.begin()+1 ; k<=m_max; ++itP,++itdP,k++){
    	  params.k=k;
    	  mIntegrator->setFunction("CLONE_PK_WD",&params);
    	  I=mIntegrator->computeIntegral(emuinf,1.);
    	  *itP=I*d2*cste1;

    	  mIntegrator->setFunctionName("CLONE_dPK_dr_WD");
    	  I=mIntegrator->computeIntegral(emuinf,1.);
    	  *itdP=(*itP)*cste2+I*d2*cste1;
          }

          // equivalent computation
          double gfn=R::gammafn(mFitness+1);
          double dg=R::digamma(mFitness+1);
          double a=pow(d2,(1-mFitness)/2.)*gfn;
          // double b=a*mFitness;
          double c=dg-0.5*log(d2);
          // double kpr;
          double kd;
          double d;
          for(k=m_max+1 ; k<=m ; k++, ++itP,++itdP){
            kd=(double)(k);
            d=a*pow(kd,-mFitness-1);
    	params.k=k;
    	mIntegrator->setFunction("CLONE_PK_WD",&params);
    	I=mIntegrator->computeIntegral(0.,emuinf);
    	*itP=cste1*(d-d2*I);

    	mIntegrator->setFunctionName("CLONE_dPK_dr_WD");
    	I=mIntegrator->computeIntegral(0.,emuinf);

      *itdP=(*itP)*cste2;
    	(*itdP)+=cste1*(d*(c-log(kd))-d2*I);
          }
        } else {
          for(itP=P.begin()+1 ; itP!= P.end(); ++itP,++itdP,k++){
    	  params.k=k;
    	  mIntegrator->setFunction("CLONE_PK_WD",&params);
    	  I=mIntegrator->computeIntegral(emuinf,1.);
    	  *itP=I*d2*cste1;

    	  mIntegrator->setFunctionName("CLONE_dPK_dr_WD");
    	  I=mIntegrator->computeIntegral(emuinf,1.);
    	  *itdP=(*itP)*cste2+I*d2*cste1;
          }
        }
      }
    }
    return List::create(_["P"]=NumericVector(P.begin(),P.end()),
  		      _["dP_dr"]=NumericVector(dP_dr.begin(),dP_dr.end()));

  }


  std::vector<double> FLAN_InhomogeneousClone::computeGeneratingFunction2(double rho,std::vector<double> Z){

    double eps =1e-8;

    std::vector<double> H(Z.size());
    std::vector<double>::iterator itH=H.begin();

    for (std::vector<double>::iterator itZ = Z.begin() ; itZ != Z.end() ; ++itZ, ++itH){
      if(fabs(*itZ) <= eps) *itH=0.;
      else if(fabs(1-(*itZ)) <= eps) *itH=1.;
      else {
        double dstar,zstar,I,enuinf,emuinf;
        emuinf=exp(-(1-2*mDeath)*mMuinf);
        enuinf=pow(emuinf,rho);
        dstar=mDeath/(1-mDeath);
        zstar=((*itZ)-dstar)/(1-(*itZ));


        MATH_Params params;
        params.rho=rho;
        params.delta=zstar;
  //       params.k=0;


        mIntegrator->setFunction("CLONE_PGF",&params);
        I=mIntegrator->computeIntegral(emuinf,1.);
  //       std::cout<<"I ="<<I<<std::endl;

        *itH=dstar+I*zstar*(1-dstar)*rho/(1-enuinf);


      }
    }

    return H;

  }


  double FLAN_InhomogeneousClone::computeGeneratingFunction1DerivativeRho(double z){

    double eps =1.e-8;

    if(fabs(z) <= eps) return 0.;
    else if(fabs(1.-z) <= eps) return 0.;
    else {
      double dstar,zstar,I1,I2,h,enuinf,emuinf;

  //     std::cout<<"rho ="<<mFitness<<std::endl;
  //     std::cout<<"delta ="<<mDeath<<std::endl;
  //     std::cout<<"muinf ="<<mMuinf<<std::endl;
  //       enuinf=pow(mEmuinf,rho);
      emuinf=exp(-(1.-2.*mDeath)*mMuinf);
      enuinf=pow(emuinf,mFitness);
      dstar=mDeath/(1-mDeath);
      zstar=(z-dstar)/(1-z);
  //     cste=(1-2*mDeath)*mMuinf*enuinf/(1-enuinf);

      MATH_Params params;
      params.rho=mFitness;
      params.delta=zstar;
  //     params.k=0;

      mIntegrator->setFunction("CLONE_PGF",&params);
      I1=mIntegrator->computeIntegral(emuinf,1.);
  //     std::cout<<"I1 ="<<I1<<std::endl;


      mIntegrator->setFunctionName("CLONE_dPGF_dr");
      I2=mIntegrator->computeIntegral(emuinf,1.);

  //     std::cout<<"I2 ="<<I2<<std::endl;

  //     h=zstar*(1-dstar)/(1-enuinf)*(mFitness*(I2-cste*I1)+I1);
      h=zstar*(1.-dstar)/(1.-enuinf)*(mFitness*I2+I1*(1.-mFitness*(1.-2.*mDeath)*mMuinf*enuinf/(1.-enuinf)));

  //     std::cout<<"dh.dr ="<<h<<std::endl;
      return h;
    }

  }
