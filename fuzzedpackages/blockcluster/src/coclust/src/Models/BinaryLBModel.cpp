/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public
 License along with this program; if not, write to the
 Free Software Foundation, Inc.,
 59 Temple Place,
 Suite 330,
 Boston, MA 02111-1307
 USA

 Contact : parmeet.bhatia@inria.fr , bhatia.parmeet@gmail.com
 */


/** @file BinaryLBModel.cpp
 *  @brief Implements concrete model class BinaryLBModel derived from ICoClustModel.
 **/

//#include <limits.h>
//#include <math.h>

#include "BinaryLBModel.h"
#ifndef RPACKAGE
using namespace cimg_library;
#endif
BinaryLBModel::BinaryLBModel( MatrixBinary const&  m_Dataij
                            , ModelParameters const& Mparam
                            , STK::Real a, STK::Real b
                            )
                            : ICoClustModel(Mparam)
                            , a_(a), b_(b)
                            , m_Dataij_(m_Dataij)
                            , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                            , m_Alphakl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_Alphaklold_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_Alphakl1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_Alphakl1old_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_Alphakltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_akl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_epsilonkl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_epsilonkltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
{
#ifdef COVERBOSE
  std::cout << "BinaryLBModel::BinaryLBModel done"<<std::endl;
  std::cout << "nbRow_="<< Mparam_.nbRow_ << std::endl;
  std::cout << "nbCol_="<< Mparam_.nbCol_ << std::endl;
  consoleOut();
#endif
}

BinaryLBModel::BinaryLBModel( MatrixBinary const& m_Dataij
                            , VectorInt    const& rowlabels
                            , VectorInt    const& collabels
                            , ModelParameters const& Mparam
                            , STK::Real a, STK::Real b
                            )
                            : ICoClustModel(Mparam,rowlabels,collabels)
                            , a_(a), b_(b)
                            , m_Dataij_(m_Dataij)
                            , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                            , m_Alphakl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_Alphaklold_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_Alphakl1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_Alphakl1old_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_Alphakltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_akl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_epsilonkl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                            , m_epsilonkltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
{
#ifdef COVERBOSE
  std::cout << "BinaryLBModel::BinaryLBModel done"<<std::endl;
  std::cout << "nbRow_="<<Mparam_.nbRow_ << std::endl;
  std::cout << "nbCol_="<< Mparam_.nbCol_ << std::endl;
  consoleOut();
#endif
}


//Compute Bernoulli log-sum for all rows
void BinaryLBModel::logSumRows(MatrixReal & m_sik)
{
  m_sik = m_Uil_*(((((m_Alphakl_+RealMin)/(((1.-m_Alphakl_))+RealMin)).log())).transpose())
        + STK::Const::VectorX (Mparam_.nbRow_)*((v_logPiek_+(( ((1.-m_Alphakl_))+RealMin).log())*v_Rl_).transpose());
}

//Compute Bernoulli log-sum for all columns
void BinaryLBModel::logSumCols(MatrixReal & m_sjl)
{
  m_sjl = m_Vjk_*( ((m_Alphakl_+RealMin)/((1.-m_Alphakl_)+RealMin)).log())
        + STK::Const::VectorX(Mparam_.nbCol_)*( (v_logRhol_+(( ((1.-m_Alphakl_))+RealMin).log()).transpose()*v_Tk_).transpose());
}


//Run CEM algorithm on data Matrix m_Uil_
bool BinaryLBModel::cemRows()
{
  //Initializations
  computeUil();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //CE-step
    if(!ceStepRows()) return false;
    //M-step
    m_Alphaklold_ = m_Alphakl_;
    mStepRows();
    //Termination check
    if((((m_Alphakl_-m_Alphaklold_).abs()/(m_Alphakl_+RealMin)).sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

//Run EM algorithm on data Matrix m_Uil_
bool BinaryLBModel::emRows()
{
  //Initializations
  computeUil();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!eStepRows()) return false;
    //M-step
    m_Alphaklold_ = m_Alphakl_;
    mStepRows();
    //Termination check
    if((((m_Alphakl_-m_Alphaklold_).abs()/(m_Alphakl_+RealMin)).sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool BinaryLBModel::semRows()
{
  computeUil();
  if(!seStepRows()) return false;
  //M-step : update row proportions and model parameters
  mStepRows();
  return true;
}

bool BinaryLBModel::GibbsRows()
{
  computeUil();
  if(!seStepRows()) return false;
  //M-step : update row proportions and model parameters
  mGibbsStepRows();
  return true;
}

// Run EM algorithm on data matrix m_Vjk_
bool BinaryLBModel::emCols()
{
  //Initializations
  computeVjk();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepCols()) return false;
     //M-step
    m_Alphaklold_ = m_Alphakl_;
    mStepCols();
    //Termination check
    if((((m_Alphakl_-m_Alphaklold_).abs()/(m_Alphakl_+RealMin)).sum())<Mparam_.epsilon_int_)
    { break;}
  }
  // Update Alpha for outer loop
  m_Alphakl1old_ = m_Alphakl1_;
  m_Alphakl1_    = m_Alphakl_;
  return true;
}

// Run CEM algorithm on data matrix m_Vjk_
bool BinaryLBModel::cemCols()
{
  //Initializations
  computeVjk();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //CE-step
    if(!ceStepCols()) return false;
    //M-step
    m_Alphaklold_ = m_Alphakl_;
    mStepCols();
    //Termination check
    if((((m_Alphakl_-m_Alphaklold_).abs()/(m_Alphakl_+RealMin)).sum())<Mparam_.epsilon_int_)
    { break;}
  }
  // Update Alpha for outer loop
  m_Alphakl1old_ = m_Alphakl1_;
  m_Alphakl1_    = m_Alphakl_;
  return true;
}

bool BinaryLBModel::semCols()
{
  computeVjk();
  if(!seStepCols()) return false;
  //M-step: Update column proportions and model parameters
  mStepCols();

  return true;
}

bool BinaryLBModel::GibbsCols()
{
  computeVjk();
  if(!seStepCols()) return false;
  //M-step: Update column proportions and model parameters
  mGibbsStepCols();

  return true;
}

void BinaryLBModel::finalizeOutput()
{
  m_akl_ = (m_Alphakl_>=.5);
  // Calculate probability for Summary Matrix
  m_epsilonkl_ = m_akl_.cast<STK::Real>().prod((1.-m_Alphakl_))
               + (1.-m_akl_.cast<STK::Real>()).prod(m_Alphakl_);
}

void BinaryLBModel::consoleOut()
{
#ifdef COVERBOSE
  std::cout <<"Model parameters:\n"
            <<"Alphakl:\n" << m_Alphakl_
            <<"\nlog(piek): "<< v_logPiek_.transpose()
            <<"\nlog(Rhol): "<< v_logRhol_.transpose() <<"\n";
  std::cout<<"likelihood: "<<likelihood() <<"\n";
  std::cout<<"ICL: "<<iclCriteriaValue() <<std::endl;
#endif
}

/* @return the number of free parmaters of the distribution of a block.*/
int BinaryLBModel::nbFreeParameters() const
{ return 2*Mparam_.nbcolclust_*Mparam_.nbrowclust_;}

STK::Real BinaryLBModel::computeLnLikelihood()
{
  likelihood_ = (v_Tk_.transpose()*
  							( m_Alphakl_.prod((m_Alphakl_+RealMin).log())
								+ ((1.-m_Alphakl_)).prod( ((1.-m_Alphakl_)+RealMin).log() )
								)*v_Rl_
			  				+ v_Tk_.dot(v_logPiek_) + v_Rl_.dot(v_logRhol_)
			          - (m_Tik_.prod( (RealMin + m_Tik_).log()) ).sum()
			          - (m_Rjl_.prod( (RealMin + m_Rjl_).log()) ).sum()
			          ); // / (Mparam_.nbRow_*nbCol_);
  return likelihood_;
}

//Compute change in Alpha and set the terminate variable accordingly
bool BinaryLBModel::initStopCriteria()
{
  return((((m_Alphakl_-m_Alphakltemp_).abs()/(m_Alphakl_+RealMin)).sum()) < Mparam_.epsilon_);
}
void BinaryLBModel::parameterStopCriteria()
{
  stopAlgo_ = ((((m_Alphakl1_-m_Alphakl1old_).abs()/(m_Alphakl1_+RealMin)).sum()) < Mparam_.epsilon_);
}

MatrixBinary const& BinaryLBModel::arrangedDataClusters()
{
  arrangedDataCluster(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

#ifndef RPACKAGE
void BinaryLBModel::displayCluster()
{
  CImg<unsigned char>  cluster(Mparam_.nbCol_,nbRow_,1,1,0);
  CImg<unsigned char>  data(Mparam_.nbCol_,nbRow_,1,1,0);
  m_ClusterDataij_ = arrangedDataClusters();

  // Assign value to images
  for ( int i = 0; i <Mparam_.nbRow_; ++i)
  {
    for ( int j = 0; j < Mparam_.nbCol_; ++j)
    {
      cluster(j,i) = (m_ClusterDataij_(i,j) == 0) ? (unsigned char)(0)
                                                  : (unsigned char)(255);
      data(j,i) = (m_Dataij_(i,j) == 0) ? (unsigned char)(0)
                                        : (unsigned char)(255);
    }
  }

  //Display data and cluster
  CImgDisplay data_disp(data,"Original Data"), cluster_disp(cluster,"Co-Clusters");
  while (!data_disp.is_closed() && !cluster_disp.is_closed()) {
    data_disp.wait();
    cluster_disp.wait();
  }

  data.save("data.jpg");
  cluster.save("cluster.jpg");
}
#endif

void BinaryLBModel::saveThetaInit()
{  m_Alphakltemp_ = m_Alphakl_;}

void BinaryLBModel::modifyTheta()
{
  m_Alphakltemp_   = m_Alphakl_;

  v_logPiektemp_   = v_logPiek_;
  v_logRholtemp_   = v_logRhol_;
  m_Rjltemp_       = m_Rjl_;
  m_Tiktemp_       = m_Tik_;

  Lmax_           = likelihood_;
}

void BinaryLBModel::copyTheta()
{
  m_Alphakl_  = m_Alphakltemp_;
  m_Alphakl1_ = m_Alphakl_;

  v_logPiek_  = v_logPiektemp_;
  v_logRhol_  = v_logRholtemp_;

  m_Tik_      = m_Tiktemp_;
  m_Rjl_      = m_Rjltemp_;
  commonFinalizeOutput();

  likelihood_ = computeLnLikelihood();
}



void BinaryLBModel::mStepFull()
{
  mSteplogRhol();
  mSteplogPiek();
  m_Alphakl_ = (m_Tik_.transpose()*m_Dataij_.cast<STK::Real>()*m_Rjl_)/(v_Tk_*v_Rl_.transpose());
}

/* compute logRhol during the m-step */
void BinaryLBModel::mSteplogRhol()
{
  if(!Mparam_.fixedproportions_)
  { v_logRhol_=((v_Rl_+a_-1)/(Mparam_.nbCol_+Mparam_.nbcolclust_*(a_-1))).log();}
}
/* compute logPiek during the m-step */
void BinaryLBModel::mSteplogPiek()
{
  if(!Mparam_.fixedproportions_)
  { v_logPiek_=((v_Tk_+a_-1)/ (Mparam_.nbRow_+Mparam_.nbrowclust_*(a_-1))).log();}
}

void BinaryLBModel::mStepRows()
{
#ifdef COVERBOSE
  std::cout << "Entering BinaryLBModel::mStepRows" << std::endl;
  std::cout << Error_msg_;
  std::cout << "v_Tk_= " << v_Tk_.transpose();
  std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif

  mSteplogPiek();
  m_Alphakl_ = (((m_Tik_.transpose())*m_Uil_)+b_-1)/((v_Tk_*(v_Rl_.transpose()))+2*(b_-1));
  m_Alphakl_ = m_Alphakl_.max(0.).min(1.);
}

void BinaryLBModel::mStepCols()
{
  mSteplogRhol();
  m_Alphakl_ = ((m_Vjk_.transpose()*m_Rjl_)+b_-1)/((v_Tk_*v_Rl_.transpose())+2*(b_-1));
  m_Alphakl_ = m_Alphakl_.max(0.).min(1.);
}

void BinaryLBModel::mGibbsStepRows()
{
  v_logPiek_=(v_Tk_+a_);
  m_Alphakl_ = (((m_Tik_.transpose())*m_Uil_)+b_);
  //generate random numbers
  VectorReal v_randgamma(Mparam_.nbrowclust_);
  STK::Real sumRng = 0.0;
  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    v_randgamma[k] = STK::Law::Gamma::rand(v_logPiek_[k],1);
    sumRng += v_randgamma[k];
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  { v_logPiek_[k] = v_randgamma[k]/sumRng;}
  v_logPiek_ = (v_logPiek_+RealMin).log();

  MatrixReal m_randgamma(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  VectorReal v_sumRng(Mparam_.nbrowclust_);
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
      m_randgamma(k,l) = STK::Law::Gamma::rand(m_Alphakl_(k,l),1);
      v_sumRng[k] += m_randgamma(k,l);
    }
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
      m_Alphakl_(k,l) = m_randgamma(k,l)/v_sumRng[k];
    }
  }

}

void BinaryLBModel::mGibbsStepCols()
{
  v_logRhol_=(v_Rl_+a_);
  m_Alphakl_ = ((m_Vjk_.transpose()*m_Rjl_)+b_);

  //generate random numbers
  VectorReal v_randgamma(Mparam_.nbcolclust_);
  STK::Real sumRng = 0.0;
  for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
    v_randgamma[l] = STK::Law::Gamma::rand(v_logRhol_[l],1);
    sumRng += v_randgamma[l];
  }

  for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
    v_logRhol_[l] = v_randgamma[l]/sumRng;
  }
  v_logRhol_ = (v_logRhol_+RealMin).log();

  MatrixReal m_randgamma(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  VectorReal v_sumRng(Mparam_.nbrowclust_);
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
      m_randgamma(k,l) = STK::Law::Gamma::rand(m_Alphakl_(k,l),1);
      v_sumRng[k] += m_randgamma(k,l);
    }
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
      m_Alphakl_(k,l) = m_randgamma(k,l)/v_sumRng[k];
    }
  }
}

STK::Real BinaryLBModel::iclCriteriaValue()
{
  STK::Real criteria = 0.0;

  criteria+= lgamma(Mparam_.nbrowclust_*a_)+lgamma(Mparam_.nbcolclust_*a_)
      -(Mparam_.nbrowclust_+Mparam_.nbcolclust_)*lgamma(a_)
      +Mparam_.nbrowclust_*Mparam_.nbcolclust_*(lgamma(2*b_)-2*lgamma(b_))
      -lgamma(Mparam_.nbRow_+Mparam_.nbrowclust_*a_)
      -lgamma(Mparam_.nbCol_+Mparam_.nbcolclust_*a_);

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    criteria+= lgamma(a_+ (v_Zi_== k).count());
  }

  for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
    criteria+= lgamma(a_+ (v_Wj_==l).count());
  }

  STK::ArrayXXi temp0(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  STK::ArrayXXi temp1(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  MatrixBinary m_tempdata = (m_Dataij_==0);
  temp0 = (m_Zik_.transpose()*m_tempdata.cast<int>()*m_Wjl_)+b_;
  temp1 = (m_Zik_.transpose()*m_Dataij_.cast<int>()*m_Wjl_)+b_;
  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
    { criteria+=lgamma(temp0(k,l))+lgamma(temp1(k,l));}
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
    { criteria-= lgamma(((v_Zi_== k).count())*((v_Wj_==l).count())+2*b_);}
  }

  return criteria;
}
