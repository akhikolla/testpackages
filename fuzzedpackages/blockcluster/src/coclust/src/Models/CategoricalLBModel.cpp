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

#include <math.h>
#include "CategoricalLBModel.h"

CategoricalLBModel::CategoricalLBModel( MatrixInt const& m_Dataij
                                      , ModelParameters const& Mparam
                                      , STK::Real a, STK::Real b
                                      )
                                      : ICoClustModel(Mparam)
                                      , a_(a), b_(b)
                                      , m_Dataij_(m_Dataij)
                                      , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                                      , v_Ui_(Mparam_.nbRow_)
                                      , v_Vj_(Mparam_.nbCol_)
                                      , r_(0)
{
  initializeStorages();
}

CategoricalLBModel::CategoricalLBModel( MatrixInt const& m_Dataij
                                      , VectorInt const & rowlabels
                                      , VectorInt const & collabels
                                      , ModelParameters const& Mparam
                                      , STK::Real a, STK::Real b
                                      )
                                      : ICoClustModel(Mparam,rowlabels,collabels)
                                      , a_(a), b_(b)
                                      , m_Dataij_(m_Dataij)
                                      , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                                      , v_Ui_(Mparam_.nbRow_)
                                      , v_Vj_(Mparam_.nbCol_)
                                      , r_(0)
{ initializeStorages();}

CategoricalLBModel::~CategoricalLBModel() {}

void CategoricalLBModel::logSumRows(MatrixReal & m_sum)
{
  m_sum = STK::Const::VectorX (Mparam_.nbRow_)*v_logPiek_.transpose();
  for (int h = 0; h < r_; ++h) {
    m_sum +=  (m3_Yhij_[h].cast<STK::Real>()*m_Rjl_*m3_logAlphahkl_[h].transpose());
  }
}
void CategoricalLBModel::logSumCols(MatrixReal & m_sum)
{
  m_sum = STK::Const::VectorX(Mparam_.nbCol_)*v_logRhol_.transpose();
  for (int h = 0; h < r_; ++h) {
    m_sum +=  (m3_Yhij_[h].transpose().cast<STK::Real>()*m_Tik_*m3_logAlphahkl_[h]);
  }
}

/* compute logRhol during the m-step */
void CategoricalLBModel::mSteplogRhol()
{
  if(!Mparam_.fixedproportions_)
  { v_logRhol_=((v_Rl_+a_-1)/(Mparam_.nbCol_+Mparam_.nbcolclust_*(a_-1))).log();}
}
/* compute logPiek during the m-step */
void CategoricalLBModel::mSteplogPiek()
{
  if(!Mparam_.fixedproportions_)
  { v_logPiek_=((v_Tk_+a_-1)/ (Mparam_.nbRow_+Mparam_.nbrowclust_*(a_-1))).log();}
}


void CategoricalLBModel::mStepRows()
{
  mSteplogPiek();
  Array2DReal m_TkbyRl = (v_Tk_*v_Rl_.transpose())+r_*(b_-1);
  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_-1)/(m_TkbyRl+RealMin);
    m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
}

void CategoricalLBModel::mStepCols()
{
  mSteplogRhol();

  Array2DReal m_TbyRkl = (v_Tk_*v_Rl_.transpose())+r_*(b_-1);
  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_-1)/(m_TbyRkl+RealMin);
    m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
}

void CategoricalLBModel::mGibbsStepRows()
{
  v_logPiek_=(v_Tk_+a_);

  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_);
    m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }

  //generate random numbers
  VectorReal v_randgamma(Mparam_.nbrowclust_);
  STK::Real sumRng = 0.0;
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    v_randgamma[k] = STK::Law::Gamma::rand(v_logPiek_[k],1);
    sumRng += v_randgamma[k];
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    v_logPiek_[k] = v_randgamma[k]/sumRng;
  }
  v_logPiek_ = (v_logPiek_+RealMin).log();

  std::vector<MatrixReal> m_randgamma;
  std::vector<VectorReal> v_sumRng(Mparam_.nbrowclust_);
  m_randgamma.resize(r_);
  v_sumRng.resize(r_);
  for (int h = 0; h < r_; ++h) {
    m_randgamma[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    v_sumRng[h].resize(Mparam_.nbrowclust_);
  }
  for (int h = 0; h < r_; ++h) {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        m_randgamma[h](k,l) = STK::Law::Gamma::rand(m3_Alphahkl_[h](k,l),1);
        v_sumRng[h][k] += m_randgamma[h](k,l);
      }
    }
  }

  for (int h = 0; h < r_; ++h) {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        m3_Alphahkl_[h](k,l) = m_randgamma[h](k,l)/v_sumRng[h][k];
      }
    }
    m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }

}

void CategoricalLBModel::mGibbsStepCols()
{
  v_logRhol_=(v_Rl_+a_);

  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_);
    m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }

  //generate random numbers
  VectorReal v_randgamma(Mparam_.nbrowclust_);
  STK::Real sumRng = 0.0;
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    v_randgamma[k] = STK::Law::Gamma::rand(v_logRhol_[k],1);
    sumRng += v_randgamma[k];
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    v_logRhol_[k] = v_randgamma[k]/sumRng;
  }
  v_logRhol_ = (v_logRhol_+RealMin).log();

  std::vector<MatrixReal> m_randgamma;
  std::vector<VectorReal> v_sumRng(Mparam_.nbcolclust_);
  m_randgamma.resize(r_);
  v_sumRng.resize(r_);
  for (int h = 0; h < r_; ++h) {
    m_randgamma[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    v_sumRng[h].resize(Mparam_.nbrowclust_);
  }

  for (int h = 0; h < r_; ++h) {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        m_randgamma[h](k,l) = STK::Law::Gamma::rand(m3_Alphahkl_[h](k,l),1);
        v_sumRng[h][k] += m_randgamma[h](k,l);
      }
    }
  }

  for (int h = 0; h < r_; ++h) {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        m3_Alphahkl_[h](k,l) = m_randgamma[h](k,l)/v_sumRng[h][k];
      }
    }
    m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
}

bool CategoricalLBModel::emRows()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!eStepRows()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    mStepRows();
    STK::Real netchange = 0.0;
    for (int h = 0; h < r_; ++h)
    {
      netchange+= ((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum();
    }
    netchange/=r_;
    //Termination check
    if (netchange<Mparam_.epsilon_int_) break;
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::cemRows()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!ceStepRows()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    mStepRows();
    STK::Real netchange = 0.0;
    for (int h = 0; h < r_; ++h)
    { netchange+= ((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum();}
    netchange/=r_;
    //Termination check
    if (netchange<Mparam_.epsilon_int_) break;
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::semRows()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}
  if(!seStepRows()) return false;
  mStepRows();
  return true;
}

bool CategoricalLBModel::emCols()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!eStepCols()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    mStepCols();
    STK::Real netchange = 0.0;
    for (int h = 0; h < r_; ++h)
    { netchange+= ((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum();}
    netchange/=r_;
    //Termination check
    if(netchange<Mparam_.epsilon_int_) break;
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::cemCols(){
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!ceStepCols()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    mStepCols();
    STK::Real netchange = 0.0;
    for (int h = 0; h < r_; ++h)
    { netchange+= ((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum();}
    netchange/=r_;
    //Termination check
    if(netchange<Mparam_.epsilon_int_) break;
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_    = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::semCols()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}
  if(!seStepCols()) return false;
  mStepCols();
  return true;
}

bool CategoricalLBModel::GibbsRows()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}
  if(!seStepRows()) return false;
  mGibbsStepRows();
  return true;
}

bool CategoricalLBModel::GibbsCols()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}
  if(!seStepCols()) return false;
  mGibbsStepCols();
  return true;
}

void CategoricalLBModel::saveThetaInit()
{
  m3_Alphahkltemp_ = m3_Alphahkl_;
}

void CategoricalLBModel::modifyTheta()
{
#ifdef COVERBOSE
  std::cout<<"Entering CategoricalLBModel::modifyTheta().\n";
#endif
  m3_Alphahkltemp_ = m3_Alphahkl_;

  v_logPiektemp_ = v_logPiek_;
  v_logRholtemp_ = v_logRhol_;
  m_Rjltemp_ = m_Rjl_;
  m_Tiktemp_ = m_Tik_;

  Lmax_ = likelihood_;
}

void CategoricalLBModel::copyTheta()
{
#ifdef COVERBOSE
  std::cout<<"Entering CategoricalLBModel::copyTheta().\n";
#endif
  m3_Alphahkl_ = m3_Alphahkltemp_;

  v_logPiek_ = v_logPiektemp_;
  v_logRhol_ = v_logRholtemp_;

  m_Tik_ = m_Tiktemp_;
  m_Rjl_ = m_Rjltemp_;

  commonFinalizeOutput();
  likelihood_ = computeLnLikelihood();
}


STK::Real CategoricalLBModel::computeLnLikelihood()
{
  Array2DReal m_Ukl = v_Tk_*v_Rl_.transpose();
  STK::Real tempsum = -(m_Ukl.prod(m_Ukl+RealMin).log()).sum();
  // Compute \sum_h \sum_{i,j,k,l} t_{ik} Y_{ij}^h r_{jl} \log(\alpha_{k,l}^h
  for (int h = 0; h < r_; ++h)
  {
    m_Ukl    = m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>()*m_Rjl_;
    // BUGFIXE (Serge Iovleff) Seems log is wrong
    //    tempsum += ( m_Ukl.prod( (m_Ukl+RealMin).log() )).sum()
    //    		 + (b_-1)*((m3_Alphahkl_[h]+RealMin).log()).sum();
    tempsum += ( m_Ukl.prod( (m3_Alphahkl_[h]+RealMin).log() )).sum();
//    		 + (b_-1)*((m3_Alphahkl_[h]+RealMin).log()).sum();
  }
  likelihood_ = tempsum
              + v_Tk_.dot(v_logPiek_) // -Mparam_.nbRow_*log(STK::Real (Mparam_.nbRow_))
              + v_Rl_.dot(v_logRhol_) // - Mparam_.nbCol_   *log(STK::Real(Mparam_.nbCol_))
              - ( m_Tik_.prod((RealMin + m_Tik_).log()) ).sum()
              - ( m_Rjl_.prod((RealMin + m_Rjl_).log()) ).sum()
			  //              + (a_-1)*(v_logPiek_.sum()+v_logRhol_.sum())
              ;
  return likelihood_;
}

/* @return the number of free parameters of the distribution of a block.*/
int CategoricalLBModel::nbFreeParameters() const
{ return Mparam_.nbcolclust_ * Mparam_.nbrowclust_ * (r_-1);}

bool CategoricalLBModel::initStopCriteria()
{
  STK::Real netchange = 0.0;
  for (int h = 0; h < r_; ++h)
  { netchange+= ((m3_Alphahkl_[h]-m3_Alphahkltemp_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum();}
  netchange/=r_;
  return (netchange<Mparam_.initepsilon_);
}

void CategoricalLBModel::parameterStopCriteria()
{
  STK::Real netchange = 0.0;
  for (int h = 0; h < r_; ++h)
  { netchange+= ((m3_Alphahkl1_[h]-m3_Alphahkl1old_[h]).abs()/(m3_Alphahkl1_[h]+RealMin)).sum();}
  netchange/=r_;
  stopAlgo_ = (netchange<Mparam_.epsilon_);
}


STK::Real CategoricalLBModel::iclCriteriaValue()
{
  STK::Real criteria = lgamma(Mparam_.nbrowclust_*a_)+lgamma(Mparam_.nbcolclust_*a_)
                     - (Mparam_.nbrowclust_+Mparam_.nbcolclust_)*lgamma(a_)
                     + Mparam_.nbrowclust_*Mparam_.nbcolclust_*(lgamma(r_*b_)-r_*lgamma(b_))
                     - lgamma(Mparam_.nbRow_+Mparam_.nbrowclust_*a_)
                     - lgamma(Mparam_.nbCol_+Mparam_.nbcolclust_*a_);

  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  { criteria += lgamma(a_+ (v_Zi_== k).count());}
  for (int l = 0; l < Mparam_.nbcolclust_; ++l)
  { criteria += lgamma(a_+ (v_Wj_==l).count());}

  STK::ArrayXXi temp(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  for (int h = 0; h < r_; ++h)
  {
    temp = ((m_Zik_.transpose()*m3_Yhij_[h].cast<int>())*m_Wjl_)+b_;
    for (int k = 0; k < Mparam_.nbrowclust_; ++k)
    {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l)
      { criteria += lgamma(temp(k,l));}
    }
  }
  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
    { criteria -= lgamma(((v_Zi_== k).count())*((v_Wj_==l).count())+r_*b_);}
  }
  return criteria;
}

void CategoricalLBModel::consoleOut()
{
#ifdef COVERBOSE
  std::cout<<"Output Model parameters\n";
  std::cout<<"\npie_k:"<<v_Piek_<<"\nrho_l:"<<v_Rhol_<<"\n";
  for (int h = 0; h < r_; ++h)
  {
    std::cout<<"Alpha_kl for category "<<h<<"\n";
    std::cout<<m3_Alphahkl_[h]<<"\n";
  }
  std::cout<<"likelihood value:"<<likelihood()<<"\n";
  std::cout<<"ICL value:"<<iclCriteriaValue()<<"\n";
  std::cout << "v_Tk_= " << v_Tk_.transpose();
  std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
}

const MatrixInt& CategoricalLBModel::arrangedDataClusters()
{
  arrangedDataCluster(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

void CategoricalLBModel::initializeStorages()
{
  int maxr =  m_Dataij_.maxElt();
  int minr =  m_Dataij_.minElt();
  r_ = maxr-minr+1;
#ifdef COVERBOSE
  std::cout<<"\nNumber of categories: "<<r_<<"\n";
  std::cout<<"Min category value: "<<minr<<"\n";
  std::cout<<"Max category value: "<<maxr<<"\n";
#endif
  //initialize data storages
  m3_Alphahkl_.resize(r_);
  m3_Alphahklold_.resize(r_);
  m3_Alphahkl1_.resize(r_);
  m3_Alphahkl1old_.resize(r_);
  m3_Alphahkltemp_.resize(r_);
  m3_logAlphahkl_.resize(r_);

  m3_Yhij_.resize(r_);
  m3_Yijh_.resize(Mparam_.nbRow_);
  m3_Yjih_.resize(Mparam_.nbCol_);

  for (int i = 0; i < Mparam_.nbRow_; ++i)
  { m3_Yijh_[i].resize(Mparam_.nbCol_,r_);}

  for (int j = 0; j < Mparam_.nbCol_; ++j)
  { m3_Yjih_[j].resize(Mparam_.nbRow_,r_);}

  for (int h = 0; h < r_; ++h)
  {
    m3_Yhij_[h] = (m_Dataij_ == minr+h);
    for (int i = 0; i < Mparam_.nbRow_; ++i)
    {
      for (int j = 0; j < Mparam_.nbCol_; ++j)
      {
        m3_Yijh_[i](j,h) = m3_Yhij_[h](i,j);
        m3_Yjih_[j](i,h) = m3_Yhij_[h](i,j);
      }
    }
    m3_Alphahkl_[h].resize(Mparam_.nbrowclust_, Mparam_.nbcolclust_);
    m3_Alphahklold_[h].resize(Mparam_.nbrowclust_, Mparam_.nbcolclust_);
    m3_Alphahkl1_[h].resize(Mparam_.nbrowclust_, Mparam_.nbcolclust_);
    m3_Alphahkl1old_[h].resize(Mparam_.nbrowclust_, Mparam_.nbcolclust_);
    m3_Alphahkltemp_[h].resize(Mparam_.nbrowclust_, Mparam_.nbcolclust_);
    m3_logAlphahkl_[h].resize(Mparam_.nbrowclust_, Mparam_.nbcolclust_);
  }
}

