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


/** @file ContinuousLBModel.cpp
 *  @brief Implements concrete model class ContinuousLBModel derived from ICoClustModel.
 **/

#include "ContinuousLBModel.h"
#include <exception>
ContinuousLBModel::ContinuousLBModel( MatrixReal const& m_Dataij, ModelParameters const& Mparam)
                                    : ICoClustModel(Mparam)
                                    , m_Dataij_(m_Dataij)
                                    , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                                    , m_Dataij2_(m_Dataij_.square())
                                    , m_Mukl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Sigma2kl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Sigma2kltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Muklold1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Muklold2_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Mukltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Vjk2_(Mparam_.nbCol_, Mparam_.nbrowclust_)
                                    , m_Uil2_(Mparam_.nbRow_, Mparam_.nbcolclust_)
{}

ContinuousLBModel::ContinuousLBModel( MatrixReal const& m_Dataij
		                                , VectorInt const & rowlabels
                                    , VectorInt const & collabels
                                    , ModelParameters const& Mparam)
                                    : ICoClustModel(Mparam,rowlabels,collabels)
                                    , m_Dataij_(m_Dataij)
                                    , m_Dataij2_(m_Dataij_.square())
                                    , m_Mukl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Sigma2kl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Sigma2kltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Muklold1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Muklold2_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Mukltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                    , m_Vjk2_(Mparam_.nbCol_, Mparam_.nbrowclust_)
                                    , m_Uil2_(Mparam_.nbRow_, Mparam_.nbcolclust_)
{}

void ContinuousLBModel::logSumRows(MatrixReal& m_sumik)
{
//  VectorReal tmp1(Mparam_.nbrowclust_);
//  tmp1 = ((0.5*((m_Sigma2kl_.log()+m_Mukl_.square()/m_Sigma2kl_)*v_Rl_)) - v_logPiek_);
  m_sumik = STK::Const::VectorX (Mparam_.nbRow_)
            * ( ( v_logPiek_ - 0.5*( (m_Sigma2kl_.log()+m_Mukl_.square()/m_Sigma2kl_)*v_Rl_))).transpose()
            - 0.5     * (m_Uil2_ * m_Sigma2kl_.inverse().transpose())
            + m_Uil_ * (m_Mukl_/m_Sigma2kl_).transpose();
}

void ContinuousLBModel::logSumCols(MatrixReal & m_sumjl)
{
//  VectorReal tmp1(Mparam_.nbcolclust_);
//  tmp1 = ((0.5*(v_Tk_.transpose()*(m_Sigma2kl_.log()+m_Mukl_.square()/m_Sigma2kl_))) - v_logRhol_);
  m_sumjl = STK::Const::VectorX(Mparam_.nbCol_)
            * ( v_logRhol_.transpose() - 0.5 * v_Tk_.transpose() * (m_Sigma2kl_.log() + m_Mukl_.square()/m_Sigma2kl_))
          - 0.5    * (m_Vjk2_ * m_Sigma2kl_.inverse())
          + m_Vjk_* (m_Mukl_/m_Sigma2kl_);
}


bool ContinuousLBModel::cemRows()
{
  //Initialization
  computeUil();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!ceStepRows()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    mStepRows();
    //Termination check
    if((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool ContinuousLBModel::emRows()
{
  //Initialization
  computeUil();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepRows()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    mStepRows();
    //Termination check
    if((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool ContinuousLBModel::semRows()
{
  //Initialization
  computeUil();
  if(!seStepRows()) return false;
  //M-step
  mStepRows();
  return true;
}
bool ContinuousLBModel::GibbsRows()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}


bool ContinuousLBModel::emCols()
{
  //Initialization
  computeVjk();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepCols()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    mStepCols();
    //Termination check
    if((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContinuousLBModel::cemCols()
{
  //Initialization
  computeVjk();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!ceStepCols()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    mStepCols();
    //Termination check
    if((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool ContinuousLBModel::semCols()
{
  //Initialization
  computeVjk();
  if(!seStepCols()) return false;
  //M-step
  mStepCols();
  return true;
}

bool ContinuousLBModel::GibbsCols()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}


bool ContinuousLBModel::initStopCriteria()
{
  return((((m_Mukl_-m_Mukltemp_)/m_Mukl_).abs().sum())<Mparam_.initepsilon_);
}

void ContinuousLBModel::parameterStopCriteria()
{
  if((((m_Mukl_-m_Muklold1_)/m_Mukl_).abs().sum())<Mparam_.epsilon_)
  {
    stopAlgo_ = true;
  }
  else
  {
    stopAlgo_ = false;
    m_Muklold1_ =m_Mukl_;
  }
}

STK::Real ContinuousLBModel::computeLnLikelihood()
{
  likelihood_ = (-0.5*(dimprod_+v_Tk_.dot(m_Sigma2kl_.log()*v_Rl_))
                 + v_Tk_.dot(v_logPiek_)
                 + v_Rl_.dot(v_logRhol_)
                 - (m_Tik_.prod( (RealMin + m_Tik_).log()) ).sum()
                 - (m_Rjl_.prod( (RealMin + m_Rjl_).log()) ).sum()
			     );
  return likelihood_;
}

/* @return the number of free parameters of the distribution of a block.*/
int ContinuousLBModel::nbFreeParameters() const
{return 2 * Mparam_.nbrowclust_ * Mparam_.nbcolclust_;}


void ContinuousLBModel::consoleOut()
{
#ifndef COVERBOSE
    std::cout<<"Output Model parameter:"<<"\nBlock Mean:\n"<<m_Mukl_<<"\nBlock Sigma:\n"<<m_Sigma2kl_<<"\npiek: "<<
        v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}


MatrixReal const& ContinuousLBModel::arrangedDataClusters()
{
  arrangedDataCluster<MatrixReal>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

void ContinuousLBModel::saveThetaInit()
{
  m_Mukltemp_     = m_Mukl_;
  m_Sigma2kltemp_ = m_Sigma2kl_;
}

void ContinuousLBModel::modifyTheta()
{
  m_Mukltemp_     = m_Mukl_;
  m_Sigma2kltemp_ = m_Sigma2kl_;

  v_logPiektemp_ = v_logPiek_;
  v_logRholtemp_ = v_logRhol_;

  m_Rjltemp_     = m_Rjl_;
  m_Tiktemp_     = m_Tik_;

  Lmax_ = likelihood_;
}

void ContinuousLBModel::copyTheta()
{
  m_Mukl_     = m_Mukltemp_;
  m_Sigma2kl_ = m_Sigma2kltemp_;
  m_Muklold1_ = m_Mukl_;

  v_logPiek_ = v_logPiektemp_;
  v_logRhol_ = v_logRholtemp_;

  m_Rjl_     = m_Rjltemp_;
  m_Tik_     = m_Tiktemp_;

  commonFinalizeOutput();
  likelihood_ = computeLnLikelihood();
}

void ContinuousLBModel::mStepFull()
{
  if(!Mparam_.fixedproportions_)
  {
    v_logRhol_=(v_Rl_/Mparam_.nbCol_).log();
    v_logPiek_=(v_Tk_/Mparam_.nbRow_).log();
  }

  MatrixReal m_trkl = (v_Tk_*v_Rl_.transpose());
  if  (Mparam_.nbRow_<Mparam_.nbCol_)
  { m_Mukl_ = ((m_Tik_.transpose()*m_Dataij_)*m_Rjl_)/m_trkl;}
  else
  { m_Mukl_ = ( m_Tik_.transpose()*(m_Dataij_*m_Rjl_))/m_trkl;}
  //m_Mukl2_    = m_Mukl_.square();
  if  (Mparam_.nbRow_<Mparam_.nbCol_)
  { m_Sigma2kl_ = ( (m_Tik_.transpose()*m_Dataij2_)*m_Rjl_)/m_trkl - m_Mukl_.square();}
  else
  { m_Sigma2kl_ = ( m_Tik_.transpose()*(m_Dataij2_*m_Rjl_))/m_trkl - m_Mukl_.square();}
}

