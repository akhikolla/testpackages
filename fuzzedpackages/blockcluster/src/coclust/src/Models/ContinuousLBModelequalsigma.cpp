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

/** @file ContinuousLBModelequalsigma.cpp
 *  @brief Implements concrete model class ContinuousLBModelequalsigma derived from ICoClustModel.
 **/

#include "ContinuousLBModelequalsigma.h"
#include <exception>
ContinuousLBModelequalsigma::ContinuousLBModelequalsigma( MatrixReal const& m_Dataij
                                                        , ModelParameters const& Mparam
                                                        )
                                                        : ICoClustModel(Mparam)
                                                        , m_Dataij_(m_Dataij)
                                                        , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                                                        , m_Dataij2_(m_Dataij_.square())
                                                        , m_Mukl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                        , Sigma2_(1.)
                                                        , Sigma2temp_(1.)
                                                        , m_Muklold1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                        , m_Muklold2_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                        , m_Mukltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                        , m_Vjk2_(Mparam_.nbCol_, Mparam_.nbrowclust_)
                                                        , m_Uil2_(Mparam_.nbRow_, Mparam_.nbcolclust_)
{}

ContinuousLBModelequalsigma::ContinuousLBModelequalsigma( MatrixReal const& m_Dataij
                                                        , VectorInt const & rowlabels
                                                        , VectorInt const & collabels
                                                        , ModelParameters const& Mparam
                                                        )
                                                        : ICoClustModel(Mparam,rowlabels,collabels)
                                                        , m_Dataij_(m_Dataij)
                                                        , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                                                        , m_Dataij2_(m_Dataij_.square())
                                                        , m_Mukl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                        , Sigma2_(1.)
                                                        , Sigma2temp_(1.)
                                                        , m_Muklold1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                        , m_Muklold2_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                        , m_Mukltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                        , m_Vjk2_(Mparam_.nbCol_, Mparam_.nbrowclust_)
                                                        , m_Uil2_(Mparam_.nbRow_, Mparam_.nbcolclust_)
{}

void ContinuousLBModelequalsigma::logSumRows(MatrixReal & m_sumik)
{
  m_sumik = STK::Const::VectorX (Mparam_.nbRow_)
            * (v_logPiek_- (0.5*(m_Mukl_.square()*v_Rl_))/Sigma2_).transpose()
          + (m_Uil_*m_Mukl_.transpose())/Sigma2_ ;
}

void ContinuousLBModelequalsigma::logSumCols(MatrixReal & m_sumjl)
{
  m_sumjl = STK::Const::VectorX(Mparam_.nbCol_)
            * (v_logRhol_.transpose() - (0.5*v_Tk_.transpose()*m_Mukl_.square())/Sigma2_)
          + (m_Vjk_*m_Mukl_)/Sigma2_ ;
}


bool ContinuousLBModelequalsigma::emRows()
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

bool ContinuousLBModelequalsigma::cemRows()
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
    if ((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}

  }
  return true;
}

bool ContinuousLBModelequalsigma::semRows()
{
  //Initialization
  computeUil();
  if(!seStepRows()) return false;
  //M-step
  mStepRows();
  return true;
}

bool ContinuousLBModelequalsigma::emCols()
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
    if((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool ContinuousLBModelequalsigma::semCols()
{
  //Initialization
  computeVjk();
  if(!seStepCols()) return false;
  //M-step
  mStepCols();
  return true;
}

bool ContinuousLBModelequalsigma::cemCols()
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

bool ContinuousLBModelequalsigma::GibbsRows()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool ContinuousLBModelequalsigma::GibbsCols()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool ContinuousLBModelequalsigma::initStopCriteria()
{
  return((((m_Mukl_-m_Mukltemp_)/m_Mukl_).abs().sum())<Mparam_.initepsilon_);
}

void ContinuousLBModelequalsigma::parameterStopCriteria()
{
  if((((m_Mukl_-m_Muklold1_)/m_Mukl_).abs().sum())<Mparam_.epsilon_)
  { stopAlgo_ = true;}
  else
  {
    stopAlgo_ = false;
    m_Muklold1_ =m_Mukl_;
  }
}

STK::Real ContinuousLBModelequalsigma::computeLnLikelihood()
{
  likelihood_ = -(dimprod_/2)*(1+std::log(Sigma2_))
              + v_Tk_.dot(v_logPiek_)
              + v_Rl_.dot(v_logRhol_)
              - (m_Tik_.prod( (RealMin + m_Tik_).log()) ).sum()
              - (m_Rjl_.prod( (RealMin + m_Rjl_).log()) ).sum();

  return likelihood_;
}

/* @return the number of free parameters of the distribution of a block.*/
int ContinuousLBModelequalsigma::nbFreeParameters() const
{
  return Mparam_.nbcolclust_*Mparam_.nbrowclust_ + 1;
}

void ContinuousLBModelequalsigma::consoleOut()
{
#ifndef COVERBOSE
    std::cout<<"Output Model parameter:"<<"\nBlock Mean:\n"<<m_Mukl_<<"\nBlock Sigma:\n"<<Sigma2_<<"\npiek: "<<
        v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}
MatrixReal const& ContinuousLBModelequalsigma::arrangedDataClusters()
{
  arrangedDataCluster<MatrixReal>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

void ContinuousLBModelequalsigma::saveThetaInit()
{
  m_Mukltemp_ = m_Mukl_;
  Sigma2temp_ = Sigma2_;
}

void ContinuousLBModelequalsigma::modifyTheta()
{
  m_Mukltemp_ = m_Mukl_;
  Sigma2temp_ = Sigma2_;

  v_logPiektemp_ = v_logPiek_;
  v_logRholtemp_ = v_logRhol_;

  m_Rjltemp_ = m_Rjl_;
  m_Tiktemp_ = m_Tik_;

  Lmax_ = likelihood_;
}

void ContinuousLBModelequalsigma::copyTheta()
{
  m_Mukl_ = m_Mukltemp_;
  Sigma2_ = Sigma2temp_;

  v_logPiek_ = v_logPiektemp_;
  v_logRhol_ = v_logRholtemp_;

  m_Rjl_ = m_Rjltemp_;
  m_Tik_ = m_Tiktemp_;

  commonFinalizeOutput();

  commonFinalizeOutput();
  likelihood_ = computeLnLikelihood();
}

void ContinuousLBModelequalsigma::mStepFull()
{
  if(!Mparam_.fixedproportions_)
  {
    v_logRhol_=(v_Rl_/Mparam_.nbCol_).log();
    v_logPiek_=(v_Tk_/Mparam_.nbRow_).log();
  }

  m_Mukl_ = (m_Tik_.transpose()*m_Dataij_.cast<STK::Real>()*m_Rjl_)/(v_Tk_*(v_Rl_.transpose()));
  //m_Mukl2_ = m_Mukl_.square();
  Sigma2_ = ((m_Tik_.transpose()*m_Dataij2_.cast<STK::Real>()*m_Rjl_).sum()-v_Tk_.transpose()*m_Mukl_.square()*v_Rl_)/dimprod_;
}


