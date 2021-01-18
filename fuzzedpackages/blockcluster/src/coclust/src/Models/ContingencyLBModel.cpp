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

/** @file ContingencyLBModel.cpp
 *  @brief Implements concrete model class ContingencyLBModel_mu_i_nu_j for Contingency Data.
 **/

#include "ContingencyLBModel.h"

ContingencyLBModel::ContingencyLBModel( MatrixReal const& m_Dataij
                                      , ModelParameters const& Mparam
                                      )
                                      : ICoClustModel(Mparam)
                                      , m_Dataij_(m_Dataij)
                                      , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                                      , DataSum_(m_Dataij_.sum())
                                      , m_Gammakl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , m_Gammaklold_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , m_Gammakl1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , m_Gammakl1old_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , m_Gammakltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , v_Ui_(Mparam_.nbRow_)
                                      , v_Vj_(Mparam_.nbCol_)
                                      , m_Ykl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
{
#ifdef COVERBOSE_CONTINGENCY
  std::cout << "ContingencyLBModel created with parameters:"<<std::endl;
  std::cout << Mparam_;
#endif
};

ContingencyLBModel::ContingencyLBModel( MatrixReal const& m_Dataij
                                      , VectorInt const& rowlabels
                                      , VectorInt const& collabels
                                      , ModelParameters const& Mparam
                                      )
                                      : ICoClustModel(Mparam,rowlabels,collabels)
                                      , m_Dataij_(m_Dataij)
                                      , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                                      , DataSum_(m_Dataij_.sum())
                                      , m_Gammakl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , m_Gammaklold_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , m_Gammakl1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , m_Gammakl1old_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , m_Gammakltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                      , v_Ui_(Mparam_.nbRow_)
                                      , v_Vj_(Mparam_.nbCol_)
                                      , m_Ykl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
{
#ifdef COVERBOSE_CONTINGENCY
  std::cout << "ContingencyLBModel created with parameters:"<<std::endl;
  std::cout << Mparam;
#endif
};

void ContingencyLBModel::logSumRows(MatrixReal & m_ik)
{
  m_ik = STK::Const::VectorX (Mparam_.nbRow_)*(v_logPiek_ - m_Gammakl_ * v_Rl_).transpose()
       + m_Uil_* m_Gammakl_.log().transpose();
}

void ContingencyLBModel::logSumCols(MatrixReal & m_jl)
{
  m_jl = STK::Const::VectorX(Mparam_.nbCol_)*(v_logRhol_.transpose() - v_Tk_.transpose() * m_Gammakl_)
       + m_Vjk_* m_Gammakl_.log();
}


bool ContingencyLBModel::emRows()
{
  //Initialization
  computeUil();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepRows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepRows();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool ContingencyLBModel::semRows()
{
  //Initialization
  computeUil();
  if(!seStepRows()) return false;
  //M-step
  mStepRows();
  return true;
}

bool ContingencyLBModel::cemRows()
{
  //Initialization
  computeUil();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!ceStepRows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepRows();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool ContingencyLBModel::emCols()
{
  //Initializations
  computeVjk();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepCols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepCols();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

bool ContingencyLBModel::semCols()
{
  //Initializations
  computeVjk();
  if(!seStepCols()) return false;
  //M-step
  mStepCols();
  return true;
}

bool ContingencyLBModel::cemCols()
{
  //Initializations
  computeVjk();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!ceStepCols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepCols();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

bool ContingencyLBModel::GibbsRows()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE_CONTINGENCY
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool ContingencyLBModel::GibbsCols()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE_CONTINGENCY
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

STK::Real ContingencyLBModel::computeLnLikelihood()
{
  likelihood_ = ( (m_Tik_.transpose()*m_Dataij_*m_Rjl_).prod((m_Gammakl_).log())
                 - m_Gammakl_.prod(v_Tk_*v_Rl_.transpose())
                ).sum()
              //- DataSum_
              + v_Tk_.dot(v_logPiek_) + v_Rl_.dot(v_logRhol_)
              - (m_Tik_.prod( (RealMin + m_Tik_).log()) ).sum()
              - (m_Rjl_.prod( (RealMin + m_Rjl_).log()) ).sum();
  return likelihood_;
}
/* @return the number of free parameters of the distribution of a block.*/
int ContingencyLBModel::nbFreeParameters() const
{
  return Mparam_.nbcolclust_ * Mparam_.nbrowclust_;
}

bool ContingencyLBModel::initStopCriteria()
{ return((((m_Gammakl_-m_Gammakltemp_)/m_Gammakl_).abs().sum())<Mparam_.initepsilon_);}

void ContingencyLBModel::parameterStopCriteria()
{ stopAlgo_ = ((((m_Gammakl1_-m_Gammakl1old_)/m_Gammakl1_).abs().sum())<Mparam_.epsilon_);}


void ContingencyLBModel::consoleOut()
{
#ifdef COVERBOSE_CONTINGENCY
  std::cout<< "Output Model parameter:\ngammakl:\n" << m_Gammakl_
           << "v_Tk_: "      << v_Tk_.transpose()
           << "v_Rl_: "      << v_Rl_.transpose()
           << "v_logPiek_: " << v_logPiek_.transpose()
           << "v_logRhol_: " << v_logRhol_.transpose();
#endif
}

MatrixReal const& ContingencyLBModel::arrangedDataClusters()
{
  arrangedDataCluster(m_ClusterDataij_, m_Dataij_);
  return m_ClusterDataij_;
}

void ContingencyLBModel::saveThetaInit()
{
  m_Gammakltemp_ = m_Gammakl_;
}
void ContingencyLBModel::modifyTheta()
{
  m_Gammakltemp_ = m_Gammakl_;

  v_logPiektemp_ = v_logPiek_;
  v_logRholtemp_ = v_logRhol_;
  m_Rjltemp_     = m_Rjl_;
  m_Tiktemp_     = m_Tik_;

  Lmax_ = likelihood_;
}

void ContingencyLBModel::copyTheta()
{
  m_Gammakl_  = m_Gammakltemp_;
  m_Gammakl1_ = m_Gammakl_;

  v_logPiek_ = v_logPiektemp_;
  v_logRhol_ = v_logRholtemp_;
  m_Rjl_     = m_Rjltemp_;
  m_Tik_     = m_Tiktemp_;

  commonFinalizeOutput();
  likelihood_ = computeLnLikelihood();
}

void ContingencyLBModel::mStepFull()
{
  if(!Mparam_.fixedproportions_)
  {
    v_logRhol_=(v_Rl_/Mparam_.nbCol_).log();
    v_logPiek_=(v_Tk_/Mparam_.nbRow_).log();
  }
  // try some optimization
  if (m_Tik_.sizeCols() < m_Rjl_.sizeCols())
  { m_Ykl_     = (m_Tik_.transpose()*m_Dataij_)*m_Rjl_;}
  else
  { m_Ykl_     = m_Tik_.transpose()*(m_Dataij_*m_Rjl_);}
  m_Gammakl_ = m_Ykl_/(v_Tk_* v_Rl_.transpose());
  //m_Gammakl_ = m_Ykl_/( STK::sumByRow(m_Ykl_)* STK::sum(m_Dataij_*m_Rjl_));
}
