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

/** @file ContingencyLBModel_mu_i_nu_j.cpp
 *  @brief Implements concrete model class ContingencyLBModel_mu_i_nu_j
 *   derived from ICoClustModel.
 **/

#include "ContingencyLBModel_mu_i_nu_j.h"
#include <iostream>

ContingencyLBModel_mu_i_nu_j::ContingencyLBModel_mu_i_nu_j( MatrixReal const& m_Dataij
                                                          , VectorReal const& v_Mui
                                                          , VectorReal const& v_Nuj
                                                          , ModelParameters const& Mparam
                                                          )
                                                          : ICoClustModel(Mparam)
                                                          , m_Dataij_(m_Dataij)
                                                          , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                                                          , v_Mui_(v_Mui),v_Nuj_(v_Nuj)
                                                          , DataSum_(m_Dataij.sum())
                                                          , m_Gammakl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , m_Gammaklold_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , m_Gammakl1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , m_Gammakl1old_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , m_Gammakltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , v_nul_(Mparam_.nbcolclust_)
                                                          , v_muk_(Mparam_.nbrowclust_)
                                                          , v_Ui_(Mparam_.nbRow_)
                                                          , v_Vj_(Mparam_.nbCol_)
                                                          , m_Ykl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
{};

ContingencyLBModel_mu_i_nu_j::ContingencyLBModel_mu_i_nu_j( MatrixReal const& m_Dataij
                                                          , VectorInt const & rowlabels
                                                          , VectorInt const & collabels
                                                          , VectorReal const& v_Mui
                                                          , VectorReal const& v_Nuj
                                                          , ModelParameters const& Mparam
                                                          )
                                                          : ICoClustModel(Mparam,rowlabels,collabels)
                                                          , m_Dataij_(m_Dataij)
                                                          , m_ClusterDataij_(Mparam_.nbRow_, Mparam_.nbCol_)
                                                          , v_Mui_(v_Mui),v_Nuj_(v_Nuj)
                                                          , DataSum_(m_Dataij.sum())
                                                          , m_Gammakl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , m_Gammaklold_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , m_Gammakl1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , m_Gammakl1old_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , m_Gammakltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                          , v_nul_(Mparam_.nbcolclust_)
                                                          , v_muk_(Mparam_.nbrowclust_)
                                                          , v_Ui_(Mparam_.nbRow_)
                                                          , v_Vj_(Mparam_.nbCol_)
                                                          , m_Ykl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
{};

void ContingencyLBModel_mu_i_nu_j::logSumRows(MatrixReal & m_ik)
{
  m_ik = STK::Const::VectorX(Mparam_.nbRow_)*v_logPiek_.transpose()
       + m_Uil_*m_Gammakl_.log().transpose()
       - v_Mui_*(m_Gammakl_* v_nul_).transpose();
}

void ContingencyLBModel_mu_i_nu_j::logSumCols(MatrixReal & m_jl)
{
  m_jl = STK::Const::VectorX(Mparam_.nbCol_)*v_logRhol_.transpose()
       + m_Vjk_*m_Gammakl_.log()
       - v_Nuj_*(v_muk_.transpose()*m_Gammakl_);
}

bool ContingencyLBModel_mu_i_nu_j::emRows()
{
  computeUil();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepRows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepRows();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_).abs()/m_Gammakl_).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::semRows()
{
  computeUil();

  if(!seStepRows()) return false;
  //M-step
  mStepRows();
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::cemRows()
{
  computeUil();

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {

    if(!ceStepRows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepRows();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_).abs()/m_Gammakl_).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::emCols()
{
  //Initialization
  computeVjk();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepCols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepCols();

    //Termination check
    if((((m_Gammakl_-m_Gammaklold_).abs()/m_Gammakl_).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }

  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::semCols()
{
  //Initialization
  computeVjk();

  if(!seStepCols()) return false;
  //M-step
  mStepCols();

  return true;
}

bool ContingencyLBModel_mu_i_nu_j::cemCols()
{
  computeVjk();

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!ceStepCols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepCols();

    //Termination check
    if((((m_Gammakl_-m_Gammaklold_).abs()/m_Gammakl_).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }

  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::GibbsRows()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool ContingencyLBModel_mu_i_nu_j::GibbsCols()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

STK::Real ContingencyLBModel_mu_i_nu_j::computeLnLikelihood()
{
  likelihood_ = (m_Ykl_.prod(m_Gammakl_.log()) ).sum() - DataSum_
		        + v_Tk_.transpose()*v_logPiek_
			      + v_Rl_.transpose()*v_logRhol_
	          - (m_Tik_.prod((RealMin + m_Tik_).log()) ).sum()
  		      - (m_Rjl_.prod((RealMin + m_Rjl_).log()) ).sum();

  return likelihood_;
}

/* @return the number of free parameters of the distribution of a block.*/
int ContingencyLBModel_mu_i_nu_j::nbFreeParameters() const
{return Mparam_.nbcolclust_*Mparam_.nbrowclust_;}

bool ContingencyLBModel_mu_i_nu_j::initStopCriteria()
{
  return(((m_Gammakl_-m_Gammakltemp_).abs()/m_Gammakl_).sum()<Mparam_.initepsilon_);
}
void ContingencyLBModel_mu_i_nu_j::parameterStopCriteria()
{
  stopAlgo_ = (((m_Gammakl1_-m_Gammakl1old_).abs()/m_Gammakl1_).sum()<Mparam_.epsilon_);
}


void ContingencyLBModel_mu_i_nu_j::consoleOut()
{
#ifndef COVERBOSE
  std::cout<<"Output Model parameter:"<<"\ngammakl:\n"<<m_Gammakl_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}

MatrixReal const& ContingencyLBModel_mu_i_nu_j::arrangedDataClusters()
{
  arrangedDataCluster<MatrixReal>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

void ContingencyLBModel_mu_i_nu_j::saveThetaInit()
{
  m_Gammakltemp_ = m_Gammakl_;
}

void ContingencyLBModel_mu_i_nu_j::modifyTheta()
{
  m_Gammakltemp_ = m_Gammakl_;
  v_logPiektemp_ = v_logPiek_;
  v_logRholtemp_ = v_logRhol_;

  m_Rjltemp_ = m_Rjl_;
  m_Tiktemp_ = m_Tik_;

  Lmax_ = likelihood_;
}

void ContingencyLBModel_mu_i_nu_j::copyTheta()
{
  m_Gammakl_  = m_Gammakltemp_;
  m_Gammakl1_ = m_Gammakl_;

  v_logPiek_ = v_logPiektemp_;
  v_logRhol_ = v_logRholtemp_;

  m_Tik_ = m_Tiktemp_;
  m_Rjl_ = m_Rjltemp_;

  commonFinalizeOutput();
  likelihood_ = computeLnLikelihood();
}


void ContingencyLBModel_mu_i_nu_j::mStepFull()
{
  if(!Mparam_.fixedproportions_)
  {
    v_logRhol_=(v_Rl_/Mparam_.nbCol_).log();
    v_logPiek_=(v_Tk_/Mparam_.nbRow_).log();
  }
  m_Ykl_     = m_Tik_.transpose()*m_Dataij_*m_Rjl_;
  m_Gammakl_ = m_Ykl_/((m_Tik_.transpose()*v_Mui_)*(v_Nuj_.transpose()*m_Rjl_));
}
