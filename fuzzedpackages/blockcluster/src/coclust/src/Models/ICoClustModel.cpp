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


/** @file ICoClustModel.cpp
 *  @brief This file only initializes the static members of ICoClustModel.
 **/
#define MINSIZE 1e-12

#include "ICoClustModel.h"

ICoClustModel::ICoClustModel( ModelParameters const& Mparam)
                            : Mparam_(Mparam)
                            , likelihood_(-STK::Arithmetic<STK::Real>::infinity())
                            , Lmax_(-STK::Arithmetic<STK::Real>::infinity())
                            , empty_cluster_(false)
                            , m_Tik_(Mparam_.nbRow_, Mparam_.nbrowclust_, 1./Mparam_.nbrowclust_)
                            , m_Rjl_(Mparam_.nbCol_, Mparam_.nbcolclust_, 1./Mparam_.nbcolclust_)
                            , m_Tiktemp_(m_Tik_)
                            , m_Rjltemp_(m_Rjl_)
                            , v_Tk_(Mparam_.nbrowclust_)
                            , v_Rl_(Mparam_.nbcolclust_)
                            , v_Piek_(Mparam_.nbrowclust_, 1./Mparam_.nbrowclust_)
                            , v_Rhol_(Mparam_.nbcolclust_, 1./Mparam_.nbcolclust_)
                            , v_logPiek_(v_Piek_.log())
                            , v_logRhol_(v_Rhol_.log())
                            , v_logPiektemp_(v_logPiek_)
                            , v_logRholtemp_(v_logRhol_)
                            , m_Zik_(Mparam_.nbRow_, Mparam_.nbrowclust_)
                            , m_Wjl_(Mparam_.nbCol_, Mparam_.nbcolclust_)
                            , v_Zi_(Mparam_.nbRow_)
                            , v_Wj_(Mparam_.nbCol_)
                            , m_Vjk_(Mparam_.nbCol_, Mparam_.nbrowclust_)
                            , m_Uil_(Mparam_.nbRow_, Mparam_.nbcolclust_)
                            , stopAlgo_(false)
							              , dimprod_(Mparam.nbRow_*Mparam.nbCol_)
{
  //
  UnknownLabelsRows_.resize(Mparam.nbRow_);
  UnknownLabelsCols_.resize(Mparam.nbCol_);
  for ( int i = UnknownLabelsRows_.begin(); i < UnknownLabelsRows_.end(); ++i)
  { UnknownLabelsRows_[i] = i;}
  for ( int j = UnknownLabelsCols_.begin(); j < UnknownLabelsCols_.end(); ++j)
  { UnknownLabelsCols_[j] =j;}
}

ICoClustModel::ICoClustModel( ModelParameters const& Mparam
                            , VectorInt const& rowlabels
                            , VectorInt const& collabels
                            )
                            : Mparam_(Mparam)
                            , likelihood_(-STK::Arithmetic<STK::Real>::infinity())
                            , empty_cluster_(false)
                            , m_Tik_(Mparam_.nbRow_, Mparam_.nbrowclust_, 1./Mparam_.nbrowclust_)
                            , m_Rjl_(Mparam_.nbCol_, Mparam_.nbcolclust_, 1./Mparam_.nbcolclust_)
                            , m_Tiktemp_(Mparam_.nbRow_, Mparam_.nbrowclust_)
                            , m_Rjltemp_(Mparam_.nbCol_, Mparam_.nbcolclust_)
                            , v_Tk_(Mparam_.nbrowclust_)
                            , v_Rl_(Mparam_.nbcolclust_)
                            , v_Piek_(Mparam_.nbrowclust_, 1./Mparam_.nbrowclust_)
                            , v_Rhol_(Mparam_.nbcolclust_, 1./Mparam_.nbcolclust_)
                            , v_logPiek_(v_Piek_.log())
                            , v_logRhol_(v_Rhol_.log())
                            , v_logPiektemp_(v_logPiek_)
                            , v_logRholtemp_(v_logRhol_)
                            , m_Zik_(Mparam_.nbRow_, Mparam_.nbrowclust_)
                            , m_Wjl_(Mparam_.nbCol_, Mparam_.nbcolclust_)
                            , v_Zi_(Mparam_.nbRow_)
                            , v_Wj_(Mparam_.nbCol_)
                            , m_Vjk_(Mparam_.nbCol_, Mparam_.nbrowclust_)
                            , m_Uil_(Mparam_.nbRow_, Mparam_.nbcolclust_)
                            , stopAlgo_(false)
							              , dimprod_ (Mparam_.nbRow_*Mparam.nbCol_)
{
  setRowLabels(rowlabels);
  setColLabels(collabels);
}


bool ICoClustModel::cemInitStep()
{
  initializeStep();
  bool flag = true;
  for ( int itr = 1; itr < Mparam_.nbinititerations_; ++itr)
  {
    saveThetaInit();

    computeUil();
    mStepRows();
    if (!ceStepRows()) { flag = false; break;}

    computeVjk();
    mStepCols();
    if (!ceStepCols()) { flag = false; break;}

    if(initStopCriteria()) { break;}
  }
  return flag;

}
bool ICoClustModel::emInitStep()
{
  initializeStep();
  bool flag = true;
  for ( int itr = 1; itr < Mparam_.nbinititerations_; ++itr)
  {
    saveThetaInit();

    computeUil();
    mStepRows();
    if (!eStepRows()) { flag = false; break;}

    computeVjk();
    mStepCols();
    if (!eStepCols()) { flag = false; break;}

    if(initStopCriteria()) { break;}
  }
  return flag;
}

bool ICoClustModel::randomInitStep()
{
  initializeStep();
  bool flag = true;
  for ( int itr = 1; itr < Mparam_.nbinititerations_; ++itr)
  {
    computeUil();
    mStepRows();
    if (!seStepRows()) { flag = false; break;}

    computeVjk();
    mStepCols();
    if (!seStepCols()) { flag = false; break;}
  }
  return flag;
}



bool ICoClustModel::eStepRows()
{
  //E-step, compute sumik = log(\psi)
  MatrixReal  m_sumik (Mparam_.nbRow_, Mparam_.nbrowclust_);
  logSumRows(m_sumik);
  m_Tik_  = (m_sumik-STK::maxByRow(m_sumik)*STK::Const::PointX(Mparam_.nbrowclust_)).exp();
  m_Tik_ /= STK::sumByRow(m_Tik_)*STK::Const::PointX(Mparam_.nbrowclust_);

  // reinitialize known labels
  for ( int i=knownLabelsRows_.begin();i< knownLabelsRows_.end();i++)
  {
    m_Tik_.row(knownLabelsRows_[i].first).setZeros();
    m_Tik_(knownLabelsRows_[i].first, knownLabelsRows_[i].second)=1;
  }
  // check empty class
  if( (empty_cluster_ = finalizeStepRows()) )
  {
    Error_msg_  = "In ICoClustModel::eStepRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
    std::cout << "v_Tk_= " << v_Tk_.transpose();
    std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
    return false;
  }
  return true;
}
bool ICoClustModel::ceStepRows()
{
  MatrixReal  m_sumik (Mparam_.nbRow_,Mparam_.nbrowclust_);
  logSumRows(m_sumik);

  for ( int i =UnknownLabelsRows_.begin(); i< UnknownLabelsRows_.end();i++)
  {
    int maxIndex;
    m_sumik.row(UnknownLabelsRows_[i]).maxElt(maxIndex);
    m_Tik_.row(UnknownLabelsRows_[i]).setZeros();
    m_Tik_(UnknownLabelsRows_[i], maxIndex)=1;
  }

  // check empty class
  if( (empty_cluster_ = finalizeStepRows()) )
  {
    Error_msg_  = "In ICoClustModel::ceStepRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
#endif
    return false;
  }
  return true;
}

bool ICoClustModel::eStepCols()
{
  //Temporary variables
  MatrixReal m_sumjl(Mparam_.nbCol_,Mparam_.nbcolclust_);
  logSumCols(m_sumjl);
  m_Rjl_ = (m_sumjl-STK::maxByRow(m_sumjl)*STK::Const::PointX(Mparam_.nbcolclust_)).exp();
  m_Rjl_ /= (sumByRow(m_Rjl_)*STK::Const::PointX(Mparam_.nbcolclust_));

  // reinitialize known labels
  for ( int j=knownLabelsCols_.begin();j< knownLabelsCols_.end();j++)
  {
    m_Rjl_.row(knownLabelsCols_[j].first).setZeros();
    m_Rjl_(knownLabelsCols_[j].first,knownLabelsCols_[j].second)=1;
  }
  // check empty class
  if( (empty_cluster_ = finalizeStepCols()) )
  {
    Error_msg_  = "In ICoClustModel::eStepCols(). Class size too small while running model.\n";
#ifdef COVERBOSE_CONTINGENCY
    std::cout << Error_msg_;
    std::cout << "v_Tk_= " << v_Tk_.transpose();
    std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
    return false;
  }
  return true;
}

bool ICoClustModel::ceStepCols()
{
  MatrixReal  m_sumjl(Mparam_.nbCol_,Mparam_.nbcolclust_);
  logSumCols(m_sumjl);

  // adjust
  for ( int j=UnknownLabelsCols_.begin();j< UnknownLabelsCols_.end();j++)
  {
    int maxIndex;
    m_sumjl.row(UnknownLabelsCols_[j]).maxElt(maxIndex);
    m_Rjl_.row(UnknownLabelsCols_[j]).setZeros();
    m_Rjl_(UnknownLabelsCols_[j], maxIndex)=1;
  }
  // check empty class
  if( (empty_cluster_ = finalizeStepCols()) )
  {
    Error_msg_  = "In ICoClustModel::ceStepCols(). Class size too small while running model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
#endif
    return false;
  }
  return true;
}

bool ICoClustModel::seStepRows()
{
  //E-step: Calculate conditional row class probabilities
  MatrixReal m_sumik (Mparam_.nbRow_,Mparam_.nbrowclust_);
  logSumRows(m_sumik);

  m_Tik_ = (m_sumik-maxByRow(m_sumik)*STK::Const::PointX(Mparam_.nbrowclust_)).exp();
  m_Tik_ /=(sumByRow(m_Tik_)*STK::Const::PointX(Mparam_.nbrowclust_));

  //
  for ( int i=knownLabelsRows_.begin();i< knownLabelsRows_.end();i++)
  {
    m_Tik_.row(knownLabelsRows_[i].first).setZeros();
    m_Tik_(knownLabelsRows_[i].first,knownLabelsRows_[i].second)=1;
  }
//  //S-step : generate class using m_Tik_ and copy back to m_Tik_
  for ( int i=UnknownLabelsRows_.begin();i< UnknownLabelsRows_.end();i++)
  {
    // simul categorical random variable
    int k = STK::Law::Categorical::rand( m_Tik_.row(UnknownLabelsRows_[i]) );
    m_Tik_.row(UnknownLabelsRows_[i]).setZeros();
    m_Tik_(UnknownLabelsRows_[i], k)=1;
  }

//  rowClassMatrixDraw();
//  m_Tik_ = m_Zik_.cast<STK::Real>();
  // check empty class
  if( (empty_cluster_ = finalizeStepRows()) )
  {
    Error_msg_  = "In ICoClustModel::seStepRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
#endif
    return false;
  }
  return true;
}

bool ICoClustModel::seStepCols()
{
  //E-step: Calculation of conditional column class probabilities
  MatrixReal m_sumjl(Mparam_.nbCol_,Mparam_.nbcolclust_);
  logSumCols(m_sumjl);
  m_Rjl_ =(m_sumjl-maxByRow(m_sumjl)*STK::Const::PointX(Mparam_.nbcolclust_)).exp();
  m_Rjl_ /=(sumByRow(m_Rjl_)*STK::Const::PointX(Mparam_.nbcolclust_));
  //
  for ( int j=knownLabelsCols_.begin();j< knownLabelsCols_.end();j++)
  {
    m_Rjl_.row(knownLabelsCols_[j].first).setZeros();
    m_Rjl_(knownLabelsCols_[j].first, knownLabelsCols_[j].second)=1;
  }
//S-step: Draw column class using m_Rjl_ and copy back to m_Rjl_
  for ( int i=UnknownLabelsCols_.begin();i< UnknownLabelsCols_.end();i++)
  {
    // simul categorical random variable
    int l = STK::Law::Categorical::rand( m_Rjl_.row(UnknownLabelsCols_[i]));
    m_Rjl_.row(UnknownLabelsCols_[i]).setZeros();
    m_Rjl_(UnknownLabelsCols_[i], l)=1;
  }

//  colClassMatrixDraw();
//  m_Rjl_ = m_Wjl_.cast<STK::Real>();
  // check empty class
  if( (empty_cluster_ = finalizeStepCols()) )
  {
    Error_msg_  = "In ICoClustModel::seStepCols(). Class size too small while running model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
#endif
    return false;
  }
  return true;
}

void ICoClustModel::initializeStep()
{
  // Calculate classification vector for rows
  m_Tik_.randUnif();
  m_Tik_ /= STK::sumByRow(m_Tik_) * STK::Const::PointX(m_Tik_.cols());
  v_Tk_ = STK::sumByCol(m_Tik_);
  m_Zik_.setZeros();
  int maxIndex;
  for ( int i = 0; i <Mparam_.nbRow_; ++i)
  {
    m_Tik_.row(i).maxElt(maxIndex);
    v_Zi_[i] = maxIndex;
    m_Zik_(i,maxIndex) = 1;
  }

  // Calculate classification vector for columns
  m_Rjl_.randUnif();
  m_Rjl_ /= STK::sumByRow(m_Rjl_) * STK::Const::PointX(m_Rjl_.cols());
  v_Rl_ = STK::sumByCol(m_Rjl_);
  m_Wjl_.setZeros();
  for ( int j = 0; j < Mparam_.nbCol_; ++j)
  {
    m_Rjl_.row(j).maxElt(maxIndex);
    v_Wj_[j] = maxIndex;
    m_Wjl_(j,maxIndex) = 1;
  }

  // Row initialization
  mSteplogPiek();
  v_Piek_ = v_logPiek_.exp();
  computeUil();
  mStepRows();
  eStepRows();
  // col intialization
  mSteplogRhol();
  computeVjk();
  mStepCols();
  eStepCols();
#ifdef COVERBOSE_CONTINGENCY
  std::cout<<"ICoClustModel::initializeStep. done\n";
  std::cout<<"v_Tk_=\n" << v_Tk_.transpose() << "\n";
  std::cout<<"v_Rl_=\n" << v_Rl_.transpose() << "\n";
  std::cout<<"v_Piek_=\n" << v_Piek_.transpose() << "\n";
  std::cout<<"v_Rhol_=\n" << v_Rhol_.transpose() << "\n";
#endif
}

/* compute the vector v_Tk_ and check if the size block is not too small
 *  @return false if the size block is under the threshold, true otherwise
 **/
bool ICoClustModel::finalizeStepRows()
{
  // compute size of the blocks by column
  v_Tk_ = STK::sumByCol(m_Tik_);
  // check empty class
  return (v_Tk_ * v_Rl_.transpose() < MINSIZE).any();
}
/* compute the vector v_Rl_ and check if the size block is not too small
 *  @return false if the size block is under the threshold, true otherwise
 **/
bool ICoClustModel::finalizeStepCols()
{
  // compute size of the blocks by column
  v_Rl_ = STK::sumByCol(m_Rjl_);
  // check empty class
  return (v_Tk_ * v_Rl_.transpose() < MINSIZE).any();
}

/*
 * Interface for calculating Stopping condition using percentage Change in Likelihood. This function will set the
 * ICoClustModel::stopAlgo_ parameter to either true or false depending on whether the change
 * in Likelihood is less than Mparam_.epsilon_ or not respectively.
 */
void ICoClustModel::likelihoodStopCriteria()
{
  STK::Real l_old = likelihood_;
  computeLnLikelihood();
  stopAlgo_ = (std::abs(likelihood_-l_old)<std::abs(likelihood_)*Mparam_.epsilon_)
            ? true : false;
}

STK::Real ICoClustModel::iclCriteriaValue()
{

  STK::Real criteria = 0.0;
  //  Error_msg_ = "ICL creteria is not yet implemented for this model.";
  criteria = likelihood_
           - std::log(Mparam_.nbRow_)*(Mparam_.nbrowclust_-1.)/2.
           - std::log(Mparam_.nbCol_)*(Mparam_.nbcolclust_-1.)/2.
           - std::log(Mparam_.nbCol_*Mparam_.nbRow_) * nbFreeParameters()/2.;
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return criteria;
}

void ICoClustModel::commonFinalizeOutput()
{
  // Calculate row and column proportions
  if(!Mparam_.fixedproportions_)
  {
    v_Piek_ = v_logPiek_.exp();
    v_Rhol_ = v_logRhol_.exp();
  }
  else
  {
    v_Piek_ = (1.0/Mparam_.nbrowclust_)*STK::Const::VectorX(Mparam_.nbrowclust_);
    v_Rhol_ = (1.0/Mparam_.nbcolclust_)*STK::Const::VectorX(Mparam_.nbcolclust_);
  }
  // compute class size
  v_Tk_ = STK::Stat::sumByCol(m_Tik_);
  v_Rl_ = STK::Stat::sumByCol(m_Rjl_);

  // Calculate classification vector for rows
  m_Zik_.setZeros();
  m_Wjl_.setZeros();
  int maxIndex;
  for ( int i = 0; i <Mparam_.nbRow_; ++i)
  {
    m_Tik_.row(i).maxElt(maxIndex);
    v_Zi_[i] = maxIndex;
    m_Zik_(i,maxIndex) = 1;
  }
  // Calculate classification vector for columns
  for ( int j = 0; j < Mparam_.nbCol_; ++j)
  {
    m_Rjl_.row(j).maxElt(maxIndex);
    v_Wj_[j] = maxIndex;
    m_Wjl_(j,maxIndex) = 1;
  }
  // check empty cluster
  empty_cluster_ = (v_Tk_ * v_Rl_.transpose() < MINSIZE).any();
}

void ICoClustModel::finalizeOutput() {}

void ICoClustModel::setRowLabels(VectorInt const& rowlabels)
{
  for ( int i = rowlabels.begin(); i < rowlabels.end(); ++i)
  {
    int cluster = rowlabels[i];
    if (cluster<0)
    {
      UnknownLabelsRows_.push_back(i);
    }
    else
    {
      knownLabelsRows_.push_back(std::pair<int,int>(i,cluster));
      v_Zi_[i] = cluster;
      m_Zik_(i,cluster) = 1;
    }
  }
}

void ICoClustModel::setColLabels(VectorInt const& collabels)
{
  for ( int j = collabels.begin(); j < collabels.end(); ++j)
  {
    int cluster = collabels[j];
    if (cluster<0) { UnknownLabelsCols_.push_back(j);}
    else
    {
      knownLabelsCols_.push_back(std::pair<int,int>(j,cluster));
      v_Wj_[j] = cluster;
      m_Wjl_(j,cluster) = 1;
    }
  }
}

#undef MINSIZE
