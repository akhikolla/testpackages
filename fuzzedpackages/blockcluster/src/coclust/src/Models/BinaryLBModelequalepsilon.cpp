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


/** @file BinaryLBModelequalepsilon.cpp
 *  @brief Implements concrete BinaryLBModelequalepsilon model class derived from ICoClustModel.
 **/

#include "BinaryLBModelequalepsilon.h"

#ifndef RPACKAGE
using namespace cimg_library;
#endif
BinaryLBModelequalepsilon::BinaryLBModelequalepsilon( MatrixBinary const&  m_Dataij
                                                    , ModelParameters const& Mparam
                                                    , STK::Real a, STK::Real b
                                                    )
                                                    : ICoClustModel(Mparam)
                                                    , a_(a), b_(b)
                                                    , m_Dataij_(m_Dataij)
                                                    , m_Xjl_(STK::sumByCol(m_Dataij_.cast<STK::Real>()).transpose()*STK::Const::Point<STK::Real>(Mparam_.nbcolclust_))
                                                    , m_Xik_(STK::sumByRow(m_Dataij_.cast<STK::Real>())*STK::Const::Point<STK::Real>(Mparam_.nbrowclust_))
                                                    , m_Tk_Rl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Ukl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , v_Ui_(Mparam_.nbRow_)
                                                    , m_Ykl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Ykl_old2_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Ykl_old1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Akl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Akltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , Epsilon_(0)
                                                    , Epsilontemp_(0)
                                                    , W1_(0),W1_old_(0)
{};

BinaryLBModelequalepsilon::BinaryLBModelequalepsilon( MatrixBinary const&  m_Dataij
                                                    , VectorInt const& rowlabels
                                                    , VectorInt const& collabels
                                                    , ModelParameters const& Mparam
                                                    , STK::Real a, STK::Real b
                                                    )
                                                    : ICoClustModel(Mparam,rowlabels,collabels)
                                                    , a_(a), b_(b)
                                                    , m_Dataij_(m_Dataij)
                                                    , m_Xjl_(STK::sumByCol(m_Dataij_.cast<STK::Real>()).transpose()*STK::Const::Point<STK::Real>(Mparam_.nbcolclust_))
                                                    , m_Xik_(STK::sumByRow(m_Dataij_.cast<STK::Real>())*STK::Const::Point<STK::Real>(Mparam_.nbrowclust_))
                                                    , m_Tk_Rl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Ukl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , v_Ui_(Mparam_.nbRow_)
                                                    , m_Ykl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Ykl_old2_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Ykl_old1_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Akl_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , m_Akltemp_(Mparam_.nbrowclust_, Mparam_.nbcolclust_)
                                                    , Epsilon_(0)
                                                    , Epsilontemp_(0)
                                                    , W1_(0),W1_old_(0)
{};

//Compute Bernoulli log-sum for all rows
void BinaryLBModelequalepsilon::logSumRows(MatrixReal & m_sum)
{
  STK::Real logepsilon = log(Epsilon_/(1.-Epsilon_));
  m_sum = STK::Const::VectorX (Mparam_.nbRow_)*(v_logPiek_+logepsilon*m_Akl_.cast<STK::Real>()*v_Rl_).transpose()
          -logepsilon*(2*m_Uil_*m_Akl_.cast<STK::Real>().transpose() + m_Xik_.cast<STK::Real>());
}

//Compute Bernoulli log-sum for all columns
void BinaryLBModelequalepsilon::logSumCols(MatrixReal & m_sum)
{
  STK::Real logepsilon = log(Epsilon_/(1-Epsilon_));
  m_sum = STK::Const::VectorX(Mparam_.nbCol_)*(v_logRhol_.transpose()+logepsilon*v_Tk_.transpose()*m_Akl_.cast<STK::Real>())
          -logepsilon*(2*m_Vjk_*m_Akl_.cast<STK::Real>() + m_Xjl_.cast<STK::Real>());
}


//Run EM algorithm on data Matrix m_Uil_
bool BinaryLBModelequalepsilon::emRows()
{
  //Initializations
  computeUil();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!eStepRows()) return false;
    //M-step
    m_Ykl_old2_ = m_Ykl_;
    mStepRows();
    if(((m_Ykl_-m_Ykl_old2_)/m_Ykl_).abs().sum()<Mparam_.epsilon_int_)
    {
      break;
    }
  }
  return true;
}

//Run CEM algorithm on data Matrix m_Uil_
bool BinaryLBModelequalepsilon::cemRows()
{
  //Initializations
  computeUil();
  MatrixReal  m_sumik (Mparam_.nbRow_,Mparam_.nbrowclust_);
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!ceStepRows()) return false;
    //M-step
    //M-step
    m_Ykl_old2_ = m_Ykl_;
    mStepRows();
    if(((m_Ykl_-m_Ykl_old2_)/m_Ykl_).abs().sum()<Mparam_.epsilon_int_)
    {
       break;
    }
  }
  return true;
}

bool BinaryLBModelequalepsilon::semRows()
{
  //Initializations
  computeUil();
  if(!seStepRows()) return false;
  //M-step
  mStepRows();
  return true;
}

// Run EM algorithm on data matrix m_Vjk_
bool BinaryLBModelequalepsilon::emCols()
{
  //Initializations
  computeVjk();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!eStepCols()) return false;
    //M-step
    mSteplogRhol();
    m_Ykl_old2_ = m_Ykl_;
    mStepCols();
    if(((m_Ykl_-m_Ykl_old2_)/m_Ykl_).abs().sum()<Mparam_.epsilon_int_)
    {
      break;
    }
  }
  return true;
}

// Run CEM algorithm on data matrix m_Vjk_
bool BinaryLBModelequalepsilon::cemCols()
{
  //Initializations
  computeVjk();
  MatrixReal  m_sumjl(Mparam_.nbCol_,Mparam_.nbcolclust_);
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //CE-step
    if(!ceStepCols()) return false;
    //M-step
    m_Ykl_old2_ = m_Ykl_;
    mStepCols();

    if(((m_Ykl_-m_Ykl_old2_)/m_Ykl_).abs().sum()<Mparam_.epsilon_int_)
    {
      break;
    }
  }

  return true;
}

bool BinaryLBModelequalepsilon::semCols()
{
  //Initializations
  computeVjk();
  if(!seStepCols()) return false;
  //M-step
  mStepCols();
  return true;
}

bool BinaryLBModelequalepsilon::GibbsRows()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool BinaryLBModelequalepsilon::GibbsCols()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

STK::Real BinaryLBModelequalepsilon::computeLnLikelihood()
{
  likelihood_ = (dimprod_*(Epsilon_*std::log(Epsilon_/(1.-Epsilon_)) + std::log(1.-Epsilon_))
              + v_Tk_.dot(v_logPiek_) + v_Rl_.dot(v_logRhol_)
              - (m_Tik_.prod((RealMin + m_Tik_).log()) ).sum()
              - (m_Rjl_.prod((RealMin + m_Rjl_).log()) ).sum()
			    );
  return likelihood_;
}
/* @return the number of free parameters of the distribution of a block.*/
int BinaryLBModelequalepsilon::nbFreeParameters() const
{ return 1;}

void BinaryLBModelequalepsilon::consoleOut()
{
#ifdef COVERBOSE
  std::cout<<"Output Model parameter:"<<"\nakl:\n"<<m_Akl_<<"\nepsilon:\n"<<Epsilon_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
  std::cout<<"ICL: "<<iclCriteriaValue()<<"\n";
#endif
}


//Compute change in Alpha and set the terminate variable accordingly
bool BinaryLBModelequalepsilon::initStopCriteria()
{
  return(std::abs((Epsilon_-Epsilontemp_)/Epsilon_) < Mparam_.initepsilon_);
}
//Compute change in Alpha and set the terminate variable accordingly
void BinaryLBModelequalepsilon::parameterStopCriteria()
{
  STK::Real relativechange = ((m_Ykl_-m_Ykl_old1_)/m_Ykl_old1_).abs().sum();
  if(relativechange<Mparam_.epsilon_)
    stopAlgo_ = true;
  else
    stopAlgo_ = false;

  // Update Ykl for outer loop
  m_Ykl_old1_ = m_Ykl_;
}

MatrixBinary const& BinaryLBModelequalepsilon::arrangedDataClusters()
{
  arrangedDataCluster<MatrixBinary>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

#ifndef RPACKAGE
void BinaryLBModelequalepsilon::displayCluster()
{
  CImg<unsigned char>  cluster(Mparam_.nbCol_,nbRow_,1,1,0);
  CImg<unsigned char>  data(Mparam_.nbCol_,nbRow_,1,1,0);

  MatrixBinary m_ClusterDataij_ = arrangedDataClusters();

  // Assign value to images
  for ( int i = 0; i <Mparam_.nbRow_; ++i) {
    for ( int j = 0; j < Mparam_.nbCol_; ++j) {
      if(m_ClusterDataij_(i,j) == 0)
        cluster(j,i) = (unsigned char)(0);
      else
        cluster(j,i) = (unsigned char)(255);

      if(m_Dataij_(i,j) == 0)
        data(j,i) = (unsigned char)(0);
      else
        data(j,i) = (unsigned char)(255);

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

void BinaryLBModelequalepsilon::saveThetaInit()
{
  m_Akltemp_ = m_Akl_;
  Epsilontemp_ = Epsilon_;
}

void BinaryLBModelequalepsilon::modifyTheta()
{
  m_Akltemp_ = m_Akl_;
  Epsilontemp_ = Epsilon_;

  v_logPiektemp_ = v_logPiek_;
  v_logRholtemp_ = v_logRhol_;
  m_Rjltemp_ = m_Rjl_;
  m_Tiktemp_ = m_Tik_;

  Lmax_ = likelihood_;
}

void BinaryLBModelequalepsilon::copyTheta()
{
  m_Akl_   = m_Akltemp_;
  Epsilon_ = Epsilontemp_;

  v_logPiek_ = v_logPiektemp_;
  v_logRhol_ = v_logRholtemp_;

  m_Tik_ = m_Tiktemp_;
  m_Rjl_ = m_Rjltemp_;
  commonFinalizeOutput();

  likelihood_ = computeLnLikelihood();
}

void BinaryLBModelequalepsilon::mStepFull()
{
  if(!Mparam_.fixedproportions_)
  {
    v_logRhol_=(v_Rl_/Mparam_.nbCol_).log();
    v_logPiek_=(v_Tk_/Mparam_.nbRow_).log();
  }

  m_Ykl_ = m_Tik_.transpose()*m_Dataij_.cast<STK::Real>()*m_Rjl_;
  m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
  m_Akl_ = (m_Ykl_>=m_Tk_Rl_);
  Epsilon_= (m_Ykl_-(v_Tk_*v_Rl_.transpose()).prod(m_Akl_.cast<STK::Real>()) ).abs().sum()/dimprod_;

}

STK::Real BinaryLBModelequalepsilon::iclCriteriaValue(){
  STK::Real criteria = 0.0;

  criteria+= lgamma(Mparam_.nbrowclust_*a_)+lgamma(Mparam_.nbcolclust_*a_)
      -(Mparam_.nbrowclust_+Mparam_.nbcolclust_)*lgamma(a_)
      +Mparam_.nbrowclust_*Mparam_.nbcolclust_*(lgamma(2*b_)-2*lgamma(b_))
      -lgamma(Mparam_.nbRow_+Mparam_.nbrowclust_*a_)
      -lgamma(Mparam_.nbCol_+Mparam_.nbcolclust_*a_);

  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    criteria+= lgamma(a_+ (v_Zi_== k).count());
  }

  for (int l = 0; l < Mparam_.nbcolclust_; ++l)
  {
    criteria+= lgamma(a_+ (v_Wj_==l).count());
  }

  STK::ArrayXXi temp0(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  STK::ArrayXXi temp1(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  MatrixBinary m_tempdata = (m_Dataij_==0);
  temp0 = (m_Zik_.transpose()*m_tempdata.cast<int>()*m_Wjl_)+b_;
  temp1 = (m_Zik_.transpose()*m_Dataij_.cast<int>()*m_Wjl_)+b_;
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
    {
      criteria+=lgamma(temp0(k,l))+lgamma(temp1(k,l));
    }
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
    {
      criteria-= lgamma(((v_Zi_== k).count())*((v_Wj_==l).count())+2*b_);
    }
  }
  return criteria;
}

