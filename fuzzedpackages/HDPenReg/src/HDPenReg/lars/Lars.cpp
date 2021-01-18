/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : quentin.grimonprez@inria.fr
*/

/*
 * Project:  MPAGenomics::
 * created on: 13 févr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Lars.cpp
 *  @brief In this file, methods associates to @c Lars.
 **/

#include "../larsRmain.h"

using namespace STK;
using namespace std;

namespace HD
{
  //Constructors
/*
  * Constructor
  * @param X matrix of data, a row=a individual
  * @param y response
  */
Lars::Lars( CArrayXX const& X,CVectorX const& y, bool intercept)
          : n_(X.sizeRows())
          , p_(X.sizeCols())
          , maxSteps_(3*min(n_,p_))
          , X_(X), y_(y)
          , muX_(Range(1,p_))
          , path_(maxSteps_)
          , isActive_(Range(1,p_), false)
          , toIgnore_(Range(1,p_), false)
          , nbActiveVariable_(0)
          , nbIgnoreVariable_(0)
          , activeVariables_(Range(1,0))
          , step_(0)
          , mu_()
          , eps_(Arithmetic<Real>::epsilon())
          , Xi_(Range(1,n_),Range(1,1))
          , qrX_(Xi_)
          , c_(Range(1,0))
          , intercept_(intercept)
          , msg_error_()
{ initialization();}

/*
 * @param X matrix of data, a row=a individual
 * @param y response
 * @param maxStep number of maximum step to do
 * @param eps epsilon (for 0)
 */
Lars::Lars( CArrayXX const& X,CVectorX const& y, int maxSteps, bool intercept, Real eps)
          : n_(X.sizeRows())
          , p_(X.sizeCols())
          , maxSteps_(maxSteps)
          , X_(X)
          , y_(y)
          , muX_(Range(1,p_))
          , path_(maxSteps_)
          , isActive_(Range(1,p_), false)
          , toIgnore_(Range(1,p_), false)
          , nbActiveVariable_(0)
          , nbIgnoreVariable_(0)
          , activeVariables_(Range(1,0))
          , step_(0)
          , mu_()
          , eps_(eps)
          , Xi_(Range(1,n_), Range(1,1), 0.)
          , qrX_(Xi_)
          , c_( Range(1,0) )
          , intercept_(intercept)
          , msg_error_()
{ initialization();}

/* initialization of algorithm
 */
void Lars::initialization()
{
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering Lars::initialization")<<endl;
#endif
  if(intercept_)
  {
    //we center y
    mu_ = y_.mean();
    y_ -= mu_;
    muX_ = meanByCol(X_);
    X_  -= Const::VectorX(X_.rows()) * muX_.transpose();
  }
  else
  {
    mu_ =0;
    muX_.zeros();
  }
#ifdef LARS_DEBUG
    print(y_,"y_","y centered");
    print(X_,"X_","X_ centered");
    print(muX_,"muX_","muX computed");
#endif
  c_ = X_.transpose()*y_;
  Xi_.reserveCols(min(n_,p_));
  Xi_.shift(1,1);

#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::initialization done")<<endl;
  print(Xi_,"Xi_","initialization done");
  print(c_,"c_","initialization done");
  print(activeVariables_,"activeVariables_","initialization done");
#endif
}


/*
 * compute Cmax the maximum correlation
 */
Real Lars::computeCmax()
{
  Real Cmax = 0.;
  for(int i = isActive_.begin(); i < isActive_.end(); i++)
  { if(!isActive_[i]) { Cmax = std::max(Cmax, std::abs(c_[i]));}}
  return Cmax;
}

/*
 * search non active variable with the greatest correlation
 * @param Cmax correlation max
 * @param newId a vector containing the index of variable to potentially add
 */
void Lars::computeAddSet(Real Cmax, vector<int>& newId ) const
{
  for(int i = isActive_.begin(); i < isActive_.end(); i++)
  {
    if(!isActive_[i])//update with only variable non active
    {
      //don't care about ignore variable
      if( !toIgnore_[i] & (abs(c_[i]) >= Cmax-eps_) )
      { newId.push_back(i);}
    }
  }
}


/*
 * update the QR decomposition of Xi
 * @param idxVar index of active variable to add
 * @param signC sign of correlation of active variable
 * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
 */
void Lars::updateR(int idxVar,VectorXi &signC, pair<bool,vector<int> > &action)
{
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering Lars::updateR")<<endl;
#endif
  //update Xi_
  Xi_.pushBackCols(1);
  Xi_.shift(1,1); // in case
  Xi_.col(Xi_.lastIdxCols()) = X_.col(idxVar);
  //update the QR decomposition
  qrX_.pushBackCol(Xi_.col(Xi_.lastIdxCols()));

#ifdef VERBOSE
    cout<<"Step "<<step_<<" : Variable "<< idxVar<<" added"<<endl;
//      cout<<std::abs(qrX_.R()( min(n_,nbActiveVariable_+1), nbActiveVariable_+1) )<<std::endl;
#endif

  //check if the variable added is not colinear with an other
  if( std::abs(qrX_.R()( min(n_,nbActiveVariable_+1), nbActiveVariable_+1) ) < eps_ )
  {
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::updateR popBackCol")<<endl;
#endif
    //we cancel the add of the variable in the qr decomposition
    qrX_.popBackCols();
#ifdef LARS_DEBUG
      cout<<"Step "<<step_<<" : Variable "<< idxVar<<" dropped (colinearity)"<<endl;
#endif
    //toIgnore_[idxVar+1]=true;//the variable is add to the ignore set
    toIgnore_[idxVar]=true;//the variable is add to the ignore set
    nbIgnoreVariable_++;
    Xi_.popBackCols(1);
    Xi_.shift(1,1); // in case
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::updateR popBack done")<<endl;
#endif
  }
  else
  {
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::updateR pushBack")<<endl;
#endif
    //update action
    action.second.push_back(idxVar);
    action.first = true;
    //add the index to the set of active variable
    activeVariables_.pushBack(1);
    activeVariables_.back() = idxVar;
    activeVariables_.shift(1); // in case

    nbActiveVariable_++;
    isActive_[idxVar] = true;
    //compute signC
    signC.pushBack(1);
    signC[nbActiveVariable_] = ( c_[idxVar] > 0 ) ? 1 : -1;
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::updateR pushBack done")<<endl;
#endif
//      cout<<"add "<<idxVar<<"  cor "<<c_[idxVar]<<" sign "<<signC[nbActiveVariable_]<<endl;
  }
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::updateR done")<<endl;
#endif
}


/*
 * compute inv(Xi'*Xi)*1 from qr decomposition
 * @param Gi1 for stock inv(Xi'*Xi)*1
 * @param signC sign of correlation of active variable
 */
void Lars::computeGi1(CVectorX &Gi1, VectorXi const& signC) const
{
  CVectorX v(Range(1,nbActiveVariable_));
  Gi1.resize(Range(1,nbActiveVariable_));
  Gi1=0;
//    signC.shift(1); // in case
  for(int i=1; i<=nbActiveVariable_; i++)
  {
    //resolve R'*v=signC(i)*e_i
    v=0.;
    v[i]=signC[i];
    v[1] /= qrX_.R()(i,i);
    for(int j=2; j<=nbActiveVariable_; j++)
    {
      for(int k=1; k<=(j-1); k++)
        v[j] -= qrX_.R()(k,j)*v[k];
      v[j]/=qrX_.R()(j,j);
    }
    //resolve R*signC*z=v
    v[nbActiveVariable_] *= signC[nbActiveVariable_]/qrX_.R()(nbActiveVariable_,nbActiveVariable_);
    for(int j=nbActiveVariable_-1; j>0; j--)
    {
      for(int k=j+1; k <= nbActiveVariable_; k++)
        v[j] -= signC[k] * qrX_.R()(j,k) * v[k];
      v[j] /= signC[j] * qrX_.R()(j,j);
    }
    Gi1+=v;
  }
}

/*
 * Compute gammahat for the update of coefficient in add case
 * @param Aa norm of the inverse of G
 * @param a X' * equiangular vector
 * @param Cmax correlation max
 * @return gammaHat a real
 */
Real Lars::computeGamHat(Real const& Aa, CVectorX const& a, Real Cmax) const
{
  Real gamHat(Cmax/Aa), gam(0);
  for(int i=isActive_.begin(); i<isActive_.end(); i++)
  {
    //only for the non active variable and non ignored variable
    if(!isActive_[i] && !toIgnore_[i])
    {
      //gamma is the min only on the positive value of gam1 and gam2
      if(Aa!=a[i]) { if( (gam = (Cmax-c_[i])/(Aa-a[i]) ) > eps_ ) {gamHat=min(gamHat,gam);}}
      if(Aa!=-a[i]){ if( (gam = (Cmax+c_[i])/(Aa+a[i]) ) > eps_ ) {gamHat=min(gamHat,gam);}}
    }
  }
  return gamHat;
}

/*
 * Compute gammaTilde for the update of coefficient in drop case
 * @param w Aa*Gi1 @see computeGi1
 * @param idxMin we stock the index (in the activeVariable vector) of the variable with the min value
 * @return gammatilde a real
 */
Real Lars::computeGamTilde(CVectorX const& w,vector<int> &idxMin) const
{
  Real gamTilde(std::numeric_limits<Real>::max()),gam(0);
  idxMin.erase(idxMin.begin(),idxMin.end());
  for(int i=1; i <= path_.lastState().size(); i++)
  {
    if(w[i]) gam = -path_.lastVarCoeff(i)/w[i];
    //we search the minimum only on positive value
    if(gam > eps_)
    {
      if(gam < gamTilde)
      {
        gamTilde = gam;
        idxMin.erase(idxMin.begin(),idxMin.end());
      }
      if(gam==gamTilde) { idxMin.push_back(i);}
    }
  }
  return gamTilde;
}

/*
 * Update the coefficient of the path
 * @param gamma gammaHat or gammaTilde
 * @param w Aa*Gi1 @see computeGi1
 * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
 * @param isAddCase true if we add a variable
 * @param dropId id top potentially drop
 */
void Lars::updateBeta(Real gamma, CVectorX const& w, pair<bool,vector<int> > action, bool isAddCase, vector<int> dropId)
{
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering Lars::updateBeta")<<endl;
  print(w,"w", "updateBeta");
#endif
  if( (action.first) && (isAddCase) )
  {
    //add situation
    path_.addCaseUpdate(gamma,w,action.second);
  }
  else
  {
    if( (action.first) && (!isAddCase) )
    {
      //add situation with a drop
      vector<int> drop(dropId.size());
      for(int i = 0; i < (int) dropId.size(); i++)
      { drop[i] = activeVariables_[dropId[i]];}
      path_.addWithDropCaseUpdate(gamma,w,action.second,drop,dropId);
    }
    else
    {
      if((!action.first) && (isAddCase))
      {
        //update after a drop situation
        path_.update(gamma,w);
      }
      else
      {
        vector<int> drop(dropId.size());
        for(int i = 0; i < (int) dropId.size(); i++)
        {  drop[i] = activeVariables_[dropId[i]];}
        //drop step after a drop step
        path_.dropAfterDropCaseUpdate(gamma,w,drop,dropId);
      }
    }
  }
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::updateBeta done")<<endl;
#endif
}


/* dropStep
 * downdate qr decomposition,  X and signC
 * @param idxVar index of active variable to drop
 * @param signC sign of correlation of active variable
 */
void Lars::dropStep(vector<int> const& idxVar, VectorXi &signC)
{
#ifdef VERBOSE
    for(int i = 0; i < (int) idxVar.size(); i++)
      cout<<"Step "<<step_+1<<" : Variable "<< activeVariables_[idxVar[i]]<<" dropped"<<endl;
#endif
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering Lars::dropStep")<<endl;
  print(qrX_.R(),"qrX_.R()", "dropStep");
  print(qrX_.Q(),"qrX_.Q()", "dropStep");
  print(Xi_,"Xi_", "dropStep");
  print(signC,"signC", "dropStep");
  print(activeVariables_,"activeVariables_", "dropStep");
#endif

    //idxVar are sort in the increasing order, we erase first the element with the greater index
  for(int i = idxVar.size()-1; i >= 0 ; i--)
  {
#ifdef LARS_DEBUG
  stk_cerr << _T("idxVar[i]=")<<idxVar[i]<<endl;
  stk_cerr << _T("qrX_.eraseCol(")<< idxVar[i] << _T(")") <<endl;
  if (idxVar[i] < qrX_.R().beginCols())
  { stk_cerr << _T("idxVar[i] < qrX_.R().beginCols()") <<endl;}
  if (qrX_.R().lastIdxCols() < idxVar[i])
  { stk_cerr << _T("qrX_.R().lastIdxCols() < idxVar[i]") <<endl;}
#endif
  //downdate R
    qrX_.eraseCol(idxVar[i]);
    //downdate Xi_
#ifdef LARS_DEBUG
  stk_cerr << _T("Xi_.eraseCols(")<< idxVar[i] << _T(")") <<endl;
#endif
    Xi_.eraseCols(idxVar[i]);
#ifdef LARS_DEBUG
  stk_cerr << _T("signC.erase(")<< idxVar[i] << _T(")") <<endl;
#endif
    signC.erase(idxVar[i]);
#ifdef LARS_DEBUG
  stk_cerr << _T("activeVariables_[idxVar[i]]=")<<activeVariables_[idxVar[i]]<<endl;
#endif
    isActive_[activeVariables_[idxVar[i]]] = false;
#ifdef LARS_DEBUG
  stk_cerr << _T("activeVariables_.erase(")<< idxVar[i] << _T(")") <<endl;
#endif
    activeVariables_.erase(idxVar[i]);
  }
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::dropStep done")<<endl;
#endif
}

/*
 * first step
 * @param Cmax correlation max
 * @param newId vector of index of active variable to add
 * @param signC sign of correlation of active variable
 * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
 * @param Aa norm of the inverse of G
 * @param Gi1 for stock inv(Xi'*Xi)*1
 * @param w Aa*Gi1
 * @param u unit vector making equal angles with the column of Xi
 * @param a X' * equiangular vector       * @param Gi1
 * @param gam the step for update coefficients
 * @return
 */
bool Lars::firstStep( Real &Cmax, vector<int> &newId
                    , VectorXi &signC
                    , pair<bool, vector<int> > &action
                    , Real &Aa
                    , CVectorX &Gi1
                    , CVectorX &w
                    , CVectorX &u
                    , CVectorX &a
                    , Real &gam)
{
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering Lars::firstStep")<<endl;
#endif
  step_++;
  //computation of correlation
  Cmax=computeCmax();
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::firstStep. computeCmax() done")<<endl;
#endif
  if( Cmax < eps_*100)
  {
    step_--;
#ifdef VERBOSE
    std::cout << "Correlation max is equal to 0.";
#endif
    msg_error_ = "Correlation max is equal to 0.";
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering Lars::firstStep. Return false")<<endl;
#endif
    return false;
  }
  newId.resize(0);
  computeAddSet(Cmax, newId);
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::firstStep. computeAddSet() done")<<endl;
#endif
  if(newId.size()==0)
  {
    step_--;
#ifdef VERBOSE
    std::cout << "No variable selected for add in the add step."<<std::endl;
#endif
    msg_error_ = "No variable selected for add in the add step.";
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering Lars::firstStep. Return false")<<endl;
#endif
    return false;
  }
  addCmax(Cmax);
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::firstStep. addCmax() done")<<endl;
#endif
  for (vector<int>::iterator it = newId.begin() ; it != newId.end(); it++)
  { firstUpdateR(*it,signC,action);}
#ifdef LARS_DEBUG
  stk_cerr << _T("firstUpdateR done")<<endl;
#endif

  //compute the inverse of G
  computeGi1(Gi1,signC);
#ifdef LARS_DEBUG
  stk_cerr << _T("computeGi1 done")<<endl;
#endif
  //compute Aa
  Aa = 1/sqrt(Gi1.sum());
  //compute w
  w = Gi1 * signC * Aa;
  //compute equiangular vector
  u = Xi_ * w;
  //computation of gamma hat
  //if the number of active variable is equal to the max number authorized, we don't search a new index
  if( nbActiveVariable_ == min(n_-1,p_-nbIgnoreVariable_) )
  {  gam=Cmax/Aa;}
  else
  {
    //computation of a
    a = X_.transpose() * u;
    //computation of gamma hat
    gam = computeGamHat(Aa,a,Cmax);
  }
  //update beta
  vector<int> vide;
  updateBeta(gam,w,action,true,vide);
#ifdef LARS_DEBUG
  stk_cerr << _T("updateBeta done")<<endl;
#endif
  //update of c
  c_ -= (X_.transpose() * u) * gam;
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::firstStep done. Return true")<<endl;
#endif
  return true;
}

/*
 * updateR only for the first step
 * @see updateR
 */
void Lars::firstUpdateR(int idxVar, VectorXi &signC, pair<bool,vector<int> > &action)
{
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering Lars::firstUpdateR")<<endl;
#endif
  //create Xi_ and qrXi_
  Xi_.col(Xi_.lastIdxCols()) = X_.col(idxVar);
  Xi_.shift(1,1); // in case
  qrX_.setData(Xi_);
  qrX_.run();
#ifdef LARS_DEBUG
  stk_cerr << _T("qrX_.run() done.")<<endl;
#endif

#ifdef VERBOSE
    cout<<"Step 1 : Variable "<< idxVar<<" added"<<endl;
#endif

  //check if the variable added is not colinear with an other
  if( std::abs(qrX_.R()( min(n_,nbActiveVariable_+1), nbActiveVariable_+1) ) < eps_ )
  {
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering popBackCols")<<endl;
#endif
    qrX_.popBackCols();
#ifdef VERBOSE
      cout<<"Step 1 : Variable "<< idxVar<<" dropped (colinearity)"<<endl;
#endif

    toIgnore_[idxVar+1] = true;//the variable is add to the ignore set
    nbIgnoreVariable_++;
    Xi_.popBackCols(1);
    Xi_.shift(1,1); // in case
#ifdef LARS_DEBUG
  stk_cerr << _T("popBackCols done.")<<endl;
#endif
  }
  else
  {
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering pushBack")<<endl;
#endif
    nbActiveVariable_++;
    //Add the idx to the active set
    activeVariables_.pushBack(1);
    activeVariables_.shift(1); // in case
    //activeVariables_.back() = idxVar;
    activeVariables_.back() = idxVar;
    isActive_[idxVar] = true;

    //compute signC
    signC.pushBack(1);
    signC.shift(1); // in case
    signC.back() = ( c_[idxVar] > 0 ) ? 1 : -1;

    action.first = true;
    action.second.push_back(idxVar);
#ifdef LARS_DEBUG
  stk_cerr << _T("push_back done.")<<endl;
#endif
  }
}

/* run lars algorithm*/
void Lars::run()
{
#ifdef VERBOSE
  cout<<"######################################"<<endl;
  cout<<"########### LARS ALGORITHM ###########"<<endl;
  cout<<"######################################"<<endl<<endl;
  Chrono::start();
#endif
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering Lars::run")<<endl;
#endif
  //initialization();
  bool isAddCase(true), continuer;
  vector<int> dropId;
  Real Aa(0),gam(0),gammaTilde(0),Cmax(0);

  CVectorX Gi1(Range(1,1));
  CVectorX w(Range(1,0));
  VectorXi signC(Range(1,0));
  CVectorX a( Range(1,p_), 0), u( Range(1,n_), 0);
  pair<bool,vector<int> > action;

  vector<int> newId;
  newId.reserve(p_);

  //Gi1.reserveCols(min(n_-1,p_));
  //w.reserveCols(min(n_-1,p_));
  signC.reserveCols(min(n_-1,p_));

  continuer=firstStep(Cmax,newId,signC,action,Aa,Gi1,w,u,a,gam);
  if (!continuer) return;
  //we stop, if we reach maxStep or if there is no more variable to add
  Real oldCmax;
  while( (step_< maxSteps_) && ( nbActiveVariable_ < min( n_-1, (p_-nbIgnoreVariable_) ) ) )
  {
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::run step_ = ") << step_<<endl;
#endif
    step_++;
    oldCmax = Cmax;
    //computation of correlation
    Cmax = computeCmax();
    if( Cmax < eps_*100)
    {
      step_--;
#ifdef LARS_DEBUG
        std::cout << "Correlation max is equal to 0."<<std::endl;
#endif
      msg_error_ = "Correlation max is equal to 0.";
      break;
    }
    //if correlation max increased, we stop, Cmax must decreased
    if( Cmax > oldCmax)
    {//          stk_cout<<Cmax<<endl;
      step_--;
#ifdef LARS_DEBUG
        std::cout << "Correlation max has increased."<<std::endl;
#endif
      msg_error_ = "Correlation max has increased";
      break;
    }
    //add case : update of QR decomposition, active set and X'*X
    if(isAddCase)
    {
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::run isAddCase")<<endl;
#endif
      newId.resize(0);
      computeAddSet(Cmax, newId);
      if(newId.size() == 0)
      {
        step_--;
#ifdef VERBOSE
          std::cout << "No variable selected for add in the add step."<<std::endl;
#endif
        msg_error_ = "No variable selected for add in the add step.";
        break;
      }
      action.second.erase(action.second.begin(),action.second.end());
      for(vector<int>::iterator it = newId.begin() ; it != newId.end(); it++)
      { updateR(*it,signC,action);}
    }
    else
    {
      action=make_pair(false,dropId);
    }
    addCmax(Cmax);
    //compute the inverse of G
    computeGi1(Gi1,signC);
    //compute Aa
    Aa = 1/sqrt(Gi1.sum());
    //compute w
    w = Gi1*signC*Aa;
    //compute equiangular vector
    u = Xi_*w;
    //computation of gamma hat
    //if the number of active variable is equal to the max number authorized, we don't search a new index
    if( nbActiveVariable_ == min(n_-1, p_-nbIgnoreVariable_) )
    {  gam = Cmax/Aa;}
    else
    {
      //computation of a
      a = X_.transpose()*u;
      //computation of gamma hat
      gam = computeGamHat(Aa,a,Cmax);
    }
    //computation of gamma tilde
    gammaTilde = computeGamTilde(w,dropId);
    if( gammaTilde < gam )
    {
      gam = gammaTilde;
      isAddCase = false;
      nbActiveVariable_ -= dropId.size();
    }
    else
    { isAddCase = true;}
    //update beta
    updateBeta(gam,w,action,isAddCase,dropId);
    //update of c_
    if( nbActiveVariable_ == min(n_-1, p_-nbIgnoreVariable_) )
    {
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::run update c_ with nbActiveVariable_ == min(n_-1, p_-nbIgnoreVariable_)")<<endl;
#endif
      c_ -= (X_.transpose()*u)*gam;
    }
    else
    {
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::run update c_")<<endl;
#endif
      c_ -= a * gam;
    }
    //drop situation
    if(!isAddCase) { dropStep(dropId,signC);}
    //path_.states(step_).printCoeff();
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::run update c_ done")<<endl;
#endif
  }
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::run while terminated")<<endl;
  stk_cerr<<endl<<"Algorithm finished"<<endl;
  stk_cerr<<"Number of steps: "<<step_<<endl;
  stk_cerr<<"Number of active variables: "<<nbActiveVariable_<<endl;
#endif

#ifdef VERBOSE
  Real t1 = Chrono::elapsed();
  cout<<endl<<"Algorithm finished in "<<t1<<"s"<<endl;
  cout<<"Number of steps: "<<step_<<endl;
  cout<<"Number of active variables: "<<nbActiveVariable_<<endl;
#endif
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::run done")<<endl;
#endif
}


/*
 * predict the path for a ratio index = l1norm/l1normmax or a specific value of lambda
 * @param X new data for predict the response
 * @param index index (lambda or fraction) where the response is estimated.
 * @param lambdaMode if TRUE, index corresponds to a value of lambda, if FALSE, index is a real between 0 and 1
 * corresponding to ratio between the l1 norm of estimates to calculate and l1 norm max of solution
 * @param predicted response (be modified)
 */
void Lars::predict(CArrayXX const& X, Real index, bool lambdaMode, CVectorX &yPred)
{
#ifdef LARS_DEBUG
stk_cerr << _T("Entering Lars::predict")<<endl;
#endif
  yPred = mu_;
  //std::cout<<"la "<<index<<"  "<<path_.lambda(0)<<"  a"<<((index == 0.) && !lambdaMode)<<"  b"<<((index>=path_.lambda(0)) && lambdaMode)<<"    c"<<(((index == 0.) && !lambdaMode) || ((index>=index-path_.lambda(0)) && lambdaMode))<<std::endl;
  //index = 0 : all coefficients are equal to 0
  //lambda >= lambda max : all coefficients equal to 0
  if( ((index == 0.) && !lambdaMode) || ((index>=path_.lambda(0)) && lambdaMode) )
    return ;

  //index = 1 : coefficients of the last step
  //lambda <= lambda min : coefficients of the last step
  if( ((index == 1.) && !lambdaMode) || ((index <= path_.lambda().back()) && lambdaMode) )
  {
    int lastStep = path_.size()-1;//stocké dans un vector index à 0
    int nbVar = path_.lastState().sizeRows();

    for(int i = yPred.begin(); i < yPred.end(); i++)
      for(int j = 1; j <= nbVar; j++)
        yPred[i] += (X(i, varIdx(lastStep,j)) -muX_[varIdx(lastStep,j)]) * coefficient(lastStep,j);
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::predict done")<<endl;
#endif
    return ;
  }

  //fraction >0 and <1 or lambda <lambda max and > lambda min
  VectorX l1norm(path_.l1norm());
  VectorX listIndex;
  int ind = 1;

  if(lambdaMode)
  {
    //std::cout<<"la1b"<<std::endl;
    int lambdasize=path_.lambda().size();
    listIndex.resize(lambdasize);
    for(int i = 0; i < lambdasize; i++)
    {  listIndex[i+1] = path_.lambda(i);}
    //stk_cout<<listIndex<<std::endl;
    //stk_cout<<l1norm<<std::endl;
    while(listIndex[ind] > index)
      ind++;
  }
  else
  {
    listIndex = l1norm;
    index *= listIndex.back();
    while(listIndex[ind] < index)
      ind++;
  }
  Real l1normNew(0.);
  if(!lambdaMode) { l1normNew = index;}
  else
  { l1normNew = l1norm[ind-1]+(index-listIndex[ind-1])/(listIndex[ind]-listIndex[ind-1])*(l1norm[ind]-l1norm[ind-1]);}
  //std::cout<<"encadre : "<<listIndex[ind-1]<<"   "<<listIndex[ind]<<"   "<<l1norm[ind-1]<<"   "<<l1norm[ind]<<std::endl;
  //std::cout<<index<<"  : "<<path_.states(ind-2).l1norm()<<"   "<<l1normNew<<"   "<<path_.states(ind-1).l1norm()<<"   "<<std::endl;
  //compute coefficient
  int nbCoeff = std::max(path_.states(ind-2).size(),path_.states(ind-1).size());
  Array2DVector< pair<int,Real> > coeff(Range(1,nbCoeff));
  //coeff.move(computeCoefficients(path_.states(ind-2),path_.states(ind-1),path_.evolution(ind-2),fraction));
  computeCoefficients(path_.states(ind-2),path_.states(ind-1),path_.evolution(ind-2),l1normNew,coeff);

  for( int i = yPred.begin(); i < yPred.end(); i++)
    for( int j = 1; j <= coeff.sizeRows(); j++)
      yPred[i] += (X(i, coeff[j].first) - muX_[coeff[j].first] ) * coeff[j].second ;
#ifdef LARS_DEBUG
  stk_cerr << _T("Lars::predict done")<<endl;
#endif
}

void Lars::computeCoefficients(PathState const& state1,PathState const& state2,pair<std::vector<int> ,std::vector<int> > const& evolution, Real const& l1norm, Array2DVector< pair<int,Real> > &coeff)
{
  //Array2DVector< pair<int,Real> > coeff(std::max(state1.size(),state2.size()));
  int maxSize = state1.size() + evolution.first.size();
  coeff.resize(Range(1,maxSize));
  if(evolution.second.size()==0)
  {//no drop variable
    int j(1);
    for(j = 1; j <= state1.size(); j++)
      coeff[j]=make_pair(state1.varIdx(j),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(j), state2.varCoeff(j)));
    //add variable case
    if(evolution.first.size()!=0)
    {
      for(int i = 0; i < (int) evolution.first.size(); i++)
        coeff[j]=make_pair( evolution.first[i], computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, 0., state2.varCoeff(j)));
    }
  }
  else
  {
    //delete variable case
    int i = 1;
    for( int j = 0; j < (int) evolution.second.size(); j++)
    {
      //while we don't meet the delete variable, variable has the same index in the two sets
      while(evolution.second[j]!=state1.varIdx(i))
      {
        coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i), state2.varCoeff(i)));
        i++;
      }
      //compute coefficient for the delete variable
      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i),0.));
      i++;
    }
    //compute coefficient for the other variable
    while(i < state1.size()+1)
    {
      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i), state2.varCoeff(i-1)));
      i++;
    }
    //drop with an add variable
    if(evolution.first.size()!=0)
    {
      for(int j = 0; j < (int) evolution.first.size(); j++)
        coeff[i]=make_pair(evolution.first[j], computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, 0., state2.varCoeff(i-1)));
    }
  }
}



}//end namespace
