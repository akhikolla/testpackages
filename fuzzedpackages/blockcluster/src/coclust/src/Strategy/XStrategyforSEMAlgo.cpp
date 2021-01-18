/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

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

 Author : Parmeet Bhatia
 Contact : bhatia.parmeet@gmail.com , serge.iovleff@stkpp.org

 */


/** @file XStrategyforSEMAlgo.cpp
 *  @brief Implements XStrategyAlgo class.
 **/

#include "XStrategyforSEMAlgo.h"


bool XStrategyforSEMAlgo::run()
{
#ifdef COVERBOSE
  std::cout<<"Running XStrategyforSEMAlgo\n";
#endif
  STK::Real Lmax = -RealMax,L1 = -RealMax,Lcurrent;
  int ntry_empty=0;
  bool non_empty = false;
  for ( int itry = 0; itry < strategyParam_.nbtry_; ++itry)
  {
    non_empty = false;
    L1 = -RealMax;
    for ( int ixem = 0; ixem < strategyParam_.nbxem_; ++ixem)
    {
      ntry_empty = -1;
      p_Model_->setEmptyCluster(true);
      while(p_Model_->isEmptyCluster() && ntry_empty<strategyParam_.nbinitmax_)
      {
        if (p_Init_->run())
        {
          for ( int itr = 0; itr < strategyParam_.nbiter_xem_; ++itr)
          {
            if (p_Algo_->run()) {
              Lcurrent = p_Model_->computeLnLikelihood();
              if(Lcurrent>=L1)
              {
                non_empty = true;
                L1 = Lcurrent;
                p_Model_->modifyTheta();
              }
            }
          }
        }
        ntry_empty++;
      }
    }
    if (!non_empty) continue;
    // start XEM iterations
    p_Model_->copyTheta();
    p_Model_->modifyTheta(); // initialize Theta max values to start values
    Lcurrent = p_Model_->computeLnLikelihood();
    if(!p_Model_->isEmptyCluster()&&Lcurrent>Lmax)
    {
      Lmax = Lcurrent;
      p_Model_->modifyTheta();
    }
    for ( int itr = 0; itr < strategyParam_.nbiter_XEM_; ++itr)
    {
      if (p_Algo_->run())
      {
        Lcurrent = p_Model_->computeLnLikelihood();
        if(Lcurrent>Lmax)
        {
          Lmax = Lcurrent;
          p_Model_->modifyTheta();
        }
      }
    }
  }
  // check if al try failed
  if(non_empty)
  {
    p_Model_->copyTheta();
    p_Model_->finalizeOutput();
#ifdef COVERBOSE
    std::cout<<"Strategy over.\n";
    p_Model_->consoleOut();
    std::cout<<"\nLmax:"<<Lmax<<"\n";
#endif
    return !p_Model_->isEmptyCluster();
  }
  return false;
}
