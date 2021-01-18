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

/** @file CoCluster.h
 *  @brief Declares class CoCluster which acts as entry point to do co-clustering.
 */

#ifndef COCLUSTER_H_
#define COCLUSTER_H_

#include "../Strategy/IStrategy.h"
class IInit;
class ICoClustModel;
class IAlgo;


/**@brief This class acts as entry point to perform co-clustering. This class provides functions for setting
 * algorithm, model and Initialization. This class also provides a CoCluster::run method that sets various pointers
 * in various classes and finally call the IAlgo::run method to perform co-clustering.
 */
class CoCluster
{
  public:
    /**Default Constructor*/
    inline CoCluster(): p_Strategy_(0), p_Model_(0), p_Init_(0), p_Algo_(0) {};
    /**Default Destructor*/
    inline ~CoCluster(){};

    /** This function perform Co-Clustering.*/
    bool run();
    /**It sets the algorithm to be run.*/
    void setStrategy(IStrategy * );
    /**It  sets the model to be run.*/
    void setModel(ICoClustModel * );
    /** It sets the Initialization method */
    void setInit(IInit * );
    /** It sets the algorithm to be run */
    void setAlgo(IAlgo * );

  private:
    /**Pointer to Strategy */
    IStrategy * p_Strategy_;
    /**Pointer to Model*/
    ICoClustModel * p_Model_;
    /**Pointer to Initialization*/
    IInit * p_Init_;
    /**Pointer to algorithm*/
    IAlgo * p_Algo_;
};

inline void CoCluster::setStrategy(IStrategy * strategy)
{ p_Strategy_ = strategy;}

inline void CoCluster::setModel(ICoClustModel * Model)
{ p_Model_ = Model;}

inline void CoCluster::setInit(IInit * Init)
{ p_Init_ = Init;}

inline void CoCluster::setAlgo(IAlgo * algo)
{ p_Algo_ = algo;}

#endif /* COCLUSTER_H_ */
