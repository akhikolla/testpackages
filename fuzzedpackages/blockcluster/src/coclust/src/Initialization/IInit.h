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

/** @file IInit.h
 *  @brief Declares abstract class IInit for Initializations.
 **/

#ifndef IINIT_H_
#define IINIT_H_

#include "../Models/ICoClustModel.h"

/** @brief This is an interface class for initialization algorithms.
 * It have only one abstract function IInit::run which when called, initialize
 * the various model parameters by calling the corresponding Initialization function
 * declared in IModel.
 * Every initialization algorithm must derive from this abstract class
 * and provide implementation for IInit::run method.
 */
class IInit
{
  protected:
    /** constructor */
    IInit():p_Model_(0) {};
    // pointer on model
    ICoClustModel * p_Model_;
  public:
    /** Interface for running initialization */
    virtual bool run() = 0;
    /** Refer IInit::p_Model_ to the instantiated model class inside CoClustermain.cpp
     */
    void setModel(ICoClustModel * model);
    /** Virtual Destructor*/
    virtual ~IInit(){};
};

inline void IInit::setModel(ICoClustModel * model)
{ p_Model_ = model;}
#endif /* IINIT_H_ */
