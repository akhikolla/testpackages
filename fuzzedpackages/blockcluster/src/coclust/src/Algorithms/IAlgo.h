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

/** @file IAlgo.h
 *  @brief Declares interface class for various algorithms.
 **/

#ifndef IALGO_H_
#define IALGO_H_
#include "../Models/ICoClustModel.h"

/** @brief This is an interface class for various algorithms available in the package.
 *
 */
class IAlgo
{
  public:
    IAlgo(){};
    virtual bool run() = 0;
    void setModel(ICoClustModel*);
    virtual ~IAlgo(){};

  protected:
    ICoClustModel * p_Model_;
};

inline void IAlgo::setModel(ICoClustModel* p_model)
{
  p_Model_ = p_model;
}

#endif /* IALGO_H_ */
