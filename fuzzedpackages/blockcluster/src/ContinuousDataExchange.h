/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2012  Parmeet Singh Bhatia

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

/*
 * Project:  RCocluster
 * created on: Feb 22, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file ContinuousCoClust.h
 *  @brief 
 **/

#ifndef ContinuousDataExchange_H_
#define ContinuousDataExchange_H_

/** @brief
 *
 */
#include "IDataExchange.h"
#include "coclust/src/Models/ContinuousLBModel.h"
#include "coclust/src/Models/ContinuousLBModelequalsigma.h"

class ContinuousDataExchange: public IDataExchange
{
  public:
    ContinuousDataExchange(){};
    virtual void dataOutput(Rcpp::S4 & obj,ICoClustModel*,bool);
    virtual void dataInput(Rcpp::S4 & obj);
    virtual void instantiateModel(ICoClustModel*& model);
    inline const MatrixReal & GetData(){return m_Dataij_;}
    ~ContinuousDataExchange(){};
  protected:
    MatrixReal m_Dataij_;
};

#endif /* ContinuousDataExchange_H_ */
