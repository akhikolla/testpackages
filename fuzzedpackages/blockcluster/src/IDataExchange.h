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

/** @file ICoClust.h
 *  @brief
 **/

#ifndef IDataExchange_H_
#define IDataExchange_H_

#include <iostream>
//#include "conversion.h"
#include "coclust/src/typedefs/typedef.h"
#include "coclust/src/Models/ICoClustModel.h"
#include "coclust/src/Algorithms/IAlgo.h"
#include "coclust/src/Initialization/IInit.h"
#include "coclust/src/Strategy/IStrategy.h"
#include "coclust/src/InputParameters/InputParameters.h"
#include <RTKpp.h>

/** @brief
 *
 */
class IDataExchange
{
  public:
    inline IDataExchange(){ initializeParamEnum();}
    void initializeParamEnum();

    virtual void dataInput(Rcpp::S4 & obj) = 0;
    virtual void dataOutput(Rcpp::S4& obj,ICoClustModel*, bool) = 0;
    virtual void instantiateModel(ICoClustModel*& model) = 0;

    void instantiateAlgo(IAlgo*& algo,IStrategy*& strat);
    void instantiateInit(IInit*& init);
    void setInput(Rcpp::S4 & obj);

    StrategyParameters& GetStrategyParameters(){return strategyParam_;}
    void setStrategyParameters(StrategyParameters& strat){strategyParam_ = strat;}
    virtual ~IDataExchange();

  protected:
    Strategy           strategy_;
    StrategyParameters strategyParam_;
    ModelParameters    Mparam_;
    VectorInt          v_rowlabels_,v_collabels_;
    //RVectorInt;

    std::map<std::string,Algorithm>      S_Algorithm;
    std::map<std::string,StopCriteria>   S_StopCriteria;
    std::map<std::string,DataType>       S_DataType;
    std::map<std::string,Initialization> S_Init;
    std::map<std::string,Model>          S_Model;
};
#endif /*IDataExchange_H_*/
