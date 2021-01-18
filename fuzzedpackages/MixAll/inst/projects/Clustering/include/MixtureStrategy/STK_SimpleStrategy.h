/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 3 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_MixtureStrategy.h
 *  @brief In this file we define the class implementing a simple strategy
 **/


#ifndef STK_SIMPLESTRATEGY_H
#define STK_SIMPLESTRATEGY_H

#include "STK_IMixtureStrategy.h"

namespace STK
{
// forward declaration
class IMixtureAlgo;

/** @ingroup Clustering
 *  helper structure encapsulating the parameters of the simple strategy
 **/
struct SimpleStrategyParam
{
  /** Constructor. Set default values */
  inline SimpleStrategyParam() : p_algo_(0) {}
  /** destructor */
  virtual ~SimpleStrategyParam();
  /** number of iterations in the Initialization */
  IMixtureAlgo* p_algo_;
};

/** @ingroup Clustering
 *  A SimpleStrategy is just nbTry long run.
 **/
class SimpleStrategy: public IMixtureStrategy
{
  public:
    /** default constructor.
     * @param p_model a reference pointer on the model to estimate
     **/
    inline SimpleStrategy( IMixtureComposer*& p_model): IMixtureStrategy(p_model), p_param_(0)
    {}
    /** copy constructor.
     *  @param strategy the strategy to copy
     **/
    inline SimpleStrategy( SimpleStrategy const& strategy): IMixtureStrategy(strategy), p_param_(0)
    {}

    /** destructor */
    inline virtual ~SimpleStrategy() { if (p_param_) delete p_param_;}
    /** clone pattern */
    inline virtual SimpleStrategy* clone() const { return new SimpleStrategy(*this);}
    /** set the parameters of the strategy
     * @param  p_param  the parameters of the strategy */
    inline void setParam(SimpleStrategyParam* p_param) { p_param_ = p_param;}

    /** run the strategy */
    virtual bool run();

  protected:
    SimpleStrategyParam* p_param_;
};


}  // namespace STK

#endif /* STK_SIMPLESTRATEGY_H */
