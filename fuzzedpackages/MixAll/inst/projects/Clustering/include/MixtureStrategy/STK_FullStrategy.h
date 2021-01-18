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

/** @file STK_FullStrategy.h
 *  @brief In this file we define the class implementing the full strategy class.
 **/


#ifndef STK_FULLSTRATEGY_H
#define STK_FULLSTRATEGY_H

#include "STK_IMixtureStrategy.h"

namespace STK
{
// forward declaration
class IMixtureAlgo;

/** @ingroup Clustering
 *  helper structure encapsulating the parameters of the Full strategy
 **/
struct FullStrategyParam
{  /** Constructor. Set default values */
    inline FullStrategyParam() : nbInitRun_(1), nbShortRun_(0), p_shortAlgo_(0) , p_longAlgo_(0)
    {}
    /** destructor */
    virtual ~FullStrategyParam();
    /** number of initialization run to perform */
    int nbInitRun_;
    /** number of short run to perform */
    int nbShortRun_;
    /** algorithm to use in short runs  */
    IMixtureAlgo* p_shortAlgo_;
    /** algorithm to use in long run  */
    IMixtureAlgo* p_longAlgo_;
};

/** @ingroup Clustering
 *  A FullStrategy is based on the following paradigm:
 *  - perform nbInitRun_ of the p_init_ initialization method, select the best initialization
 *  - perform nbShortRun of the shortAlgo with a small number of iterations and
 *   a high tolerance,
 *  - pick the best model obtained,
 *  - on this best model perform a long run.
 *  This strategy is used in Rmixmod R package.
 **/
class FullStrategy: public IMixtureStrategy
{
  public:
    /** default constructor.
     * @param p_model a reference pointer on the model to estimate
     **/
    inline FullStrategy( IMixtureComposer*& p_model): IMixtureStrategy(p_model), p_param_()
    {}
    /** copy constructor.
     *  @param strategy the strategy to copy
     **/
    inline FullStrategy( FullStrategy const& strategy): IMixtureStrategy(strategy), p_param_(0)
    {}
    /** destructor */
    inline virtual ~FullStrategy() { if (p_param_) delete p_param_;}
    /** clone pattern */
    inline virtual FullStrategy* clone() const { return new FullStrategy(*this);}
    /** set the parameters of the strategy
     * @param  p_param  the parameters of the Xem strategy
     **/
    void setParam(FullStrategyParam * p_param) { p_param_ = p_param;}

    /** run the strategy */
    virtual bool run();

  protected:
    FullStrategyParam* p_param_;
    /** Perform the Initialization step
     *  Initialize nbInitRun_ (should be  > 0) model and select the best model among them.
     *  @param p_bestModel a pointer initialized to
     **/
    bool initStep(IMixtureComposer*& p_bestModel);
};

}  // namespace STK

#endif /* STK_FULLSTRATEGY_H */
