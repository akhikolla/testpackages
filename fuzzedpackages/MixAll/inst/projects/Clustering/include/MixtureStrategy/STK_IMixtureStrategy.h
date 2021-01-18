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

/** @file STK_IMixtureStrategy.h
 *  @brief In this file we define the interface base class for strategies to use in order
 *  to estimate a mixture model.
 **/


#ifndef STK_IMIXTURESTRATEGY_H
#define STK_IMIXTURESTRATEGY_H

namespace STK
{
// forward declarations
class IMixtureComposer;
class IMixtureInit;

/** @ingroup Clustering
 *  Interface base class for all the strategies */
class IMixtureStrategy: public IRunnerBase
{
  public:
    /** default constructor.
     *  @param p_model the model to estimate
     **/
    IMixtureStrategy( IMixtureComposer*& p_model);
    /** copy constructor
     *  @param strategy the strategy to copy
     **/
    IMixtureStrategy( IMixtureStrategy const& strategy);
    /** destructor */
    virtual ~IMixtureStrategy();

    /** @return get number of tries of strategy */
    inline int nbTry() const { return nbTry_;}

    /** set the number of tries of each strategies.
     * @param nbTry the number of tries to set */
    inline void setNbTry(int nbTry) { nbTry_ = nbTry;}
    /** set the initialization method to use
     * @param  p_init the initialization method to use */
    inline void setMixtureInit(IMixtureInit* p_init) { p_init_ = p_init;}

  protected:
    /** number of tries of each strategies (1 by default) */
    int nbTry_;
    /** reference on the main model */
    IMixtureComposer*& p_model_;
    /** initialization method */
    IMixtureInit* p_init_;
    /** store a model in p_model_ if it is better.
     * @param p_otherModel the model to store
     **/
    void storeModel(IMixtureComposer*& p_otherModel);
};


}  // namespace STK

#endif /* STK_IMIXTURESTRATEGY_H */
