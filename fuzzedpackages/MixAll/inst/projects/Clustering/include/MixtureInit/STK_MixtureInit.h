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
 * created on: 24 août 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_MixtureInit.h
 *  @brief In this file we define the initialization methods.
 **/


#ifndef STK_MIXTUREINIT_H
#define STK_MIXTUREINIT_H

#include "STK_IMixtureInit.h"

namespace STK
{
/** @ingroup Clustering
 *  Implementation of the random initialization. This class will initialize the
 *  parameter by calling the randomInit() method of the model. */
class RandomInit: public IMixtureInit
{
  public:
    /** default constructor */
    inline RandomInit(): IMixtureInit() {}
    /** copy constructor
     * @param init the initialization to copy
     **/
    inline RandomInit(RandomInit const& init): IMixtureInit(init) {}
    /** destructor */
    inline virtual ~RandomInit(){}
    /** clone pattern */
    inline virtual RandomInit* clone() const { return new RandomInit(*this);}
    /** run the initialization by calling the randomInit method of the model.
     * @return @c true if no error occur, @c false otherwise*/
    virtual bool run();
};

/** @ingroup Clustering
 *  Initialization by simulating a realization of the class labels zi accordingly
 *  to the initial proportions.
 **/
class ClassInit: public IMixtureInit
{
  public:
    /** default constructor */
    inline ClassInit(): IMixtureInit() {}
    /** copy constructor
     * @param init the initialization to copy
     **/
    inline ClassInit(ClassInit const& init): IMixtureInit(init) {}
    /** destructor */
    inline virtual ~ClassInit(){}
    /** clone pattern */
    inline virtual ClassInit* clone() const { return new ClassInit(*this);}
    /** run the initialization by calling the randomClassInit method of the model.
     * @return @c true if no error occur, @c false otherwise*/
    virtual bool run();
};

/** @ingroup Clustering
 *  Initialization by simulating the tik accordingly to the initial
 *  proportions. */
class FuzzyInit: public IMixtureInit
{
  public:
    /** default constructor */
    inline FuzzyInit(): IMixtureInit() {}
    /** copy constructor
     *   @param init the initialization to copy
     **/
    inline FuzzyInit(FuzzyInit const& init): IMixtureInit(init) {}
    /** destructor */
    inline virtual ~FuzzyInit(){}
    /** clone pattern */
    inline virtual FuzzyInit* clone() const { return new FuzzyInit(*this);}
    /** run the algorithm on the model calling E step and M step.
     * @return @c true if no error occur, @c false otherwise*/
    virtual bool run();
};

} // namespace STK

#endif /* STK_MIXTUREINIT_H */
