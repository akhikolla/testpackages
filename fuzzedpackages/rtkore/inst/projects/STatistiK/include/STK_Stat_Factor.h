/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::STatistiK
 * Purpose:  Compute factors of a set of variables.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_Factor.h
 *  @brief In this file we define and implement the Factor class.
 **/

#ifndef STK_STAT_FACTOR_H
#define STK_STAT_FACTOR_H

#include <Sdk/include/STK_IRunner.h>
#include <Sdk/include/STK_Macros.h>

namespace STK
{


namespace Stat
{
/** @ingroup StatDesc
 *  @brief Computation of the Factors of a 1D Container.
 *
 *  The class @c Factor is a factory class for computing the factors of a vector.
 *  The values can be of any type. the Coding is performed from the previous type
 *  in integer. The mapping is stored and preserved in a map array.
 **/
template <class Array>
class Factor: public IRunnerWithData<Array>
{
  public:
    typedef IRunnerWithData<Array> Base;
    typedef typename hidden::Traits<Array>::Row RowVector;
    typedef typename hidden::Traits<Array>::Col ColVector;
    typedef typename Array::Type Type;

    typedef std::map<Type, int> EncodingMap;
    typedef std::map<int, Type> DecodingMap;

    using Base::p_data_;

    /** Default Constructor. */
    Factor();
    /** Constructor.
     *  @param data a reference on the data set
     **/
    Factor( Array const& data);
    /** Constructor.
     *  @param p_data a pointer on the data set
     **/
    Factor( Array const* p_data);
    /** copy constructor.
     *  @param f the Factor to copy
     **/
    Factor( Factor const& f);
    /** virtual destructor.*/
    inline virtual ~Factor() {}

    /** clone pattern */
    inline virtual Factor* clone() const { return new Factor(*this);}

    /** @return vector with the factors encoded as integers */
    inline CVectorXi const& asInteger() const { return asInteger_;}
    /** @return vector with the levels */
    inline Array2DVector<Type> const& levels() const {return levels_;}
    /** @return vector with the counts of each level */
    inline VectorXi const& counts() const {return counts_;}
    /** @return the value of the first level */
    inline int const& firstLevel() const { return firstLevel_;}
    /** @return number of levels */
    inline int const& nbLevels() const { return nbLevels_;}
    /** @return encoding maps factor from Type -> int */
    inline EncodingMap const& encoder() const { return encoder_;}
    /** @return decoding maps factor from int -> Type */
    inline DecodingMap const& decoder() const { return decoder_;}

    /** set the value of the first level */
    inline void setFirstLevel(int firstLevel) { firstLevel_ = firstLevel;}

    /** run the estimation of the Factor statistics. **/
    virtual bool run();

  protected:
    /** vector with the levels in an integer format*/
    CVectorXi asInteger_;
    /** first level */
    int firstLevel_;
    /** Number of levels of each variables */
    int nbLevels_;
    /** vector with the levels */
    Array2DVector<Type> levels_;
    /** Array with the counts of each factor */
    VectorXi counts_;
    /** encoder of the levels */
    EncodingMap encoder_;
    /** decoder of the levels */
    DecodingMap decoder_;

    /** udpating method in case we set a new data set */
    virtual void update();
};

template <class Array>
Factor<Array>::Factor(): Base(), asInteger_(), firstLevel_(baseIdx), nbLevels_()
                       , levels_(), counts_(), encoder_() {}

template <class Array>
Factor<Array>::Factor( Array const& data): Base(data)
                                         , asInteger_(p_data_->range())
                                         , firstLevel_(baseIdx)
                                         , nbLevels_(0)
                                         , levels_()
                                         , counts_()
                                         , encoder_()
                                         , decoder_()
{}

/* Constructor.
 *  @param p_data a pointer on the data set
 **/
template <class Array>
Factor<Array>::Factor( Array const* p_data): Base(p_data)
                                           , asInteger_()
                                           , firstLevel_(baseIdx)
                                           , nbLevels_(0)
                                           , levels_()
                                           , counts_()
                                           , encoder_()
                                           , decoder_()
{
  if (p_data_)
  {
    asInteger_.resize(p_data_->range());
    nbLevels_= 0;
  }
}

template <class Array>
Factor<Array>::Factor( Factor const& f): Base(f), asInteger_(f.asInteger_)
                                       , firstLevel_(f.firstLevel_), nbLevels_(f.nbLevels_)
                                       , levels_(f.levels_), counts_(f.counts_)
                                       , encoder_(f.encoder_)
                                       , decoder_(f.decoder_)
{}

template <class Array>
void Factor<Array>::update()
{
   // if there is no data there is nothing to update
   if (p_data_)
   {
     asInteger_.resize(p_data_->rows());
     firstLevel_ = baseIdx;
     nbLevels_=0;
     levels_.clear();
     counts_.clear();
     encoder_.clear();
     decoder_.clear();
   }
}

template <class Array>
bool Factor<Array>::run()
{
  if (!p_data_)
  { this->msg_error_ = STKERROR_NO_ARG(FactorArray::run,data is not set);
    return false;
  }
  try
  {
    for (int i=p_data_->begin(); i< p_data_->end(); ++i)
    {
      // find coding
      Type idData = p_data_->elt(i);
      typename EncodingMap::const_iterator it = encoder_.find(idData);
      if (it != encoder_.end()) // levels already exist, just update the levels array
      { asInteger_[i] = it->second;
        counts_[it->second]++;  // add one to this level
      }
      else // find a new level to add
      {
        // create a new level and set it
        int lev = firstLevel_ + nbLevels_;
        asInteger_[i] = lev;
        encoder_.insert(std::pair<Type, int>(idData, lev));
        decoder_.insert(std::pair<int, Type>(lev, idData));
        levels_.push_back(idData);
        counts_.push_back(1); // start counting for this new level
        nbLevels_++;
      }
    }
  }
  catch (Exception const& error)
  {
    this->msg_error_ += _T("Error in Factor::run():\nWhat: ");
    this->msg_error_ += error.error();
    return false;
  }
  // no error
  return true;
}

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_FACTOR_H */
