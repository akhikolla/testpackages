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
 * Project:  stkpp::Model
 * created on: 22 juil. 2011
 * Purpose: define the class IUnivStatModel.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_ModelBernoulli_pj.h
 *  @brief In this file we define the class ModelBernoulli_pj.
 **/

#ifndef STK_MODELBERNOULLI_PJ_H
#define STK_MODELBERNOULLI_PJ_H

#include <cmath>

#include "STK_Model_Util.h"
#include "STK_IMultiStatModel.h"
#include <STatistiK/include/STK_Law_Bernoulli.h>

namespace STK
{
// forward declaration
template <class Data, class WColVector = CVectorX> class ModelBernoulli_pj;
class Bernoulli_pjParameters;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the StatModelTraits for the ModelBernoulli_pj
 *  model
 **/
template<class Data_, class WColVector_>
struct StatModelTraits< ModelBernoulli_pj<Data_, WColVector_> >
{
  /** Type of the container storing the data */
  typedef DataBridge<Data_> Data;
  /** Type of the array storing the weights of the data */
  typedef WColVector_ WColVector;
  /** Type of the data in the container */
  typedef typename Data::Type Type;
  /** Type of the parameters of the ModelBernoulli_pj */
  typedef Bernoulli_pjParameters Parameters;
};

} // hidden

/** @ingroup StatModels
 *  Structure encapsulating the parameters of a Joint Bernoulli model.
 */
class Bernoulli_pjParameters
{
  public:
    /** default constructor */
    Bernoulli_pjParameters() {}
    /** constructor with fixed size */
    Bernoulli_pjParameters( Range const& size)
                          : prob_(size, 0.5)
                          , lnProb_(size, -Const::_LN2_)
                          , ln1mProb_(size, -Const::_LN2_)
    {}
    /** copy constructor. @param param the parameters to copy. */
    Bernoulli_pjParameters( Bernoulli_pjParameters const& param)
                          : prob_(param.prob_)
                          , lnProb_(param.lnProb_)
                          , ln1mProb_(param.ln1mProb_)
    {}
    /** destructor */
    ~Bernoulli_pjParameters() {}
    /** @return the range of the data set */
    inline Range const& range() const { return range_;}
    /** @return the probability of success of the jth law */
    inline CPointX const& prob() const { return prob_;}
    /** @return the probability of success of the jth law */
    inline CPointX const& lnProb() const { return lnProb_;}
    /** @return the probability of success of the jth law */
    inline CPointX const& ln1mProb() const { return ln1mProb_;}

    /** set the probability of success of the jth law */
    inline void setProb(int j, Real const& prob)
    { prob_[j] = prob;
      if (prob>0) { lnProb_[j] = std::log(prob);}
      else        { lnProb_[j] = -Arithmetic<Real>::infinity();}
      if (prob<1) { ln1mProb_[j] = std::log(1.-prob);}
      else        { ln1mProb_[j] = -Arithmetic<Real>::infinity();}
    }
    /** resize the parameters only if the range is modified, otherwise, stay
     *  with the current values.
     *  @param range the range of the parameters (= range of the variables of the model)
     **/
    inline void resize(Range const& range)
    {
      if (range != range_)
      {
        prob_.resize(range); lnProb_.resize(range); ln1mProb_.resize(range);
        range_ = range;
      }
    }

    CPointX prob_;
    CPointX lnProb_;
    CPointX ln1mProb_;
    Range range_;
};

/** @ingroup StatModels
 * A joint Bernoulli model is a statistical model of the form:
 * following form
 * \f[
 *     f(\mathbf{x}_i = x|\theta) =
 *     \prod_{j=1}^p p_{j}^{x_i^j} (1-p_{j})^{1-x_i^j},
 *      \quad x_i^j\in\{0,1\}, \quad j=1,\ldots,p, \quad i=1,\ldots,n.
 * \f]
 *
 **/
template <class Data_, class WColVector_>
class ModelBernoulli_pj: public IMultiStatModel< ModelBernoulli_pj<Data_, WColVector_> >
{
  public:
    /** Type of the container storing the data */
    typedef DataBridge<Data_> Data;
    typedef typename hidden::Traits<Data_>::Row RowVector;
    /** Type of the array storing the weights of the data */
    typedef WColVector_ WColVector;
    /** Type of the data in the container */
    typedef typename Data::Type Type;
    /** Base class */
    typedef IMultiStatModel< ModelBernoulli_pj<Data_, WColVector_> > Base;
    using Base::p_data;
    using Base::param;

    /** default constructor. */
    ModelBernoulli_pj(): Base() {}
    /** Constructor with data set. */
    ModelBernoulli_pj( Data const& data): Base(data) {}
    /** Constructor with a ptr on the data set. */
    ModelBernoulli_pj( Data const* p_data): Base(p_data) {}
    /** Copy constructor. */
    ModelBernoulli_pj( ModelBernoulli_pj const& model): Base(model) {}
    /** destructor */
    ~ModelBernoulli_pj() {}

    /** @return the vector of the probabilities of the observations */
    inline CPointX const& prob() const { return param().prob();}
    /** vector of the log probabilities of the observations */
    inline CPointX const& lnProb() const { return param().lnProb();}
    /** vector of the log probabilities of the reversed observations */
    inline CPointX const& ln1mProb() const { return param().ln1mProb();}

    /** compute the number of free parameters */
    inline int computeNbFreeParameters() const { return p_data()->dataij().sizeCols();}
    /** compute the log Likelihood of an observation. */
    Real computeLnLikelihood( RowVector const& rowData) const;
    /** compute the parameters */
    void computeParameters();
    /** compute the weighted parameters
     *  @param weights the weights of the samples
     **/
    void computeParameters(WColVector const& weights);
    /** Write the parameters on the output stream os */
    void writeParametersImpl(ostream& os) const;
};

/* compute the log Likelihood of an observation. */
template<class Data_, class WColVector_>
Real ModelBernoulli_pj<Data_, WColVector_>::computeLnLikelihood( RowVector const& rowData) const
{
  Real sum =0.;
  for (Integer j= rowData.begin(); j < rowData.end(); ++j)
  {
    sum += rowData[j] * lnProb()[j] + (1-rowData[j] * ln1mProb()[j] );
  }
  return sum;
}

/* compute the parameters */
template<class Data_, class WColVector_>
void ModelBernoulli_pj<Data_, WColVector_>::computeParameters()
{
  for (int j=p_data()->dataij().beginCols(); j < p_data()->dataij().endCols(); ++j)
  {
    Real sum=0.;
    int nbObs=p_data()->dataij().sizeRows();

    for (int i=p_data()->dataij().beginRows(); i<p_data()->dataij().endRows(); ++i)
    { (p_data()->dataij().elt(i,j) == binaryNA_) ? --nbObs : sum += p_data()->dataij().elt(i,j);}
    if (nbObs != 0) { param().prob_[j] = sum/nbObs;}
               else { param().prob_[j] = Arithmetic<Real>::NA();}
  }
}
/* compute the weighted parameters */
template<class Data_, class WColVector_>
void ModelBernoulli_pj<Data_, WColVector_>::computeParameters( WColVector const& weights)
{
  // compute
  for (int j=p_data()->dataij().beginCols(); j < p_data()->dataij().endCols(); ++j)
  {
    Real sum=0., wsum = 0.;
    for (int i=p_data()->dataij().beginRows(); i<p_data()->dataij().endRows(); ++i)
    { if (p_data()->dataij().elt(i,j) != binaryNA_)
      { sum  += weights[i]*p_data()->dataij().elt(i,j);
        wsum += weights[i];
      }
    }
    if (wsum != 0) { param().prob_[j] = sum/wsum;}
              else { param().prob_[j] = Arithmetic<Real>::NA();}
  }
}

/* Write the parameters on the output stream os */
template<class Data_, class WColVector_>
void ModelBernoulli_pj<Data_, WColVector_>::writeParametersImpl(ostream& os) const
{
  os << _T("prob = ") << prob();
  os << _T("lnProb = ") << lnProb();
  os << _T("ln1mProb = ") << ln1mProb();
}

} // namespace STK

#endif /* STK_MODELBERNOULLI_PJ_H */
