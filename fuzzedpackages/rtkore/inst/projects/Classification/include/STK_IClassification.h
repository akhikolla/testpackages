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
 * Project:  stkpp::Classif
 * created on: 23 juin 2010
 * Purpose:  Interface base class for classification methods.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IClassification.h
 *  @brief In this file we define the Interface base class IClassification.
 **/

#ifndef STK_ICLASSIFICATION_H
#define STK_ICLASSIFICATION_H

#include <Sdk/include/STK_IRunner.h>
#include <Sdk/include/STK_Macros.h>

namespace STK
{

/** @ingroup Classif
 * @brief Interface base class for Regression methods.
 *
 * In machine learning and statistics, classification is the problem of
 * identifying to which of a set of categories (sub-populations) a new
 * observation belongs, on the basis of a training set of data containing
 * observations (or instances) whose category membership is known. An example
 * would be assigning a given email into "spam" or "non-spam" classes or
 * assigning a diagnosis to a given patient as described by observed
 * characteristics of the patient (gender, blood pressure, presence or absence
 * of certain symptoms, etc.). Classification is an example of pattern
 * recognition.
 *
 * In the terminology of machine learning, classification is considered an
 * instance of supervised learning, i.e. learning where a training set of
 * correctly identified observations is available. The corresponding
 * unsupervised procedure is known as clustering, and involves grouping data
 * into categories based on some measure of inherent similarity or distance.
 *
 * Often, the individual observations are analyzed into a set of quantifiable
 * properties, known variously as explanatory variables or features. These
 * properties may variously be categorical (e.g. "A", "B", "AB" or "O", for
 * blood type), ordinal (e.g. "large", "medium" or "small"), integer-valued
 * (e.g. the number of occurrences of a particular word in an email) or
 * real-valued (e.g. a measurement of blood pressure). Other classifiers work
 * by comparing observations to previous observations by means of a similarity
 * or distance function.
 *
 * An algorithm that implements classification, especially in a concrete
 * implementation, is known as a classifier. The term "classifier" sometimes
 * also refers to the mathematical function, implemented by a classification
 * algorithm, that maps input data to a category.
 *
 *
 * In this interface, the pure virtual function to implement are
 * @code
 *   virtual bool estimationStep() =0;
 *   virtual bool estimationStep( Weights_ const& weights) =0;
 *   virtual bool predictionStep() =0;
 *   virtual int computeNbFreeParameter() const =0;
 *   virtual YArray_ extrapolate(XArray_ const& x) const =0;
 * @endcode
 * The virtual function
 * @code
 *   virtual void initializeStep();
 * @endcode
 * can be overloaded.
 *
 * The default behavior of the @c run methods is
 * @code
 *    if (!initializeStep()) { flag = false;}
 *    if (!estimationStep()) { flag = false;}
 *    if (!finalizeStep()) { flag = false;}
 *    nbFreeParameter_ = computeNbFreeParameter();
 *    if (!predictionStep()) { flag = false;}
 * @endcode
 * and it can be overloaded in derived class.
 */
template <class YArray_, class XArray_, class Weights_>
class IClassification: public IRunnerSupervised<YArray_, XArray_, Weights_>
{
  protected:
    typedef IRunnerSupervised<YArray_, XArray_, Weights_> Base;
    using Base::p_x_;
    using Base::p_y_;
    /** Default constructor. Initialize the data members. */
    IClassification();
    /** Constructor. Initialize the data members.
     *  @param p_y,p_x pointer array with the observed output and output
     **/
    IClassification( YArray_ const* p_y, XArray_ const* p_x);
    /** Constructor. Initialize the data members.
     *  @param y,x arrays with the observed output and input
     **/
    IClassification( YArray_ const& y, XArray_ const& x);

  public:
    /** virtual destructor. */
    virtual ~IClassification() {}

    /**  @return number of class */
    inline int nbClass() const { return nbClass_;}
    /**  @return number of parameters of the classification function */
    inline int nbFreeParameter() const { return nbFreeParameter_;}

    /** run the computations. Default Implementation. */
    virtual bool run()
    {
      if (!p_x_ || !p_y_)
      { this->msg_error_ = STKERROR_NO_ARG(IClassification::run,missing data sets);
        return false;
      }
      bool flag=true;
      // perform any initialization step needed before the classification step
      if (!initializeStep()) { flag = false;}
      // compute the classification
      if (!estimationStep()) { flag = false;}

      // perform any post-operation needed after the classification step
      if (!finalizeStep()) { flag = false;}
      // Compute the number of parameter of the classification function.
      nbFreeParameter_ = computeNbFreeParameter();
      // create array of the predicted value and compute prediction
      if (!predictionStep()) { flag = false;}
      // return the result of the computations
      this->hasRun_ = true;
      return flag;
    }
    /** run the weighted computations.
     *  @param weights weights of the samples
     **/
    virtual bool run( Weights_ const& weights)
    {

      if (!p_x_ || !p_y_)
      { this->msg_error_ = STKERROR_NO_ARG(IClassification::run,missing data sets);
        return false;
      }
      bool flag=true;
      // perform any pre-operation needed before the classification step
      if (!initializeStep()) { flag = false;}
      // compute weighted classification
      if (!estimationStep(weights)) { flag = false;}

      // perform any post-operation needed after the classification step
      if (!finalizeStep()) { flag = false;}
      // Compute the number of parameter of the classification function.
      nbFreeParameter_ = computeNbFreeParameter();
      // create array of the predicted value and compute prediction
      if (!predictionStep()) { flag = false;}

      // return the result of the computations
      this->hasRun_ = true;
      return flag;
    }

  protected:
    /** number of class */
    int nbClass_;
    /** number of parameter of the classification method. */
    int nbFreeParameter_;

    /** @brief perform any computation needed before the call of the classification
     *  method. Default implementation is do nothing.
     */
    virtual bool initializeStep() {return true;}
    /** Compute the predicted outputs by the classification function and store the
     *  result in the p_predicted_ array. Default implementation is do nothing.
     **/
    virtual bool predictionStep() {return true;};
    /** @brief perform any computation needed after the call of the classification
     *  method. Default implementation is do nothing.
     */
    virtual bool finalizeStep() {return true;}

  private:
    /** compute the classification function. */
    virtual bool estimationStep() =0;
    /** compute the weighted classification function.
     * @param weights the weights of the samples
     **/
    virtual bool estimationStep(Weights_ const& weights) =0;
    /** Compute the number of parameter of the classification function.
     * @return the number of parameter of the classification function
     **/
    virtual int computeNbFreeParameter() const =0;
};

/* Default constructor. Initialize the data members. */
template <class YArray_, class XArray_, class Weights_>
IClassification<YArray_,XArray_,Weights_>::IClassification(): Base(), nbClass_(0), nbFreeParameter_(0) {}

/* Constructor. Initialize the data members.
 *  @param p_y,p_x pointer array with the observed output and output
 **/
template <class YArray_, class XArray_, class Weights_>
IClassification<YArray_,XArray_,Weights_>::IClassification( YArray_ const* p_y, XArray_ const* p_x)
                                                          : Base(p_y, p_x)
                                                          , nbClass_(0)
                                                          , nbFreeParameter_(0)
{}

/* Constructor. Initialize the data members.
 * @param y,x arrays with the observed output and input
 **/
template <class YArray_, class XArray_, class Weights_>
IClassification<YArray_,XArray_,Weights_>::IClassification( YArray_ const& y, XArray_ const& x)
                                                          : Base(y, x)
                                                          , nbClass_(0)
                                                          , nbFreeParameter_(0)
{}

} // namespace STK

#endif /* STK_ICLASSIFICATION_H */
