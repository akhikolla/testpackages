/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016 Serge Iovleff

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
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IMixtureBridge.h
 *  @brief In this file we define the interface for all the bridge classes
 *  between the mixtures and the composer.
 **/

#ifndef STK_IMIXTUREBRIDGE_H
#define STK_IMIXTUREBRIDGE_H

#include "STK_IMixture.h"

#include <DManager/include/STK_DataBridge.h>

namespace STK
{

/** @ingroup Clustering
 *  @brief Interface base class for the bridges of the STK++ mixture.
 *
 *  This interface produce a default implementation to all the virtual
 *  functions required by the STK::IMixture interface by delegating
 *  to the derived bridged Mixture the computations.
 *
 *  The virtual methods to implement in derived class are
 *  @code
 *    virtual Derived* create() const;
 *    virtual Derived* clone() const;
 *  @endcode
 *
 *  Pseudo virtual method to implement in derived class are
 *  @code$
 *    safeValue(j); // return a safe value for column j
 *  @endcode
 *  The members of this interface are the following
 *  @code
 *    Mixture mixture_;           // concrete class deriving from STK::IMixtureDensity
 *    MissingIndexes v_missing_;  // array with the indexes (i,j) of the missing values
 *    Data* p_dataij_;
 *  @endcode
 *
 */
template<class Derived>
class IMixtureBridge: public IMixture, public IRecursiveTemplate<Derived>
{
  public:
    typedef typename hidden::MixtureBridgeTraits<Derived>::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Data Data;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Type Type;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Parameters Parameters;
    typedef typename hidden::MixtureBridgeTraits<Derived>::MissingIndexes MissingIndexes;
    typedef typename hidden::MixtureBridgeTraits<Derived>::MissingValues MissingValues;
    typedef typename MissingIndexes::const_iterator ConstIterator;

    enum
    {
      // class of mixture
      idMixtureClass_ = hidden::MixtureBridgeTraits<Derived>::idMixtureClass_
    };


  protected:
    /** default constructor.
     *  @param p_data pointer on the DataBridge wrapping the data set
     *  @param idData id name of the data
     *  @param nbCluster number of cluster
     **/
    IMixtureBridge( Data* p_data, String const& idData, int nbCluster);
    /** copy constructor
     *  @param bridge the IMixtureBridge to copy
     **/
    IMixtureBridge( IMixtureBridge const& bridge);
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    IMixtureBridge( Mixture const& mixture, String const& idData, int nbCluster);
    /** destructor */
    virtual ~IMixtureBridge() {}

  public:
    // getters
    /** @return the mixture */
    inline Mixture const& mixture() const { return mixture_;}
    /** @return coordinates of the missing values in the data set */
    inline MissingIndexes const& v_missing() const { return v_missing_;}
    /** @return the pointer on the data set */
    inline Data* const p_dataij() const { return p_dataij_;}
    /** @return number of the missing values in data set */
    inline int nbMissing() const { return v_missing_.size();}
    /** get the (imputed) missing values of a data set.
     *  @note In C++11, it will be possible to use a tuple rather that this pair of pair...
     *  @param data the array to return with the missing values
     **/
    void getMissingValues( MissingValues& data) const;

    // getter and setter of the parameters
    /** get the parameters of the model */
    inline void getParameters( Parameters& param) const { param = mixture_.param_;}
    /** set the parameters of the model */
    inline void setParameters( Parameters const& param) { mixture_.param_ = param;}

    // getter and setter of the parameters using an array
    /** This function is used in order to get the current values of the parameters
     *  in an array.
     *  @param param array that will store the parameters of the mixture.
     */
    template<class Array>
    inline void getParameters(Array& param) const { mixture_.getParameters(param);}
    /** This function is used in order to set the current values of the
     *  parameters to the parameters using an array.
     *  @param param the array/expression with the parameters of the mixture to
     *  store in the Parameters.
     **/
    template<class Array>
    inline void setParameters(ExprBase<Array> const& param)
    { mixture_.param_.setParameters(param.asDerived());}

    // start default implementation of virtual method inherited from IMixture
    /** @brief Initialize the mixture model before its use by the composer.
     *
     *  The parameters values are set to their default values if the mixture_
     *  is newly created. if IMixtureBridge::initializeStep is used during a
     *  cloning, mixture class have to take care of the existing values of the
     *  parameters.
     **/
    inline virtual void initializeStep()
    {
      if (!p_composer())
        STKRUNTIME_ERROR_NO_ARG(IMixtureBridge::initializeStep,composer is not set);
      if (!mixture_.initializeStep()) throw Clust::initializeStepFail_;
    }
     /** @brief This function must be defined to return the component probability distribution
      *  function (PDF) for corresponding sample i and cluster k.
     *   @param i,k sample and cluster numbers
     *   @return the log-component probability
     **/
    inline virtual Real lnComponentProbability(int i, int k)
    { return mixture_.lnComponentProbability(i, k);}
    /** @brief This function is equivalent to Mstep and must be defined to update
     *  parameters.
     **/
    inline virtual void paramUpdateStep()
    { if (!mixture_.run( p_tik(), p_tk())) throw Clust::mStepFail_;}
    /** @brief This function should be used in order to initialize randomly the
     *  parameters of the mixture.
     **/
    inline virtual void randomInit() { mixture_.randomInit( p_tik(), p_tk());};
    /** This function must return the number of free parameters.
     *  @return Number of free parameters
     **/
    inline virtual int nbFreeParameter() const { return mixture_.computeNbFreeParameters();}
    /** This function must return the number of missing value in data set identified by idData_.
     *  @return Number of missing values
     */
    virtual int nbMissingValues() const { return v_missing_.size();}
    /** @brief This function should be used to store any intermediate results
     *  during various iterations after the burn-in period.
     *  @param iteration Provides the iteration number beginning after the burn-in
     *  period.
     **/
    inline virtual void storeIntermediateResults(int iteration)
    { mixture_.param_.updateStatistics();}
    /**@brief This step can be used to signal to the mixtures that they must
     * release the stored results. This is usually called if the estimation
     * process failed.
     **/
    inline virtual void releaseIntermediateResults()
    { mixture_.param_.releaseStatistics();}
    /** @brief set the parameters of the model.
     *  This step should be used to set the initialize the statistics to the
     *  current value of the parameters.
     **/
    inline virtual void setParametersStep()
    { mixture_.param_.setStatistics();}
    /** @brief This step can be used by developer to finalize any thing. It will
     *  be called only once after we finish running the estimation algorithm.
     */
    inline virtual void finalizeStep() { mixture_.finalizeStep();}
    /** @brief This function should be used for imputation of data.
     *  The default implementation (in base class IMixture) is to do nothing.
     */
    inline virtual void imputationStep();
    /** @brief This function must be defined for simulation of all the latent
     * variables and/or missing data excluding class labels. The class labels
     * will be simulated by the framework itself because to do so we have to
     * take into account all the mixture laws.
     */
    inline virtual void samplingStep();
    /** This function can be used to write summary of parameters to the output stream.
     *  @param os Stream where you want to write the summary of parameters.
     */
    inline virtual void writeParameters(ostream& os) const
    { mixture_.writeParameters(p_tik(), os);}

  protected:
    /** utility function for lookup the data set and find missing values coordinates.
     *  @return the number of missing values
     **/
   virtual std::vector< std::pair<int,int> >::size_type findMissing();

    /** @brief This function will be used once for imputation of missing data
     *  at the initialization step (@see initializeStep).
     *
     *  By default missing values in a column will be replaced by the mean or map (maximum a
     *  posteriori) value of the column.
     *  Derived class as to implement the safeValue(i,j) method.
     *  @note This method is virtual as derived class may want to implement a more
     *  efficient algorithm.
     **/
    virtual void removeMissing();

    /** The Mixture to bridge with the composer */
    Mixture mixture_;
    /** vector with the coordinates of the missing values */
    MissingIndexes v_missing_;
    /** pointer on the data set */
    Data* p_dataij_;
};

/* default constructor.
 *  @param p_data pointer on the DataBridge that will be used by the bridge.
 *  @param idData id name of the mixture model
 *  @param nbCluster number of cluster
 **/
template< class Derived>
IMixtureBridge<Derived>::IMixtureBridge( Data* p_dataij, String const& idData, int nbCluster)
                                       : IMixture(idData)
                                       , mixture_(nbCluster)
                                       , v_missing_()
                                       , p_dataij_(p_dataij)
{// find coordinates of the missing values
  findMissing();
}
/* copy constructor */
template< class Derived>
IMixtureBridge<Derived>::IMixtureBridge( IMixtureBridge const& bridge)
                                       : IMixture(bridge)
                                       , mixture_(bridge.mixture_)
                                       , v_missing_(bridge.v_missing_)
                                       , p_dataij_(bridge.p_dataij_)
{}
// implementation
template< class Derived>
void IMixtureBridge<Derived>::imputationStep()
{
  for(ConstIterator it = v_missing().begin(); it!= v_missing().end(); ++it)
  { p_dataij_->elt(it->first, it->second) = mixture_.impute(it->first, it->second, p_tik()->row(it->first) );}
}
// implementation
template< class Derived>
void IMixtureBridge<Derived>::samplingStep()
{
  for(ConstIterator it = v_missing().begin(); it!= v_missing().end(); ++it)
  { p_dataij_->elt(it->first, it->second) = mixture_.sample(it->first, it->second, p_tik()->row(it->first));}
}

// implementation
/* protected constructor to use in order to create a bridge.
 *  @param mixture the mixture to copy
 *  @param idData id name of the mixture
 *  @param nbCluster number of cluster
 **/
template< class Derived>
IMixtureBridge<Derived>::IMixtureBridge( Mixture const& mixture, String const& idData, int nbCluster)
                                       : IMixture(idData)
                                       , mixture_(mixture)
                                       , v_missing_()
                                       , p_dataij_(0)
{}

template< class Derived>
std::vector< std::pair<int,int> >::size_type IMixtureBridge<Derived>::findMissing()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering IMixtureBridge<Derived>::findMissing()\n");
#endif
  if (p_dataij_)
  {
    for (int j=p_dataij_->beginCols(); j< p_dataij_->endCols(); ++j)
    {
      for (int i=p_dataij_->beginRows(); i< p_dataij_->endRows(); ++i)
      {
        if (Arithmetic<Type>::isNA(p_dataij_->elt(i,j)))
        {
          v_missing_.push_back(std::pair<int,int>(i,j));
        }
      }
    }
  }
 return v_missing_.size();
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("IMixtureBridge<Derived>::findMissing() terminated, nbMiss= ") << v_missing_.size() << _T("\n");
#endif
}


/* This function will be used for imputation of the missing data
 *  at the initialization (initializationStep). By default missing values
 *  in a column will be replaced by the mean or map (maximum a
 *  posteriori) value of the column.
 **/
template< class Derived>
void IMixtureBridge<Derived>::removeMissing()
{
  if (p_dataij_)
  {
    Type value = Type();
    int j, old_j = Arithmetic<int>::NA();
    for(ConstIterator it = v_missing().begin(); it!= v_missing().end(); ++it)
    {
      j = it->second; // get column
      if (j != old_j)
      {
        old_j = j;
        value =  this->asDerived().safeValue(j);
      }
      p_dataij_->elt(it->first, j) = value;
    }
  }
}

/* get the (imputed) missing values of a data set.
 *  @note In C++11, it will be possible to use a tuple rather that this pair of pair...
 *  @param data the array to return with the missing values
 **/
template< class Derived>
void IMixtureBridge<Derived>::getMissingValues( MissingValues& data) const
{
  data.resize(v_missing_.size());
  for(size_t i = 0; i< v_missing_.size(); ++i)
  {
    data[i].first  = v_missing_[i];
    data[i].second =  p_dataij_->elt(v_missing_[i].first, v_missing_[i].second);
  }
}

} // namespace STK

#endif /* IMIXTUREBRIDGE_H */
