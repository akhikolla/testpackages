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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 14 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/


#ifndef STK_IMIXTURE_H
#define STK_IMIXTURE_H

/**@file STK_IMixture.h
 * @brief define the main interface for linking specific mixture model to the
 * composer.
 */

#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>
#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayNumber.h>

namespace STK
{

class IMixtureStatModel;

/** @ingroup Clustering
 *  @brief Interface base class for all the mixture models that will be
 *  processed by the composer.
 *
 *  Any mixture that can be used and add to a mixed model have to implement
 *  the pure virtual methods defined in this interface. With each mixture is
 *  associated a data set identified by an idData string.
 *
 *  Using an idData is not mandatory but is useful in case data set is stored
 *  in a DataHandler.
 *
 *  @sa STK::IDataHandler
 */
class IMixture
{
  public:
    /**Constructor with identification character
     * @param idData Identification string of the data associated to this mixture.
     * @note The Id is provided by the framework if the associated data set is in
     * a STK::IDataHandler struct.
     */
    IMixture( String const& idData = String());
    /**copy constructor.
     * @warning The pointer on the composer is not copied and is set to 0: it has
     * to be set again.
     * @param mixture the mixture to copy */
    IMixture( IMixture const& mixture);
    /** Virtual destructor. */
    virtual ~IMixture();

    /** @return the Id data of the mixture */
    inline String const& idData() const { return idData_;}
    /** @return A constant pointer on the composer. */
    inline IMixtureStatModel const* const p_composer() const { return p_composer_;}

    /** set the mixture composer to the mixture
     *  @param p_composer the composer to set
     **/
    void setMixtureModel( IMixtureStatModel const* p_composer);

    /**This is a standard clone function in usual sense. It must be defined to
     * provide new object of your class with values of various parameters equal
     * to the values of calling object. In other words, this is equivalent to
     * polymorphic copy constructor.
     * @return New instance of class as that of calling object.
     */
    virtual IMixture* clone() const  = 0;
    /**This is a standard create function in usual sense. It must be defined to
     * provide new object of your class with correct behavior.
     * In other words, this is equivalent to virtual constructor.
     * @return New instance of class as that of calling object.
     */
    virtual IMixture* create() const  = 0;
    /** @brief This function must be defined in derived class for initialization
     *  of the mixture parameters.
     *  This method should create any container needed by the model, resize
     *  them and initialize them.
     *  Since this method can be used when create is called, its main
     *  purpose should be to reset the mixture parameters, while leaving the
     *  data unchanged.
     *  This function will be called once the model is created and data is set.
     */
    virtual void initializeStep() = 0;
    /** @brief This function should be used in order to initialize randomly the
     *  parameters of the mixture.
     */
    virtual void randomInit()  = 0;
    /** @brief This function is equivalent to mStep and must be defined to update
     *  parameters.
     */
    virtual void paramUpdateStep() = 0;
    /** This function must be defined to return the posterior probability (PDF)
     * for corresponding sample_num and Cluster_num.
     * @param sample_num Sample number
     * @param Cluster_num Cluster number
     * @return the value of component probability in log scale
     */
    virtual Real lnComponentProbability(int sample_num, int Cluster_num) = 0;
    /** This function must return the number of free parameters.
     *  @return Number of free parameters
     */
    virtual int nbFreeParameter() const = 0;
    /** This function must return the number of missing value in data set identified by idData_.
     *  @return Number of missing values
     */
    virtual int nbMissingValues() const = 0;
    /** @brief This function should be used for Imputation of data.
     *  The default implementation (in the base class) is to do nothing.
     */
    virtual void imputationStep() {/**Do nothing by default*/}
    /** @brief This function must be defined for simulation of all the latent
     * variables and/or missing data excluding class labels. The class labels
     * will be simulated by the framework itself because to do so we have to
     * take into account all the mixture laws.
     */
    virtual void samplingStep() {/**Do nothing by default*/};
    /** @brief This function should be used to store any intermediate results
     * during various iterations after the burn-in period.
     * @param iteration Provides the iteration number beginning after the burn-in
     * period.
     */
    virtual void storeIntermediateResults(int iteration) {/**Do nothing by default*/}
    /**@brief This step can be used to signal to the mixtures that they must
     * release the stored results. This is usually called if the estimation
     * process failed.
     **/
    virtual void releaseIntermediateResults() {/**Do nothing by default*/}
    /** @brief set the parameters of the model.
     *  This function should be used to set the parameters computed using the
     *  intermediate results. This method will be called after the long-run and
     *  before the finalize step.
     **/
    virtual void setParametersStep() {/**Do nothing by default*/}
    /** @brief This step can be used by developer to finalize any thing. It will
     *  be called only once after we finish running the estimation algorithm.
     */
    virtual void finalizeStep() {/**Do nothing by default*/}
    /** This function can be used to write summary of parameters on to the output stream.
     *  @param out Stream where you want to write the summary of parameters.
     */
    inline virtual void writeParameters(std::ostream& out) const
    {
#ifdef STK_MIXTURE_DEBUG
    stk_cout<< _T("You need to override this method in your mixture!");
#endif
    }

  protected:
    /** This function can be used in derived classes to get number of samples.
     *  @return Number of samples.
     */
    int nbSample() const;
    /** This function can be used in derived classes to get number of classes.
     *  @return Number of classes.
     */
    int nbCluster() const;
    /** This function can be used in derived classes to get estimated number
     *  of individuals from the framework.
     *  @return Pointer to the number of individuals.
     */
    CPointX const* p_pk() const;
    /** This function can be used in derived classes to get proportions from the framework.
     *  @return Pointer to proportions.
     */
    CPointX const* p_tk() const;
    /** This function can be used in derived classes to get posterior probabilities from the framework.
     *  @return Pointer to tik.
     */
    CArrayXX const* p_tik() const;
    /** This function can be used in derived classes to get class labels from the framework.
     *  @return Pointer to zi.
     */
    CVectorXi const* p_zi() const;

  private:
    /** pointer on the main composer model */
    const IMixtureStatModel* p_composer_;
    /** Id name of the mixture */
    String idData_;
};

} // namespace STK

#endif /* STK_IMIXTURE_H */
