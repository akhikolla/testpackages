/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2015  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 10 May 2016
 * Author:   Iovleff, serge.iovleff@stkpp.org
 **/

/** @file ILauncherBase.h
 *  @brief In this file we define the interface base class for launchers.
 **/


#ifndef STK_ILAUNCHERBASE_H
#define STK_ILAUNCHERBASE_H

#include "RDataHandler.h"

namespace STK
{

/** The ILauncherBase allow to create the composer or learner for estimate, learn
 *  or predict a mixture model with less effort.
 **/
class ILauncherBase: public IRunnerBase
{
  public:
    /** constructor with a list of component.
     *  @param model a reference on the current model
     **/
    ILauncherBase( Rcpp::S4 model);
    /** destructor. */
    virtual ~ILauncherBase();
    /** @return the model */
    inline Rcpp::S4 const& s4_model() const { return s4_model_;}

  protected:
    /** model from the R side */
    Rcpp::S4 s4_model_;

    /** data handler */
    RDataHandler handler_;
    /** kernel handler */
    KernelHandler kerHandler_;

    /** diagonal Gaussian mixture models manager */
    DiagGaussianMixtureManager<RDataHandler> diagGaussianManager_;
    /** Poisson mixture models manager */
    PoissonMixtureManager<RDataHandler> poissonManager_;
    /** gamma mixture models manager */
    GammaMixtureManager<RDataHandler> gammaManager_;
    /** categorical mixture models manager */
    CategoricalMixtureManager<RDataHandler> categoricalManager_;
    /** categorical mixture models manager */
    KernelMixtureManager kernelManager_;

    /** set parameters @c param to component identified by idData to the mixture model */
    bool setParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param);
    /** set parameters @c param to diagGaussian component identified by idData to the mixture model */
    void setDiagGaussianParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param);
    /** set parameters @c param to Poisson component identified by idData to the mixture model */
    void setPoissonParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param);
    /** set parameters @c param to gamma component identified by idData to the mixture model */
    void setGammaParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param);
    /** set parameters @c param to Categorical component identified by idData to the mixture model */
    void setCategoricalParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param);
    /** set parameters @c param to KMM component identified by idData to the mixture model */
    void setKmmParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param);

    /** set the parameters of a component given by its ID using defined Mixture manager
     *  @param p_model model with the parameters to get
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     **/
    void setParametersToComponent(IMixtureStatModel* p_model, String const& idData, Rcpp::S4 s4_component);
    /** set the parameters to a component given by its ID
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     **/
    void setParametersToComponent( IMixtureStatModel* p_model
                                 , KernelMixtureManager const& manager
                                 , String const& idData
                                 , Rcpp::S4 s4_component
                                 );
    /** set diagonal Gaussian parameters and data to S4 component
     *  @param p_model model with the parameters to get
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setDiagGaussianParametersToComponent( IMixtureStatModel* p_model
                                             , String const& idData
                                             , Rcpp::S4 s4_component
                                             );
    /** set Poisson parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setPoissonParametersToComponent( IMixtureStatModel* p_model
                                        , String const& idData
                                        , Rcpp::S4 s4_component
                                        );
    /** set gamma parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     *  */
    void setGammaParametersToComponent( IMixtureStatModel* p_model
                                      , String const& idData
                                      , Rcpp::S4 s4_component
                                      );
    /** set categorical parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setCategoricalParametersToComponent( IMixtureStatModel* p_model
                                            , String const& idData
                                            , Rcpp::S4 s4_component
                                            );
    /** set kernel parameters
     *  @param p_model model with the parameters to get
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setKernelParametersToComponent( IMixtureStatModel* p_model
                                       , String const& idData
                                       , Rcpp::S4 s4_component
                                       );

    /** set (estimated) missing values to a R matrix
     *  @param p_model model with the missing values
     *  @param idData ID of the data
     *  @param data data set with missing values to set
     **/
    template< typename Type>
    void setMissingValuesToMatrix( IMixtureStatModel* p_model
                                 , String const& idData
                                 , RMatrix<Type>& data);
    /** set (estimated) missing values by Diagonal Gaussian model to a R matrix
     *  @param p_model model with the missing values
     *  @param idData ID of the data
     *  @param data data set with missing values to set
     **/
    void setDiagGaussianMissingValuesToMatrix( IMixtureStatModel* p_model
                                             , String const& idData
                                             , RMatrix<Real>& data
                                             );
    /** set (estimated) missing values by Poisson model to a R matrix
     *  @param p_model model with the missing values
     *  @param idData ID of the data
     *  @param data data set with missing values to set
     **/
    void setPoissonMissingValuesToMatrix( IMixtureStatModel* p_model
                                        , String const& idData
                                        , RMatrix<Integer>& data
                                        );
    /** set (estimated) missing values by Gamma model to a R matrix
     *  @param p_model model with the missing values
     *  @param idData ID of the data
     *  @param data data set with missing values to set
     **/
    void setGammaMissingValuesToMatrix( IMixtureStatModel* p_model
                                      , String const& idData
                                      , RMatrix<Real>& data
                                      );
    /** set (estimated) missing values by Categorical model to a R matrix
     *  @param p_model model with the missing values
     *  @param idData ID of the data
     *  @param data data set with missing values to set
     **/
    void setCategoricalMissingValuesToMatrix( IMixtureStatModel* p_model
                                            , String const& idData
                                            , RMatrix<Integer>& data
                                            );

    /** get parameters of a component given by its ID
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     **/
    ArrayXX getParameters( String const& idData, Rcpp::S4 s4_component);
    /** get diagonal Gaussian parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getDiagGaussianParameters( String const& idData, Rcpp::S4 s4_component);
    /** get Poisson parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getPoissonParameters( String const& idData, Rcpp::S4 s4_component);
    /** get gamma parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     *  */
    ArrayXX getGammaParameters( String const& idData, Rcpp::S4 s4_component);
    /** get categorical parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getCategoricalParameters( String const& idData, Rcpp::S4 s4_component);
    /** get kernel parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getKernelParameters( String const& idData, Rcpp::S4 s4_component);

};

} // namespace STK

#endif /* STK_ILAUNCHERBASE_H */
