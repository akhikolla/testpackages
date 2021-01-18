/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 15 may 2016
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_ILauncherBase.cpp
 *  @brief In this file we implement the ILauncherBase which
 *  construct properly a mixture model.
 **/


#include "../inst/projects/MixAll/ILauncherBase.h"

using namespace Rcpp;

namespace STK
{
/* facade design pattern.
 * The ILauncherBase allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ILauncherBase::ILauncherBase( Rcpp::S4 model)
                            : IRunnerBase()
                            , s4_model_(model)
                            , handler_()
                            , diagGaussianManager_(handler_)
                            , poissonManager_(handler_)
                            , gammaManager_(handler_)
                            , categoricalManager_(handler_)
                            , kernelManager_(kerHandler_)
{}
/* destructor. */
ILauncherBase::~ILauncherBase() {}

/* set model parameters with param */
bool ILauncherBase::setParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param)
{
    // get idModel from idData
    String idModel;
    if(!handler_.getIdModelName(idData, idModel)) { return false;}
    // get mixture id
    Clust::Mixture mix = Clust::stringToMixture(idModel);
    if (mix == Clust::unknown_mixture_) { return false;}

    // get mixture class (gamma, Gaussian, Poisson, etc.) and parameters for this class of mixture
    switch (Clust::mixtureToMixtureClass(mix))
    {
      case Clust::DiagGaussian_:
        setDiagGaussianParameters(p_model, idData, param);
        return true;
        break;
      case Clust::Poisson_:
        setPoissonParameters(p_model, idData, param);
        return true;
        break;
      case Clust::Gamma_:
        setGammaParameters(p_model, idData, param);
        return true;
        break;
      case Clust::Categorical_:
        setCategoricalParameters(p_model, idData, param);
        return true;
        break;
      case Clust::Kmm_:
        setKmmParameters(p_model, idData, param);
        return true;
        break;
      case Clust::unknown_mixture_class_:
        return false;
        break;
      default:
        return false;
        break;
    }
    return true; // avoid compiler warning
}

// setters
/* set the diagonal Gaussian parameters */
void ILauncherBase::setDiagGaussianParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param)
{ p_model->setParameters(diagGaussianManager_,idData, param);}
/* set the diagonal Gaussian parameters */
void ILauncherBase::setPoissonParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param)
{ p_model->setParameters(poissonManager_,idData, param);}
/* set the gamma parameters */
void ILauncherBase::setGammaParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param)
{ p_model->setParameters(gammaManager_,idData, param);}
/* set the Categorical parameters */
void ILauncherBase::setCategoricalParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param)
{ p_model->setParameters(categoricalManager_,idData, param);}
/* set the Categorical parameters */
void ILauncherBase::setKmmParameters( IMixtureStatModel* p_model, String const& idData, ArrayXX const& param)
{ p_model->setParameters(kernelManager_,idData, param);}

/* fill the s4_component with the parameters */
void ILauncherBase::setParametersToComponent(IMixtureStatModel* p_model, String const& idData, Rcpp::S4 s4_component)
{
  String rModelName = s4_component.slot("modelName");
  switch (Clust::mixtureToMixtureClass(Clust::stringToMixture(rModelName)))
  {
    case Clust::DiagGaussian_:
      setDiagGaussianParametersToComponent(p_model, idData, s4_component);
      break;
    case Clust::Poisson_:
      setPoissonParametersToComponent(p_model, idData, s4_component);
      break;
    case Clust::Gamma_:
      setGammaParametersToComponent(p_model, idData, s4_component);
      break;
    case Clust::Categorical_:
      setCategoricalParametersToComponent(p_model, idData, s4_component);
      break;
    case Clust::Kmm_:
      setKernelParametersToComponent(p_model, idData, s4_component);
      break;
    case Clust::unknown_mixture_class_:
      break;
    default:
      break;
  }
}

/* fill the s4_component with the parameters */
/* get the parameters of a component given by its ID
 *  @param p_model model with the parameters to get
 *  @param manager mixture manager
 *  @param idData ID of the data
 *  @param s4_component component storing the parameters
 **/
void ILauncherBase::setParametersToComponent( IMixtureStatModel* p_model
                                            , KernelMixtureManager const& manager
                                            , String const& idData
                                            , Rcpp::S4 s4_component
                                            )
{
  String rModelName = s4_component.slot("modelName");
  switch (Clust::mixtureToMixtureClass(Clust::stringToMixture(rModelName)))
  {
    case Clust::DiagGaussian_:
      break;
    case Clust::Poisson_:
      break;
    case Clust::Gamma_:
      break;
    case Clust::Categorical_:
      break;
    case Clust::Kmm_:
      setKernelParametersToComponent(p_model, idData, s4_component);
      break;
    case Clust::unknown_mixture_class_:
      break;
    default:
      break;
  }
}


/* get the diagonal Gaussian parameters */
void ILauncherBase::setDiagGaussianParametersToComponent( IMixtureStatModel* p_model
                                                        , String const& idData
                                                        , Rcpp::S4 s4_component)
{
  // get parameters
  ArrayXX params;
  p_model->getParameters(diagGaussianManager_ ,idData, params);
  // get dimensions
  int K = params.sizeRows()/2, nbVariable = params.sizeCols();
  // get results
  ArrayXX mean(K, nbVariable), sigma(K, nbVariable);
  for (int k=0; k<K; ++k)
  {
    mean.row(k)  = params.row(2*k);
    sigma.row(k) = params.row(2*k+1);
  }
  // save results in s4_model
  s4_component.slot("mean")  = Rcpp::wrap(mean);
  s4_component.slot("sigma") = Rcpp::wrap(sigma);
  // get data
  RMatrix<Real> data = (SEXP)s4_component.slot("data");
  setDiagGaussianMissingValuesToMatrix(p_model, idData, data);
    //= Rcpp::wrap(manager.getData(idData));
}

/* get the diagonal Gaussian parameters */
void ILauncherBase::setPoissonParametersToComponent( IMixtureStatModel* p_model
                                                   , String const& idData
                                                   , Rcpp::S4 s4_component)
{
  // get parameters
  ArrayXX params;
  p_model->getParameters(poissonManager_, idData, params);
  // save results in s4_model
  s4_component.slot("lambda")  = Rcpp::wrap(params);
  // get data
  RMatrix<Integer> data = (SEXP)s4_component.slot("data");
  setPoissonMissingValuesToMatrix(p_model, idData, data);
  //s4_component.slot("data") =Rcpp::wrap( manager.getData(idData));
}

/* get the gamma parameters */
void ILauncherBase::setGammaParametersToComponent( IMixtureStatModel* p_model
                                                 , String const& idData
                                                 , Rcpp::S4 s4_component)
{
  // get parameters
  ArrayXX params;
  p_model->getParameters(gammaManager_, idData, params);
  // get dimensions
  int K = params.sizeRows()/2, nbVariable = params.sizeCols();
  // get results
  ArrayXX shape(K, nbVariable), scale(K, nbVariable);
  for (int k=0; k<K; ++k)
  {
    shape.row(k) = params.row(2*k);
    scale.row(k) = params.row(2*k+1);
  }
  // save results in s4_model
  s4_component.slot("shape") = Rcpp::wrap(shape);
  s4_component.slot("scale") = Rcpp::wrap(scale);
  // get data
  RMatrix<Real> data = (SEXP)s4_component.slot("data");
  setGammaMissingValuesToMatrix(p_model, idData, data);
  //s4_component.slot("data") = Rcpp::wrap(manager.getData(idData));
}

/* get the diagonal Categorical parameters */
void ILauncherBase::setCategoricalParametersToComponent( IMixtureStatModel* p_model
                                                       , String const& idData
                                                       , Rcpp::S4 s4_component)
{
  // get parameters
  ArrayXX params;
  p_model->getParameters(categoricalManager_, idData, params);
  params.shift(0,0);
  // save results in s4_model
  s4_component.slot("plkj") = Rcpp::wrap(params);
  // get data
  RMatrix<Integer> data = (SEXP)s4_component.slot("data");
  setCategoricalMissingValuesToMatrix(p_model, idData, data);
}

/* get the kernel parameters */
void ILauncherBase::setKernelParametersToComponent( IMixtureStatModel* p_model
                                                  , String const& idData
                                                  , Rcpp::S4 s4_component)
{
  // get parameters
  ArrayXX param;
  p_model->getParameters(kernelManager_, idData, param);
  // save results in s4_model
  s4_component.slot("sigma2") = Rcpp::wrap(param.col(0));
  s4_component.slot("dim")    = Rcpp::wrap(param.col(1));
}


/* get the parameters of a component given by its ID
 *  @param s4_component the component with the parameters to get
 *  @param idData ID of the data
 **/
 ArrayXX ILauncherBase::getParameters( String const& idData, Rcpp::S4 s4_component)
{
  String rModelName = s4_component.slot("modelName");
  switch (Clust::mixtureToMixtureClass(Clust::stringToMixture(rModelName)))
  {
    case Clust::DiagGaussian_:
      return getDiagGaussianParameters(idData, s4_component);
      break;
    case Clust::Poisson_:
      return getPoissonParameters(idData, s4_component);
      break;
    case Clust::Gamma_:
      return getGammaParameters(idData, s4_component);
      break;
    case Clust::Categorical_:
      return getCategoricalParameters(idData, s4_component);
      break;
    case Clust::Kmm_:
      return getKernelParameters(idData, s4_component);
      break;
    case Clust::unknown_mixture_class_:
      break;
    default:
      break;
  }
  return ArrayXX();
}


/* utility function allowing to set parameters when model is DiagGaussian */
ArrayXX ILauncherBase::getDiagGaussianParameters( String const& idData, Rcpp::S4 s4_component)
{
  // get mean and sigma
  RMatrix<double> mean((SEXP)s4_component.slot("mean"));
  RMatrix<double> sigma((SEXP)s4_component.slot("sigma"));
  int nbCluster  = s4_model_.slot("nbCluster");
  ArrayXX params(nbCluster*2, mean.cols());
  // get results
  for (int k=0; k<nbCluster; ++k)
  {
    params.row(2*k)   = mean.row(k);
    params.row(2*k+1) = sigma.row(k);
  }
  return params;
}

/* utility function allowing to set parameters when model is DiagGaussian */
ArrayXX ILauncherBase::getGammaParameters( String const& idData, Rcpp::S4 s4_component)
{
  // get mean and sigma
  RMatrix<double> shape((SEXP)s4_component.slot("shape"));
  RMatrix<double> scale((SEXP)s4_component.slot("scale"));
  int nbCluster  = s4_model_.slot("nbCluster");
  ArrayXX params(nbCluster*2, shape.cols());
  // get results
  for (int k=0; k<nbCluster; ++k)
  {
    params.row(2*k)   = shape.row(k);
    params.row(2*k+1) = scale.row(k);
  }
  return params;
}

/* utility function allowing to set parameters when model is DiagGaussian */
ArrayXX ILauncherBase::getCategoricalParameters( String const& idData, Rcpp::S4 s4_component)
{
  // get plkj array
  RMatrix<double> plkj((SEXP)s4_component.slot("plkj"));
  ArrayXX params = plkj;
  return params;
}

/* utility function allowing to set parameters when model is Poisson */
ArrayXX ILauncherBase::getPoissonParameters( String const& idData, Rcpp::S4 s4_component)
{
  // get mean and sigma
  RMatrix<double> lambda((SEXP)s4_component.slot("lambda"));
  ArrayXX params = lambda;
  return params;
}

/* utility function allowing to set parameters when model is DiagGaussian */
ArrayXX ILauncherBase::getKernelParameters( String const& idData, Rcpp::S4 s4_component)
{
  // get mean and sigma
  RMatrix<double> sigma2((SEXP)s4_component.slot("sigma2"));
  RMatrix<double> dim((SEXP)s4_component.slot("dim"));
  int nbCluster  = s4_model_.slot("nbCluster");
  ArrayXX params(nbCluster*2, sigma2.cols());
  // get results
  for (int k=0; k<nbCluster; ++k)
  {
    params.row(2*k)   = sigma2.row(k);
    params.row(2*k+1) = dim.row(k);
  }
  return params;
}

/* set (estimated) missing values to a R matrix
 *  @param idData ID of the data
 *  @param data data set with missing values to set
 **/
void ILauncherBase::setDiagGaussianMissingValuesToMatrix( IMixtureStatModel* p_model
                                                        , String const& idData
                                                        , RMatrix<Real>& data
                                                        )
{
  typedef typename hidden::MixtureManagerTraits< DiagGaussianMixtureManager<RDataHandler> >::MissingValues Missingvalues;
  Missingvalues missing;
  p_model->getMissingValues(diagGaussianManager_, idData, missing);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ILauncherBase::setDiagGaussianMissingValuesToMatrix. missing loaded\n");
#endif
  for(size_t m = 0; m< missing.size(); ++m)
  {
    int i = missing[m].first.first, j = missing[m].first.second;
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("i=") << i << _T(", j=") << j << _T(", value=") << missing[m].second << _T("\n");
#endif
    data(i, j) = missing[m].second;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ILauncherBase::setDiagGaussianMissingValuesToMatrix. missing values set\n");
#endif
}

/* set (estimated) missing values to a R matrix
 *  @param idData ID of the data
 *  @param data data set with missing values to set
 **/
void ILauncherBase::setPoissonMissingValuesToMatrix( IMixtureStatModel* p_model
                                                   , String const& idData
                                                   , RMatrix<Integer>& data
                                                   )
{
  typedef typename hidden::MixtureManagerTraits< PoissonMixtureManager<RDataHandler> >::MissingValues Missingvalues;
  Missingvalues missing;
  p_model->getMissingValues(poissonManager_, idData, missing);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ILauncherBase::setPoissonMissingValuesToMatrix. missing loaded\n");
#endif
  for(size_t m = 0; m< missing.size(); ++m)
  {
    int i = missing[m].first.first, j = missing[m].first.second;
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("i=") << i << _T(", j=") << j << _T(", value=") << missing[m].second << _T("\n");
#endif
    data(i, j) = missing[m].second;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ILauncherBase::setPoissonMissingValuesToMatrix. missing values set\n");
#endif
}

/* set (estimated) missing values to a R matrix
 *  @param idData ID of the data
 *  @param data data set with missing values to set
 **/
void ILauncherBase::setGammaMissingValuesToMatrix( IMixtureStatModel* p_model
                                                 , String const& idData
                                                 , RMatrix<Real>& data
                                                 )
{
  typedef typename hidden::MixtureManagerTraits< GammaMixtureManager<RDataHandler> >::MissingValues Missingvalues;
  Missingvalues missing;
  p_model->getMissingValues(gammaManager_, idData, missing);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ILauncherBase::setGammaMissingValuesToMatrix. missing loaded\n");
#endif
  for(size_t m = 0; m< missing.size(); ++m)
  {
    int i = missing[m].first.first, j = missing[m].first.second;
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("i=") << i << _T(", j=") << j << _T(", value=") << missing[m].second << _T("\n");
#endif
    data(i, j) = missing[m].second;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ILauncherBase::setGammaMissingValuesToMatrix. missing values set\n");
#endif
}

/* set (estimated) missing values to a R matrix
 *  @param idData ID of the data
 *  @param data data set with missing values to set
 **/
void ILauncherBase::setCategoricalMissingValuesToMatrix( IMixtureStatModel* p_model
                                                       , String const& idData
                                                       , RMatrix<Integer>& data
                                                       )
{
  typedef typename hidden::MixtureManagerTraits< CategoricalMixtureManager<RDataHandler> >::MissingValues Missingvalues;
  Missingvalues missing;
  p_model->getMissingValues(categoricalManager_, idData, missing);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ILauncherBase::setCategoricalMissingValuesToMatrix. missing loaded\n");
#endif
  for(size_t m = 0; m< missing.size(); ++m)
  {
    int i = missing[m].first.first, j = missing[m].first.second;
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("i=") << i << _T(", j=") << j << _T(", value=") << missing[m].second << _T("\n");
#endif
    data(i, j) = missing[m].second;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ILauncherBase::setCategoricalMissingValuesToMatrix. missing values set\n");
#endif
}

} // namespace STK

