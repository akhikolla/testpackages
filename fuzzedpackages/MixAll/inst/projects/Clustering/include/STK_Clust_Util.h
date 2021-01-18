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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 2 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_Clust_Util.h
 *  @brief In this file we define the enum, constants and utilities functions
 *  of the Clustering project.
 **/


#ifndef STK_CLUST_UTIL_H
#define STK_CLUST_UTIL_H

#include <STKernel.h>

namespace STK
{

// forward declaration
class IMixtureAlgo;
class IMixtureAlgoPredict;
class IMixtureAlgoLearn;
class IMixtureInit;
class IMixtureStrategy;
class IMixtureComposer;
class IMixtureCriterion;

/** @ingroup Clustering
 *  @brief struct storing the parameters of the mixture.
 *  Parameters of a mixture model have an unique Id defined
 *  in STK::Clust::Mixture enumeration.
 **/
template <int Id> struct ModelParameters;

/** @ingroup Clustering
 *  @brief struct storing the parameters of the matrix valued mixture.
 *  Parameters of a matrix mixture model have two Id defined
 *  in STK::Clust::HDCovarianceModel enumeration.
 **/
template <class Array_> struct HDMatrixModelParameters;

namespace hidden
{
/** @ingroup hidden
 *  MixtureBridgeTraits struct for bridged mixtures
 *  The traits struct MixtureBridgeTraits must be specialized for any
 *  Mixture deriving from the Interface IMixtureBridge.
 **/
template<class Derived> struct MixtureBridgeTraits;
/**  @ingroup hidden
 *  Main class for the mixtures traits policy.
 *  The traits struct MixtureTraits must be specialized for any
 *  Mixture deriving from the Interface IMixtureDensity.
 **/
template <class Mixture> struct MixtureTraits;

/** @ingroup hidden
 *  Main class for the mixture managers traits policy.
 *  The traits struct MixtureManagerTraits must be specialized for any
 *  STK::Clust::MixtureClass manager.
 *  @sa CategoricalMixtureManager, DiagGaussianMixtureManager
 **/
template <class Manager> struct MixtureManagerTraits;

} // namespace hidden


namespace Clust
{

/** @ingroup Clustering
 *  @brief initialization type.
 *  There is four ways to initialize the mixture model:
 *  - using random values for the parameters
 *  - using random class for the sampling
 *  - using random probabilities class for the sampling
 *  - using parameters values
 **/
enum initType
{
  noInit_ = -1,         ///< no initialization
  randomInit_ = -2,     ///< DEPRECATED
  randomParamInit_ = 0, ///< initialize randomly the parameters
  randomClassInit_ = 1, ///< initialize randomly the class labels
  randomFuzzyInit_ = 2, ///< initialize randomly the partnership class probabilities
  valueParamInit_ = 3   ///< initialize parameters using given values
};

/** @ingroup Clustering
 *  Convert a String to a initType. The recognized strings are
 * <table>
 * <tr> <th> Initialization   </th></tr>
 * <tr> <td> "randomInit" (DEPECATED) </td></tr>
 * <tr> <td> "randomParamInit"</td></tr>
 * <tr> <td> "randomClassInit"</td></tr>
 * <tr> <td> "randomFuzzyInit"</td></tr>
 * <tr> <td> "random"         </td></tr>
 * <tr> <td> "class"          </td></tr>
 * <tr> <td> "fuzzy"          </td></tr>
 * </table>
 *  @param type the type of initialization wanted
 *  @return the initType corresponding (default is randomClassInit)
 *  @note if the string is not found in the list above,the type Clust::randomClassInit_
 *  is returned.
 **/
initType stringToInit( String const& type);

/** @ingroup Clustering
 *  Estimation algorithms
 **/
enum algoType
{
  emAlgo_ = 0,
  cemAlgo_ = 1,
  semAlgo_ = 2,
  semiSemAlgo_ = 3
};


/** @ingroup Clustering
 *  Convert a String to an algoType. The recognized strings are
 * <table>
 * <tr> <th> Algorithm     </th></tr>
 * <tr> <td> "emAlgo"      </td></tr>
 * <tr> <td> "cemAlgo"     </td></tr>
 * <tr> <td> "semAlgo"     </td></tr>
 * <tr> <td> "semiSemAlgo" </td></tr>
 * <tr> <td> "em"          </td></tr>
 * <tr> <td> "cem"         </td></tr>
 * <tr> <td> "sem"         </td></tr>
 * <tr> <td> "semiSem"     </td></tr>
 * </table>
 *  @param type the type of algorithm wanted
 *  @return the algoType corresponding (default is emAlgo)
 *  @note The capitalized letters have no effect and if the string is not found
 *  in the list above, the type Clust::emAlgo_ is returned.
 **/
algoType stringToAlgo( String const& type);

/** @ingroup Clustering
 *  Learning estimation algorithms
 **/
enum algoPredictType
{
  emPredictAlgo_,
  semiSEMPredictAlgo_
};

/** @ingroup Clustering
 *  Convert a String to an algoPredictType. The recognized strings are
 * <table>
 * <tr> <th> Algorithm     </th></tr>
 * <tr> <td> "em"          </td></tr>
 * <tr> <td> "semiSem"         </td></tr>
 * </table>
 *  @param type the type of algorithm wanted
 *  @return the algoPredictType corresponding (default is em)
 *  @note The capitalized letters have no effect and if the string is not found
 *  in the list above, the type Clust::emPredictAlgo_ is returned.
 **/
algoPredictType stringToPredictAlgo( String const& type);

/** @ingroup Clustering
 *  Learning estimation algorithms
 **/
enum algoLearnType
{
  imputeAlgo_,
  simulAlgo_
};

/** @ingroup Clustering
 *  Convert a String to an algoLearnType. The recognized strings are
 * <table>
 * <tr> <th> Algorithm     </th></tr>
 * <tr> <td> "imputeAlgo"  </td></tr>
 * <tr> <td> "simulAlgo"   </td></tr>
 * <tr> <td> "impute"      </td></tr>
 * <tr> <td> "simul"       </td></tr>
 * </table>
 *  @param type the type of algorithm wanted
 *  @return the algoType corresponding (default is emAlgo)
 *  @note The capitalized letters have no effect and if the string is not found
 *  in the list above,the type Clust::emAlgo_ is returned.
 **/
algoLearnType stringToLearnAlgo( String const& type);

/** @ingroup Clustering
 *  strategy of estimation
 **/
enum strategyType
{
  simpleStrategy_ = 0,
  XemStrategy_    = 1, // not implemented
  SemStrategy_    = 2,
  FullStrategy_   = 3
};

/** @ingroup Clustering
 *  type of criterion to use in order to select the mixture model
 **/
enum criterionType
{
  aic_  = 0,
  bic_  = 1,
  icl_  = 2,
  ml_   = 3
};

/** @ingroup Clustering
 *  Convert a String to an criterionType. The recognized strings are
 * <table>
 * <tr> <th> Criterion </th></tr>
 * <tr> <td> "AIC"     </td></tr>
 * <tr> <td> "BIC"     </td></tr>
 * <tr> <td> "ICL      </td></tr>
 * <tr> <td> "ML"      </td></tr>
 * </table>
 *  @param type the type of criterion wanted
 *  @return the criterionType corresponding
 *  @note The capitalized letters have no effect and if the string is not found
 *  in the list above,the type Clust::bic_ is returned.
 **/
criterionType stringToCriterion( String const& type);

/** @ingroup Clustering
 *  Specific exceptions allowing to handle the erroros that can occur in the
 *  estimation process.
 **/
enum exceptions
{
  randomInitFail_,
  randomParamInitFail_,
  randomClassInitFail_,
  randomFuzzyInitFail_,
  estimFail_,
  initializeStepFail_,
  mStepFail_,
  eStepFail_,
  mapStepFail_,
  cStepFail_,
  sStepFail_
};

/** @ingroup Clustering
 *  convert a Clust::exceptions to a String.
 *  @param type the type of exception that occur
 *  @return the string associated to this exception.
 **/
String exceptionToString( exceptions const& type);

/** @ingroup Clustering
 *  Give the state of the model.
 **/
enum modelState
{
  modelCreated_ =0,         ///< the model has been created but is not initialized
  modelInitialized_ =1,     ///< the model is initialized and its parameters are initialized to default values
  modelParamInitialized_=2, ///< The parameters of the model have been initialized
  shortRun_,                ///< A short run has been done
  longRun_,                 ///< A long run has been done
  modelFinalized_           ///< the model is finalized
};

/** @ingroup Clustering
 *  list of the parsimonious covariance models that can be used
 **/
enum ParsimoniousCovarianceModel
{
  Covariance_EII_ =100,
  Covariance_VII_,
  Covariance_EEI_,
  Covariance_VEI_,
  Covariance_EVI_,
  Covariance_VVI_,
  Covariance_EEE_,
  Covariance_VEE_,
  Covariance_EVE_,
  Covariance_VVE_,
  Covariance_EEV_,
  Covariance_VEV_,
  Covariance_EVV_,
  Covariance_VVV_  // =114
};
/** @ingroup Clustering
 *  list of the HD covariance models that can be used
 **/
enum HDCovarianceModel
{
  HDCovariance_AjkBkQkDk_ =120,
  HDCovariance_AjkBkQkD_,
  HDCovariance_AjkBkQDk_,
  HDCovariance_AjkBkQD_,
  HDCovariance_AjkBQkDk_,
  HDCovariance_AjkBQkD_,
  HDCovariance_AjkBQDk_,
  HDCovariance_AjkBQD_,
  HDCovariance_AkBkQkDk_,
  HDCovariance_AkBkQkD_,
  HDCovariance_AkBkQDk_,
  HDCovariance_AkBkQD_,
  HDCovariance_AkBQkDk_,
  HDCovariance_AkBQkD_,
  HDCovariance_AkBQDk_,
  HDCovariance_AkBQD_,
  HDCovariance_AjBkQkDk_,
  HDCovariance_AjBkQkD_,
  HDCovariance_AjBkQDk_,
  HDCovariance_AjBkQD_,
  HDCovariance_AjBQkDk_,
  HDCovariance_AjBQkD_,
  HDCovariance_AjBQDk_,
  HDCovariance_AjBQD_,
  HDCovariance_ABkQkDk_,
  HDCovariance_ABkQkD_,
  HDCovariance_ABkQDk_,
  HDCovariance_ABkQD_,
  HDCovariance_ABQkDk_,
  HDCovariance_ABQkD_,
  HDCovariance_ABQDk_,
  HDCovariance_ABQD_ // =151
};

} // namespace Clust

namespace hidden
{
template<int Id> struct HDCovarianceChooser;

template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjkBkQkDk_>
{
  const bool isAj_ = true;
  const bool isAk_ = true;
  const bool isBk_ = true;
  const bool isQk_ = true;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjkBkQkD_>
{
  const bool isAj_ = true;
  const bool isAk_ = true;
  const bool isBk_ = true;
  const bool isQk_ = true;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjkBkQDk_>
{
  const bool isAj_ = true;
  const bool isAk_ = true;
  const bool isBk_ = true;
  const bool isQk_ = false;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjkBkQD_>
{
  const bool isAj_ = true;
  const bool isAk_ = true;
  const bool isBk_ = true;
  const bool isQk_ = false;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjkBQkDk_>
{
  const bool isAj_ = true;
  const bool isAk_ = true;
  const bool isBk_ = false;
  const bool isQk_ = true;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjkBQkD_>
{
  const bool isAj_ = true;
  const bool isAk_ = true;
  const bool isBk_ = false;
  const bool isQk_ = true;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjkBQDk_>
{
  const bool isAj_ = true;
  const bool isAk_ = true;
  const bool isBk_ = false;
  const bool isQk_ = false;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjkBQD_>
{
  const bool isAj_ = true;
  const bool isAk_ = true;
  const bool isBk_ = false;
  const bool isQk_ = false;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AkBkQkDk_>
{
  const bool isAj_ = false;
  const bool isAk_ = true;
  const bool isBk_ = true;
  const bool isQk_ = true;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AkBkQkD_>
{
  const bool isAj_ = false;
  const bool isAk_ = true;
  const bool isBk_ = true;
  const bool isQk_ = true;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AkBkQDk_>
{
  const bool isAj_ = false;
  const bool isAk_ = true;
  const bool isBk_ = true;
  const bool isQk_ = false;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AkBkQD_>
{
  const bool isAj_ = false;
  const bool isAk_ = true;
  const bool isBk_ = true;
  const bool isQk_ = false;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AkBQkDk_>
{
  const bool isAj_ = false;
  const bool isAk_ = true;
  const bool isBk_ = false;
  const bool isQk_ = true;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AkBQkD_>
{
  const bool isAj_ = false;
  const bool isAk_ = true;
  const bool isBk_ = false;
  const bool isQk_ = true;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AkBQDk_>
{
  const bool isAj_ = false;
  const bool isAk_ = true;
  const bool isBk_ = false;
  const bool isQk_ = false;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AkBQD_>
{
  const bool isAj_ = false;
  const bool isAk_ = true;
  const bool isBk_ = false;
  const bool isQk_ = false;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjBkQkDk_>
{
  const bool isAj_ = true;
  const bool isAk_ = false;
  const bool isBk_ = true;
  const bool isQk_ = true;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjBkQkD_>
{
  const bool isAj_ = true;
  const bool isAk_ = false;
  const bool isBk_ = true;
  const bool isQk_ = true;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjBkQDk_>
{
  const bool isAj_ = true;
  const bool isAk_ = false;
  const bool isBk_ = true;
  const bool isQk_ = false;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjBkQD_>
{
  const bool isAj_ = true;
  const bool isAk_ = false;
  const bool isBk_ = true;
  const bool isQk_ = false;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjBQkDk_>
{
  const bool isAj_ = true;
  const bool isAk_ = false;
  const bool isBk_ = false;
  const bool isQk_ = true;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_AjBQkD_>
{
  const bool isAj_ = true;
  const bool isAk_ = false;
  const bool isBk_ = false;
  const bool isQk_ = true;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_ABkQkDk_>
{
  const bool isAj_ = false;
  const bool isAk_ = false;
  const bool isBk_ = true;
  const bool isQk_ = true;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_ABkQkD_>
{
  const bool isAj_ = false;
  const bool isAk_ = false;
  const bool isBk_ = true;
  const bool isQk_ = true;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_ABkQDk_>
{
  const bool isAj_ = false;
  const bool isAk_ = false;
  const bool isBk_ = true;
  const bool isQk_ = false;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_ABkQD_>
{
  const bool isAj_ = false;
  const bool isAk_ = false;
  const bool isBk_ = true;
  const bool isQk_ = false;
  const bool isDk_ = false;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_ABQkDk_>
{
  const bool isAj_ = false;
  const bool isAk_ = false;
  const bool isBk_ = false;
  const bool isQk_ = true;
  const bool isDk_ = true;
};
template<>
struct HDCovarianceChooser<Clust::HDCovariance_ABQkD_>
{
  const bool isAj_ = false;
  const bool isAk_ = false;
  const bool isBk_ = false;
  const bool isQk_ = true;
  const bool isDk_ = false;
};

}  // namespace hidden

namespace Clust
{
/** @ingroup Clustering
 * list of the mixtures that can be used by the composer
 **/
enum Mixture
{
  Gamma_ajk_bjk_ =0,
  Gamma_ajk_bk_,
  Gamma_ajk_bj_,
  Gamma_ajk_b_,
  Gamma_ak_bjk_,
  Gamma_ak_bk_,
  Gamma_ak_bj_,
  Gamma_ak_b_,
  Gamma_aj_bjk_,
  Gamma_aj_bk_,
  Gamma_a_bjk_,
  Gamma_a_bk_, // = 11
  Gaussian_sjk_ =20,
  Gaussian_sk_,
  Gaussian_sj_,
  Gaussian_s_,
  Gaussian_sjsk_, // = 24
  Categorical_pjk_ =40,
  Categorical_pk_,
  Poisson_ljk_ = 60,
  Poisson_lk_,
  Poisson_ljlk_,
  Kmm_sk_ = 80,
  Kmm_s_,
  HDGaussian_AjkBkQkDk_ =120,
  HDGaussian_AjkBkQkD_,
  HDGaussian_AjkBkQDk_,
  HDGaussian_AjkBkQD_,
  HDGaussian_AjkBQkDk_,
  HDGaussian_AjkBQkD_,
  HDGaussian_AjkBQDk_,
  HDGaussian_AjkBQD_,
  HDGaussian_AkBkQkDk_,
  HDGaussian_AkBkQkD_,
  HDGaussian_AkBkQDk_,
  HDGaussian_AkBkQD_,
  HDGaussian_AkBQkDk_,
  HDGaussian_AkBQkD_,
  HDGaussian_AkBQDk_,
  HDGaussian_AkBQD_,
  HDGaussian_AjBkQkDk_,
  HDGaussian_AjBkQkD_,
  HDGaussian_AjBkQDk_,
  HDGaussian_AjBkQD_,
  HDGaussian_AjBQkDk_,
  HDGaussian_AjBQkD_,
  HDGaussian_ABkQkDk_,
  HDGaussian_ABkQkD_,
  HDGaussian_ABkQDk_,
  HDGaussian_ABkQD_,
  HDGaussian_ABQkDk_,
  HDGaussian_ABQkD_,
  HDGaussian_ABQD_,   // =148
  unknown_mixture_ = -1
};

/** @ingroup Clustering
 *  list of the class of mixture implemented in stkpp
 **/
enum MixtureClass
{
  Gamma_,
  DiagGaussian_,
  Categorical_,
  Poisson_,
  Kmm_,
  Matrix_,
  HDGaussian_,
  HDMatrixGaussian_,
  unknown_mixture_class_ = -1
};

/** @ingroup Clustering
 *  Convert a String to a Mixture. The recognized strings are
 * <table >
 * <tr> <th> Model             </th> </tr>
 * <tr> <td> "Gamma_Ajkbjk"   </td></tr>
 * <tr> <td> "Gamma_Ajkbk"    </td></tr>
 * <tr> <td> "Gamma_Ajkbj"    </td></tr>
 * <tr> <td> "Gamma_Ajkb"     </td></tr>
 * <tr> <td> "Gamma_ak_bjk"    </td></tr>
 * <tr> <td> "Gamma_ak_bk"     </td></tr>
 * <tr> <td> "Gamma_ak_bj"     </td></tr>
 * <tr> <td> "Gamma_ak_b"      </td></tr>
 * <tr> <td> "Gamma_aj_bjk"    </td></tr>
 * <tr> <td> "Gamma_aj_bk"     </td></tr>
 * <tr> <td> "Gamma_a_bjk"     </td></tr>
 * <tr> <td> "Gamma_a_bk"      </td></tr>
 * <tr> <td> "Gaussian_sjk"    </td></tr>
 * <tr> <td> "Gaussian_sk"     </td></tr>
 * <tr> <td> "Gaussian_sj"     </td></tr>
 * <tr> <td> "Gaussian_s"      </td></tr>
 * <tr> <td> "Gaussian_sjsk"   </td></tr>
 * <tr> <td> "HDGaussian_AjkBkQkDk" </td></tr>
 * <tr> <td> "HDGaussian_AjkBkQkD"  </td></tr>
 * <tr> <td> "HDGaussian_AjkBkQdk"  </td></tr>
 * <tr> <td> "HDGaussian_AjkBkQd"   </td></tr>
 * <tr> <td> "HDGaussian_AjkBQkDk"  </td></tr>
 * <tr> <td> "HDGaussian_AjkBQkD"   </td></tr>
 * <tr> <td> "HDGaussian_AjkBQdk"   </td></tr>
 * <tr> <td> "HDGaussian_AjkBQd"    </td></tr>
 * <tr> <td> "HDGaussian_AkBkQkDk"  </td></tr>
 * <tr> <td> "HDGaussian_AkBkQkD"   </td></tr>
 * <tr> <td> "HDGaussian_AkBkQdk"   </td></tr>
 * <tr> <td> "HDGaussian_AkBkQd"    </td></tr>
 * <tr> <td> "HDGaussian_AkBQkDk"   </td></tr>
 * <tr> <td> "HDGaussian_AkBQkD"    </td></tr>
 * <tr> <td> "HDGaussian_AkBQdk"    </td></tr>
 * <tr> <td> "HDGaussian_AkBQd"     </td></tr>
 * <tr> <td> "HDGaussian_AjBkQkDk"  </td></tr>
 * <tr> <td> "HDGaussian_AjBkQkD"   </td></tr>
 * <tr> <td> "HDGaussian_AjBkQdk"   </td></tr>
 * <tr> <td> "HDGaussian_AjBkQd"    </td></tr>
 * <tr> <td> "HDGaussian_AjBQkDk"   </td></tr>
 * <tr> <td> "HDGaussian_AjBQkD"    </td></tr>
 * <tr> <td> "HDGaussian_ABkQkDk"   </td></tr>
 * <tr> <td> "HDGaussian_ABkQkD"    </td></tr>
 * <tr> <td> "HDGaussian_ABkQdk"    </td></tr>
 * <tr> <td> "HDGaussian_ABkQd,"    </td></tr>
 * <tr> <td> "HDGaussian_ABQkDk"    </td></tr>
 * <tr> <td> "HDGaussian_ABQkD"     </td></tr>
 * <tr> <td> "HDGaussian_ABQD"      </td></tr>
 * <tr> <td> "Categorical_pjk"      </td></tr>
 * <tr> <td> "Categorical_pk"       </td></tr>
 * <tr> <td> "Poisson_ljk"          </td></tr>
 * <tr> <td> "Poisson_lk"           </td></tr>
 * <tr> <td> "Poisson_ljlk"         </td></tr>
 * <tr> <td> "Kmm_sk"               </td></tr>
 * <tr> <td> "Kmm_s"                </td></tr>
 * </table>
 *  @param type the String we want to convert
 *  @return the Mixture represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_mixture_ type is returned.
 **/
Mixture stringToMixture( String const& type);

/** @ingroup Clustering
 *  convert a string to a Mixture and specify if the model is with free proportions
 *  or fixed proportions. The recognized strings are
 * <table border >
 * <tr> <th> Free proportions     </th><th> Fixed Proportions   </th> </tr>
 * <tr> <td> "Gamma_pk_ajk_bjk"   </td><td> "Gamma_p_ajk_bjk"   </td> </tr>
 * <tr> <td> "Gamma_pk_ajk_bk"    </td><td> "Gamma_p_ajk_bk"    </td> </tr>
 * <tr> <td> "Gamma_pk_ajk_bj"    </td><td> "Gamma_p_ajk_bj"    </td> </tr>
 * <tr> <td> "Gamma_pk_ajk_b"     </td><td> "Gamma_p_ajk_b"     </td> </tr>
 * <tr> <td> "Gamma_pk_ak_bjk"    </td><td> "Gamma_p_ak_bjk"    </td> </tr>
 * <tr> <td> "Gamma_pk_ak_bk"     </td><td> "Gamma_p_ak_bk"     </td> </tr>
 * <tr> <td> "Gamma_pk_ak_bj"     </td><td> "Gamma_p_ak_bj"     </td> </tr>
 * <tr> <td> "Gamma_pk_ak_b"      </td><td> "Gamma_p_ak_b"      </td> </tr>
 * <tr> <td> "Gamma_pk_aj_bjk"    </td><td> "Gamma_p_aj_bjk"    </td> </tr>
 * <tr> <td> "Gamma_pk_aj_bk"     </td><td> "Gamma_p_aj_bk"     </td> </tr>
 * <tr> <td> "Gamma_pk_a_bjk"     </td><td> "Gamma_p_a_bjk"     </td> </tr>
 * <tr> <td> "Gamma_pk_a_bk"      </td><td> "Gamma_p_a_bk"      </td> </tr>
 * <tr> <td> "Gaussian_pk_sjk"    </td><td> "Gaussian_p_sjk"    </td> </tr>
 * <tr> <td> "Gaussian_pk_sk"     </td><td> "Gaussian_p_sk"     </td> </tr>
 * <tr> <td> "Gaussian_pk_sj"     </td><td> "Gaussian_p_sj"     </td> </tr>
 * <tr> <td> "Gaussian_pk_s"      </td><td> "Gaussian_p_s"      </td> </tr>
 * <tr> <td> "Gaussian_pk_sjsk"   </td><td> "Gaussian_p_sjsk"   </td> </tr>
 * <tr> <td> "HDGaussian_pk_AjkBkQkDk" </td><td> "HDGaussian_p_AjkBkQkDk" </td></tr>
 * <tr> <td> "HDGaussian_pk_AjkBkQkD"  </td><td> "HDGaussian_p_AjkBkQkD"  </td></tr>
 * <tr> <td> "HDGaussian_pk_AjkBkQdk"  </td><td> "HDGaussian_p_AjkBkQdk"  </td></tr>
 * <tr> <td> "HDGaussian_pk_AjkBkQd"   </td><td> "HDGaussian_p_AjkBkQd"   </td></tr>
 * <tr> <td> "HDGaussian_pk_AjkBQkDk"  </td><td> "HDGaussian_p_AjkBQkDk"  </td></tr>
 * <tr> <td> "HDGaussian_pk_AjkBQkD"   </td><td> "HDGaussian_p_AjkBQkD"   </td></tr>
 * <tr> <td> "HDGaussian_pk_AjkBQdk"   </td><td> "HDGaussian_p_AjkBQdk"   </td></tr>
 * <tr> <td> "HDGaussian_pk_AjkBQd"    </td><td> "HDGaussian_p_AjkBQd"    </td></tr>
 * <tr> <td> "HDGaussian_pk_AkBkQkDk"  </td><td> "HDGaussian_p_AkBkQkDk"  </td></tr>
 * <tr> <td> "HDGaussian_pk_AkBkQkD"   </td><td> "HDGaussian_p_AkBkQkD"   </td></tr>
 * <tr> <td> "HDGaussian_pk_AkBkQdk"   </td><td> "HDGaussian_p_AkBkQdk"   </td></tr>
 * <tr> <td> "HDGaussian_pk_AkBkQd"    </td><td> "HDGaussian_p_AkBkQd"    </td></tr>
 * <tr> <td> "HDGaussian_pk_AkBQkDk"   </td><td> "HDGaussian_p_AkBQkDk"   </td></tr>
 * <tr> <td> "HDGaussian_pk_AkBQkD"    </td><td> "HDGaussian_p_AkBQkD"    </td></tr>
 * <tr> <td> "HDGaussian_pk_AkBQdk"    </td><td> "HDGaussian_p_AkBQdk"    </td></tr>
 * <tr> <td> "HDGaussian_pk_AkBQd"     </td><td> "HDGaussian_p_AkBQd"     </td></tr>
 * <tr> <td> "HDGaussian_pk_AjBkQkDk"  </td><td> "HDGaussian_p_AjBkQkDk"  </td></tr>
 * <tr> <td> "HDGaussian_pk_AjBkQkD"   </td><td> "HDGaussian_p_AjBkQkD"   </td></tr>
 * <tr> <td> "HDGaussian_pk_AjBkQdk"   </td><td> "HDGaussian_p_AjBkQdk"   </td></tr>
 * <tr> <td> "HDGaussian_pk_AjBkQd"    </td><td> "HDGaussian_p_AjBkQd"    </td></tr>
 * <tr> <td> "HDGaussian_pk_AjBQkDk"   </td><td> "HDGaussian_p_AjBQkDk"   </td></tr>
 * <tr> <td> "HDGaussian_pk_AjBQkD"    </td><td> "HDGaussian_p_AjBQkD"    </td></tr>
 * <tr> <td> "HDGaussian_pk_ABkQkDk"   </td><td> "HDGaussian_p_ABkQkDk"   </td></tr>
 * <tr> <td> "HDGaussian_pk_ABkQkD"    </td><td> "HDGaussian_p_ABkQkD"    </td></tr>
 * <tr> <td> "HDGaussian_pk_ABkQdk"    </td><td> "HDGaussian_p_ABkQdk"    </td></tr>
 * <tr> <td> "HDGaussian_pk_ABkQd,"    </td><td> "HDGaussian_p_ABkQd,"    </td></tr>
 * <tr> <td> "HDGaussian_pk_ABQkDk"    </td><td> "HDGaussian_p_ABQkDk"    </td></tr>
 * <tr> <td> "HDGaussian_pk_ABQkD"     </td><td> "HDGaussian_p_ABQkD"     </td></tr>
 * <tr> <td> "HDGaussian_pk_ABQD"      </td><td> "HDGaussian_p_ABQD"      </td></tr>
 * <tr> <td> "Categorical_pk_pjk"     </td><td> "Categorical_p_pjk" </td> </tr>
 * <tr> <td> "Categorical_pk_pk"      </td><td> "Categorical_p_pk"  </td> </tr>
 * <tr> <td> "Poisson_pk_ljk"         </td><td> "Poisson_p_ljk"     </td> </tr>
 * <tr> <td> "Poisson_pk_lk"          </td><td> "Poisson_p_lk"      </td> </tr>
 * <tr> <td> "Poisson_pk_ljlk"        </td><td> "Poisson_p_ljlk"    </td> </tr>
 * <tr> <td> "Kmm_pk_sk"              </td><td> "Kmm_p_sk"          </td> </tr>
 * <tr> <td> "Kmm_pk_s"               </td><td> "Kmm_p_s"           </td> </tr>
 * </table>
 *  @param type the String we want to convert
 *  @param[out] freeProp @c true if the model have free proportions, @c false otherwise.
 *  @return the Mixture represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_mixture_ type is returned.
 **/
Mixture stringToMixture( String const& type, bool& freeProp);

/** @ingroup Clustering
 *  convert a Mixture to a String.
 *  @param type the type of Mixture we want to convert
 *  @return the string associated to this type.
 **/
String mixtureToString( Mixture const& type);

/** @ingroup Clustering
 *  convert a Mixture to a string specifying if the model is with free
 *  proportions.
 *  @sa stringToMixture
 *  @param type the Mixture we want to convert
 *  @param freeProp @c true if the model have free proportions, @c false otherwise.
 *  @return the string represented by the Mixture @c type.
 **/
String mixtureToString(Mixture type, bool freeProp);

/** @ingroup Clustering
 *  convert a Mixture to a MixtureClass.
 *  @param type the type of Mixture
 *  @return the MixtureClass associated to this Mixture.
 **/
MixtureClass mixtureToMixtureClass( Mixture const& type);

/** @ingroup Clustering
 * Default number of try in an estimation strategy */
const int defaultNbTry = 5;

/** @ingroup Clustering
 * Default algorithm type in short run */
const Clust::initType defaultInitType = randomFuzzyInit_;
/** @ingroup Clustering
 * Default number of initializations to perform */
const int defaultNbInit = 5;
/** @ingroup Clustering
 * Default algorithm type in initialization */
const Clust::algoType defaultAlgoInInit = emAlgo_;
/** @ingroup Clustering
 * Default number of iteration in an initialization algorithm */
const int defaultNbIterMaxInInit = 20;
/**  @ingroup Clustering
 * Default epsilon in the short runs (used in strategy) */
const Real defaultEpsilonInInit = 1e-02;

/** @ingroup Clustering
 * Default algorithm type in short run */
const Clust::algoType defaultAlgoShortRun = emAlgo_;
/** @ingroup Clustering
 * Default number of iterations in the short runs (used in FullStrategy) */
const int defaultMaxIterShortRun = 200;
/** @ingroup Clustering
 *  Default epsilon in the short runs (used in strategy) */
const Real defaultEpsilonShortRun = 1e-04;

/** @ingroup Clustering
 * Default algorithm type in long run */
const Clust::algoType defaultAlgoLongRun = emAlgo_;
/**  @ingroup Clustering
 * Default number of iterations in the long run (used in FullStrategy) */
const int defaultMaxIterLongRun = 1000;
/**  @ingroup Clustering
 * Default epsilon in the long run (used in strategy) */
const Real defaultEpsilonLongRun = 1e-08;

/** @ingroup Clustering
 *  @param criterion selection criterion to use
 *  @return a pointer on the class computing the criterion
 **/
IMixtureCriterion* createCriterion( Clust::criterionType criterion);

/** @return a pointer on the class computing the criterion
 *  @param criterion string with the criterion name
 **/
STK::IMixtureCriterion* createCriterion( String const& criterion);

/** @ingroup Clustering
 *  utility function for creating an estimation algorithm.
 *  @param algo the algorithm to create
 *  @param nbIterMax,epsilon the maximal number of iteration and the tolerance of the algorithm
 **/
IMixtureAlgo* createAlgo( Clust::algoType algo, int nbIterMax, Real epsilon);

/** @ingroup Clustering
 *  utility function for creating a learning algorithm.
 *  @param algo the algorithm to create
 *  @param nbIterMax,epsilon the maximal number of iteration and the tolerance of the algorithm
 **/
IMixtureAlgoLearn* createLearnAlgo(Clust::algoLearnType algo, int nbIterMax, Real epsilon);

/** @ingroup Clustering
 *  utility function for creating a predicting algorithm.
 *  @param algo the algorithm to create
 *  @param nbIterBurn,nbIterLong,epsilon number of iteration of the burning and estimation steps
 *  and tolerance of the algorithm
 **/
IMixtureAlgoPredict* createPredictAlgo(Clust::algoPredictType algo, int nbIterBurn, int nbIterLong, Real epsilon);

/** @ingroup Clustering
 *  Utility function for creating a model initializer.
 *  @param init the kind of initializer to create
 *  @param nbInits the number of initialization to try
 *  @param algo the kind of algorithm to add to the initializer
 *  @param nbIterMax,epsilon the maximal number of iteration and the tolerance of the initialization algorithm
 **/
IMixtureInit* createInit( Clust::initType init = defaultInitType
                        , int nbInits          = defaultNbInit
                        , Clust::algoType algo = defaultAlgoInInit
                        , int nbIterMax        = defaultNbIterMaxInInit
                        , Real epsilon         = defaultEpsilonInInit);
/** @ingroup Clustering
 *  utility function for creating a a short Run algorithm.
 *  @param algo the algorithm to create
 *  @param nbIterMax the maximal number of iteration of the algorithm
 *  @param epsilon the tolerance of the algorithm
 **/
inline IMixtureAlgo* createShortRunAlgo( Clust::algoType algo = defaultAlgoShortRun
                                       , int nbIterMax        = defaultMaxIterShortRun
                                       , Real epsilon         = defaultEpsilonShortRun)
{ return createAlgo(algo, nbIterMax, epsilon);}
/** @ingroup Clustering
 *  utility function for creating a long Run algorithm.
 *  @param algo the algorithm to create
 *  @param nbIterMax the maximal number of iteration of the algorithm
 *  @param epsilon the tolerance of the algorithm
 **/
inline IMixtureAlgo* createLongRunAlgo( Clust::algoType algo = defaultAlgoLongRun
                                      , int nbIterMax        = defaultMaxIterLongRun
                                      , Real epsilon         = defaultEpsilonLongRun)
{ return createAlgo(algo, nbIterMax, epsilon);}

/** @ingroup Clustering
 *  Utility function for creating a SimpleStrategy.
 *  @param p_composer the composer to which we want to apply a the strategy.
 *  @param nbTry the number of tries.
 *  @param p_init the initializer to use.
 *  @param algo the algorithm to use in the long run.
 *  @return an instance of the SimpleStrategy
 **/
IMixtureStrategy* createSimpleStrategy( IMixtureComposer*& p_composer
                                      , int nbTry
                                      , IMixtureInit* const& p_init
                                      , IMixtureAlgo* const& algo);

/** @ingroup Clustering
 *  Utility function for creating a FullStrategy.
 *  @param p_composer the composer to which we want to apply a the strategy.
 *  @param nbTry the maximal number of tries.
 *  @param nbInitRun the number of initialization to perform.
 *  @param p_init the initializer to use.
 *  @param nbShortRun the number of shortRun.
 *  @param shortRunAlgo the algorithm to use in the short run.
 *  @param longRunAlgo the algorithm to use in the long run.
 *  @return an instance of the FullStrategy
 **/
IMixtureStrategy* createFullStrategy( IMixtureComposer*& p_composer
                                    , int nbTry, int nbInitRun
                                    , IMixtureInit* const& p_init
                                    , int nbShortRun, IMixtureAlgo* const& shortRunAlgo
                                    , IMixtureAlgo* const& longRunAlgo);
}  // namespace Clust

}  // namespace STK

#endif /* STK_CLUST_UTIL_H */
