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
 * Purpose:  Primary include file for STatistiK project.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STatistiK.h
 *  @brief This file include all the header files of the project STatistiK.
 *
 * @defgroup STatistiK STatistiK (Statistical tools).
 * @brief The StatistiK project contains the main tools for computing the usual
 * statistics.
 *
 * It is divided in two sub-projects:
 * - @ref Laws
 * - @ref StatDesc
 **/

/** @ingroup STatistiK
 *  @defgroup Laws  The probabilities laws sub-project.
 *  In this sub-project, we compute and simulate the usual probabilities laws:
 *  normal law, binomial law, Cauchy law,...
 **/

/** @ingroup STatistiK
 *  @defgroup StatDesc The descriptive statistics sub-project.
 *  In this sub-project, we compute the usual descriptive statistics of variables.
 **/

/** @ingroup STatistiK
 *  @defgroup Kernel Positive kernels.
 *  In this sub-project, we compute the usual positive kernels used by rkhs methods.
 **/

/** @ingroup STatistiK
 *  @namespace STK::Law
 *  @brief This is the namespace enclosing the classes handling the usual
 *  probabilities Laws.
 **/

/** @ingroup STatistiK
 *  @namespace STK::MultiLaw
 *  @brief This is the namespace enclosing the classes handling the multivariate
 *  probabilities Laws.
 **/

/** @ingroup STatistiK
 *  @namespace STK::Stat
 *  @brief this is the namespace for the statistical treatment.
 *  The namespace Stat enclose the methods and classes  for usual statistical
 *  treatment of the variable like mean, variance, covariance, ...
 **/

/** @ingroup STatistiK
 *  @namespace STK::Kernel
 *  @brief this is the namespace for the Kernels.
 *  The namespace Kernel enclose the methods and classes for usual positive
 *  kernels needed by rkhs methods.
 **/

#ifndef STATISTIK_H
#define STATISTIK_H

// random number generators
#include <STatistiK/include/MersenneTwister.h>
#include <STatistiK/include/STK_RandBase.h>

// namespace Law
// probabilities laws
#include <STatistiK/include/STK_Law_Util.h>
#include <STatistiK/include/STK_Law_Functors.h>

#include <STatistiK/include/STK_Law_Bernoulli.h>
#include <STatistiK/include/STK_Law_Beta.h>
#include <STatistiK/include/STK_Law_Binomial.h>
#include <STatistiK/include/STK_Law_Categorical.h>
#include <STatistiK/include/STK_Law_Cauchy.h>
#include <STatistiK/include/STK_Law_ChiSquared.h>
#include <STatistiK/include/STK_Law_Exponential.h>
#include <STatistiK/include/STK_Law_FisherSnedecor.h>
#include <STatistiK/include/STK_Law_Gamma.h>
#include <STatistiK/include/STK_Law_Geometric.h>
#include <STatistiK/include/STK_Law_HyperGeometric.h>
#include <STatistiK/include/STK_Law_Logistic.h>
#include <STatistiK/include/STK_Law_LogNormal.h>
#include <STatistiK/include/STK_Law_NegativeBinomial.h>
#include <STatistiK/include/STK_Law_Normal.h>
#include <STatistiK/include/STK_Law_Poisson.h>
#include <STatistiK/include/STK_Law_Student.h>
#include <STatistiK/include/STK_Law_Uniform.h>
#include <STatistiK/include/STK_Law_UniformDiscrete.h>
#include <STatistiK/include/STK_Law_Weibull.h>

#include <STatistiK/include/STK_MultiLaw_Normal.h>
#include <STatistiK/include/STK_MultiLaw_JointBernoulli.h>
#include <STatistiK/include/STK_MultiLaw_JointCauchy.h>
#include <STatistiK/include/STK_MultiLaw_JointNormal.h>
#include <STatistiK/include/STK_MultiLaw_JointGamma.h>

#include <STatistiK/include/STK_Stat_Functors.h>
#include <STatistiK/include/STK_Stat_UnivariateReal.h>

// bivariate Statistics
#include <STatistiK/include/STK_Stat_Bivariate.h>
#include <STatistiK/include/STK_Stat_Covariance.h>

// Multivariate Statistics
#include <STatistiK/include/STK_Stat_Multivariate.h>
#include <STatistiK/include/STK_Stat_MultivariateReal.h>

// perform the usual Computations on categorical variables
#include <STatistiK/include/STK_Stat_Factor.h>
#include <STatistiK/include/STK_Stat_MultiFactor.h>
#include <STatistiK/include/STK_Stat_ConfusionMatrix.h>
#include <STatistiK/include/STK_Stat_Transform.h>

// Kernels
#include <STatistiK/include/STK_Kernel_Util.h>

#include <STatistiK/include/STK_Kernel_Gaussian.h>
#include <STatistiK/include/STK_Kernel_Laplace.h>
#include <STatistiK/include/STK_Kernel_Linear.h>
#include <STatistiK/include/STK_Kernel_Polynomial.h>
#include <STatistiK/include/STK_Kernel_RationalQuadratic.h>
#include <STatistiK/include/STK_Kernel_Hamming.h>

#endif /*STATISTIK_H*/

