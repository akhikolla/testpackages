/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015  Serge Iovleff

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
 * Project:  rtkore
 * created on: 13 octobre 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file fastRand.cpp
 *  @brief In this file we test the binding with the R random number generators.
 **/

#include "RTKpp.h"

using namespace STK;

RcppExport SEXP fastBetaRand( SEXP n, SEXP alpha, SEXP beta)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::Beta law(Rcpp::as<double>(alpha), Rcpp::as<double>(beta));
  tab.rand(law);
  return tab.vector();

  END_RCPP;
}

RcppExport SEXP fastBinomialRand( SEXP n, SEXP nb, SEXP prob)
{
  BEGIN_RCPP;
  STK::RVector<int> tab(Rcpp::as<int>(n));
  STK::Law::Binomial law(Rcpp::as<int>(nb), Rcpp::as<double>(prob));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastCategoricalRand( SEXP n, SEXP pr)
{
  BEGIN_RCPP;
  STK::RVector<int> tab(Rcpp::as<int>(n));
  STK::RVector<Real> prob(pr);
  STK::Law::Categorical law(prob);
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastCauchyRand( SEXP n, SEXP mu, SEXP scale)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::Cauchy law(Rcpp::as<double>(mu), Rcpp::as<double>(scale));
  tab.rand(law);
  return tab.vector();

  END_RCPP;
}

RcppExport SEXP fastChiSquaredRand( SEXP n, SEXP df)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::ChiSquared law(Rcpp::as<int>(df));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastExponentialRand( SEXP n, SEXP lambda)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::Exponential law(Rcpp::as<Real>(lambda));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastFisherSnedecorRand( SEXP n, SEXP df1, SEXP df2)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::FisherSnedecor law(Rcpp::as<int>(df1), Rcpp::as<int>(df2));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastGammaRand( SEXP n, SEXP shape, SEXP scale)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::Gamma law(Rcpp::as<Real>(shape), Rcpp::as<Real>(scale));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastGeometricRand( SEXP n, SEXP prob)
{
  BEGIN_RCPP;
  STK::RVector<int> tab(Rcpp::as<int>(n));
  STK::Law::Geometric law(Rcpp::as<Real>(prob));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastHyperGeometricRand( SEXP n,  SEXP nbSuccesses, SEXP nbFailures, SEXP nbDraws)
{
  BEGIN_RCPP;
  STK::RVector<int> tab(Rcpp::as<int>(n));
  STK::Law::HyperGeometric law(Rcpp::as<int>(nbSuccesses), Rcpp::as<int>(nbFailures), Rcpp::as<int>(nbDraws));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastLogisticRand( SEXP n, SEXP mu, SEXP scale)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::Logistic law(Rcpp::as<Real>(mu), Rcpp::as<Real>(scale));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastLogNormalRand( SEXP n, SEXP mu, SEXP sigma)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::LogNormal law(Rcpp::as<Real>(mu), Rcpp::as<Real>(sigma));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastNegativeBinomialRand( SEXP n, SEXP size, SEXP prob)
{
  BEGIN_RCPP;
  STK::RVector<int> tab(Rcpp::as<int>(n));
  STK::Law::NegativeBinomial law(Rcpp::as<int>(size), Rcpp::as<Real>(prob));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastNormalRand( SEXP n, SEXP mu, SEXP sigma)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::Normal law(Rcpp::as<Real>(mu), Rcpp::as<Real>(sigma));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastPoissonRand( SEXP n, SEXP lambda)
{
  BEGIN_RCPP;
  STK::RVector<int> tab(Rcpp::as<int>(n));
  STK::Law::Poisson law(Rcpp::as<Real>(lambda));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastStudentRand( SEXP n, SEXP df)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::Student law(Rcpp::as<int>(df));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastUniformRand( SEXP n, SEXP a, SEXP b)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::Uniform law(Rcpp::as<double>(a), Rcpp::as<double>(b));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastUniformDiscreteRand( SEXP n, SEXP a, SEXP b)
{
  BEGIN_RCPP;
  STK::RVector<int> tab(Rcpp::as<int>(n));
  STK::Law::UniformDiscrete law(Rcpp::as<int>(a), Rcpp::as<int>(b));
  tab.rand(law);
  return tab.vector();
  END_RCPP;
}

RcppExport SEXP fastWeibullRand( SEXP n, SEXP k, SEXP lambda)
{
  BEGIN_RCPP;
  STK::RVector<double> tab(Rcpp::as<int>(n));
  STK::Law::Weibull law(Rcpp::as<double>(k), Rcpp::as<double>(lambda));
  tab.rand(law);
  return tab.vector();

  END_RCPP;
}





