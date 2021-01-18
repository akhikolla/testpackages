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
 * created on: 3 mars 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file rtkore.h
 *  @brief In this file .
 **/


/** @file rtkore.h
 *  @brief define the entry C/C++ entry functions for R
 **/
#ifndef RTKORE_H
#define RTKORE_H

#ifdef __cplusplus
extern "C"
{
#endif

SEXP stk_version(SEXP robj);
SEXP fastBetaRand( SEXP n, SEXP alpha, SEXP beta);
SEXP fastBinomialRand( SEXP n, SEXP nb, SEXP prob);
SEXP fastCategoricalRand( SEXP n, SEXP pr);
SEXP fastCauchyRand( SEXP n, SEXP mu, SEXP scale);
SEXP fastChiSquaredRand( SEXP n, SEXP df);
SEXP fastExponentialRand( SEXP n, SEXP lambda);
SEXP fastFisherSnedecorRand( SEXP n, SEXP df1, SEXP df2);
SEXP fastGammaRand( SEXP n, SEXP shape, SEXP scale);
SEXP fastGeometricRand( SEXP n, SEXP prob);
SEXP fastHyperGeometricRand( SEXP n,  SEXP nbSuccesses, SEXP nbFailures, SEXP nbDraws);
SEXP fastLogisticRand( SEXP n, SEXP mu, SEXP scale);
SEXP fastLogNormalRand( SEXP n, SEXP mu, SEXP sigma);
SEXP fastNegativeBinomialRand( SEXP n, SEXP size, SEXP prob);
SEXP fastNormalRand( SEXP n, SEXP mu, SEXP sigma);
SEXP fastPoissonRand( SEXP n, SEXP lambda);
SEXP fastStudentRand( SEXP n, SEXP df);
SEXP fastUniformRand( SEXP n, SEXP a, SEXP b);
SEXP fastUniformDiscreteRand( SEXP n, SEXP a, SEXP b);
SEXP fastWeibullRand( SEXP n, SEXP k, SEXP lambda);

#ifdef __cplusplus
extern "C"
}
#endif /*__cplusplus */

#endif /* RTKORE_H*/





