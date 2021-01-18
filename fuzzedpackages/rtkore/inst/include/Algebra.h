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
 * Project:  Algebra
 * Purpose:  Primary include file for Algebra sub-project.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file Algebra.h
 *  @brief This file include all the header files of the project Algebra.
 **/

/** @defgroup Algebra Algebra
 *  @brief The Algebra project provides structures, tools and methods of the
 *  usual algebra techniques.
 *
 * The Algebra project proposes some set of template function for computing
 * dot product,  weighted dot product, vector norm, weighted vector norm
 * and so on... It implement also linear algebra methods for the ArrayXX
 * and the ArraySquareX classes:
 *    @li The Qr decomposition of an arbitrary matrix of Real, @sa Qr
 *    @li The svd decomposition of an arbitrary matrix of real, @sa Svd
 *    @li An Eigenvalue decomposition for symmetric (square) matrices and a
 *    generalized inverse method for such matrices, @sa SymEigen.
*    @li methods for solving the linear least square problem.
 *
 * It proposes also some set of method for performing
 * @li Givens rotations on a matrix
 * @li GramSchmidt orthogonalization of the columns of a matrix/array
 * @li Householder rotations of an array/matrix.
 * @li Cholesky decomposition.
 *
 **/

/** @ingroup Algebra
 *  @namespace STK::lapack
 *  @brief namespace enclosing the wrappers of the lapack routines.
 **/


#ifndef Algebra_H
#define Algebra_H

/* Utilities Algebra methods. */
#include <Algebra/include/STK_transpose.h>
#include <Algebra/include/STK_Givens.h>
#include <Algebra/include/STK_GramSchmidt.h>
#include <Algebra/include/STK_Householder.h>
#include <Algebra/include/STK_Cholesky.h>

/* Algebra methods */
#include <Algebra/include/STK_CG.h>
#include <Algebra/include/STK_InvertMatrix.h>

/* built-in */
#include <Algebra/include/STK_Qr.h>
#include <Algebra/include/STK_Svd.h>
#include <Algebra/include/STK_SymEigen.h>
#include <Algebra/include/STK_MultiLeastSquare.h>

/* lapack */
#include <Algebra/include/STK_lapack_Qr.h>
#include <Algebra/include/STK_lapack_Svd.h>
#include <Algebra/include/STK_lapack_SymEigen.h>
#include <Algebra/include/STK_lapack_MultiLeastSquare.h>

#endif /*Algebra_H*/

