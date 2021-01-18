/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CMAESBOUND_H

#define _CMAESBOUND_H 1

#include "define.h"
#include "domain.h"
#include "newutil.h"

namespace realea {

class IEvalReal {
   public:
	virtual double eval(Real *sol, unsigned n)=0;
	virtual ~IEvalReal(void) {}
};


namespace internal {




/**
 *  @class CMAESBound
 *
 *  @brief it is responsable of the CMAES behaviour among the boundaries
 */
class CMAESBound {
public:
    /***
     * Constructor.
     *
     * @param domain domain search
     */
    CMAESBound(IEvalReal *eval, DomainRealPtr domain);
    
 /**
 * Init the parameter values
 *
 * @param lambda lambda value.
 * @param mueff.
 * @param sigma.
 * @param C covariance matrix.
 */
void setParam(int lambda, double mueff, ColumnVector &sigma, MyMatrix &C);

/**
 * Given an individual it return two fitness, the real and the used for obtain the mean of best lambda.  
 *
 * (Note: if domain.isBound() is false it does not do anything).
 *
 * If it is avanzed, apply the fitness criterion only when the mean is out of the bounds
 * It it is not avanzed, only return the solution clipped in the bounds, and fitness == fitness_sel
 *
 * @param 
 * @param arx current solutions to check
 * @param xmean mean 
 *
 * @param arxvalid Resulting fitness.
 * @param fitness vector of obtained Fitness 
 * @param fitness_sel vector of fitness to use into the selection processs
 *
 */
void evalSols(ColumnVector &xmean, MyMatrix &arx, MyMatrix &arxvalid, RowVector &fitness_raw, RowVector &fitness_sel);

private:
    IEvalReal *m_eval;

    DomainRealPtr m_domain; /** Domain */

    /**
     * @var True means that it admits points out the bounderies, and only it act when the prope mean
     * 	is out. False mean that it act when a point. 
     * Si está a true sigue el criterio de CMAES avanzado (admite puntos
     * fuera de los bordes, y sólo reacciona cuando la media se encuentra fuera
     * de los mismos).
     */
    bool m_avanzed;

    unsigned m_ndim; /***< Ndim */

    bool m_isactive; 

   /**
    * @var Sigma
    */
   ColumnVector m_sigma;
   /**
    * @var Lambda
    */
   int m_lambda;

   /**
    * @var mueff
    */
   double m_mueff;

   /**
    * @var Diagonal de la Matrix de Covarianza
    */
   ColumnVector m_diagC;
private:
    /**
     * Vector bounds border
     */
    ColumnVector m_scale;

    /**
     * Vector boundaries weights 
     */
    ColumnVector m_weights;

    /**
     * @var BOundaries Histrodigram
     */
    queue<Real> m_dfithist;
    
    /**
     * Actual size 
     */
    unsigned int m_dfithist_size;

    /// delta fit for setting weights
    int m_validfitval;
    bool m_iniphase;

    // The num of application
    unsigned int m_numapplied;
};

}

}

#endif
