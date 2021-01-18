/***************************************************************************
 *   AstonGeostats, algorithms for low-rank geostatistical models          *
 *                                                                         *
 *   Copyright (C) Ben Ingram, 2008                                        *
 *                                                                         *
 *   Ben Ingram, IngramBR@Aston.ac.uk                                      *
 *   Neural Computing Research Group,                                      *
 *   Aston University,                                                     *
 *   Aston Street, Aston Triangle,                                         *
 *   Birmingham. B4 7ET.                                                   *
 *   United Kingdom                                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef GAUSSIANPROCESS_H_
#define GAUSSIANPROCESS_H_

#include "ForwardModel.h"
#include "Optimisable.h"
#include "CovarianceFunction.h"
#include "psgp_common.h"

class GaussianProcess : public ForwardModel, public Optimisable
{
public:
	GaussianProcess(int Inputs, int Outputs, mat& Xdata, vec& ydata, CovarianceFunction& cf);
	virtual ~GaussianProcess();

	void   makePredictions(vec& Mean, vec& Variance, const mat& Xpred, const mat& C) const;
	void   makePredictions(vec& Mean, vec& Variance, const mat& Xpred, CovarianceFunction &cf) const;
	void   makePredictions(vec& Mean, vec& Variance, const mat& Xpred) const;
	double loglikelihood() const;

	vec    getParametersVector() const;
	void   setParametersVector(const vec p);

	double objective() const;
	vec    gradient() const;

	void   estimateParameters();

	
	
	
private:

	mat    computeCholesky(const mat& iM) const;
	mat    computeInverseFromCholesky(const mat& C) const;

	vec    getGradientVector() const;

	CovarianceFunction& covFunc;
	mat& Locations;
	vec& Observations;

};

#endif /*GAUSSIANPROCESS_H_*/
