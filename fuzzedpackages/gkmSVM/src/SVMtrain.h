/* SVMtrain.h
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include "global.h"
class CSVMtrain
{
public:
	CSVMtrain(void);

	void train(double **kernel, int npos, int nneg, double *lambdas); // calculates labdas using iterative method 

    double evaluateObjFunc(double **kernel, int npos, int nneg, double *lambdas); // returns the quadratic obj func value. This is not used by train. it is added to be used to optimize parameters used to define kernel
	~CSVMtrain(void);
	int niter20; 
private:
	void initLambdas(double *x, int n);

};

