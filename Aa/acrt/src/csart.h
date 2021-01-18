/*   Copyright (C) 2016  David Preinerstorfer
#    david.preinerstorfer@econ.au.dk
#
#    This file is a part of acrt.
#
#    acrt is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/
*/


#ifndef CSART_H_
#define CSART_H_

Eigen::MatrixXd premult(Eigen::VectorXd partial, int dim);

Eigen::VectorXd  ctest(Eigen::MatrixXd umat,
                Eigen::MatrixXd Rbmat,
                Eigen::MatrixXd Wmat, 
                Eigen::MatrixXd Bmat,
                int cores); 

Eigen::VectorXd  ctestE(Eigen::MatrixXd umat,
                Eigen::MatrixXd Rbmat,
                Eigen::MatrixXd Wmat, 
                Eigen::MatrixXd Bmat,
                int cores);
 
#endif
