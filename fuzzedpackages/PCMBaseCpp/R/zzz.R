# Copyright 2018 Venelin Mitov
#
# This file is part of PCMBaseCpp.
#
# PCMBase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMBase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMBase.  If not, see <http://www.gnu.org/licenses/>.

#' @import Rcpp
#' @import methods
#' @useDynLib PCMBaseCpp
loadModule( "PCMBaseCpp__Tree", TRUE )
loadModule( "PCMBaseCpp__OrderedTree", TRUE )
loadModule( "PCMBaseCpp__QuadraticPolyBM", TRUE )
loadModule( "PCMBaseCpp__QuadraticPolyBM1D", TRUE )
loadModule( "PCMBaseCpp__QuadraticPolyDOU", TRUE )
loadModule( "PCMBaseCpp__QuadraticPolyMixedGaussian", TRUE )
loadModule( "PCMBaseCpp__QuadraticPolyMixedGaussian1D", TRUE )
loadModule( "PCMBaseCpp__QuadraticPolyOU", TRUE )
loadModule( "PCMBaseCpp__QuadraticPolyOU1D", TRUE )
loadModule( "PCMBaseCpp__QuadraticPolyJOU", TRUE )
loadModule( "PCMBaseCpp__QuadraticPolyWhite", TRUE )
