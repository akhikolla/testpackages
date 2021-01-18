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

#' Check if the PCMBaseCpp version correpsonds to a dev release
#' @description This function is used during unit-testing, to disable some unit-
#' tests which run extremely long or are consistently failing on some systems.
#' @return a logical
#' @export
PCMBaseCppIsADevRelease <- function() {
  TRUE
}

