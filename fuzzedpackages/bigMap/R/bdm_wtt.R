# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# -----------------------------------------------------------------------------
# +++ Watertrack transform
# -----------------------------------------------------------------------------

wtt.get <- function(pakde)
{
	wtt <- wtt_cpp(pakde$x, pakde$y, pakde$z)
	# clustering information
	wtt$P <- wtt$P +1;
	wtt$M <- wtt$M +1;
	wtt$C <- wtt$C +1
	wtt$s <- length(unique(wtt$C))
	# message out
	cat('clusters:', wtt$s, '\n')
	return(wtt)
}
