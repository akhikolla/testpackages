# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

.onLoad <- function(libname, pkgname) { # nocov start
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
} # nocov end

.onAttach <- function(...) {
  ssMousetrackLib <- dirname(system.file(package = "ssMousetrack"))
  pkgdesc <- suppressWarnings(utils::packageDescription("ssMousetrack", lib.loc = ssMousetrackLib))
  if (length(pkgdesc) > 1) {
    builddate <- gsub(';.*$', '', pkgdesc$Packaged)
    #packageStartupMessage(paste("ssMousetrack (Version ", pkgdesc$Version, ", packaged: ", builddate, ")", sep = ""))
    pkgVs <- paste("ssMousetrack (Version ", pkgdesc$Version, ", packaged: ", builddate, ")", sep = "")
    msg <- paste0("\n",
                  pkgVs,
                  "\n",
                  "Please cite our work! Type citation(\"ssMousetrack\") for details.\n",
                  "Submit your suggestions and bug-reports at: https://github.com/antcalcagni/ssMousetrack/issues\n"
    )
    packageStartupMessage(msg)
  }
  packageStartupMessage("")
}

