# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
#                    Julien Papaix <julien.papaix@inrae.fr>
#                    Jean-François Rey <jean-francois.rey@inrae.fr>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,i
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#' @title Class LandsepiParams
#' @name LandsepiParams
#' @description Landsepi simulation parameters
#' @details An object of class LandsepiParams that can be created by calling \code{\link{createSimulParams}}
#' @aliases LandsepiParams-class
#'
#' @slot Landscape a landscape as sf object.See \code{\link{loadLandscape}}, \code{\link{loadLandscape}}
#' @slot Croptypes a dataframe with three columns named 'croptypeID' for croptype index,
#' 'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar within the croptype.
#' See \code{\link{loadCroptypes}}, \code{\link{setCroptypes}} and See \code{\link{allocateCroptypeCultivars}}
#' @slot Cultivars a dataframe of parameters associated with each host genotype (i.e. cultivars, lines)
#' when cultivated in pure crops.See \code{\link{loadCultivar}} and \code{\link{setCultivars}}
#' @slot CultivarsGenes a list containing, for each host genotype, the indices of carried resistance genes.
#' See \code{\link{allocateCultivarGenes}}
#' @slot Genes a data.frame of parameters associated with each resistance gene and with the evolution of
#' each corresponding pathogenicity gene. See \code{\link{loadGene}} and \code{\link{setGenes}}
#' @slot Pathogen a list of pathogen aggressiveness parameters on a susceptible host
#' for a pathogen genotype not adapted to resistance. See \code{\link{loadPathogen}} and \code{\link{setPathogen}}
#' @slot PI0 initial probability for the first host (whose index is 0) to be infectious (i.e. state I)
#' at the beginning of the simulation. Must be between 0 and 1. See \code{\link{setInoculum}}
#' @slot DispHost a vectorized matrix giving the probability of host dispersal
#' from any field of the landscape to any other field. See \code{\link{loadDispersalHost}} and \code{\link{setDispersalHost}}
#' @slot DispPatho a vectorized matrix giving the probability of pathogen dispersal
#' from any field of the landscape to any other field. See \code{\link{loadDispersalPathogen}} and \code{\link{setDispersalPathogen}}
#' @slot OutputDir the directory for simulation outputs 
#' @slot OutputGPKG the name of the output GPKG file containing parameters of the deployment strategy
#' @slot Outputs a list of outputs parameters. See \code{\link{setOutputs}}
#' @slot TimeParam a list of time parameters. See \code{\link{setTime}}
#' @slot Seed an integer used as seed value (for random number generator). See \code{\link{setTime}}
#' @import sf
#' @exportClass LandsepiParams

setClass(
  Class = "LandsepiParams",
  slots = list(
    Landscape = "sf",
    Croptypes = "data.frame",
    Cultivars = "data.frame",
    CultivarsGenes = "data.frame",
    Genes = "data.frame",
    Pathogen = "list",
    PI0 = "numeric",
    DispHost = "vector",
    DispPatho = "vector",
    OutputDir = "character",
    OutputGPKG = "character",
    Outputs = "list",
    TimeParam = "list", # Nyear Nstep
    Seed = "numeric"
  )
)

setValidity("LandsepiParams", function(object) {
  if (is.null(object@OutputDir) || length(object@OutputDir) == 0) { # check is valid dir .dir.file
    "Need an output directory"
  } else {
    TRUE
  }
})
