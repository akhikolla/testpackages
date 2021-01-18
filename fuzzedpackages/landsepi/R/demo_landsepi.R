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

#' @title Package demonstration
#' @name demo_landsepi
#' @description run a simulation demonstration with landsepi
#' @param seed an interger used as seed for Random Number Generator (default 12345).
#' @param strat a string specifying the deployment strategy: "MO" for mosaic of resistant 
#' cultivars, "MI" for intra-fied mixtures, "RO" for cultivar rotations, and "PY" for 
#' resistance gene pyramiding in a cultivar.
#' @param Nyears number of cropping seasons (years) to simulate.
#' @param nTSpY number of time-steps (days) per cropping season.
#' @param videoMP4 a logical indicating if a video must be generated (TRUE, default) or not (FALSE).
#' @details In these examples, 2 completely efficient resistance sources (typical of major 
#' resistance genes) are deployed in the landscape according 
#' to one of the following strategies:\itemize{
#' \item Mosaic: 3 pure crops (S + R1 + R2) with very high spatial aggregation.
#' \item Mixture: 1 pure susceptible crop + 1 mixture of two resistant cultivars, with high aggregation.
#' \item Rotation: 1 susceptible pure crop + 2 resistant crops in alternation every 2 years
#' , with moderate aggregation.
#' \item Pyramiding: 1 susceptible crop + 1 pyramided cultivar in a fragmented landscape (low aggregation).
#' }
#' @return A set of text files, graphics and a video showing epidemic dynamics.
#' @seealso \link{runSimul}, \link{runShinyApp}
#' @examples 
#' \dontrun{
#' ## Run demonstrations (in 20-year simulations) for different deployment strategies:
#' demo_landsepi(strat = "MO") ## for a mosaic of cultivars
#' demo_landsepi(strat = "MI") ## for a mixture of cultivars
#' demo_landsepi(strat = "RO") ## for a rotation of cultivars
#' demo_landsepi(strat = "PY") ## for a pyramid of resistance genes
#' }
#' @include Methods-LandsepiParams.R RcppExports.R
#' @export
demo_landsepi <- function(seed = 12345, strat = "MO", Nyears = 20, nTSpY = 120, videoMP4 = TRUE) {
  # seed = 1; strat="MO"; Nyears = 5; nTSpY = 120; videoMP4 = FALSE  ## for debugging
  initPath <- getwd()

  ## Simulation parameters
  simul_params <- createSimulParams(outputDir = getwd())
  
  ## Seed
  simul_params <- setSeed(simul_params, seed)
  
  ## Time parameters
  simul_params <- setTime(simul_params, Nyears, nTSpY)

  ## Pathogen parameters
  basic_patho_param <- loadPathogen("rust")
  # basic_patho_param <- list(infection_rate = 0.4, latent_period_exp = 10, latent_period_var = 9
  #                           , propagule_prod_rate = 3.125, infectious_period_exp = 24, infectious_period_var = 105
  #                           , survival_prob = 1e-4, repro_sex_prob = 0
  #                           , sigmoid_kappa = 5.333, sigmoid_sigma = 3, sigmoid_plateau = 1)
  simul_params <- setPathogen(simul_params, basic_patho_param)


  ## Initial conditions
  simul_params <- setInoculum(simul_params, 5e-4)


  ## Landscape parameters
  id_landscape <- 1
  landscape <- loadLandscape(id_landscape)
  # landscape <- get("landscapeTEST1")
  simul_params <- setLandscape(simul_params, landscape)


  ## Dispersal parameters
  disp_patho <- loadDispersalPathogen(id_landscape)
  # disp_patho <- get("dispP_1")
  simul_params <- setDispersalPathogen(simul_params, disp_patho)


  ## Genes
  gene1 <- loadGene(name = "MG 1", type = "majorGene")
  gene2 <- loadGene(name = "MG 2", type = "majorGene")
  if (strat == "PY") {
    gene1$mutation_prob <- 1E-4
    gene2$mutation_prob <- 1E-4
  }
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  # genes <- data.frame(geneName =               c("MG1", "MG2"),
  #                     efficiency =             c(1.0  , 1.0  ),
  #                     time_to_activ_exp =      c(0.0  , 0.0  ),
  #                     time_to_activ_var =      c(0.0  , 0.0  ),
  #                     mutation_prob =          c(1E-7 , 1E-7 ),
  #                     Nlevels_aggressiveness = c(2    , 2    ),
  #                     fitness_cost =           c(0.5  , 0.5  ),
  #                     tradeoff_strength =      c(1.0  , 1.0  ),
  #                     target_trait =           c("IR" , "IR"  ),
  #                     stringsAsFactors = FALSE)
  simul_params <- setGenes(simul_params, genes)


  ## Cultivars
  cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
  if (strat == "PY") { ## 1 susceptible cultivar + 1 resistant cultivar
    cultivar2 <- loadCultivar(name = "Resistant", type = "growingHost")
    cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
  } else { ## 1 susceptible cultivar + 2 resistant cultivars
    cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
    cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
    cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  }
  # cultivars <- data.frame(cultivarName = c("Susceptible", "Resistant1", "Resistant2"),
  #                         initial_density =   rep(0.1, 3),
  #                         max_density =       rep(2.0, 3),
  #                         growth_rate =       rep(0.1, 3),
  #                         reproduction_rate = rep(0.0, 3),
  #                         death_rate =        rep(0.0, 3),
  #                         yield_H =           rep(2.5, 3),
  #                         yield_L =           rep(0.0, 3),
  #                         yield_I =           rep(0.0, 3),
  #                         yield_R =           rep(0.0, 3),
  #                         production_cost =   rep(225, 3),
  #                         market_value =      rep(200, 3),
  #                         stringsAsFactors = FALSE)
  simul_params <- setCultivars(simul_params, cultivars)


  ## Allocate genes to cultivars
  if (strat == "PY") { ## 2 genes for 1 cultivar
    simul_params <- allocateCultivarGenes(simul_params, "Resistant", c("MG 1", "MG 2"))
  } else { ## 2 genes for 2 cultivars
    simul_params <- allocateCultivarGenes(simul_params, "Resistant1", c("MG 1"))
    simul_params <- allocateCultivarGenes(simul_params, "Resistant2", c("MG 2"))
  }


  ## Allocate cultivars to croptypes
  switch(strat,
    "MO" = { ## 3 pure crops, very high aggregation
      croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Resistant crop 1", "Resistant crop 2"))
      croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
      croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
      croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
      simul_params <- setCroptypes(simul_params, croptypes)


      rotation_sequence <- croptypes$croptypeID ## No rotation -> 1 rotation_sequence element
      rotation_period <- 0 ## same croptypes every years
      prop <- c(1 / 3, 1 / 3, 1 / 3) ## croptypes proportions
      aggreg <- 10 ## aggregated landscape
    },
    "MI" = { ## 1 pure crop + 1 balanced mixture, high aggregation
      croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Mixture"))
      croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
      croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Resistant1", "Resistant2"))
      # croptypes <- data.frame(croptypeID = c(0, 1), croptypeName = c("Susceptible crop", "Mixture")
      #                         , Susceptible = c(1.0, 0), Resistant1 = c(0, 0.5), Resistant2 = c(0, 0.5)
      #                         , stringsAsFactors = FALSE)
      simul_params <- setCroptypes(simul_params, croptypes)


      rotation_sequence <- croptypes$croptypeID
      rotation_period <- 0
      prop <- c(1 / 2, 1 / 2)
      aggreg <- 0.75
    },
    "RO" = { ## 1 pure crop + 2 pure crops in alternation every 2 years, moderate aggregation
      croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Resistant crop 1", "Resistant crop 2"))
      croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
      croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
      croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
      # croptypes <- data.frame(croptypeID = c(0, 1, 2), croptypeName = c("Susceptible crop","Resistant crop 1","Resistant crop 2")
      #                         , Susceptible = c(1.0, 0, 0), Resistant1 = c(0, 1.0, 0), Resistant2 = c(0, 0, 1.0)
      #                         , stringsAsFactors = FALSE)
      simul_params <- setCroptypes(simul_params, croptypes)


      rotation_sequence <- list(
        c(croptypes$croptypeID[1], croptypes$croptypeID[2]),
        c(croptypes$croptypeID[1], croptypes$croptypeID[3])
      )
      rotation_period <- 2 ## rotation every 2 years
      prop <- list(c(1 / 2, 1 / 2))
      aggreg <- 0.25
    },
    "PY" = { ## 2 pure crops, fragmented landscape (low aggregation)
      croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Pyramid"))
      croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
      croptypes <- allocateCroptypeCultivars(croptypes, "Pyramid", "Resistant")
      # croptypes <- data.frame(croptypeID = c(0, 1), croptypeName = c("Susceptible crop", "Pyramid")
      #                         , Susceptible = c(1.0, 0), Resistant = c(0, 1.0)
      #                         , stringsAsFactors = FALSE)
      simul_params <- setCroptypes(simul_params, croptypes)


      rotation_sequence <- croptypes$croptypeID
      rotation_period <- 0
      prop <- c(1 / 2, 1 / 2)
      aggreg <- 0.07
    },
    {
      stop('Unknown strategy. Possible values for argument "strat" are: "MO", "MI"", "RO", "PY"')
    }
  )


  ## Allocate croptypes to landscape
  simul_params <- allocateLandscapeCroptypes(simul_params,
    rotation_period = rotation_period,
    rotation_sequence = rotation_sequence,
    rotation_realloc = FALSE,
    prop = prop,
    aggreg = aggreg
  )

  ## list of outputs to be generated
  outputlist <- loadOutputs()
  simul_params <- setOutputs(simul_params, outputlist)

  ## Check simulation parameters
  checkSimulParams(simul_params)

  ## Save deployment strategy into GPKG file
  simul_params <- saveDeploymentStrategy(simul_params)

  ## Run simulation
  runSimul(simul_params, graphic = TRUE, videoMP4 = videoMP4)

  ## Set back the path to the initial repository
  setwd(initPath)
}
