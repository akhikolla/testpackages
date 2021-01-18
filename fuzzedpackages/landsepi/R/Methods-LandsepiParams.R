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


# Default data.frame column names
.croptypesColNames <- c("croptypeID", "croptypeName")
.cultivarsColNames <- c("cultivarName", "initial_density", "max_density", "growth_rate"
                        , "reproduction_rate", "death_rate", "yield_H", "yield_L", "yield_I"
                        , "yield_R", "production_cost", "market_value")
.cultivarsGenesColNames <- c()
.geneColNames <- c("geneName", "efficiency", "time_to_activ_exp", "time_to_activ_var"
                   , "mutation_prob", "Nlevels_aggressiveness", "fitness_cost"
                   , "tradeoff_strength", "target_trait")


#' @title LandsepiParams
#' @description Creates and initialises a LandespiParams object with default parameters.
#' @param .Object a LandsepiParam object. 
#' @param Landscape a landscape as sf object.
#' @param Croptypes a dataframe with three columns named 'croptypeID' for croptype index,
#' 'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar 
#' within the croptype.
#' @param Cultivars a dataframe of parameters associated with each host genotype 
#' (i.e. cultivars, lines)
#' when cultivated in pure crops.
#' @param CultivarsGenes a list containing, for each host genotype, the indices of 
#' carried resistance genes.
#' @param Genes a data.frame of parameters associated with each resistance gene and with 
#' the evolution of each corresponding pathogenicity gene.
#' @param Pathogen a list of pathogen aggressiveness parameters on a susceptible host
#' for a pathogen genotype not adapted to resistance.
#' @param PI0 initial probability for the first host (whose index is 0) to be infectious 
#' (i.e. state I) at the beginning of the simulation. Must be between 0 and 1.
#' @param DispHost a vectorized matrix giving the probability of host dispersal
#' from any field of the landscape to any other field
#' @param DispPatho a vectorized matrix giving the probability of pathogen dispersal
#' from any field of the landscape to any other field.
#' @param OutputDir the directory for simulation outputs 
#' @param OutputGPKG the name of the output GPKG file containing parameters of the 
#' deployment strategy
#' @param Outputs a list of outputs parameters.
#' @param TimeParam a list of time parameters.
#' @param Seed an integer used as seed value (for random number generator).
#' @param ... more options
#' @rdname initialize-methods
# @aliases LandsepiParams-method
#' @include Class-LandsepiParams.R
#' @include GPKGTools.R tools.R
setMethod(
  "initialize", "LandsepiParams",
  function(.Object,
           Landscape = st_sf(st_sfc()),
           Croptypes = data.frame(),
           Cultivars = data.frame(matrix(ncol = length(.cultivarsColNames)
                                         , nrow = 0, dimnames = list(NULL, .cultivarsColNames))),
           CultivarsGenes = data.frame(),
           Genes = data.frame(matrix(ncol = length(.geneColNames)
                                     , nrow = 0, dimnames = list(NULL, .geneColNames))),
           Pathogen = list(
             name = "no pathogen",
             survival_prob = 0,
             repro_sex_prob = 0,
             infection_rate = 0,
             propagule_prod_rate = 0,
             latent_period_exp = 0,
             latent_period_var = 0,
             infectious_period_exp = 0,
             infectious_period_var = 0,
             sigmoid_kappa = 0,
             sigmoid_sigma = 0,
             sigmoid_plateau = 1
           ),
           PI0 = 0,
           DispHost = vector(),
           DispPatho = vector(),
           OutputDir = normalizePath(character(getwd())),
           OutputGPKG = "landsepi_landscape.gpkg",
           Outputs = list(epid_outputs = "", evol_outputs = ""
                          , thres_breakdown = NA, GLAnoDis = NA, audpc100S = NA),
           TimeParam = list(),
           Seed = NULL,
           ...) {
    # .Object <- callNextMethod(...)
    .Object@Landscape <- Landscape
    .Object@Croptypes <- Croptypes
    .Object@Cultivars <- Cultivars
    .Object@CultivarsGenes <- CultivarsGenes
    .Object@Genes <- Genes
    .Object@Pathogen <- Pathogen
    .Object@PI0 <- PI0
    .Object@DispHost <- DispHost
    .Object@DispPatho <- DispPatho
    .Object@OutputDir <- OutputDir
    .Object@OutputGPKG <- OutputGPKG
    .Object@Outputs <- Outputs
    .Object@TimeParam <- TimeParam
    .Object@Seed <- Seed

    validObject(.Object)
    .Object
  }
)

#' @name print
#' @title print
#' @description Prints a LandespiParams object.
#' @param x a LandsepiParams object
#' @param ... print options
#' @rdname print-methods
#' @aliases print,LandsepiParams-method
#' @export
setMethod("print", "LandsepiParams", function(x, ...) {
  message("## LandsepiParams values :")
  message("### Landscape")
  print(x@Landscape)
  if (nrow(x@Landscape) != 0) {
    plot(st_geometry(x@Landscape))
  } else {
    message("Nothings to plot")
  }
  message("### Croptypes")
  message(x@Croptypes)
  message("### Cultivars")
  message(x@Cultivars)
  message("### CultivarsGenes")
  message(x@CultivarsGenes)
  message("### Genes")
  message(x@Genes)

  message("### Inoculum : ", x@PI0)
  message("### Nyears : ", x@TimeParam$Nyears)
  message("### nTSpY Number of step by year : ", x@TimeParam$nTSpY)
  message("### Output Directory : ", x@OutputDir)
  message("### Output GPKG : ", x@OutputGPKG)
  message("### Seed : ", x@Seed)
  message("### Outputs")
  message(x@Outputs)
})


#' @name summary
#' @title summary
#' @description Prints the summary of a LandespiParams object.
#' @param object a LandsepiParams object.
#' @rdname summary-methods
#' @aliases summary,LandsepiParams-method
#' @export
setMethod("summary", "LandsepiParams", function(object) {
  message("## LandsepiParam Object slots:\n")

  message("### Landscape : ")
  if (nrow(object@Landscape) == 0) {
    message("\tnot set (see setLandscape method)")
  } else {
    summary(object@Landscape)
  }

  message("### Croptypes (proportions of Cultivars in each croptype) : ")
  if (nrow(object@Croptypes) == 0) {
    message("\tnot set (see setCroptypes method)")
  } else {
    summary(object@Croptypes)
  }

  message("### Cultivars (cultivars parameters) : ")
  if (nrow(object@Croptypes) == 0) {
    message("\tnot set (see setCultivars method)")
  } else {
    summary(object@Cultivars)
  }

  message("### CultivarsGenes (List of Genes by Cultivars) : ")
  if (nrow(object@Croptypes) == 0) {
    message("\tnot set (see setCultivarsGenes method)")
  } else {
    summary(object@CultivarsGenes)
  }

  message("### Genes (Genes parameters) : ")
  if (nrow(object@Genes) == 0) {
    message("\tnot set (see setGenes method)")
  } else {
    summary(object@Genes)
  }

  message("### Pathogen parameters : ")
  if (length(object@Pathogen) == 0) {
    message("\tnot set (see loadPathogen and setPathogen methods)")
  } else {
    summary(object@Pathogen)
  }

  message("### Pathogen Dispersal Matrix (as vector) : ")
  if (length(object@DispPatho) == 0) {
    message("\tnot set (see loadDispersalPathogen and setDispersalPathogen methods)")
  } else {
    summary(object@DispPathogen)
  }

  message("### Host Dispersal Matrix (as vector) : ")
  if (length(object@DispHost) == 0) {
    message("\tnot set (see loadDispersalHost and setDispersalHost methods)")
  } else {
    summary(object@DispHost)
  }

  message("### Inoculum : ", object@PI0)
  message("### Nyears : ", object@TimeParam$Nyears)
  message("### nTSpY Number of step by year : ", object@TimeParam$nTSpY)
  message("### Output Directory : ", object@OutputDir)
  message("### Output GPKG : ", object@OutputGPKG)
  message("### Seed : ", object@Seed)
  message("### Outputs")
  message(object@Outputs)
})


#' @name show
#' @title show
#' @description Shows a LandespiParams object.
#' @param object a LandsepiParams object
#' @rdname show-methods
#' @aliases show,LandsepiParams-method
#' @export
setMethod("show", "LandsepiParams", function(object) {
  print(object)
})



#' @name checkSimulParams
#' @title Check simulation parameters
#' @description Checks validity of a LandsepiParams object.
#' @param params a LandsepiParams Object.
#' @return TRUE if OK for simulation, FALSE otherwise
#' @export
checkSimulParams <- function(params) {
  validity <- TRUE
  validity <- validity && checkLandscape(params)
  validity <- validity && checkCultivars(params)
  validity <- validity && checkCroptypes(params)
  validity <- validity && checkGenes(params)
  validity <- validity && checkCultivarsGenes(params)
  validity <- validity && checkPathogen(params)
  validity <- validity && checkDispersalHost(params)
  validity <- validity && checkDispersalPathogen(params)
  validity <- validity && checkInoculum(params)
  validity <- validity && checkTime(params)
  validity <- validity && checkOutputs(params)
  
  return(validity)
}



#' @name createSimulParams
#' @title Create a LandsepiParams object.
#' @description Creates a default object of class LandsepiParams.
#' @param outputDir ouput directory for simulation (default: current directory)
#' @details Create a default object of class LandsepiParams used to store all 
#' simulation parameters. It also creates a subdirectory in \code{outputDir} 
#' using the date; this directory will contain all simulation outputs.
#' @return a LandsepiParams object initialised with the following context:
#' \itemize{
#' \item random seed
#' \item all pathogen parameters fixed at 0
#' \item no between-field dispersal (neither pathogen nor host)
#' \item no pathogen introduction
#' \item no resistance gene
#' \item no output to generate.
#' }
#' @examples \dontrun{
#' createSimulParams()
#' }
#' @export
createSimulParams <- function(outputDir = "./") {

  ## Avoid subdirectory creation
  if (length(grep("simul_landsepi_", normalizePath(outputDir))) != 0) {
    outputDir <- dirname(normalizePath(outputDir))
  }

  ## create a subdirectory with time
  timeSimul <- paste(strsplit(as.character(Sys.time()), " ")[[1]], collapse = "_")
  nameDir <- paste(outputDir, "/simul_landsepi_", gsub(":", "-", timeSimul), sep = "")
  dir.create(nameDir)

  message("Created output directory : ", normalizePath(nameDir))

  lp <- new("LandsepiParams",
    OutputDir = normalizePath(nameDir),
    Seed = setSeedValue()
  )

  return(lp)
}


#' @name loadSimulParams
#' @title Load simulation parameters
#' @description Loads a GPKG file from the output of a landsepi simulation.
#' @details See \code{\link{saveDeploymentStrategy}}.
#' @param inputGPKG name of the GPKG file.
#' @return a LandsepiParams object.
#' @export
loadSimulParams <- function(inputGPKG = "") {
  lp <- new("LandsepiParams",
    OutputDir = normalizePath(dirname(inputGPKG)),
    OutputGPKG = basename(inputGPKG)
    # Seed = setSeedValue(seed),
    # TimeParam = list(Nyears = Nyears, nTSpY = nTSpY)
  )

  lp <- setLandscape(lp, st_read(dsn = inputGPKG, layer = "croptypeID"))
  lp <- setCroptypes(lp, CroptypeBDD2Params(inputGPKG))
  lp <- setGenes(lp, GeneBDD2Params(inputGPKG))
  lp <- setCultivars(lp, CultivarBDD2Params(inputGPKG))
  lp@CultivarsGenes <- CultivarGeneBDD2Params(inputGPKG)

  ## TODO get all parameters from GPKG and parameters.txt if exist
  ## TODO: doesn't seem to work with croptypes, cultivars and cultivarGenes

  message("not implemented yet for DispHost, DispPatho, PI0, pathogen, seed, time, outputs...")

  return(lp)
}


#' @name saveDeploymentStrategy
#' @title Save landscape and deployment strategy 
#' @description Generates a GPKG file containing the landscape and all parameters of 
#' the deployment strategy
#' @details The function generates a GPKG file in the simulation path. 
#'  The GPKG file contains all input parameters needed to restore the landscape (sf object) 
#'  and deployment strategy (croptypes, cultivars and genes).
#' @param params a LandsepiParams Object.
#' @param outputGPKG name of the GPKG output (default: "landsepi_landscape.gpkg") to be generated.
#' @param overwrite a boolean specifying if existing files can be overwritten (TRUE) or not 
#' (FALSE, default).
#' @return an updated LandsepiParams object.
#' @examples
#' \dontrun{
#' ## Initialisation
#' simul_params <- createSimulParams(outputDir = getwd())
#' ## Time parameters
#' simul_params <- setTime(simul_params, Nyears = 10, nTSpY = 120)
#' ## Landscape
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' ## Genes
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' ## Cultivars
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' ## Allocate genes to cultivars
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant1", c("MG 1"))
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant2", c("MG 2"))
#' ## Allocate cultivars to croptypes
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
#' , "Resistant crop 1"
#' , "Resistant crop 2"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
#' simul_params <- setCroptypes(simul_params, croptypes)
#' ## Allocate croptypes to landscape        
#' rotation_sequence <- croptypes$croptypeID ## No rotation -> 1 rotation_sequence element
#' rotation_period <- 0 ## same croptypes every years
#' prop <- c(1 / 3, 1 / 3, 1 / 3) ## croptypes proportions
#' aggreg <- 10 ## aggregated landscape
#' simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = rotation_period,
#' rotation_sequence = rotation_sequence,
#' rotation_realloc = FALSE, prop = prop, aggreg = aggreg)
#' ## Save into a GPKG file
#' simul_params <- saveDeploymentStrategy(simul_params)
#' }
#' @export
saveDeploymentStrategy <- function(params, outputGPKG = "landsepi_landscape.gpkg", overwrite = FALSE) {
  
  params@OutputGPKG <- outputGPKG
  
  if (!dir.exists(params@OutputDir)) {
    warning("Directory ", params@OutputDir, " does not exist")
    return(params)
  }

  # setwd(params@OutputDir)
  # message("Move to ", params@OutputDir, " directory for simulation")

  if (file.exists(paste0(params@OutputDir, "/", params@OutputGPKG))) {
    if (overwrite == FALSE) {
      warning(params@OutputGPKG, " already exists, can't overwrite it.")
      warning("use overwrite = TRUE to allow replacement of existing files")
      return(params)
    }
    else {
      message("Will overwrite existing files in ", params@OutputDir)
      file.remove(paste0(params@OutputDir, "/", params@OutputGPKG))
    }
  }

  # try to add one more fields (year_Nyears+1), if missing
  if (length(grep("^year_", colnames(params@Landscape))) == params@TimeParam$Nyears) {
    params@Landscape[, paste0("year_", params@TimeParam$Nyears + 1)] <- as.data.frame(
      params@Landscape[, paste0("year_", params@TimeParam$Nyears)])[, 1]
    message("Add one more year of simulation (only for simulation model constraints)")
  }
  ## create gpkg file
  message("Create ", paste0(params@OutputDir, "/", params@OutputGPKG), " file")
  ## save only years_ and area ?
  land <- params@Landscape[, grep("^year_", colnames(params@Landscape))]
  gpkgFile <- createLandscapeGPKG(land, paste0(params@OutputDir, "/", params@OutputGPKG))
  ## add data tables
  GPKGAddTables(gpkgFile)
  ## Fill data tables
  if (nrow(params@Cultivars) > 0) GPKGAddInputData(gpkgFile, table = "Cultivar"
                                                   , data = params2CultivarBDD(params)
                                                   , deleteExistingData = TRUE)
  if (nrow(params@Croptypes) > 0) GPKGAddInputData(gpkgFile, table = "CultivarList"
                                                   , data = params2CroptypeBDD(params)
                                                   , deleteExistingData = TRUE)
  if (nrow(params@Genes) > 0) GPKGAddInputData(gpkgFile, table = "Gene"
                                               , data = params2GeneBDD(params)
                                               , deleteExistingData = TRUE)
  if (nrow(params@CultivarsGenes) > 0) GPKGAddInputData(gpkgFile, table = "GeneList"
                                                        , data = params2GeneListBDD(params)
                                                        , deleteExistingData = TRUE)

  return(params)
}


#' @name runSimul
#' @title Run a simulation
#' @description Runs a simulation with landsepi, 
#' a stochastic, spatially-explicit, demo-genetic model simulating the spread and evolution
#' of a pathogen in a heterogeneous landscape and generating a wide range of epidemiological, 
#' evolutionary and economic outputs.
#' @param params a LandsepiParams Object containing all simulation parameters. Must be initialised 
#' with \code{\link{createSimulParams}} and updated using \code{set*()} methods 
#' (see vignettes for details).
#' @param graphic a logical indicating if graphics must be generated (TRUE, default) 
#' or not (FALSE).
#' @param writeTXT a logical indicating if outputs must be written in text files (TRUE, default) 
#' or not (FALSE).
#' @param videoMP4 a logical indicating if a video must be generated (TRUE) or not (FALSE, default).
#' Works only if graphic=TRUE and audpc is computed.
#' @param keepRawResults a logical indicating if binary files must be kept after the end of 
#' the simulation (default=FALSE). Careful, many files may be generated if keepRawResults=TRUE.
#' @details See ?landsepi for details on the model, assumptions and outputs, and our vignettes 
#' for tutorials (\code{browseVignettes("landsepi")}). The function runs the model simulation using 
#' a LandsepiParams object.
#' Briefly, the model is stochastic, spatially explicit (the basic spatial unit is an 
#' individual field), based on a SEIR (‘susceptible-exposed-infectious-removed’, renamed HLIR 
#' for 'healthy-latent-infectious-removed' to avoid confusions
#'  with 'susceptible host') structure with a discrete time step. It simulates the spread and
#'  evolution of a pathogen in a heterogeneous cropping landscape, across cropping seasons split 
#'  by host harvests which impose potential bottlenecks to the pathogen. A wide array of 
#'  resistance deployment strategies can be simulated and evaluated using several possible 
#'  outputs to assess the epidemiological, evolutionary and economic performance
#'  of deployment strategies.
#' 
#' @return A list containing all required outputs.
#' A set of text files, graphics and a video showing epidemic dynamics can be generated.
#' If keepRawResults=TRUE, a set of binary files is generated for every year of simulation and
#' every compartment: \itemize{
#'  \item H: healthy hosts,
#'  \item Hjuv: juvenile healthy hosts,
#'  \item L: latently infected hosts,
#'  \item I: infectious hosts,
#'  \item R: removed hosts,
#'  \item P: propagules.}
#' Each file indicates for every time-step the number of individuals in each field, and when 
#' appropriate for each host and pathogen genotypes.
#' @examples \dontrun{
#' ## Here is an example of simulation of a mosaic of three cultivars (S + R1 + R2). See our 
#' ## tutorials for more examples.
#' ## Initialisation
#' simul_params <- createSimulParams(outputDir = getwd())
#' ## Seed & Time parameters
#' simul_params <- setSeed(simul_params, seed = 1)
#' simul_params <- setTime(simul_params, Nyears = 10, nTSpY = 120)
#' ## Pathogen & inoculum parameters
#' simul_params <- setPathogen(simul_params, loadPathogen("rust"))
#' simul_params <- setInoculum(simul_params, 5e-4)
#' ## Landscape & dispersal
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' simul_params <- setDispersalPathogen(simul_params, loadDispersalPathogen(1))
#' ## Genes
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' ## Cultivars
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' ## Allocate genes to cultivars
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant1", c("MG 1"))
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant2", c("MG 2"))
#' ## Allocate cultivars to croptypes
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
#' , "Resistant crop 1"
#' , "Resistant crop 2"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
#' simul_params <- setCroptypes(simul_params, croptypes)
#' ## Allocate croptypes to landscape        
#' rotation_sequence <- croptypes$croptypeID ## No rotation -> 1 rotation_sequence element
#' rotation_period <- 0 ## same croptypes every years
#' prop <- c(1 / 3, 1 / 3, 1 / 3) ## croptypes proportions
#' aggreg <- 10 ## aggregated landscape
#' simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = rotation_period,
#' rotation_sequence = rotation_sequence,
#' rotation_realloc = FALSE, prop = prop, aggreg = aggreg)
#' ## list of outputs to be generated
#' simul_params <- setOutputs(simul_params, loadOutputs())
#' ## Check simulation parameters
#' checkSimulParams(simul_params)
#' ## Save deployment strategy into GPKG file
#' simul_params <- saveDeploymentStrategy(simul_params)
#' ## Run simulation
#' runSimul(simul_params)
#' }
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018).
#' Assessing the durability andefficiency of landscape-based strategies to deploy plant 
#' resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
#' @export
runSimul <- function(params, graphic=TRUE, writeTXT=TRUE, videoMP4=FALSE, keepRawResults=FALSE) {

  ### !!!!! BE CAREFUL !!!!! ###
  ### croptypes, cultivars, genes and cultivarsGenes have to be ordored all in the same way
  ### ID have to match row index and col index

  initPath <- getwd()
  setwd(params@OutputDir)

  cultivars_genes_list <- lapply(1:nrow(params@Cultivars), FUN = function(i) {
    return(which(params@CultivarsGenes[i, ] == 1) - 1)
  })

  cdf <- as.data.frame(params@Landscape)
  ncol <- length(grep("^year_", colnames(cdf)) %in% colnames(cdf))
  ## TODO: use value of Nyears in previous line?
  rotation <- as.matrix(cdf[, grep("^year_", colnames(cdf))], ncol = ncol)
  croptypes_cultivars <- params2CroptypeBDD(params)[, c(2, 3, 4)]

  ## Run the simulation
  outputs <- simul_landsepi(
    seed = params@Seed,
    time_param = params@TimeParam,
    croptype_names = params@Croptypes$croptypeName,
    cultivars = params@Cultivars,
    cultivars_genes_list = cultivars_genes_list,
    genes = params@Genes,
    landscape = as_Spatial(st_geometry(params@Landscape)),
    area = as.vector(params@Landscape$area[, 1]),
    rotation = rotation,
    croptypes_cultivars_prop = croptypes_cultivars,
    basic_patho_param = params@Pathogen,
    disp_patho = params@DispPatho,
    disp_host = params@DispHost,
    pI0 = params@PI0,
    epid_outputs = params@Outputs[["epid_outputs"]],
    evol_outputs = params@Outputs[["evol_outputs"]],
    thres_breakdown = params@Outputs[["thres_breakdown"]],
    GLAnoDis = params@Outputs[["GLAnoDis"]],
    audpc100S = params@Outputs[["audpc100S"]],
    graphic = graphic,
    writeTXT = writeTXT,
    videoMP4 = videoMP4,
    keepRawResults = keepRawResults
  )

  setwd(initPath)
  return(outputs)
}


#' @name setSeed
#' @title Set the seed
#' @description Updates a LandsepiParams object with a seed value for random number generator
#' @param params a LandsepiParams Object.
#' @param seed an integer used as seed value (for random number generator).
#' @return a LandsepiParams object.
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setSeed(simul_params, 100)
#' simul_params@Seed
#' }
#' @export
setSeed <- function(params, seed) {
  params@Seed <- setSeedValue(seed)
  return(params)
}


#' @name setTime
#' @title Set time parameters
#' @description Updates a LandsepiParams object with time parameters : Nyears and nTSpY
#' @param params a LandsepiParams Object.
#' @param Nyears an integer giving the number of cropping seasons (e.g. years) to simulate.
#' @param nTSpY an integer giving the number of time steps per cropping season (e.g. days).
#' @return a LandsepiParams object.
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setTime(simul_params, Nyears=10, nTSpY=120)
#' simul_params@TimeParam
#' }
#' @export
setTime <- function(params, Nyears, nTSpY) {
  land <- params@Landscape
  
  if (nrow(land) > 0){
    st_geometry(land) <- NULL
    ldf <- as.data.frame(land)
    
    if (length(grep("^year_", colnames(ldf))) < Nyears) {
      message("Landscape croptypes affectation by year have to be regenerated")
    }
  }
  
  params@TimeParam <- list(Nyears = Nyears, nTSpY = nTSpY)
  checkTime(params)
  return(params)
}

#' @name checkTime
#' @title Check time
#' @description Checks time parameters validity
#' @param params a LandsepiParams Object.
#' @return a boolean TRUE if times are setted.
checkTime <- function(params) {
  validity <- TRUE
  if(is.null(params@TimeParam$Nyears) || !is.numeric(params@TimeParam$Nyears) 
     || length(params@TimeParam$Nyears) != 1 ) {
    warning("Invalid nb of years of simulation: use setTime()")
    validity <- FALSE
  } else if (!is.wholenumber(params@TimeParam$Nyears) 
             || !is.strict.positive(params@TimeParam$Nyears) ){
    warning("Nb of years of simulation must be a whole number > 0")
    validity <- FALSE
  }
  
  
  if(is.null(params@TimeParam$nTSpY) || !is.numeric(params@TimeParam$nTSpY) 
     || length(params@TimeParam$nTSpY) > 1 ) {
    warning("Invalid nb of steps per year: use setTime()")
    validity <- FALSE
  } else if (!is.wholenumber(params@TimeParam$nTSpY) 
             || !is.strict.positive(params@TimeParam$nTSpY)){
    warning("Nb of steps per year must be a whole number > 0")
    validity <- FALSE
  }
  
  return(validity)
}


#' @name setLansdcape
#' @title Set the landscape
#' @description Updates a LandsepiParams object with a sp or sf object as landscape.
#' @details The landscape should be a sp or sf object. Built-in landscape are available using 
#' \code{\link{loadLandscape}}. 
#' See our tutorial (vignettes) for details on how to use your own landscape.
#' If the landscape contains only polygons, croptypes can be allocated later using 
#' \code{\link{allocateLandscapeCroptypes}}.
#' Otherwise the landscape has to contain a data.frame specifying for every year, the index 
#' of the croptype cultivated in each polygon.
#' Each features has a field identified by "year_XX" (XX <- seq(1:Nyears+1)) and containing 
#' the croptype ID.
#'
#' | Features/fields | year_1 | year_2 | ... year_Nyears+1 |
#' |---------------- | ------ | ------ | ----------------- |
#' | polygons1       | 13     | 10     | 13                |
#' | polygonsX       | 2      | 1      | 2                 |
#'
#' @param params a LandsepiParams Object.
#' @param land a landscape as sp or sf object
#' @return a LandsepiParams object.
#' @seealso \link{loadLandscape}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' simul_params@Landscape
#' }
#' @importFrom sf st_as_sf
#' @export
setLandscape <- function(params, land) {
  if (class(land)[1] != "sf") {
    params@Landscape <- st_as_sf(land)
  } else {
    params@Landscape <- land
  }

  params@Landscape$area <- data.frame(area = st_area(params@Landscape))
  
  ## Initialise host and pathogen dispersal with diagonal matrices
  if (length(params@DispHost) == 0){
    disp_host <- loadDispersalHost(params, type = "no")
    params <- setDispersalHost(params, disp_host)
  }
  if (length(params@DispPatho) == 0){
    disp_patho <- loadDispersalHost(params, type = "no") # (same function)
    params <- setDispersalPathogen(params, disp_patho)
  }
  return(params)
}


#' @name loadLandscape
#' @title Load a landscape
#' @description Loads one of the five built-in landscapes simulated using a T-tesselation algorithm 
#' and composed of 155, 154, 152, 153 and 156 fields, respectively.
#' Each landscape is identified by a numeric from 1 to 5.
#' @param id a landscape ID between 1 to 5 (default = 1)
#' @return a landscape in sp format
#' @seealso \link{landscapeTEST}, \link{setLandscape}
#' @examples
#' land <- loadLandscape(1)
#' length(land)
#' @export
loadLandscape <- function(id = 1) {
  if (id >= 1 && id <= 5) {
    land <- get(paste0("landscapeTEST", id))
  }
  else {
    stop("Indices of available landscapes are 1 to 5")
  }

  return(land)
}


#' @name checkLandscape
#' @title Check the landscape
#' @description Checks landscape validity
#' @param params a LandsepiParams Object.
#' @return TRUE if Ok, FALSE otherwise
checkLandscape <- function(params) {
  ret <- TRUE
  ## TODO : check bbox, proj4string, epsg and geometry as POLYGON

  land <- params@Landscape
  st_geometry(land) <- NULL
  ldf <- as.data.frame(land)

  # check CroptypeID present in Landscape
  if (sum(!unique(as.integer(unlist(ldf[, grep("^year_", colnames(ldf))]))) 
          %in% params@Croptypes$croptypeID) != 0) {
    warning("Croptypes undef in Landscape")
    warning("Croptypes ID : ", params@Croptypes$croptypeID)
    warning("Croptypes ID in Landscape : ", unique(unlist(ldf)))
    ret <- FALSE
  }

  # layer croptypeID need Nyears + 1 fields
  if (length(grep("^year_", colnames(ldf))) != params@TimeParam$Nyears + 1) {
    warning("Landscape Fields 'year_X' numbers [", length(grep("^year_", colnames(ldf)))
            , "] differ from Nyears ", params@TimeParam$Nyears, " +1")
    ret <- FALSE
  }

  return(ret)
}


#' @name setDispersalPathogen
#' @title Set pathogen dispersal
#' @description Updates a LandsepiParams object with a pathogen dispersal matrix.
#' @details See tutorial (vignettes) on how to 
#' use your own landscape and compute your own pathogen dispersal kernel. 
#' The disersal matrix a square matrix whose size is the number of fields in the landscape and whose elements are, 
#' for each line i and each column i' the probability that propagules migrate 
#' from field i to field i'.  
#' @param params a LandsepiParams Object.
#' @param mat a square matrix giving the probability of pathogen dispersal
#' from any field of the landscape to any other field. 
#' It can be generated manually, or, alternatively, via \code{\link{loadDispersalPathogen}}. 
#' The size of the matrix must match the number of fields in the landscape, and lines of the matrix must sum 
#' to 1.
#' @return a LandsepiParam object.
#' @seealso \link{loadDispersalPathogen}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' d <- loadDispersalPathogen(1)
#' simul_params <- setDispersalPathogen(simul_params, d)
#' simul_params@DispPatho
#' }
#' @export
setDispersalPathogen <- function(params, mat) {
  if (class(mat)[1] == "matrix") {
    params@DispPatho <- as.vector(mat)
  } else {
    params@DispPatho <- mat
  }

  checkDispersalPathogen(params)

  return(params)
}


#' @name loadDispersalPathogen
#' @title Load a pathogen dispersal matrix
#' @description It loads one of the five built-in vectorised dispersal matrices of rust fungi associated with the 
#' five built-in landscapes. Landscape and DispersalPathogen ID must be the same.
#' @param id a matrix ID between 1 to 5 (must match the ID of the landscape loaded with 
#' \code{\link{loadLandscape}}).
#' @details *landsepi* includes built-in dispersal matrices to represent rust dispersal in the 
#' five built-in landscapes. These have been computed from a power-law dispersal kernel: 
#' \eqn{ g(d) = ((b-2)*(b-1) / (2*pi*a^2)) * (1 +  d/a)^(-b) }
#'  with a=40 the scale parameter and b=7 a parameter related to the width of the dispersal kernel. 
#'  The expected mean dispersal distance is given by 2*a/(b-3)=20 m.
#' @return a vectorised dispersal matrix.
#' @seealso \link{dispP}, \link{setDispersalPathogen}
#' @examples
#' d <- loadDispersalPathogen(1)
#' d
#' @export
loadDispersalPathogen <- function(id = 1) {
  if (id >= 1 && id <= 5) {
    disp <- get(paste0("dispP_", id))
  }
  else {
    warning("Indices of available pathogen dispersal matrices are 1 to 5")
    disp <- numeric()
  }

  return(disp)
}



#' @name checkDispersalPathogen
#' @title Check pathogen dispersal
#' @description Checks pathogen dispersal validity
#' @param params a LandsepiParams Object.
#' @return a boolean TRUE if OK, FALSE otherwise
checkDispersalPathogen <- function(params) {
  ret <- TRUE
  if (length(params@DispPatho) != nrow(params@Landscape) * nrow(params@Landscape)) {
    warning("Size of pathogen dispersal is not landscape features^2")
    ret <- FALSE
  }

  if (sum(params@DispPatho > 1) != 0 || sum(params@DispPatho < 0) != 0) {
    warning("Probabilities of pathogen dispersal must to be in [0,1]")
    warning(params@DispPatho[which(params@DispPatho > 1)])
    warning(params@DispPatho[which(params@DispPatho < 0)])
    ret <- FALSE
  }

  return(ret)
}

#' @name setDispersalHost
#' @title Set host dispersal
#' @description Updates a LandsepiParams object with a host dispersal matrix.
#' @details the dispersal matrix gives the probability for a host individual in a field i (row)
#' to migrate to field j (column) through dispersal. 
#' If the host is a cultivated plant: seeds are harvested and do not disperse. 
#' Thus the dispersal matrix is the identity matrix.
#' @param params a LandsepiParams Object.
#' @param mat a square matrix giving the probability of host dispersal
#' from any field of the landscape to any other field. 
#' It can be generated manually, or, alternatively, via \code{\link{loadDispersalHost}}.
#' The size of the matrix must match the number of fields in the landscape.
#' @return a LandsepiParam object.
#' @seealso \link{loadDispersalHost}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' d <- loadDispersalHost(simul_params)
#' simul_params <- setDispersalHost(simul_params, d)
#' simul_params@DispHost
#' }
#' @export
setDispersalHost <- function(params, mat) {
  if (class(mat)[1] == "matrix") {
    params@DispHost <- as.vector(mat)
  } else {
    params@DispHost <- mat
  }

  checkDispersalHost(params)

  return(params)
}



#' @name loadDispersalHost
#' @title Load a host dispersal matrix
#' @description It loads a vectorised diagonal matrix to simulate no host dispersal.
#' @details as the size of the matrix depends on the number of fields in the landscape, 
#' the landscape must be defined before calling \code{loadDispersalHost}. 
#' @param params a LandsepiParams Object.
#' @param type a character string specifying the type of dispersal ("no" for no dispersal)
#' @return a vectorised dispersal matrix.
#' @seealso \link{setDispersalHost}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' d <- loadDispersalHost(simul_params)
#' d
#' }
#' @export
loadDispersalHost <- function(params, type = "no") {
  if (nrow(params@Landscape) == 0) {
    warning("Lanscape has to be set before loading host dispersal matrix")
  }
  else {
    disp <- switch(type,
                   "no" = diag(1, nrow(params@Landscape), nrow(params@Landscape)),
                   diag(1, nrow(params@Landscape), nrow(params@Landscape))
    )
  }
  
  return(disp)
}


#' @name checkDispersalHost
#' @title Check host dispersal
#' @description Checks host dispersal matrix validity.
#' @param params a LandsepiParams Object.
#' @return a boolean TRUE if OK, FALSE otherwise
checkDispersalHost <- function(params) {
  ret <- TRUE
  if (length(params@DispHost) != nrow(params@Landscape) * nrow(params@Landscape)) {
    warning("Size of vector of host dispersal is not landscape features^2")
    ret <- FALSE
  }

  if (sum(params@DispHost > 1) != 0 || sum(params@DispHost < 0) != 0) {
    warning("Host dispersal probabilities must to be in [0,1]")
    warning(params@DispHost[which(params@DispHost > 1)])
    warning(params@DispHost[which(params@DispHost < 0)])
    ret <- FALSE
  }

  return(ret)
}



#' @name allocateLandscapeCroptypes
#' @title Allocate croptypes to the landscape
#' @description Updates the landscape of a LandsepiParams object with croptype allocation in 
#' every field of the landscape and every year of simulation. Allocation is based on an algorithm 
#' which controls croptype proportions (in surface) and spatio-temporal aggregation.
#' @param params a LandsepiParams Object.
#' @param rotation_period number of years before rotation of the landscape. There is no rotation 
#' if rotation_period=0 or rotation_period=Nyears.
#' @param rotation_sequence a list, each element of the list contains indices of croptypes that 
#' are cultivated during a period given by "rotation_period". There is no change in cultivated 
#' croptypes if the list contains only one element (e.g. only one vector c(0,1,2), indicating 
#' cultivation of croptypes 0, 1 and 2).
#' @param rotation_realloc a logical indicating if a new random allocation of croptypes is 
#' performed when the landscape is rotated (FALSE=static allocation, TRUE=dynamic allocation). 
#' Note that if rotation_realloc=FALSE, all elements of the list "rotation_sequence" must have 
#' the same length, and only the first element of the lists "prop" and "aggreg" will be used.
#' @param prop a list of the same size as "rotation_sequence", each element of the list contains 
#' a vector of the proportions (in surface) associated with the croptypes in "rotation_sequence". 
#' A single vector can be given instead of a list if all elements of "rotation_sequence" are 
#' associated with the same proportions.
#' @param aggreg a list of the same size as "rotation_sequence", each element of the list is a 
#' single double indicating the degree of
#' aggregation of the landscape. This double must greater or equal 0; the greater its value, 
#' the higher the degree of spatial aggregation (roughly, aggreg between 0 and 0.1 for fragmented 
#' landscapes, between 0.1 and 0.5 for balanced landscapes, between 0.5 and 3 for aggregated 
#' landscapes, and above 3 for highly aggregated landscapes). A single double can be given 
#' instead of a list if all elements of "rotation_sequence" are associated with the same level 
#' of aggregation.
#' @param algo the algorithm used for the computation of the variance-covariance matrix 
#' of the multivariate normal distribution: "exp" for exponential function, "periodic" 
#' for periodic function, "random" for random draw (see details of function multiN). 
#' If algo="random", the parameter aggreg is not used.
#' @param graphic a logical indicating if graphics must be generated (TRUE) or not (FALSE).
#' @details An algorithm based on latent Gaussian fields is used to allocate two different croptypes 
#' across the simulated landscapes (e.g. a susceptible and a resistant cultivar, denoted as 
#' SC and RC, respectively). This algorithm allows the control of the proportions of each croptype 
#' in terms of surface coverage, and their level of spatial aggregation. 
#' A random vector of values is drawn from a multivariate normal distribution with expectation 0 
#' and a variance-covariance matrix which depends on the pairwise distances between
#' the centroids of the fields. Next, the croptypes are allocated to different fields 
#' depending on whether each value drawn from the multivariate normal distribution is above 
#' or below a threshold. The proportion of each cultivar in the landscape is controlled by the value
#' of this threshold. To allocate more than two croptypes, \code{AgriLand} uses sequentially 
#' this algorithm. For instance, the allocation of three croptypes (e.g. SC, RC1 and RC2) 
#' is performed as follows:
#' \enumerate{
#' \item the allocation algorithm is run once to segregate the fields where the susceptible 
#' cultivar is grown, and
#' \item the two resistant cultivars (RC1 and RC2) are assigned to the remaining candidate 
#' fields by re-running the allocation algorithm.
#' }
#' @return a LandsepiParams object with Landscape updated with the layer "croptypeID". 
#' It contains croptype allocation in every field of the landscape for all years of simulation.
#' @examples
#' \dontrun{
#' ## Initialisation
#' simul_params <- createSimulParams(outputDir = getwd())
#' ## Time parameters
#' simul_params <- setTime(simul_params, Nyears = 10, nTSpY = 120)
#' ## Landscape
#' simul_params <- setLandscape(simul_params, loadLandscape(1))
#' ## Cultivars
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' ## Allocate cultivars to croptypes
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
#' , "Resistant crop 1"
#' , "Resistant crop 2"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
#' simul_params <- setCroptypes(simul_params, croptypes)
#' ## Allocate croptypes to landscape        
#' rotation_sequence <- croptypes$croptypeID ## No rotation -> 1 rotation_sequence element
#' rotation_period <- 0 ## same croptypes every years
#' prop <- c(1 / 3, 1 / 3, 1 / 3) ## croptypes proportions
#' aggreg <- 10 ## aggregated landscape
#' simul_params <- allocateLandscapeCroptypes(simul_params, rotation_period = rotation_period,
#' rotation_sequence = rotation_sequence,
#' rotation_realloc = FALSE, prop = prop, aggreg = aggreg)
#' simul_params@Landscape
#' }
#' @export
allocateLandscapeCroptypes <- function(params, rotation_period, rotation_sequence
                                       , rotation_realloc = FALSE
                                       , prop, aggreg, algo = "periodic", graphic = TRUE) {
  ### TODO Check validity of params slot before Agriland

  orig_landscape <- params@Landscape

  croptypeSP <- AgriLand(as_Spatial(params@Landscape),
    Nyears = params@TimeParam$Nyears,
    rotation_period = rotation_period,
    rotation_sequence = rotation_sequence,
    rotation_realloc = rotation_realloc,
    prop = prop,
    aggreg = aggreg,
    algo = algo,
    croptype_names = params@Croptypes$croptypeName,
    graphic = graphic,
    outputDir = params@OutputDir
  )
  params@Landscape <- st_as_sf(croptypeSP)
  params@Landscape$area <- data.frame(area = st_area(params@Landscape))
  if (length(orig_landscape$ID) != 0 && length(orig_landscape$Name) != 0) {
    params@Landscape$ID <- orig_landscape$ID
    params@Landscape$NAME <- orig_landscape$NAME
  }

  return(params)
}


#' @name loadPathogen
#' @title Load pathogen parameters
#' @description Loads default pathogen parameters for a specific disease
#' @details Available diseases:
#' * "rust"
#' @param disease a disease name (default: "rust")
#' @return a list of pathogen aggressiveness parameters on a susceptible host
#' for a pathogen genotype not adapted to resistance
#' @seealso \link{setPathogen}
#' @examples
#' basic_patho_params <- loadPathogen()
#' basic_patho_params
#' @export
loadPathogen <- function(disease = "rust") {
  patho <- switch(disease,
    "rust" = list(
      name = "rust",
      survival_prob = 1e-4,
      repro_sex_prob = 0,
      infection_rate = 0.4,
      propagule_prod_rate = 3.125,
      latent_period_exp = 10,
      latent_period_var = 9,
      infectious_period_exp = 24,
      infectious_period_var = 105,
      sigmoid_kappa = 5.333,
      sigmoid_sigma = 3,
      sigmoid_plateau = 1
    ),
    list()
  )

  if (length(patho) == 0) {
    warning('Unknown type of disease: "', disease
            , '". Currently the only possible type is: "rust"')
  }

  return(patho)
}


#' @name setPathogen
#' @title Set the pathogen
#' @description Updates a LandsepiParams object with pathogen parameters
#' @details a set of parameters representative of rust fungi can be loaded via 
#' \code{\link{loadPathogen}}.
#' @param params a LandsepiParams Object.
#' @param patho_params a list of pathogen aggressiveness parameters on a susceptible host
#' for a pathogen genotype not adapted to resistance: \itemize{
#' \item infection_rate = maximal expected infection rate of a propagule on a healthy host,
#' \item propagule_prod_rate = maximal expected effective propagule production rate of an 
#' infectious host per timestep,
#' \item latent_period_exp = minimal expected duration of the latent period,
#' \item latent_period_var = variance of the latent period duration,
#' \item infectious_period_exp = maximal expected duration of the infectious period,
#' \item infectious_period_var = variance of the infectious period duration,
#' \item survival_prob = probability for a propagule to survive the off-season,
#' \item repro_sex_prob = probability for an infectious host to reproduce via sex rather 
#' than via cloning,
#' \item sigmoid_kappa = kappa parameter of the sigmoid contamination function,
#' \item sigmoid_sigma = sigma parameter of the sigmoid contamination function,
#' \item sigmoid_plateau = plateau parameter of the sigmoid contamination function.
#' }
#' It can be generated manually, or, alternatively, via \code{\link{loadPathogen}}.
#' @return a LandsepiParams object
#' @seealso \link{loadPathogen}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setPathogen(simul_params, loadPathogen())
#' simul_params@Pathogen
#' }
#' @export
setPathogen <- function(params, patho_params) {
  params@Pathogen <- patho_params
  checkPathogen(params)

  return(params)
}


#' @name checkPathogen
#' @title Check pathogen
#' @description Checks papthogen validity
#' @param params a LandsepiParams Object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkPathogen <- function(params) {
  
  ret <- TRUE
  if (length(params@Pathogen) == 0 ||
      sum( sapply(params@Pathogen, length) != rep(1,length(params@Pathogen)) ) > 0 ){
    warning("Invalid parameters for Pathogen, use setPathogen()")
    ret <- FALSE
    return(ret)
  }

  if ( !is.numeric(params@Pathogen$infection_rate) ||
       !is.in.01(params@Pathogen$infection_rate) ){
    warning("Infection rate must be between 0 and 1")
    ret <- FALSE
  }
  if ( !is.numeric(params@Pathogen$propagule_prod_rate) ||
       !is.positive(params@Pathogen$propagule_prod_rate) ){
    warning("Propagule production rate must be >= 0")
    ret <- FALSE
  }
  if (!is.numeric(params@Pathogen$survival_prob) ||
      !is.in.01(params@Pathogen$survival_prob) ){
    warning("Survival probability must be between 0 and 1")
    ret <- FALSE
  }
  if (!is.numeric(params@Pathogen$repro_sex_prob) || 
      params@Pathogen$repro_sex_prob < 0 || 
      params@Pathogen$repro_sex_prob > 1){
    warning("Probability of sexual reproduction must be between 0 and 1")
    ret <- FALSE
  }
  if (!is.numeric(params@Pathogen$sigmoid_plateau) || 
      !is.in.01(params@Pathogen$sigmoid_plateau) ){
    warning("sigmoid_plateau must be between 0 and 1")
    ret <- FALSE
  }
  if (!is.numeric(params@Pathogen$latent_period_exp) ||
      !is.numeric(params@Pathogen$latent_period_var) ||
      !is.numeric(params@Pathogen$infectious_period_exp) ||
      !is.numeric(params@Pathogen$infectious_period_var) ||
      !is.numeric(params@Pathogen$sigmoid_kappa) ||
      !is.numeric(params@Pathogen$sigmoid_sigma) ||
      
      !is.positive(params@Pathogen$latent_period_exp) || 
      !is.positive(params@Pathogen$infectious_period_exp) ||
      !is.positive(params@Pathogen$latent_period_var) ||
      !is.positive(params@Pathogen$infectious_period_var) ||
      !is.positive(params@Pathogen$sigmoid_sigma) || 
      !is.positive(params@Pathogen$sigmoid_kappa) ){
    warning("Latent period, infectious period and sigmoid parameters must be >= 0")
    ret <- FALSE
  }

  return(ret)
}



#' @name loadCroptypes
#' @title Load Croptypes
#' @description Creates a data.frame containing croptype parameters and filled with 0
#' @param params a LandsepiParams Object.
#' @param croptypeIDs a vector of indices of croptypes (must match with croptypes in the landscape)
#' @param names a vector containing the names of all croptypes
#' @details Croptypes need to be later updated with \code{\link{allocateCroptypeCultivars}}.
#' If neither croptypeIDs nor names are given, it will automatically generate
#' 1 croptype per cultivar.
#' @return a data.frame with croptype parameters
#' @seealso \link{setCroptypes}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Mixture"))
#' croptypes
#' }
#' @export
loadCroptypes <- function(params, croptypeIDs = NULL, names = NULL) {
  cultivar_names <- params@Cultivars$cultivarName
  Ncultivars <- length(cultivar_names)

  if (is.null(croptypeIDs) & is.null(names)) {
    croptypeIDs <- 0:(Ncultivars - 1)
  }

  if (is.null(croptypeIDs)) {
    Ncroptypes <- length(names)
    croptypeIDs <- 0:(Ncroptypes - 1)
  } else if (is.null(names)) {
    Ncroptypes <- length(croptypeIDs)
    names <- paste("Crop", 1:Ncroptypes)
  }


  Ncroptypes <- length(croptypeIDs)
  cropt <- list(
    croptypeID = croptypeIDs,
    croptypeName = names
  )
  prop_tmp <- rep(0, Ncultivars)
  names(prop_tmp) <- c(cultivar_names)
  cropt <- c(cropt, prop_tmp)

  cropt <- as.data.frame(cropt, stringsAsFactors = FALSE)
  return(cropt)
}


#' @name allocateCroptypeCultivars
#' @title Allocate cultivars to one croptype
#' @description Updates a given croptype by allocating cultivars composing it.
#' @param croptypes a dataframe containing all croptypes, initialised via 
#' \code{\link{loadCroptypes}}
#' @param croptypeName the name of the croptype to be allocated
#' @param cultivarsInCroptype name of cultivars composing the croptype
#' @param prop vector of proportions of each cultivar in the croptype. Default to 
#' balanced proportions.
#' @return a croptype data.frame updated for the concerned croptype.
#' @seealso \link{setCroptypes}, \link{setCultivars} 
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Mixture"))
#' croptypes
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Resistant1", "Resistant2"))
#' croptypes
#' }
#' @export
allocateCroptypeCultivars <- function(croptypes, croptypeName, cultivarsInCroptype, prop = NULL) {
  n <- length(cultivarsInCroptype) ## number of cultivars composing the croptype
  if (is.null(prop)) {
    prop <- rep(1 / n, n)
  } else if (length(prop) == 1) {
    prop <- rep(prop, n)
  }

  for (k in 1:n) {
    croptypes[croptypes$croptypeName == croptypeName, cultivarsInCroptype[k]] <- prop[k]
  }

  return(croptypes)
}



#' @name setCroptypes
#' @title Set croptypes
#' @description Updates a LandsepiParams object with croptypes and their composition with regard 
#' to cultivar proportions
#' @details
#' The data.frame for cultivar allocations into croptypes must take this format (example):
#'
#' | croptypeID | croptypeName  | cultivarName1 | cultivarName2 | ... |
#' | ---------- | ------------- | ------------- | ------------- | --- |
#' | 0          |  "cropt1"     |  1            | 0             | ... |
#' | 1          |  "cropt2"     |  0.5          | 0.5           | ... |
#'
#' croptypeIDs have to match values from landscape "croptypeID" layer with feature year_X. 
#' Cultivars names have to match cultivar names in the cultivars data.frame.
#'
#' @param params a LandsepiParams Object.
#' @param dfCroptypes a data.frame containing cultivar proportions in each croptype (see details). 
#' It can be generated manually, or initialised with \code{\link{loadCroptypes}} and later 
#' updated with \code{\link{allocateCroptypeCultivars}}.
#' @return a LandsepiParams object
#' @seealso \link{loadCroptypes}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#' cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop", "Mixture"))
#' croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#' croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Resistant1", "Resistant2"))
#' simul_params <- setCroptypes(simul_params, croptypes)
#' simul_params@Croptypes
#' }
#' @export
setCroptypes <- function(params, dfCroptypes) {

  # no croptypeID and croptypeName
  if (is.null(dfCroptypes$croptypeID) && is.null(dfCroptypes$croptypeName)) {
    warning("Can't find croptype ID or Name in the data.frame")
    return(params)
  }
  else {
    # Try to add Name
    if (!is.null(dfCroptypes$croptypeID)) {
      rownames(dfCroptypes) <- dfCroptypes$croptypeID
      if (length(params@Landscape$Name) != 0 && length(params@Landscape$ID) != 0) {
        id_name <- as.data.frame(params@Landscape[, c("ID", "Name")], stringsAsFactors = FALSE)
        id_name$geometry <- NULL
        id_name <- unique(id_name)
        id_name <- id_name[which(dfCroptypes$croptypeID == id_name$ID), "Name"]
        dfCroptypes <- data.frame(croptypeName = id_name, dfCroptypes, stringsAsFactors = FALSE)
      }
    }
    else {
      # Try to add ID
      if (length(params@Landscape$ID) != 0 && length(params@Landscape$Name) != 0) {
        id_name <- as.data.frame(params@Landscape[, c("ID", "Name")], stringsAsFactors = FALSE)
        id_name$geometry <- NULL
        id_name <- unique(id_name)
        id_name <- id_name[which(dfCroptypes$croptypeName == id_name$Name), "ID"]
        dfCroptypes <- data.frame(croptypeID = id_name, dfCroptypes, stringsAsFactors = FALSE)
        rownames(dfCroptypes) <- dfCroptypes$croptypeID
      }
      else {
        warning("Can't retrieve croptypeID from croptypeName")
      }
    }
  }

  params@Croptypes <- data.frame(dfCroptypes, stringsAsFactors = FALSE)

  checkCroptypes(params)

  return(params)
}


#' @name checkCroptypes
#' @title Check croptypes
#' @description checks croptypes validity
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkCroptypes <- function(params) {
  ret <- TRUE

  if (nrow(params@Croptypes) == 0) {
    message("Croptypes data.frame undef")
    ret <- FALSE
  }

  # check cultivars proportion by croptypes
  lcultivars <- as.matrix(params@Croptypes[, -which(.croptypesColNames 
                                                    %in% colnames(params@Croptypes))])
  ret_tmp <- apply(params@Croptypes[, -which("croptypeName" == colnames(params@Croptypes))],
    MARGIN = 1,
    FUN = function(l) {
      if (sum(as.numeric(l[-1])) != 0 && sum(as.numeric(l[-1])) != 1.0) {
        message("Croptypes ", l[1], " have a proportion of cultivars not egal to 1")
        return(FALSE)
      }
      else {
        return(TRUE)
      }
    }
  )
  if (sum(!ret_tmp) > 0) ret <- FALSE

  if (nrow(params@Landscape) > 0) {
    lc <- as.data.frame(params@Landscape)[, grep("^year_", colnames(params@Landscape))]
    if (length(lc) > 0 &
        sum(!params@Croptypes$croptypeID %in% unique(unlist(lc))) != 0) {
      ret <- FALSE
      message("croptypeID from Croptypes not found in the landscape")
    }
  }

  if (nrow(params@Cultivars) > 0) {
    if ((ncol(params@Croptypes) - length(which(.croptypesColNames 
                                               %in% colnames(params@Croptypes)))) 
        > nrow(params@Cultivars)) {
      message("Croptypes have more Cultivars than those defined in Cultivars data.frame")
      ret <- FALSE
    }
    if (length(which(colnames(params@Croptypes)[-which(.croptypesColNames 
                                                       %in% colnames(params@Croptypes))] 
                     %in% rownames(params@Cultivars)))
    != (ncol(params@Croptypes) - length(which(.croptypesColNames 
                                              %in% colnames(params@Croptypes))))) {
      message("Cultivars in Croptypes data.frame are undef in Cultivars data.frame")
      ret <- FALSE
    }
  }

  return(ret)
}


#' @name loadCultivar
#' @title Load a cultivar
#' @description create a data.frame containing cultivar parameters depending of his type
#' @param name a character string (without space) specifying the cultivar name.
#' @param type the cultivar type: "growingHost", "nongrowingHost" or "nonhost" 
#' (default = "nonhost").
#' @details 
#' * "growingHost" is adapted to situations where the infection unit is a piece of leaf 
#' (e.g. where a fungal lesion can develop); the number of available infection units 
#' increasing during the season due to plant growth. 
#' * "nongrowingHost" corresponds to situations where the infection unit is the whole plant 
#' (e.g. for viral systemic infection); thus the number of infection units is constant. 
#' * "nonCrop" is not planted, does not cost anything and does not yield anything 
#' (e.g. forest, fallow).
#' @return a dataframe of parameters associated with each host genotype 
#' (i.e. cultivars, lines) when cultivated in pure crops.
#' @seealso \link{setCultivars}
#' @examples
#' c1 <- loadCultivar("winterWheat", type = "growingHost")
#' c1
#' c2 <- loadCultivar("forest", type = "nonhost")
#' c2
#' @export
loadCultivar <- function(name, type = "growingHost") {
  culti <- switch(type,
    "growingHost" = list(
      "cultivarName" = name,
      "initial_density" = 0.1,
      "max_density" = 2.0,
      "growth_rate" = 0.1,
      "reproduction_rate" = 0.0,
      "death_rate" = 0.0,
      "yield_H" = 2.5,
      "yield_L" = 0.0,
      "yield_I" = 0.0,
      "yield_R" = 0.0,
      "production_cost" = 225,
      "market_value" = 200
    ),
    "nongrowingHost" = list(
      "cultivarName" = name,
      "initial_density" = 2.0,
      "max_density" = 2.0,
      "growth_rate" = 0.0,
      "reproduction_rate" = 0.0,
      "death_rate" = 0.0,
      "yield_H" = 2.5,
      "yield_L" = 0.0,
      "yield_I" = 0.0,
      "yield_R" = 0.0,
      "production_cost" = 225,
      "market_value" = 200
    ),
    "nonCrop" = list(
      "cultivarName" = name,
      "initial_density" = 0.0,
      "max_density" = 2.0,
      "growth_rate" = 0.0,
      "reproduction_rate" = 0.0,
      "death_rate" = 0.0,
      "yield_H" = 0.0,
      "yield_L" = 0.0,
      "yield_I" = 0.0,
      "yield_R" = 0.0,
      "production_cost" = 0,
      "market_value" = 0
    ),
    list()
  )

  culti <- as.data.frame(culti, stringsAsFactors = FALSE)
  if (length(culti) == 0) {
    warning('Unknown type of host: "', type
            , '". Possible types are: "growingHost", "nongrowingHost", "nonCrop"')
  } else {
    # To be sure of the columns names
    colnames(culti) <- .cultivarsColNames
  }

  return(culti)
}


#' @name setCultivars
#' @title Set cultivars
#' @description Updates a LandsepiParams object with cultivars parameters
#' @param params a landsepiParams object.
#' @param dfCultivars a data.frame defining the cultivars (see details). It can be generated 
#' manually or, alternatively, via \code{\link{loadCultivar}}.
#'
#' @details dfCultivars is a dataframe of parameters associated with each host genotype 
#' (i.e. cultivars, lines) when cultivated in pure crops. Columns of the dataframe are:\itemize{
#' \item cultivarName: cultivar names (cannot accept space),
#' \item initial_density: host densities (per square meter) at the beginning of the cropping season,
#' \item max_density: maximum host densities (per square meter) at the end of the cropping season,
#' \item growth rate: host growth rates,
#' \item reproduction rate: host reproduction rates,
#' \item death rate: host death rates,
#' \item yield_H: yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status H,
#' \item yield_L: yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status L,
#' \item yield_I: yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status I,
#' \item yield_R: yield (in weight or volume units / ha / cropping season) 
#' associated with hosts in sanitary status R,
#' \item production_cost = overall production costs (in monetary units / ha / cropping season)
#' including planting costs, amortisation, labour etc.,
#' \item market_value = market values of the productions (in monetary units / weight or volume unit).
#' }
#' 
#' The data.frame must be defined as follow (example):
#' 
#' | cultivarName | initial_density | max_density | growth_rate | reproduction_rate | death_rate | yield_H | yield_L | yield_I |yield_R | production_cost | market_value |
#' | ------------ | --------------- | ----------- | ----------- | ----------------- | ---------- | ------- | ------- | ------- | ------ | --------------- | ------------ |
#' | Susceptible  | 0.1             |  2.0        | 0.1         | 0.0               | 0.0        | 2.5     | 0.0     | 0.0     | 0.0    | 225             | 200          |
#' | Resistant1   | 0.1             |  2.0        | 0.1         | 0.0               | 0.0        | 2.5     | 0.0     | 0.0     | 0.0    | 225             | 200          |
#' | Resistant2   | 0.1             |  2.0        | 0.1         | 0.0               | 0.0        | 2.5     | 0.0     | 0.0     | 0.0    | 225             | 200          |
#'
#' @return a LandsepiParams object
#' @seealso \link{loadCultivar}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' simul_params@Cultivars
#' }
#' @export
setCultivars <- function(params, dfCultivars) {
  if (!is.null(dfCultivars$cultivarName)) {
    rownames(dfCultivars) <- dfCultivars$cultivarName
  }

  params@Cultivars <- data.frame(dfCultivars[, .cultivarsColNames], stringsAsFactors = FALSE)

  checkCultivars(params)

  return(params)
}


#' @name checkCultivars
#' @title Check cultivars
#' @description check cultivars validity
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkCultivars <- function(params) {

  ret <- TRUE

  if (is.null(params@Cultivars) || nrow(params@Cultivars) == 0) {
    warning("Cultivars is NULL, use setCultivars()")
    ret <- FALSE
    return(ret)
  }

  if (sum(.cultivarsColNames %in% colnames(params@Cultivars)) != length(.cultivarsColNames)) {
    warning("Missing columns in Cultivars data.frame : ", .cultivarsColNames)
    ret <- FALSE
    return(ret)
  }
  
  if (!is.character(params@Cultivars$cultivarName) ||
      sum(grepl(" ", params@Cultivars$cultivarName)) > 0){
    warning("Cultivar names must be character strings without spaces")
    ret <- FALSE
  }
  
  if (!is.numeric(params@Cultivars$growth_rate) ||
      !is.numeric(params@Cultivars$reproduction_rate) ||
      !is.numeric(params@Cultivars$death_rate) ||
      
      sum(!is.in.01(params@Cultivars$growth_rate) > 0) || 
      sum(!is.in.01(params@Cultivars$reproduction_rate) > 0) || 
      sum(!is.in.01(params@Cultivars$death_rate) > 0) ){
    warning("growth, reproduction and death rates must be between 0 and 1")
    ret <- FALSE
  } 
  
  if(!is.numeric(params@Cultivars$initial_density) ||
     !is.numeric(params@Cultivars$yield_H) ||
     !is.numeric(params@Cultivars$yield_L) ||
     !is.numeric(params@Cultivars$yield_I) ||
     !is.numeric(params@Cultivars$yield_R) ||
     !is.numeric(params@Cultivars$production_cost) ||
     !is.numeric(params@Cultivars$market_value) ||
     
     sum( !is.positive(params@Cultivars$initial_density) ) > 0 ||
     sum( !is.positive(params@Cultivars$yield_H) ) > 0 ||
     sum( !is.positive(params@Cultivars$yield_L) ) > 0 ||
     sum( !is.positive(params@Cultivars$yield_I) ) > 0 ||
     sum( !is.positive(params@Cultivars$yield_R) ) > 0 ||
     sum( !is.positive( params@Cultivars$production_cost) ) > 0 ||
     sum( !is.positive(params@Cultivars$market_value) ) > 0 ){
    warning("initial_density, yield, production_cost and market_value must be >= 0")
    ret <- FALSE
  }
  
  if(!is.numeric(params@Cultivars$max_density) ||
     sum( !is.strict.positive(params@Cultivars$max_density) ) > 0 ||
     sum( params@Cultivars$max_density < params@Cultivars$initial_density ) > 0 ){
    warning("Maximal density must be strictly positive and greater or equal to initial density")
    ret <- FALSE
  }
  

  if (ncol(params@Croptypes) > 0) {
    if (nrow(params@Cultivars) 
        < (length(params@Croptypes[, -which(.croptypesColNames 
                                            %in% colnames(params@Croptypes))]) - 1)) {
      warning("Cultivars number is less than thoses defined in Croptypes data.frame")
      ret <- FALSE
    }
    if (length(which(rownames(params@Cultivars) 
                     %in% colnames(params@Croptypes)[
                       -which(.croptypesColNames %in% colnames(params@Croptypes))])) 
        != (ncol(params@Croptypes) - length(which(.croptypesColNames
                                                  %in% colnames(params@Croptypes))))) {
      warning("Cultivars are undef in Cultivars data.frame compared to croptypes")
      ret <- FALSE
    }
  }

  return(ret)
}


#' @name loadGene
#' @title Load a gene
#' @description Creates a data.frame containing parameters of a gene depending of his type
#' @param name name of the gene
#' @param type type of the gene: "majorGene", "APR", "QTL" or "immunity" (default = "majorGene")
#' @details 
#' * "majorGene" means a completely efficient gene that can be broken down via a single 
#' pathogen mutation
#' * "APR" means a major gene that is active only after a delay of 30 days after planting
#' * "QTL" means a partial resistance (50% efficiency) that requires several pathogen mutations 
#' to be completely eroded
#' * "immunity" means a completely efficient resistance that the pathogen has no way to adapt 
#' (i.e. the cultivar is nonhost).  
#' 
#' For different scenarios, the data.frame can be manually updated later.
#' @return a data.frame with gene parameters
#' @seealso \link{setGenes}
#' @examples
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene1
#' gene2 <- loadGene(name = "Lr34", type = "APR")
#' gene2
#' @export
loadGene <- function(name, type = "majorGene") {
  gene <- switch(type,
    "majorGene" = list(
      "geneName" = name,
      "efficiency" = 1.0,
      "time_to_activ_exp" = 0.0,
      "time_to_activ_var" = 0.0,
      "mutation_prob" = 0.0000001,
      "Nlevels_aggressiveness" = 2,
      "fitness_cost" = 0.5,
      "tradeoff_strength" = 1.0,
      "target_trait" = "IR"
    ),
    "APR" = list(
      "geneName" = name,
      "efficiency" = 1.0,
      "time_to_activ_exp" = 30.0,
      "time_to_activ_var" = 30.0,
      "mutation_prob" = 0.0000001,
      "Nlevels_aggressiveness" = 2,
      "fitness_cost" = 0.5,
      "tradeoff_strength" = 1.0,
      "target_trait" = "IR"
    ),
    "QTL" = list(
      "geneName" = name,
      "efficiency" = 0.5,
      "time_to_activ_exp" = 0.0,
      "time_to_activ_var" = 0.0,
      "mutation_prob" = 0.0001,
      "Nlevels_aggressiveness" = 6,
      "fitness_cost" = 0.5,
      "tradeoff_strength" = 1.0,
      "target_trait" = "IR"
    ),
    "immunity" = list(
      "geneName" = name,
      "efficiency" = 1.0,
      "time_to_activ_exp" = 0.0,
      "time_to_activ_var" = 0.0,
      "mutation_prob" = 0,
      "Nlevels_aggressiveness" = 1,
      "fitness_cost" = 0,
      "tradeoff_strength" = 1,
      "target_trait" = "IR"
    ),
    list()
  )

  gene <- as.data.frame(gene, stringsAsFactors = FALSE)
  if (length(gene) == 0) {
    warning('Unknown type of gene: "', type
            , '". Possible types are: "majorGene", "APR", "QTL", "immunity')
  } else {
    # To be sure of the columns names
    colnames(gene) <- .geneColNames
  }
  return(gene)
}


#' @name setGenes
#' @title Set genes
#' @description Updates a LandsepiParams object with parameters associated with resistance genes
#' and pathogen adaptation.
#' @details dfGenes is a data.frame of parameters associated with each resistance gene and 
#' with the evolution of each corresponding pathogenicity gene. Columns of the dataframe are:
#' \itemize{
#' \item geneName: names of resistance genes,
#' \item target_trait: aggressiveness components (IR, LAT, IP, or PR) targeted by resistance genes,
#' \item efficiency: resistance gene efficiencies, i.e. the percentage of reduction of the targeted 
#' aggressiveness component (IR, 1/LAT, IP and PR),
#' \item time_to_activ_exp: expected delays to resistance activation (for APRs),
#' \item time_to_activ_var: variances of the delay to resistance activation (for APRs),
#' \item mutation_prob: mutation probabilities for pathogenicity genes (each of them 
#' corresponding to a resistance gene),
#' \item Nlevels_aggressiveness: number of adaptation levels related to each resistance gene 
#' (i.e. 1 + number of required mutations for a pathogenicity gene to fully adapt to the 
#' corresponding resistance gene),
#' \item fitness_cost: fitness penalties paid by pathogen genotypes 
#' fully adapted to the considered resistance genes on host that do not carry these genes, 
#' \item tradeoff_strength: strengths of the trade-off relationships between the 
#' level of aggressiveness on hosts that do and do not carry the resistance genes.
#' }
#' 
#' The data.frame must be defined as follow (example):
#'
#' | geneName | efficiency | time_to_activ_exp | time_to_activ_var | mutation_prob | Nlevels_agressiveness | fitness_cost | tradeoff_strength | target_trait |
#' | -------- | ---------- | ----------------- | ----------------- | ------------- | --------------------- | ------------ | ----------------- | ------------ |
#' | MG1      |  1         |  0                | 0                 | 1e-07         | 2                     | 0.5          | 1                 | IR           |
#' | QTL1     | 0.5        |  0                | 0                 | 0.0001        | 10                    | 0.74         | 1                 | LAT          |
#'
#' @param params a LandsepiParams object
#' @param dfGenes a data.frame containing gene parameters. It can be defined manually, or, 
#' alternatively, with \code{\link{loadGene}}.
#' @return a LandsepiParams object.
#' @seealso \link{loadGene}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' simul_params@Genes
#' }
#' @export
setGenes <- function(params, dfGenes) {
  if (!is.null(dfGenes$geneName)) {
    rownames(dfGenes) <- dfGenes$geneName
  }
  params@Genes <- dfGenes[, .geneColNames]

  checkGenes(params)

  return(params)
}


#' @name checkGenes
#' @title Check genes
#' @description checks Genes data.frame validity
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkGenes <- function(params) {
  ret <- TRUE

  if (nrow(params@Genes) == 0) {
    warning("Simulation with no resistance gene")
    return(ret)
  }

  if (is.null(params@Genes$geneName)) warning("missing 'geneName' column into genes data.frame")

  if (sum(!.geneColNames %in% colnames(params@Genes)) > 0) {
    warning("Genes data.frame column(s) missing")
    warning("Genes colnames are ", .geneColNames)
    ret <- FALSE
    return(ret)
  }
  
  validTraits <- c("IR","LAT","PR","IP")
  if(!is.character(params@Genes$target_trait) ||
     is.na( sum(match(params@Genes$target_trait, validTraits)) ) ){
    warning( "Error: valid target traits are:", paste(validTraits, collapse = ", ") )
    ret <- FALSE
  }

  if (!is.numeric(params@Genes$efficiency) ||
      !is.numeric(params@Genes$mutation_prob) ||
      !is.numeric(params@Genes$fitness_cost) ||
      
      sum(!is.in.01(params@Genes$efficiency) > 0) || 
      sum(!is.in.01(params@Genes$mutation_prob) > 0) || 
      sum(!is.in.01(params@Genes$fitness_cost) > 0) ){
    warning("efficiencies, mutation probabilities and fitness costs must be between 0 and 1")
    ret <- FALSE
  } 
  
  if(!is.numeric(params@Genes$time_to_activ_exp) ||
     !is.numeric(params@Genes$time_to_activ_var) ||
     
     sum(!is.positive(params@Genes$time_to_activ_exp) > 0) ||
     sum(!is.positive(params@Genes$time_to_activ_var) > 0) ){
    warning("Expectation and variance of the times to resistance activation must be >= 0")
    ret <- FALSE
  }
  
  if(!is.numeric(params@Genes$tradeoff_strength) ||
     sum(!is.strict.positive(params@Genes$tradeoff_strength) > 0) ){
    warning("tradeoff strengths must be > 0")
    ret <- FALSE
  }
  
  if(!is.numeric(params@Genes$Nlevels_aggressiveness) ||
     sum( !is.wholenumber(params@Genes$Nlevels_aggressiveness) > 0) ||
     sum( !is.strict.positive(params@Genes$Nlevels_aggressiveness) > 0) ){
    warning("Number of levels of aggressiveness must be a whole number >= 1")
    ret <- FALSE
  }
  
  return(ret)
}


#' @name allocateCultivarGenes
#' @title Allocate genes to a cultivar
#' @description Updates a LandsepiParams object with, for a given cultivar, the list of genes 
#' it carries
#' @param params a LandsepiParams object.
#' @param cultivarName the name of the cultivar to be allocated.
#' @param listGenesNames the names of the genes the cultivar carries
#' @return a LandsepiParams object
#' @seealso \link{setGenes}, \link{setCultivars}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' gene1 <- loadGene(name = "MG 1", type = "majorGene")
#' gene2 <- loadGene(name = "MG 2", type = "majorGene")
#' genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#' simul_params <- setGenes(simul_params, genes)
#' cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#' cultivar2 <- loadCultivar(name = "Resistant", type = "growingHost")
#' cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
#' simul_params <- setCultivars(simul_params, cultivars)
#' simul_params <- allocateCultivarGenes(simul_params, "Resistant", c("MG 1", "MG 2"))
#' simul_params@CultivarsGenes
#' }
#' @export
allocateCultivarGenes <- function(params, cultivarName, listGenesNames = c("")) {
  if (length(params@CultivarsGenes) == 0 
      || nrow(params@CultivarsGenes) != nrow(params@Cultivars) 
      || nrow(params@Genes) != ncol(params@CultivarsGenes)) {
    params@CultivarsGenes <- data.frame(matrix(rep(0, nrow(params@Genes) * nrow(params@Cultivars))
                                               , nrow = nrow(params@Cultivars))
                                        , row.names = rownames(params@Cultivars))
    colnames(params@CultivarsGenes) <- params@Genes$geneName
  }

  if (cultivarName %in% params@Cultivars$cultivarName &&
    sum(listGenesNames %in% params@Genes$geneName) == length(listGenesNames)) {
    params@CultivarsGenes[which(rownames(params@CultivarsGenes) 
                                == cultivarName), listGenesNames] <- 1
    params@CultivarsGenes[which(rownames(params@CultivarsGenes) 
                                != cultivarName), listGenesNames] <- 0
  } else {
    stop("Can't find cultivarName or geneName from data.frame")
  }
  return(params)
}


#' @name resetCultivarsGenes
#' @title Reset cultivars genes
#' @description Resets the lists of genes carried by all cultivars
#' @param params a LandsepiParams object.
#' @return a LandsepiParams object
resetCultivarsGenes <- function(params) {
  ## (used in shinyApp)
  params@CultivarsGenes <- data.frame()
  return(params)
}


#' @name checkCultivarsGenes
#' @title Check cultivars genes
#' @description Checks CultivarsGene data.frame validity
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkCultivarsGenes <- function(params) {
  ret <- TRUE
  if (length(params@CultivarsGenes) > 0 & 
      ( nrow(params@CultivarsGenes) != nrow(params@Cultivars) || 
       nrow(params@Genes) != ncol(params@CultivarsGenes)) ) {
    warning("Cultivars Genes undef (some genes are not allocated to any cultivar)")
    ret <- FALSE
  }
  return(ret)
}

#' @name setInoculum
#' @title Set inoculum
#' @description Updates a LandsepiParams with the initial probability for the first host 
#' (whose index is 0) to be infectious (i.e. state I) at the beginning of the simulation.
#' @param params a LandsepiParams object.
#' @param val a numeric value (default = 5e-4). Must be between 0 and 1.
#' @return a LandsepiParams object
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setInoculum(simul_params, 1E-3)
#' simul_params@PI0
#' }
#' @export
setInoculum <- function(params, val = 5e-4) {
  params@PI0 <- val
  checkInoculum(params)
  
  return(params)
}

#' @name checkInoculum
#' @title Check inoculum
#' @description Checks inoculum validity.
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkInoculum <- function(params) {
  ret <- TRUE
  
  if( !is.numeric(params@PI0) ||
      length(params@PI0) != 1 ||
      !is.in.01(params@PI0) ) {
    warning("Invalid inoculum value: must be in [0,1]")
    ret <- FALSE
  }
  return(ret)
}

#' @name loadOutputs
#' @title Load outputs
#' @description Creates an output list
#' @param epid_outputs a character string (or a vector of character strings if several outputs 
#' are to be computed) specifying the type of epidemiological and economic outputs to generate 
#' (see details):\itemize{
#' \item "audpc" : Area Under Disease Progress Curve (average proportion of diseased hosts relative 
#' to the carryng capacity) 
#' \item "gla_abs" : Absolute Green Leaf Area (average number of healthy hosts per square meter)
#' \item "gla_rel" : Relative Green Leaf Area (average proportion of healthy hosts relative to the 
#' total number of existing hosts)
#' \item "eco_cost" : total crop costs (in weight or volume units per ha) 
#' \item "eco_product" : total crop production (in monetary units per ha) 
#' \item "eco_benefit" : total crop benefits (in monetary units per ha) 
#' \item "eco_grossmargin" : Gross Margin (benefits - costs, in monetary units per ha) 
#' \item "HLIR_dynamics", "H_dynamics", "L_dynamics", "IR_dynamics", "HLI_dynamics", etc.: 
#' Epidemic dynamics related to the specified sanitary status (H, L, I or R and all their 
#' combinations). Graphics only, works only if graphic=TRUE.
#' \item "all" : compute all these outputs (default)
#' \item "" : none of these outputs will be generated.
#' }
#' @param evol_outputs a character string (or a vector of character strings if several 
#' outputs are to be computed) specifying the type of evolutionary outputs to generate :\itemize{
#' \item "evol_patho": Dynamics of pathogen genotype frequencies
#' \item "evol_aggr": Evolution of pathogen aggressiveness
#' \item "durability": Durability of resistance genes
#' \item "all": compute all these outputs (default)
#' \item "": none of these outputs will be generated.
#' }
#' @return a list of outputs and parameters for output generation
#' @seealso \link{setOutputs}
#' @examples 
#' outputList <- loadOutputs(epid_outputs = "audpc", evol_outputs = "durability")
#' outputList
#' @export
loadOutputs <- function(epid_outputs = "all", evol_outputs = "all"){
  outputList <- list(epid_outputs = epid_outputs
                     , evol_outputs = evol_outputs
                     , thres_breakdown = 50000
                     , GLAnoDis = 1.48315
                     , audpc100S = 0.38)
  return(outputList)
}


#' @name setOutputs
#' @title Set outputs
#' @description Updates a LandsepiParams object with a list of output parameters.
#' @param params a LandsepiParams object.
#' @param output_list a list of outputs to be generated and parameters for output generation. 
#' It can be generated manually or, alternatively, via \code{\link{loadOutputs}}. This list 
#' is composed of:\itemize{
#' \item epid_outputs = epidemiological outputs to compute (see details)
#' \item evol_outputs = evolutionary outputs to compute (see details)
#' \item thres_breakdown = an integer (or vector of integers) giving the threshold 
#' (i.e. number of infections) above which a pathogen genotype is unlikely to go extinct, 
#' used to characterise the time to invasion of resistant hosts (several values are computed 
#' if several thresholds are given in a vector).
#' \item GLAnoDis = the absolute Green Leaf Area in absence of disease (used to compute 
#' economic outputs).
#' \item audpc100S = the audpc in a fully susceptible landscape (used as reference value 
#' for graphics and video).
#' }
#' @details "epid_outputs" is a character string (or a vector of character strings if several 
#' outputs are to be computed) specifying the type of epidemiological and economic outputs 
#' to generate:  
#' \itemize{
#' \item "audpc" : Area Under Disease Progress Curve (average proportion of diseased hosts relative 
#' to the carryng capacity) 
#' \item "gla_abs" : Absolute Green Leaf Area (average number of healthy hosts per square meter)
#' \item "gla_rel" : Relative Green Leaf Area (average proportion of healthy hosts relative to the 
#' total number of existing hosts)
#' \item "eco_cost" : total crop costs (in weight or volume units per ha) 
#' \item "eco_product" : total crop production(in monetary units per ha) 
#' \item "eco_benefit" : total crop benefits (in monetary units per ha) 
#' \item "eco_grossmargin" : Gross Margin (benefits - costs, in monetary units per ha) 
#' \item "HLIR_dynamics", "H_dynamics", "L_dynamics", "IR_dynamics", "HLI_dynamics", etc.: 
#' Epidemic dynamics related to the specified sanitary status (H, L, I or R and all their 
#' combinations). Graphics only, works only if graphic=TRUE.
#' \item "all" : compute all these outputs (default)
#' \item "" : none of these outputs will be generated.
#' }
#' "evol_outputs" is a character string (or a vector of character strings if several outputs 
#' are to be computed) specifying the type of evolutionary outputs to generate :\itemize{
#' \item "evol_patho": Dynamics of pathogen genotype frequencies
#' \item "evol_aggr": Evolution of pathogen aggressiveness
#' \item "durability": Durability of resistance genes
#' \item "all": compute all these outputs (default)
#' \item "": none of these outputs will be generated.
#' }
#' 
#' @return a LandsepiParams object.
#' @seealso \link{loadOutputs}
#' @examples
#' \dontrun{
#' simul_params <- createSimulParams()
#' simul_params <- setOutputs(simul_params, loadOutputs())
#' simul_params@Outputs
#' }
#' @export
setOutputs <- function(params, output_list){
  params@Outputs <- output_list
  
  checkOutputs(params) 
  
  return(params)
}

#' @name checkOutputs
#' @title Check outputs
#' @description Checks outputs validity.
#' @param params a LandsepiParams object.
#' @return a boolean, TRUE if OK, FALSE otherwise
checkOutputs <- function(params) {
  ret <- TRUE
  
  if( !is.character(params@Outputs$epid_outputs) ||
      !is.character(params@Outputs$evol_outputs)) {
    warning("Invalid epidemiological or evolutionary outputs")
    ret <- FALSE
  }
  
  if (!is.na(params@Outputs$audpc100S)){
    if( !is.numeric(params@Outputs$audpc100S) ||
        !is.in.01(params@Outputs$audpc100S, exclude0 = TRUE) ) {
      warning("AUDPC in a fully susceptible landscape must be in ]0;1]")
      ret <- FALSE
    }
  }
  if (!is.na(params@Outputs$GLAnoDis)){
    if( !is.numeric(params@Outputs$GLAnoDis) || 
        !is.strict.positive(params@Outputs$GLAnoDis) ) {
      warning("GLA in absence of disease must be > 0")
      ret <- FALSE
    }
  }
  
  if (!is.na(params@Outputs$thres_breakdown)){
    if( !is.wholenumber(params@Outputs$thres_breakdown) || 
        !is.strict.positive(params@Outputs$thres_breakdown) ) {
      warning("Threshold for resistance breakdown must be a whole number > 0")
      ret <- FALSE
    }
  }
  return(ret)
}

###### PRIVATE ######

#' params2CroptypeBDD
#' @description Converts a LandsepiParams object to a value compatible with BDD croptype Table
#' @param params a LandsepiParams object.
#' @return a data.frame BDD compatible
params2CroptypeBDD <- function(params) {
  if (is.null(params@Croptypes$croptypeName)) {
    croptypes <- params@Croptypes
  } else {
    croptypes <- params@Croptypes[, -which(colnames(params@Croptypes) == "croptypeName")]
  }

  colnames(croptypes) <- 
    c("croptypeID", which(rownames(params@Cultivars) 
                          %in% colnames(croptypes)[-which(.croptypesColNames 
                                                          %in% colnames(croptypes))]))
  res <- data.frame()
  for (i in 1:nrow(croptypes)) {
    for (culti in which(croptypes[i, -1] > 0)) { ## without croptypeID and croptypeName index
      # message(paste0(croptypes[i,1]," ",culti," ",croptypes[i,culti+1]))
      res <- rbind(res, c(croptypes[i, "croptypeID"]
                          , culti - 1, croptypes[i, which(colnames(croptypes) == culti)]))
    }
  }
  res <- cbind(0:(nrow(res) - 1), res)
  colnames(res) <- c("rowid", "croptypeID", "cultivarID", "proportion")
  return(res)
}

#' CroptypeBDD2Params
#' @description Converts BDD table to Croptype LandspeParams object
#' @param inputGPKG a GPKG filename
#' @return a data.frame LandsepiParams@@Croptypes compatible
CroptypeBDD2Params <- function(inputGPKG) {
  return(getGPKGCroptypes(inputGPKG))
}

#' params2CultivarBDD
#' @description Converts a LandsepiParams object to a value compatible with BDD cultivar Table
#' @param params a LandsepiParam object.
#' @return a data.frame BDD compatible
params2CultivarBDD <- function(params) {
  cultivars <- data.frame(params@Cultivars, stringsAsFactors = FALSE)
  cultivars <- data.frame("cultivarID" = 0:(nrow(cultivars) - 1), cultivars
                          , stringsAsFactors = FALSE)
  return(cultivars)
}

#' CultivarBDD2Params
#' @description Converts BDD table Cultivar to a Cultivar LandsepiParams object
#' @param inputGPKG a GPKG filename
#' @return a data.frame LandsepiParams@@Cultivars compatible
CultivarBDD2Params <- function(inputGPKG) {
  return(getGPKGCultivars(inputGPKG))
}

#' params2GeneBDD
#' @description Converts a LandsepiParams object to a value compatible with BDD Gene Table
#' @param params a LandsepiParam object.
#' @return a data.frame BDD compatible
params2GeneBDD <- function(params) {
  gene <- data.frame(params@Genes, stringsAsFactors = FALSE)
  gene <- data.frame("geneID" = 0:(nrow(gene) - 1), gene, stringsAsFactors = TRUE)
  return(gene)
}

#' GeneBDD2Params
#' @description Converts BDD table Gene to LandsepiParams object
#' @param inputGPKG a LandsepiParams
#' @return a data.frame LandsepiParams@@Genes compatible
GeneBDD2Params <- function(inputGPKG) {
  return(getGPKGGenes(inputGPKG))
}

#' params2GeneListBDD
#' @description Converts a LandsepiParams object to a value compatible with BDD cultivarsGenes Table
#' @param params a LandsepiParam object.
#' @return a data.frame BDD compatible
params2GeneListBDD <- function(params) {
  cgList <- params@CultivarsGenes
  res <- data.frame(stringsAsFactors = FALSE)
  for (i in 1:nrow(cgList)) {
    culti <- rownames(cgList)[i]
    for (gene in which(cgList[i, ] > 0)) {
      res <- rbind(res, c(which(rownames(params@Cultivars) == culti) - 1
                          , which(rownames(params@Genes) == colnames(cgList)[gene]) - 1))
    }
  }
  res <- cbind(0:(nrow(res) - 1), res)
  colnames(res) <- c("rowid", "cultivarID", "geneID")
  return(res)
}

#' CultivarGeneBDD2Params
#' @description Converts BDD table to LandsepiParams object
#' @param inputGPKG a GPKG filename
#' @return a data.frame LandsepiParams@@CultivarsGenes compatible
CultivarGeneBDD2Params <- function(inputGPKG) {
  return(getGPKGCultivarsGenes(inputGPKG))
}
