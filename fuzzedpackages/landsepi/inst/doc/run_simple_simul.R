## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- results="hide", message=FALSE-------------------------------------------
library(landsepi)

## ---- results="hide", message="FALSE"-----------------------------------------
simul_params <- createSimulParams(outputDir = getwd())

## -----------------------------------------------------------------------------
simul_params@Seed
simul_params <- setSeed(simul_params, seed = 1)
simul_params@Seed

## -----------------------------------------------------------------------------
simul_params <- setTime(simul_params, Nyears = 6, nTSpY = 120)
simul_params@TimeParam

## -----------------------------------------------------------------------------
basic_patho_param <- loadPathogen(disease = "rust")
basic_patho_param

## -----------------------------------------------------------------------------
basic_patho_param <- loadPathogen("rust")
basic_patho_param$infection_rate <- 0.5
basic_patho_param

## -----------------------------------------------------------------------------
basic_patho_param <- list(infection_rate = 0.4
                          , latent_period_exp = 10
                          , latent_period_var = 9
                          , propagule_prod_rate = 3.125
                          , infectious_period_exp = 24
                          , infectious_period_var = 105
                          , survival_prob = 1e-4
                          , repro_sex_prob = 0
                          , sigmoid_kappa = 5.333, sigmoid_sigma = 3, sigmoid_plateau = 1)

## -----------------------------------------------------------------------------
simul_params <- setPathogen(simul_params, patho_params = basic_patho_param)
simul_params@Pathogen

## -----------------------------------------------------------------------------
simul_params <- setInoculum(simul_params, val = 5e-4)
simul_params@PI0

## -----------------------------------------------------------------------------
landscape <- loadLandscape(id = 1)
length(landscape)
plot(landscape, main = "Landscape structure")

## -----------------------------------------------------------------------------
disp_patho <- loadDispersalPathogen(id = 1)
head(disp_patho)
length(landscape)^2 == length(disp_patho)

## -----------------------------------------------------------------------------
simul_params <- setLandscape(simul_params, land = landscape)
simul_params <- setDispersalPathogen(simul_params, mat = disp_patho)

## -----------------------------------------------------------------------------
cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
cultivar4 <- loadCultivar(name = "Resistant3", type = "nongrowingHost")
cultivar5 <- loadCultivar(name = "Forest", type = "nonCrop")
cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3, cultivar4, cultivar5)
                        , stringsAsFactors = FALSE)
cultivars

## -----------------------------------------------------------------------------
cultivars[cultivars$cultivarName == "Susceptible", "growth_rate"] <- 0.2
cultivars

## -----------------------------------------------------------------------------
cultivars_new <- data.frame(cultivarName = c("Susceptible", "Resistant"),
                            initial_density =   c(0.1, 0.2),
                            max_density =       c(2.0, 3.0),
                            growth_rate =       c(0.1, 0.2),
                            reproduction_rate = c(0.0, 0.0),
                            death_rate =        c(0.0, 0.0),
                            yield_H =           c(2.5, 2.0),
                            yield_L =           c(0.0, 0.0),
                            yield_I =           c(0.0, 0.0),
                            yield_R =           c(0.0, 0.0),
                            production_cost =   c(225, 300),
                            market_value =      c(200, 150),
                            stringsAsFactors = FALSE)
cultivars_new

## -----------------------------------------------------------------------------
gene1 <- loadGene(name = "MG 1", type = "majorGene")
gene2 <- loadGene(name = "Lr34", type = "APR")
gene3 <- loadGene(name = "gene 3", type = "QTL")
gene4 <- loadGene(name = "nonhost resistance", type = "immunity")
genes <- data.frame(rbind(gene1, gene2, gene3, gene4), stringsAsFactors = FALSE)
genes

## -----------------------------------------------------------------------------
genes[genes$geneName == "MG 1", "mutation_prob"] <- 1e-3
genes

## -----------------------------------------------------------------------------
genes_new <- data.frame(geneName =               c("MG1", "MG2"),
                        efficiency =             c(1.0  , 0.8  ),
                        time_to_activ_exp =      c(0.0  , 0.0  ),
                        time_to_activ_var =      c(0.0  , 0.0  ),
                        mutation_prob =          c(1E-7 , 1E-4),
                        Nlevels_aggressiveness = c(2    , 2    ),
                        fitness_cost =           c(0.50 , 0.75 ),
                        tradeoff_strength =      c(1.0  , 1.0  ),
                        target_trait =           c("IR" , "LAT"),
                        stringsAsFactors = FALSE)
genes_new

## -----------------------------------------------------------------------------
simul_params <- setGenes(simul_params, dfGenes = genes)
simul_params <- setCultivars(simul_params, dfCultivars = cultivars)
simul_params@Genes
simul_params@Cultivars

## -----------------------------------------------------------------------------
simul_params <- allocateCultivarGenes(simul_params
                                      , cultivarName = "Resistant1"
                                      , listGenesNames = c("MG 1"))
simul_params <- allocateCultivarGenes(simul_params
                                      , cultivarName = "Resistant2"
                                      , listGenesNames = c("Lr34", "gene 3"))
simul_params <- allocateCultivarGenes(simul_params
                                      , cultivarName = "Resistant3"
                                      , listGenesNames = c("nonhost resistance"))
simul_params@Cultivars

## -----------------------------------------------------------------------------
croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
                                                   , "Pure resistant crop"
                                                   , "Mixture"
                                                   , "Other"))
croptypes

## -----------------------------------------------------------------------------
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Susceptible crop"
                                       , cultivarsInCroptype = "Susceptible")
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Pure resistant crop"
                                       , cultivarsInCroptype = "Resistant1")
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Mixture"
                                       , cultivarsInCroptype = c("Resistant2","Resistant3")
                                       , prop = c(0.4, 0.6))
croptypes <- allocateCroptypeCultivars(croptypes
                                       , croptypeName = "Other"
                                       , cultivarsInCroptype = "Forest")
croptypes

## -----------------------------------------------------------------------------
simul_params <- setCroptypes(simul_params, dfCroptypes = croptypes)
simul_params@Croptypes

## -----------------------------------------------------------------------------
croptypes <- data.frame(croptypeID = c(0, 1, 2, 3)
                        , croptypeName = c("Susceptible crop"
                                           , "Pure resistant crop"
                                           , "Mixture"
                                           , "Other")
                        , Susceptible = c(1,0,0  ,0)
                        , Resistant1  = c(0,1,0  ,0)
                        , Resistant2  = c(0,0,0.5,0)
                        , Resistant3  = c(0,0,0.5,0)
                        , Forest      = c(0,0,0  ,1)
                        , stringsAsFactors = FALSE)
simul_params <- setCroptypes(simul_params, croptypes)

## -----------------------------------------------------------------------------
# croptypeIDs cultivated in each element of the rotation sequence:
rotation_sequence <- list(c(0,1,3), c(0,2,3))
rotation_period <- 2  # number of years before rotation of the landscape
prop <- list(rep(1/3, 3), rep(1/3, 3))  # proportion (in surface) of each croptype
aggreg <- 1    # level of spatial aggregation
simul_params <- allocateLandscapeCroptypes(simul_params
                                           , rotation_period = rotation_period
                                           , rotation_sequence = rotation_sequence
                                           , prop = prop
                                           , aggreg = aggreg
                                           , graphic = FALSE)
# plot(simul_params@Landscape)

## -----------------------------------------------------------------------------
outputlist <- loadOutputs(epid_outputs = "all", evol_outputs = "all")
outputlist

## -----------------------------------------------------------------------------
simul_params <- setOutputs(simul_params, outputlist)

## ---- eval=FALSE--------------------------------------------------------------
#  checkSimulParams(simul_params)
#  simul_params <- saveDeploymentStrategy(simul_params)

## ---- eval=FALSE--------------------------------------------------------------
#  runSimul(simul_params, graphic = TRUE, videoMP4 = FALSE)

## ---- include=FALSE-----------------------------------------------------------
system(paste("rm -rf ", simul_params@OutputDir))

