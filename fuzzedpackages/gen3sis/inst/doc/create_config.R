## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(gen3sis)

## ----eval=FALSE---------------------------------------------------------------
#  # get data path
#  datapath <- system.file(file.path("extdata", "EmptyConfig"), package="gen3sis")
#  # set config_empty.R file path
#  config_file_path <- file.path(tempdir(), "config_empty.R")

## ----eval=FALSE, message=FALSE------------------------------------------------
#  #writes out a config skeleton
#  write_config_skeleton(config_file_path)

## ----eval=FALSE---------------------------------------------------------------
#  # Version: 1.0
#  # Author:
#  # Date:
#  # Landscape:
#  # Publications:
#  # Description:

## ----eval=FALSE---------------------------------------------------------------
#  # set the random seed for the simulation.
#  random_seed = NA
#  # set the starting time step or leave NA to use the earliest/highest time-step.
#  start_time = NA
#  # set the end time step or leave as NA to use the latest/lowest time-step (0).
#  end_time = NA
#  # maximum total number of species in the simulation before it is aborted.
#  max_number_of_species = 25000
#  # maximum number of species within one cell before the simulation is aborted.
#  max_number_of_coexisting_species = 2500
#  # a list of traits to include with each species.
#  # a "dispersal" trait is implicitly added in any case.
#  trait_names = c("dispersal")
#  # ranges to scale the input environments with:
#  # not listed variable:         no scaling takes place
#  # listed, set to NA:           env. variable will be scaled from [min, max] to [0, 1]
#  # listed with a given range r: env. variable will be scaled from [r1, r2] to [0, 1]
#  environmental_ranges = list( )

## ----eval=FALSE---------------------------------------------------------------
#  end_of_timestep_observer = function(data, vars, config){
#    # the list of all species can be found in data$all_species.
#    # the current landscape can be found in data$landscape.
#  }

## ----eval=FALSE---------------------------------------------------------------
#  # the initial abundance of a newly colonized cell, both during setup and later when
#  # colonizing a cell during the dispersal.
#  initial_abundance = 1
#  # place species in the landscape:
#  create_ancestor_species <- function(landscape, config) {
#   stop("create the initial species here")
#  }

## ----eval=FALSE---------------------------------------------------------------
#  ### Dispersal
#  # the maximum range to consider when calculating the distances from local distance inputs.
#  max_dispersal <- Inf
#  # returns n dispersal values.
#  get_dispersal_values <- function(n, species, landscape, config) {
#    stop("calculate dispersal values here")
#  }
#  ### Speciation
#  # threshold for genetic distance after which a speciation event takes place.
#  divergence_threshold = NULL
#  # factor by which the divergence is increased between geographically isolated population.
#  # can also be a matrix between the different population clusters.
#  get_divergence_factor <- function(species, cluster_indices, landscape, config) {
#    stop("calculate divergence factor here")
#  }
#  ### Evolution
#  # mutate the traits of a species and return the new traits matrix.
#  apply_evolution <- function(species, cluster_indices, landscape, config) {
#    stop("mutate species traits here")
#  }
#  ### Ecology
#  # called for every cell with all occurring species, this function calculates abundances and/or who survives for each sites
#  # returns a vector of abundances.
#  # set the abundance to 0 for every species supposed to die.
#  apply_ecology <- function(abundance, traits, environment, config) {
#    stop("calculate species abundances and deaths here")
#  }

## ----eval=FALSE---------------------------------------------------------------
#  # Version: 1.0
#  # Author: Oskar Hagen
#  # Date: 26.10.2020
#  # Landscape: SouthAmerica
#  # Publications: R-package gen3sis
#  # Description: Example config used at the introduction vignette and similar to case study global configs in Hagen et al. 2020.
#  # O. Hagen, B. Flück, F. Fopp, J.S. Cabral, F. Hartig, M. Pontarp, T.F. Rangel, L. Pellissier. gen3sis: the general engine for eco-evolutionary simulations on the origins of biodiversity.

## ----eval=FALSE---------------------------------------------------------------
#  random_seed = 6
#  start_time = 40
#  end_time = NA
#  max_number_of_species = 50000
#  max_number_of_coexisting_species = 10000
#  trait_names = c("temp",  "dispersal")
#  environmental_ranges = list("temp" = c(-45, 55), "area"=c(2361.5, 12923.4),
#      "arid"=c(1,0.5))

## ----eval=FALSE---------------------------------------------------------------
#  end_of_timestep_observer = function(data, vars, config){
#    save_species()
#    plot_richness(data$all_species, data$landscape)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  end_of_timestep_observer = function(data, vars, config){
#      par(mfrow=c(2,3))
#      plot_raster_single(data$landscape$environment[,"temp"], data$landscape, "temp", NA)
#      plot_raster_single(data$landscape$environment[,"arid"], data$landscape, "arid", NA)
#      plot_raster_single(data$landscape$environment[,"area"], data$landscape, "area", NA)
#      plot_richness(data$all_species, data$landscape)
#      plot_species_presence(data$all_species[[1]], data$landscape)
#      plot(0,type='n',axes=FALSE,ann=FALSE)
#      mtext("STATUS",1)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  initial_abundance = 1
#  create_ancestor_species <- function(landscape, config) {
#    range <- c(-95, -24, -68, 13)
#    co <- landscape$coordinates
#    selection <- co[, "x"] >= range[1] &
#      co[, "x"] <= range[2] &
#      co[, "y"] >= range[3] &
#      co[, "y"] <= range[4]
#    new_species <- list()
#    for(i in 1:10){
#      initial_cells <- rownames(co)[selection]
#      initial_cells <- sample(initial_cells, 1)
#      new_species[[i]] <- create_species(initial_cells, config)
#      #set local adaptation to max optimal temp equals local temp
#      new_species[[i]]$traits[ , "temp"] <- landscape$environment[initial_cells,"temp"]
#      new_species[[i]]$traits[ , "dispersal"] <- 1
#      plot_species_presence(landscape, species=new_species[[i]])
#    }
#    return(new_species)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  get_dispersal_values <- function(n, species, landscape, config) {
#    values <- rweibull(n, shape = 1.5, scale = 133)
#    return(values)
#  }

## ----eval=T, fig.width=6, fig.height=3.2--------------------------------------
n <- 100
hist(rweibull(n, shape = 1.5, scale = 133), col="black")

## ----eval=FALSE---------------------------------------------------------------
#  divergence_threshold = 2
#  get_divergence_factor <- function(species, cluster_indices, landscape, config) {
#    return(1)
#    }

## ----eval=FALSE---------------------------------------------------------------
#  # mutate the traits of a species and return the new traits matrix
#  apply_evolution <- function(species, cluster_indices, landscape, config) {
#    trait_evolutionary_power <- 0.001
#    traits <- species[["traits"]]
#    cells <- rownames(traits)
#    #homogenize trait based on abundance
#    for(cluster_index in unique(cluster_indices)){
#      cells_cluster <- cells[which(cluster_indices == cluster_index)]
#      mean_abd <- mean(species$abundance[cells_cluster])
#      weight_abd <- species$abundance[cells_cluster]/mean_abd
#      traits[cells_cluster, "temp"] <- mean(traits[cells_cluster, "temp"]*weight_abd)
#    }
#    #mutations
#    mutation_deltas <-rnorm(length(traits[, "temp"]), mean=0, sd=trait_evolutionary_power)
#    traits[, "temp"] <- traits[, "temp"] + mutation_deltas
#    return(traits)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  apply_ecology <- function(abundance, traits, landscape, config) {
#    abundance_scale = 10
#    abundance_threshold = 1
#    #abundance treashold
#    survive <- abundance>=abundance_threshold
#    abundance[!survive] <- 0
#    abundance <- (( 1-abs( traits[, "temp"] - landscape[, "temp"]))*abundance_scale)*as.numeric(survive)
#    #abundance thhreashold
#    abundance[abundance<abundance_threshold] <- 0
#    k <- ((landscape[,"area"]*(landscape[,"arid"]+0.1)*(landscape[,"temp"]+0.1))*abundance_scale^2)
#    total_ab <- sum(abundance)
#    subtract <- total_ab-k
#    if (subtract > 0) {
#      # print(paste("should:", k, "is:", total_ab, "DIFF:", round(subtract,0) ))
#      while (total_ab>k){
#        alive <- abundance>0
#        loose <- sample(1:length(abundance[alive]),1)
#        abundance[alive][loose] <- abundance[alive][loose]-1
#        total_ab <- sum(abundance)
#      }
#      #set negative abundances to zero
#      abundance[!alive] <- 0
#    }
#    return(abundance)
#  }

## ----eval=T-------------------------------------------------------------------
# get data path
datapath <- system.file(file.path("extdata", "SouthAmerica"), package="gen3sis")
# creates config object from config file
config_object <- create_input_config(file.path(datapath, "config/config_southamerica.R"))

# modify the initialization function
config_object$gen3sis$initialization$create_ancestor_species <- function(landscape, config) {
  range <- c(-95, -24, -68, 13)
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] &
    co[, "x"] <= range[2] &
    co[, "y"] >= range[3] &
    co[, "y"] <= range[4]

  initial_cells <- rownames(co)[selection]
  new_species <- create_species(initial_cells, config)
  #set local adaptation to max optimal temp equals local temp
  new_species$traits[ , "temp"] <- landscape$environment[initial_cells,"temp"]
  new_species$traits[ , "dispersal"] <- 1
  new_species$abundance <- new_species$abundance*10
  return(list(new_species))
}

## ----eval=T-------------------------------------------------------------------
# modify the ecology function
config_object$gen3sis$ecology$apply_ecology <- function(abundance, traits, landscape, config, abundance_scale = 10, abundance_threshold = 8) {
  #abundance threshold
  abundance <- as.numeric(!abundance<abundance_threshold)
  abundance <- (( 1-abs( traits[, "temp"] - landscape[, "temp"]))*abundance_scale)*abundance
  #abundance threshold
  abundance[abundance<abundance_threshold] <- 0
  return(abundance)
}

## ----eval=T, fig.width=7, fig.height=3----------------------------------------
config_object_old <- create_input_config(file.path(datapath, "config/config_southamerica.R"))

# define observer
config_object_old$gen3sis$general$end_of_timestep_observer <- function(data, vars, config){
  plot_ranges(data$all_species, data$landscape, disturb=0, max_sps=10)
}

config_object$gen3sis$general$end_of_timestep_observer <- config_object_old$gen3sis$general$end_of_timestep_observer 

# define only 3 time-steps
config_object_old$gen3sis$general$end_time <- 38
config_object$gen3sis$general$end_time <- 38


## -----------------------------------------------------------------------------
verify_config(config_object_old)
verify_config(config_object)

## ----eval=FALSE, fig.width=7, fig.height=3------------------------------------
#  sim_old <- run_simulation(config = config_object_old,
#                 landscape = file.path(datapath, "landscape"),
#                 output_directory=tempdir(),
#                 call_observer = 4,
#                 verbose=0) # no progress printed

## ----eval=FALSE, fig.width=7, fig.height=3------------------------------------
#  sim_new <- run_simulation(config = config_object,
#                            landscape = file.path(datapath, "landscape"),
#                            verbose=0, # no progress printed
#                            output_directory=tempdir())

