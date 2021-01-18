## ----include = FALSE----------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  tidy=T,
  fig.align='center',
  tidy.opts = list(width.cutoff=80),
  results='hold'
)

## ----setup, message=F, eval=FALSE, results='hide'-----------------------------
#  #first we load the required packages
#  library(gen3sis)
#  library(raster)
#  
#  #next we set the working directory to the downloaded data
#  datapath <- "data_from_project_gen3sis_simulations"

## ----echo=FALSE, message=F, results='hide'------------------------------------
#knitr::opts_knit$set(root.dir = '../inst/extdata/')
# getwd()
library(gen3sis)
library(raster)

## ----eval=FALSE, message=F----------------------------------------------------
#  temperature_brick <- brick('InputRasters/SouthAmerica/temperature_rasters.grd')
#  aridity_brick <- brick('InputRasters/SouthAmerica/aridity_rasters.grd')
#  area_brick <- brick('InputRasters/SouthAmerica/area_rasters.grd')

## ----eval=FALSE---------------------------------------------------------------
#  landscapes_list <- list(temp=NULL, arid=NULL, area=NULL)
#  for(i in 1:nlayers(temperature_brick)){
#    landscapes_list$temp <- c(landscapes_list$temp, temperature_brick[[i]])
#    landscapes_list$arid <- c(landscapes_list$arid, aridity_brick[[i]])
#    landscapes_list$area <- c(landscapes_list$area, area_brick[[i]])
#  }

## ----eval=FALSE---------------------------------------------------------------
#  cost_function_null <- function(source, habitable_src, dest, habitable_dest) {
#      return(1/1000)
#  }

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  landscapes_list_t0 <- list(temp=NULL)
#  landscapes_list_t0$temp <- c(landscapes_list_t0$temp, temperature_brick[[1]])
#  
#  create_input_landscape(landscapes = landscapes_list_t0,
#                                 cost_function = cost_function_null,
#                                 output_directory =  tempdir(), # a directory name to save the files in
#                                 directions = 8, # all surrounding sites from a focus site
#                                 calculate_full_distance_matrices = TRUE,  # full distance matrix
#                                 overwrite_output = T,
#                                 crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#  )
#  

## ---- eval=T, echo=FALSE, fig.width=4, fig.height=7, fig.cap='This figure shows the connection costs from one site in the middle of South America to all other sites. To travel to the site of Antarctica that is indicated with the arrow, the travelling cost is 3790. The distance matrix was computed using the very simple cost function that has been introduced before and is not adding any penalty.', fig.align='center'----
knitr::include_graphics("../inst/extdata/SouthAmerica/images/const_cost.png")

#dist_matrix_null_t0 <- readRDS(file.path(datapath, 'CostFunctionExamples/cost_function_null/distances_full/distances_full_0.rds'))
#landscapes_null <- readRDS(file.path(datapath, 'CostFunctionExamples/cost_function_null/landscapes.rds'))

#landscapes_null <- readRDS('inst/extdata/SouthAmerica/landscape/landscapes.rds')
#landscapes_null_t0 <- na.omit(landscapes_null$temp[, c('x', 'y', '0')])

#dist_null_t0_mat <- cbind(landscapes_null_t0[,c('x', 'y')], cost=as.numeric(dist_matrix_null_t0[,1300]))
#dist_null_t0 <- rasterFromXYZ(dist_null_t0_mat)

#maxcost <- 6200
#mincost <- 0
#cost_breaks <- seq(mincost, maxcost, by=20)
#cost_colors <- rev(gray(seq(0.03, 0.9, length.out=length(cost_breaks)-1)))

#par(mar=c(1,2,1,2))
#layout(matrix(c(1,1,1,1,2), ncol=1))

#par(mar=c(1,1,1,2))
#image(dist_null_t0, col=cost_colors, breaks=cost_breaks)
#plot(dist_null_t0, col=cost_colors, breaks=cost_breaks, axes=F, box=F, legend.args = list(text = 'connection cost', side = 2, 
#         font = 2, line = 0.5, cex = 0.8))
#arrows(-61.5, -30.5, -61.5, -64.5, lwd=3, col='red')
#text(-55, -40, labels=paste('connection cost =', round(dist_null_t0_mat$cost[dist_null_t0_mat$x==(-61.5) & dist_null_t0_mat$y==(-64.5)], 0)), cex=1.5, font=2)

#plot.new()
#legend_df <- as.data.frame(cbind(seq(0, length(cost_breaks)-1, length.out=(length(cost_breaks))), rep(0.25, (length(cost_breaks))), cost_breaks))
#legend_image <- rasterFromXYZ(legend_df, res=0.01)
#plot(legend_image, legend.only=T, col=cost_colors, horizontal=T, smallplot=c(0.2, 0.8, 0.45, 0.6), 
#   axis.args=list(at=seq(mincost, maxcost, 500),labels=seq(mincost, maxcost, 500)), legend.args=list(text='connection cost'))

## ----eval=FALSE---------------------------------------------------------------
#  cost_function_water <- function(source, habitable_src, dest, habitable_dest) {
#    if(!all(habitable_src, habitable_dest)) {
#      return(2/1000)
#    } else {
#      return(1/1000)
#    }
#  }

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  create_input_landscape(landscapes = landscapes_list_t0,
#                                 cost_function = cost_function_water,
#                                 output_directory = tempdir(),# a directory name to save the files in
#                                 directions = 8, # all surrounding sites from a focus site
#                                 calculate_full_distance_matrices = TRUE,  # full distance matrix
#                                 overwrite_output = T,
#                                 crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#  )
#  

## ---- eval=TRUE, echo=FALSE, fig.width=4, fig.height=7, fig.align='center', fig.cap='By using a cost function which penalises the crossing of water, the cost of travelling from our cell in the middle of South America to Antarctica increases to 5204.'----
knitr::include_graphics("../inst/extdata/SouthAmerica/images/var_cost.png")
# dist_matrix_water_t0 <- readRDS('inst/extdata/SouthAmerica/landscape/distances_full/cost_function_water/distances_full/distances_full_0.rds')
# # C:\VITAL LOCAL\Meus Documentos\ETH PhD\Code\R\package\Gen3sis\inst\extdata\SouthAmerica\landscape\distances_full\cost_function_null\distances_full
# landscapes_water <- readRDS('inst/extdata/SouthAmerica/landscape/andscapes.rds')
# landscapes_water_t0 <- na.omit(landscapes_water$temp[, c('x', 'y', '0')])
# dist_water_t0_mat <- cbind(landscapes_water_t0[,c('x', 'y')], cost=as.numeric(dist_matrix_water_t0[,1300]))
# dist_water_t0 <- rasterFromXYZ(dist_water_t0_mat)
# 
# maxcost <- 6200
# mincost <- 0
# cost_breaks <- seq(mincost, maxcost, by=20)
# cost_colors <- rev(gray(seq(0.03, 0.9, length.out=length(cost_breaks)-1)))
# 
# par(mar=c(1,2,1,2))
# layout(matrix(c(1,1,1,1,2), ncol=1))
# 
# #par(mar=c(1,1,1,2))
# image(dist_water_t0, col=cost_colors, breaks=cost_breaks)
# #plot(dist_water_t0, col=cost_colors, breaks=cost_breaks, axes=F, box=F, legend.args = list(text = 'connection cost', side = 2, 
# #         font = 2, line = 0.5, cex = 0.8))
# arrows(-61.5, -30.5, -61.5, -64.5, lwd=3, col='red')
# text(-55, -40, labels=paste('connection cost =', round(dist_water_t0_mat$cost[dist_water_t0_mat$x==(-61.5) & dist_water_t0_mat$y==(-64.5)], 0)), cex=1.5, font=2)
# 
# plot.new()
# legend_df <- as.data.frame(cbind(seq(0, length(cost_breaks)-1, length.out=(length(cost_breaks))), rep(0.25, (length(cost_breaks))), cost_breaks))
# legend_image <- rasterFromXYZ(legend_df, res=0.01)
# plot(legend_image, legend.only=T, col=cost_colors, horizontal=T, smallplot=c(0.2, 0.8, 0.45, 0.6), 
#    axis.args=list(at=seq(mincost, maxcost, 500),labels=seq(mincost, maxcost, 500)), legend.args=list(text='connection cost'))

## ----eval=FALSE---------------------------------------------------------------
#  create_input_landscape(landscapes = landscapes_list,
#               cost_function = cost_function_water,
#               directions=8,
#               output_directory = file.path(tempdir(), "SouthAmerica"),
#               timesteps = paste0(seq(65, 0, by=-1), "Ma"),
#               crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
#               calculate_full_distance_matrices = F)

