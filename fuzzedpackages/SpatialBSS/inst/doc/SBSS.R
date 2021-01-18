## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----sim_coords---------------------------------------------------------------
coords <- runif(1000 * 2) * 10
dim(coords) <- c(1000, 2)

## ----sim_field----------------------------------------------------------------
if (requireNamespace("RandomFields", quietly = TRUE)) {
  mix_mat <- matrix(rnorm(9), 3, 3)
  
  RandomFields::RFoptions(spConform = FALSE)
  field_1 <- RandomFields::RFsimulate(model = RandomFields::RMexp(scale = 0.5), 
                                      x = coords)
  field_2 <- RandomFields::RFsimulate(model = RandomFields::RMexp(scale = 1), 
                                      x = coords)
  field_3 <- RandomFields::RFsimulate(model = RandomFields::RMspheric(scale = 2), 
                                      x = coords)
  field <- tcrossprod(cbind(field_1, field_2, field_3), mix_mat)
} else {
  stop('The package RandomFields is needed to run this example.')
}

## -----------------------------------------------------------------------------
kernel_type <- 'ring'
kernel_parameters <- c(0, 1.5, 1.5, 3, 3, 4.5, 4.5, 6)

## ----sbss_func----------------------------------------------------------------
library('SpatialBSS')
sbss_res <- sbss(x = field, coords = coords, 
                 kernel_type = kernel_type, 
                 kernel_parameters = kernel_parameters)

## ----sbss_plot----------------------------------------------------------------
plot(sbss_res, colorkey = TRUE, as.table = TRUE, cex = 1)

## ----sbss_predict-------------------------------------------------------------
predict(sbss_res, p = 2, n_grid = 50, colorkey = TRUE, as.table = TRUE, cex = 1)

## ----sbss_sp------------------------------------------------------------------
field_sp <- sp::SpatialPointsDataFrame(coords = coords, data = data.frame(field))
res_sbss_sp <- sbss(x = field_sp, kernel_type = kernel_type, 
                    kernel_parameters = kernel_parameters)

## ----sbss_sf------------------------------------------------------------------
if (requireNamespace('sf', quietly = TRUE)) {
  field_sf <- sf::st_as_sf(data.frame(coords = coords, field), 
                           coords = c(1,2))
  res_sbss_sf <- sbss(x = field_sf, kernel_type = kernel_type,
                      kernel_parameters = kernel_parameters)
} else {
  stop('Please install the package sf to run the example code.')
}


## ----sbss_func_ldiff----------------------------------------------------------
sbss_res_lcov <- sbss(x = field, coords = coords, 
                 kernel_type = kernel_type, lcov = 'ldiff',
                 kernel_parameters = kernel_parameters)

## ----k_mat--------------------------------------------------------------------
ring_kernel_matrices <- spatial_kernel_matrix(coords, kernel_type, kernel_parameters)

## ----sbss_k_list--------------------------------------------------------------
sbss_k <- sbss(x = field, kernel_list = ring_kernel_matrices)

## ----lcov_mat-----------------------------------------------------------------
local_cov <- local_covariance_matrix(field, kernel_list =  ring_kernel_matrices, 
                                     whitening = TRUE)

## ----ldiff_mat----------------------------------------------------------------
local_diff <- local_covariance_matrix(field, kernel_list =  ring_kernel_matrices, 
                                      lcov = 'ldiff', whitening = TRUE)

