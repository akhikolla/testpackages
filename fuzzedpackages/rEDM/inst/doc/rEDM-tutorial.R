## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, tidy = TRUE,
                      fig.pos = 'h', fig.align = 'center',
                      fig.width = 4, fig.height = 3.3)
knitr::opts_knit$set(global.par = TRUE, progress = FALSE)
options(digits = 2)
par.orig = par()
par(mar = c(4, 4, 2, 1), mgp = c(2.5, 1, 0))

## ----CRAN installation instructions, eval = FALSE-----------------------------
#  install.packages("rEDM")

## ----GitHub installation instructions, eval = FALSE---------------------------
#  devtools::install_github("SugiharaLab/rEDM")

## ---- LorenzProjection, echo = FALSE, fig.cap = "Time Series Projection from the Lorenz Attractor.", out.width = 250, fig.align = 'center'----
knitr::include_graphics("Lorenz_Projection.png")

## ----fig_attractor_reconstruction, echo = FALSE, fig.cap = "Attractor Reconstruction from 3 Lagged Coordinates", out.width = 250, fig.align = 'center'----
knitr::include_graphics("Lorenz_Reconstruct.png")

## ----load package-------------------------------------------------------------
library(rEDM)
str(TentMap)

## ----simplex on tentmap-------------------------------------------------------
simplex_out <- Simplex( dataFrame = TentMap, lib = "1 100", pred = "201 500",
columns = 'TentMap', target = 'TentMap', E = 3 )
simplex_out[ c(1:2, 300:301), ]

## ----tentmap simplex stats----------------------------------------------------
ComputeError( simplex_out $ Observations, simplex_out $ Predictions )

## ----Embed_Dim_TentMap, fig.cap="TentMap data prediction skill vs. embedding dimension."----
rho_E <- EmbedDimension( dataFrame = TentMap, lib = "1 100", pred = "201 500",
columns = 'TentMap', target = 'TentMap' )

## ----Prediction_interval_on_TentMap, fig.cap = "Tent map first differences simplex prediction skill as a function of forecast interval."----
rho_Tp <- PredictInterval( dataFrame = TentMap, lib = "1 100",
pred = "201 500", target = 'TentMap', columns = 'TentMap', E = 2 )

## ----Predict_nonlinear_on_TentMap, fig.cap = "Tent map first differences S-map prediction skill as a function of S-map localisation parameter."----
rho_theta <- PredictNonlinear( dataFrame = TentMapNoise, lib = "1 100",
pred = "201 500", target = 'TentMap', columns = 'TentMap', E = 2 )

## ----Simplex_TentMap----------------------------------------------------------
tentMapPredict <- Simplex( dataFrame = TentMap, lib = "1 100",
pred = "201 500", target = 'TentMap', columns = 'TentMap', E = 2 )

ComputeError( tentMapPredict $ Observations, tentMapPredict $ Predictions ) $ rho

## ----Tentmap_Noise_S-map------------------------------------------------------
smap = SMap( dataFrame = TentMapNoise, lib = "1 100", pred = "201 500",
target = 'TentMap', columns = 'TentMap', E = 2, theta = 3 )

## ----S-map_tentMap_results----------------------------------------------------
head( cbind( smap $ predictions, smap $ coefficients ), 2 )
tail( cbind( smap $ predictions, smap $ coefficients ), 2 )

## ----block_3sp_data-----------------------------------------------------------
head( block_3sp, 3 )

## ----3_Species_simplex, fig.cap="Simplex Tp=1 forecast of $x_t$."-------------
smplx_3species = Simplex( dataFrame = block_3sp, lib = "1 100",
pred = "101 190", E = 3, columns = "x_t x_t-1 z_t", target = "x_t",
embedded = TRUE )

## ----plot_3_species, fig.cap="Scatter plot of simplex forecast of $x_t$ vs. observations."----
err = ComputeError( smplx_3species $ Observations,
                    smplx_3species $ Predictions )
plot( smplx_3species $ Observations, smplx_3species $ Predictions,
pch = 19, cex = 0.5, xlab = 'Observations', ylab = 'Predictions',
main = '3 Species x_t' )
abline( a = 0, b = 1, lty = 2, col = 'blue' )
text( -1, 1, paste( capture.output( cbind( err ) ), collapse = '\n' ) )

## ----S-map_Lorenz96-----------------------------------------------------------
smap_Lorenz <- SMap( dataFrame = Lorenz5D, lib = "1 500", pred = "601 900",
E = 4, theta = 3, columns = "V1 V2 V3 V4", target = "V1", embedded = TRUE )

## ----set precision, echo = FALSE----------------------------------------------
options(digits=4)

## ----S-map_Lorenz_results, warning = FALSE------------------------------------
head( cbind( smap_Lorenz $ predictions, smap_Lorenz $ coefficients[,2:6] ), 3 )

## ---- echo = FALSE------------------------------------------------------------
par(mfrow = c(4, 1), mar = c(2, 4, 0.5, 1), oma = c(0, 0, 0, 0),
    mgp = c(2.1, 1, 0) )

## ----smap_Lorenz_plot, out.width = 300, fig.cap = "S-map prediction and coefficients of Lorenz'96 5-D system."----
predictions = smap_Lorenz $ predictions
coefficients = smap_Lorenz $ coefficients
Time = predictions $ Time

plot( Time, predictions $ Observations, type = "l", col = "blue",
ylab = "V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3 )
lines( Time, predictions $ Predictions, lwd = 2, col = 'red' )
legend("topright", legend = c("observed", "predicted"),
fill = c( 'blue', 'red' ), bty = "n", cex = 1.3 )

plot( Time, coefficients[, 6], type = "l", col = "brown",
ylab = "∂V4/∂V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3 )
plot( Time, coefficients[, 5], type = "l", col = "darkgreen",
ylab = "∂V3/∂V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3 )
plot( Time, coefficients[, 4], type = "l", col = "blue",
ylab = "∂V2/∂V1", xlab = "", lwd = 2, cex.lab = 1.3, cex.axis = 1.3 )

## ---- echo = FALSE------------------------------------------------------------
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), mgp = c(2.5, 1, 0))

## ----multiview----------------------------------------------------------------
Mview = Multiview( dataFrame = block_3sp, lib = "1 100", pred = "101 190", E = 3, columns = "x_t y_t z_t", target = "x_t" )

## ----multiview results--------------------------------------------------------
Mview $ View[ which( Mview $ View $ rho > 0.91 ),  ]

## ----fig_cross_map, echo = FALSE, fig.cap = "Cross Mapping Between Reconstructions of the Lorenz Attractor"----
knitr::include_graphics("CrossMap.png")

## ----sardine_anchovy_ccm, fig.cap = "Convergent cross mapping of Newport sea surface temperature with anchovy landings."----
cmap <- CCM( dataFrame = sardine_anchovy_sst, E = 3, Tp = 0, columns = "anchovy", target = "np_sst", libSizes = "10 70 5", sample = 100, showPlot = TRUE )

## ----Thrips data--------------------------------------------------------------
head(Thrips, 2)

## ----thrips plot, echo = FALSE, fig.width = 3, fig.height = 3.5, fig.cap = "Thrips abundance and environmental variables."----
par(mfrow = c(4, 1), mar = c(2, 4, 0, 1), mgp = c(2., 0.5, 0))
time_dec <- Thrips $ Year + (Thrips $ Month)/12
plot(time_dec, Thrips $ Thrips_imaginis, type = "l", lwd = 2, col = "darkgreen", ylab = "Thrips", xlab = "")
plot(time_dec, Thrips $ maxT_degC, type = "l", lwd = 2, col = "red", ylab = "maxT (C)", xlab = "")
plot(time_dec, Thrips $ Rain_mm, type = "l", lwd = 2, col = "blue", ylab = "Rain (mm)", xlab = "")
plot(time_dec, Thrips $ Season, type = "l", lwd = 2, col = "magenta", ylab = "Season", xlab = "" )
mtext("Year", side = 1, outer = TRUE, line = 1)

## ---- echo = FALSE------------------------------------------------------------
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), mgp = c(2.5, 1, 0))

## ----univariate thrips, fig.cap = "Simplex embedding dimension for Thrips abundance.", fig.width = 4, fig.height = 3.----
rho_E <- EmbedDimension( dataFrame = Thrips, columns = "Thrips_imaginis", target = "Thrips_imaginis", lib = "1 72", pred = "1 72", showPlot = TRUE )

## ----smap for thrips, fig.cap = "SMap localisation parameter for Thrips abundance.", fig.width = 4, fig.height = 3.----
E = 8
rho_theta_e3 = PredictNonlinear(dataFrame = Thrips, columns = "Thrips_imaginis", target = "Thrips_imaginis", lib = "1 73", pred = "1 73", E = E )

## ----compute ccm matrix for thrips, results='hold', cache = TRUE--------------
vars       = colnames(Thrips[3:6])
var_pairs  = combn( vars, 2 )  # Combinations of vars, 2 at a time
libSize    = paste( NROW(Thrips) - E, NROW(Thrips) - E, 10, collapse = " " )
ccm_matrix = array(NA, dim = c(length(vars), length(vars)), 
                   dimnames = list(vars, vars))

for( i in 1:ncol( var_pairs ) ) {
   ccm_out = CCM(dataFrame = Thrips,
                 columns = var_pairs[1,i], target = var_pairs[2,i],
                 libSizes = libSize, Tp = 0, E = E, sample = 100)
                 
   outVars = names( ccm_out )

   var_out = unlist( strsplit( outVars[2], ':' ) )
   ccm_matrix[ var_out[2], var_out[1] ] = ccm_out[1,2]
                 
   var_out = unlist( strsplit( outVars[3], ':' ) )
   ccm_matrix[ var_out[2], var_out[1] ] = ccm_out[1,3]
}

## ----compute corr matrix for thrips, cache = TRUE-----------------------------
corr_matrix <- array(NA, dim = c(length(vars), length(vars)), 
                     dimnames = list(vars, vars))

for(ccm_from in vars) {
    for(ccm_to in vars[vars != ccm_from]) {
        ccf_out <- ccf(Thrips[,ccm_from], Thrips[,ccm_to], 
                       type = "correlation", lag.max = 6, plot = FALSE)$acf
        corr_matrix[ccm_from, ccm_to] <- max(abs(ccf_out))
    }
}

## ----xmap vs. corr matrix for thrips------------------------------------------
ccm_matrix
corr_matrix

## ----thrips setup, echo = FALSE-----------------------------------------------
par( mfrow = c( 1, 3 ), mar = c( 3.5, 1.3, 2.3, 0.2 ), mgp = c(1.5, 0.5, 0) )

## ----ccm thrips, fig.cap = "Thrips cross mapped to climatic variables. Vertical axis is cross map prediction skill (rho).", fig.width = 7.5----
thrips_xmap_maxT <- CCM(dataFrame = Thrips, E = E, Tp = 0,
                        columns = "Thrips_imaginis", target = "maxT_degC",
                        libSizes = "13 73 3", sample = 300, showPlot = TRUE)
abline(h = corr_matrix['Thrips_imaginis', 'maxT_degC'], col = "black", lty = 2)

thrips_xmap_maxT <- CCM(dataFrame = Thrips, E = E, Tp = 0,
                        columns = "Thrips_imaginis", target = "Rain_mm",
                        libSizes = "13 73 3", sample = 300, showPlot = TRUE)
abline(h = corr_matrix['Thrips_imaginis', 'Rain_mm'], col = "black", lty = 2)

thrips_xmap_maxT <- CCM(dataFrame = Thrips, E = E, Tp = 0,
                        columns = "Thrips_imaginis", target = "Season",
                        libSizes = "13 73 3", sample = 300, showPlot = TRUE)
abline(h = corr_matrix['Thrips_imaginis', 'Season'], col = "black", lty = 2)

## ----thrips setdown, echo = FALSE---------------------------------------------
par(mar = c(4, 4, 2, 1), mgp = c(2.5, 1, 0))

## ----seasonal surrogates thrips, cache = TRUE---------------------------------
# Create matrix with temperature and rain surrogates (1000 time series vectors)
surr_maxT = SurrogateData( Thrips $ maxT_degC, method = "seasonal", T_period = 12, num_surr = 1000, alpha = 3)
surr_rain = SurrogateData( Thrips $ Rain_mm, method = "seasonal", T_period = 12, num_surr = 1000, alpha = 3)

# Rain cannot be negative
surr_rain = apply( surr_rain, 2, function(x){ i = which(x<0); x[i] = 0; x } )

# data.frame to hold CCM rho values between Thrips abundance and variable
rho_surr <- data.frame( maxT = numeric(1000), Rain = numeric(1000) )

# data.frames with time, Thrips, and 1000 surrogate climate variables for CCM()
maxT_data = as.data.frame( cbind( seq(1:nrow(Thrips)), Thrips $ Thrips_imaginis, surr_maxT ) )
names( maxT_data ) = c( 'time', 'Thrips_imaginis', paste( 'T', as.character( seq( 1,1000 ) ), sep = '' ) )

rain_data = as.data.frame( cbind( seq(1:nrow(Thrips)), Thrips $ Thrips_imaginis, surr_rain ) )
names( rain_data ) = c( 'time', 'Thrips_imaginis', paste( 'R', as.character( seq( 1,1000 ) ), sep = '' ) )

# Cross mapping
for (i in 1:1000) {
  targetCol = paste( 'T', i, sep = '' ) # as in maxT_data
  
  ccm_out = CCM( dataFrame = maxT_data, E = E, Tp = 0,
                 columns = 'Thrips_imaginis', target = targetCol,
                 libSizes = "73 73 5", sample = 1 )
                 
  col = paste( 'Thrips_imaginis', ':', targetCol, sep = '' )
   
  rho_surr $ maxT[i] = ccm_out[1,col]
}

for (i in 1:1000) {
  targetCol = paste( 'R', i, sep = '' ) # as in rain_data
  
  ccm_out = CCM( dataFrame = rain_data, E = E, Tp = 0,
                 columns = 'Thrips_imaginis', target = targetCol,
                 libSizes = "73 73 5", sample = 1 )
                 
  col = paste( 'Thrips_imaginis', ':', targetCol, sep = '' )
   
  rho_surr $ Rain[i] = ccm_out[1,col]
}

## ----significance of randomization test---------------------------------------
1 - ecdf( rho_surr $ maxT )( ccm_matrix["maxT_degC", "Thrips_imaginis" ] )
1 - ecdf( rho_surr $ Rain )( ccm_matrix["Rain_mm",   "Thrips_imaginis" ] )

## ----reset par, include = FALSE-----------------------------------------------
par( par.orig )

## ----parameter-table, echo = FALSE--------------------------------------------
paramTable = read.csv( "ParameterTable.csv", as.is = TRUE )
knitr::kable( paramTable, caption = '' )

