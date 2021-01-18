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

## ----tidy=T, out.width='70%', echo=F, fig.cap='*gen3sis*: **gen**eral **e**ngine for **e**co-**e**volutionary **si**mulation**s**', fig.margin=F----
knitr::include_graphics("../inst/logo/gen3sis_logo.png")

## ----setup, message=F---------------------------------------------------------
library(gen3sis)

## ----eval=TRUE----------------------------------------------------------------
print(paste("gen3sis version:", packageVersion("gen3sis")))

## ----eval=TRUE, echo=FALSE----------------------------------------------------
datapath <- system.file(file.path("extdata", "SouthAmerica"), package="gen3sis")

## ----eval=T, fig.width=7, fig.height=4, message=F, echo=F, warning=T, results='hide', fig.cap='This figure shows the temperature and aridity of the lanscape used in this vignette at 65, 30 and 0Ma.'----
library(raster)
# read landscapes.rds
landscapes <- readRDS(file.path(datapath, "landscape/landscapes.rds"))
# create rasters for each timestep
temp_65 <- rasterFromXYZ(landscapes$temp[, c(1,2,68)])
temp_30 <- rasterFromXYZ(landscapes$temp[, c(1,2,33)])
temp_0 <- rasterFromXYZ(landscapes$temp[, c(1,2,3)])   

# TEMPERATURE
par(mar=c(1,2,3,1), oma=c(0,0,3,0))
layout(matrix(c(1,1,1,4,
                2,2,2,4,
                3,3,3,4), ncol=3))
maxtemp <- 33
mintemp <- -27
temp_breaks <- seq(mintemp, maxtemp, by=1)
temp_colors <- rev(heat.colors(length(temp_breaks)-1))
image(temp_65, col=temp_colors, breaks=temp_breaks, main='65 Ma')
image(temp_30, col=temp_colors, breaks=temp_breaks, main='30 Ma')
title(main=expression('Temperature [\u00B0C]'), line=1, outer=T, cex.main=2)
image(temp_0, col=temp_colors, breaks=temp_breaks, main='0 Ma')

plot.new()
legend_df <- as.data.frame(cbind(seq(0, length(temp_breaks)-1, length.out=(length(temp_breaks))), rep(0.25, (length(temp_breaks))), temp_breaks))
legend_image <- rasterFromXYZ(legend_df, res=0.01)
plot(legend_image, legend.only=T, col=temp_colors, horizontal=T, smallplot=c(0.2, 0.8, 0.55, 0.7), 
   axis.args=list(at=seq(mintemp, maxtemp, 5),labels=seq(mintemp, maxtemp, 5)))

# ARIDITY
arid_300 <- rasterFromXYZ(landscapes$arid[, c(1,2,68)])
arid_65 <- rasterFromXYZ(landscapes$arid[, c(1,2,33)])
arid_0 <- rasterFromXYZ(landscapes$arid[, c(1,2,3)])   

par(mar=c(1,2,3,1), oma=c(0,0,3,0))
layout(matrix(c(1,1,1,4,
                2,2,2,4,
                3,3,3,4), ncol=3))
maxarid <- 1
minarid <- 0
arid_breaks <- seq(minarid, maxarid, by=0.01)
arid_colors <- colorRampPalette(c('grey95', 'peru'))(length(arid_breaks)-1)
image(arid_300, col=arid_colors, breaks=arid_breaks, main='65 Ma')
image(arid_65, col=arid_colors, breaks=arid_breaks, main='30 Ma')
title(main=expression('Aridity [index]'), line=1, outer=T, cex.main=2)
image(arid_0, col=arid_colors, breaks=arid_breaks, main='0 Ma')

plot.new()
legend_df <- as.data.frame(cbind(seq(0, length(arid_breaks)-1, length.out=(length(arid_breaks))), rep(0.25, (length(arid_breaks))), arid_breaks))
legend_image <- rasterFromXYZ(legend_df, res=0.01)
plot(legend_image, legend.only=T, col=arid_colors, horizontal=T, smallplot=c(0.2, 0.8, 0.55, 0.7), 
axis.args=list(at=seq(minarid, maxarid, 0.2),labels=seq(minarid, maxarid, 0.2)))

# # AREA
# area_300 <- rasterFromXYZ(landscapes$area[, c(1,2,68)])
# area_65 <- rasterFromXYZ(landscapes$area[, c(1,2,38)])
# area_0 <- rasterFromXYZ(landscapes$area[, c(1,2,3)])   
# 
# par(mar=c(1,2,3,1), oma=c(0,0,3,0))
# layout(matrix(c(1,1,1,4,
#                 2,2,2,4,
#                 3,3,3,4), ncol=3))
# maxarea <- 197000
# minarea <- 151000
# area_breaks <- seq(minarea, maxarea, by=1000)
# area_colors <- rev(gray(seq(0.05, 0.9, length.out=length(area_breaks)-1)))
# image(area_300, col=area_colors, breaks=area_breaks, main='65 Ma')
# image(area_65, col=area_colors, breaks=area_breaks, main='30 Ma')
# title(main=expression('Area [' ~km^2~ ']'), line=1, outer=T, cex.main=2)
# image(area_0, col=area_colors, breaks=area_breaks, main='0 Ma')
#  
# plot.new()
# legend_df <- as.data.frame(cbind(seq(0, length(area_breaks)-1, length.out=(length(area_breaks))), rep(0.25, (length(area_breaks))), area_breaks))
# legend_image <- rasterFromXYZ(legend_df, res=0.01)
# plot(legend_image, legend.only=T, col=area_colors, horizontal=T, smallplot=c(0.2, 0.8, 0.55, 0.7), 
# axis.args=list(at=seq(minarea, maxarea, 7500),labels=seq(minarea, maxarea, 7500)))

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  #we set verbose to 0 to avoid a large console outputs of how the simulation is developing
#  sim <- run_simulation(config = file.path(datapath, "config/config_southamerica.R"),
#                 landscape = file.path(datapath, "landscape"),
#                 output_directory = tempdir(),
#                 call_observer = 1,
#                 verbose=0)

## ----eval=T, echo=F-----------------------------------------------------------
sim <- readRDS(file.path(datapath, "output/config_southamerica/sgen3sis.rds"))

## ----eval=F, echo=T, fig.width=7, fig.height=2--------------------------------
#  timesteps <- c(40, 20, 0)
#  par(mfrow=c(1,3))
#  for(i in timesteps){
#    landscape_i <- readRDS(file.path(datapath, paste0('output/config_southamerica/landscapes/landscape_t_', i ,'.rds')))
#    species_i <- readRDS(file.path(datapath, paste0('output/config_southamerica/species/species_t_', i ,'.rds')))
#    plot_richness(species_i, landscape_i)
#  }

## ----eval=T, echo=F, fig.width=7, fig.height=2--------------------------------
knitr::include_graphics("../inst/extdata/SouthAmerica/images/introduction_richness.png")

## ----eval=TRUE, fig.width=7, fig.height=7-------------------------------------
plot_summary(sim)

## ----eval=TRUE, fig.width=8, fig.height=5, fig.retina=2-----------------------
landscapes <- readRDS(file.path(datapath, "landscape", "landscapes.rds")) #get the input landscape
landscape_t0 <- as.data.frame(cbind(landscapes$temp[, 1:2], temp=landscapes$temp[, 3], arid=landscapes$arid[,3], area=landscapes$area[,3])) #get landscape at last time-step
landscape_t0 <- cbind(landscape_t0, rich=sim$summary$`richness-final`[,3]) #add richness to the dataframe
landscape_t0 <- na.omit(landscape_t0)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
glm.uni <- glm(rich ~ poly(temp, 2), data=landscape_t0, family=poisson)
cor(landscape_t0$temp, landscape_t0$rich)

## ----eval=TRUE, fig.width=7, fig.height=6, fig.retina=2, message=F, results='hide'----
# prepare data with temperature and predicted richness from our model
data_plot <- data.frame(cbind(landscape_t0$temp, predict(glm.uni,type = "response")))
# sort data for plotting and ommit NA's
data_plot <- na.omit(data_plot[order(data_plot[,1], decreasing = FALSE),])
# get the number of observations
n <- paste0('observations (n = ', length(landscape_t0$rich), ')')
# plot model curve
plot(data_plot[,1],data_plot[,2], xlab="Temperature [\u00B0C]", ylab=expression(paste(alpha," richness")), frame.plot=F, type="l", col='red', lwd=2, xlim=c(min(landscape_t0$temp), max(landscape_t0$temp)), ylim=c(min(landscape_t0$rich), max(landscape_t0$rich)))
# add observed points
points(landscape_t0$temp, landscape_t0$rich, col=rgb(0.5,0.5,0.5, alpha=0.4), pch=16) 
# add legend
legend(-20,30, col=c(rgb(0.5,0.5,0.5,0.4), 'red'), legend=c(n, 'model fit'), pch=c(16, NA), lty=c(NA, 1), lwd=c(NA, 2), bty='n')

## ----eval=TRUE, fig.width=7, fig.height=6, fig.retina=2, results='hide', echo=F----
glm.uni <- glm(rich ~ poly(arid, 2), data=landscape_t0, family=poisson)
cor(landscape_t0$arid, landscape_t0$rich)

#Plot the response curve
data_plot <- data.frame(cbind(landscape_t0$arid, predict(glm.uni,type = "response")))
sorted <- na.omit(data_plot[order(data_plot[,1], decreasing = FALSE),])
data_plot <- data.frame(cbind(sorted[,1],sorted[,2]))
plot(data_plot[,1],data_plot[,2], xlab="Aridity", ylab=expression(paste(alpha," richness")), frame.plot=F, type="l", col='blue', 
     lwd=2, xlim=c(min(landscape_t0$arid), max(landscape_t0$arid)), ylim=c(min(landscape_t0$rich), max(landscape_t0$rich)))
points(landscape_t0$arid, landscape_t0$rich, col=rgb(0.5,0.5,0.5, alpha=0.4), pch=16) 
legend(0.4,30, col=c(rgb(0.5,0.5,0.5,0.4), 'blue'), legend=c(n, 'model fit'), pch=c(16, NA), lty=c(NA, 1), lwd=c(NA, 2), bty='n')

## ----eval=TRUE----------------------------------------------------------------
richness_model <- glm(rich ~ poly(temp, 2) + poly(arid, 2), family='quasipoisson', data=landscape_t0)
summary(richness_model)$coefficients

