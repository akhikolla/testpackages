## ---- echo=FALSE--------------------------------------------------------------
library(chillR)

## -----------------------------------------------------------------------------
library(chillR)
library(ggplot2)
data(KA_weather)
data(KA_bloom)
hourtemps <- stack_hourly_temps(KA_weather, latitude=50.4)

## -----------------------------------------------------------------------------
yc <- 40
zc <- 190
iSeason <- genSeason(hourtemps, years=c(2009))
res <- PhenoFlex(temp=hourtemps$hourtemps$Temp[iSeason[[1]]],
                 times=c(1: length(hourtemps$hourtemps$Temp[iSeason[[1]]])),
                 zc=zc, stopatzc=TRUE, yc=yc, basic_output=FALSE)

## ---- fig.width = 6, fig.height=4, fig.cap="Chill accumulation over time. The dashed line respresents $y_c$, the critical amount of chill units for ecodormancy to be broken."----
DBreakDay <- res$bloomindex
seasontemps<-hourtemps$hourtemps[iSeason[[1]],]
seasontemps[,"x"]<-res$x
seasontemps[,"y"]<-res$y
seasontemps[,"z"]<-res$z
seasontemps<-add_date(seasontemps)

ggplot(data=seasontemps[1:DBreakDay,],aes(x=Date,y=y)) +
  geom_line(col="blue",lwd=1.5) +
  theme_bw(base_size=20) +
  geom_hline(yintercept=yc,lty=2) +
  labs(title="Chill (y) accumulation")

## ---- fig.width = 6, fig.height=4, fig.cap="Heat accumulation over time. The dashed line respresents the $z_c$, the critical amount of heat units for ecodormancy to be broken."----
ggplot(data=seasontemps[1:DBreakDay,],aes(x=Date,y=z)) +
  geom_line(col="red",lwd=1.5) +
  theme_bw(base_size=20) +
  geom_hline(yintercept=zc,lty=2) +
  labs(title="Heat (z) accumulation")


## -----------------------------------------------------------------------------
SeasonList <- genSeasonList(hourtemps$hourtemps, years=c(2003:2008))

## -----------------------------------------------------------------------------
par <-   c(40, 190, 0.5, 25, 3372.8,  9900.3, 6319.5, 5.939917e13,  4, 36,  4,  1.60)
upper <- c(41, 200, 1.0, 30, 4000.0, 10000.0, 7000.0,       6.e13, 10, 40, 10, 50.00)
lower <- c(38, 180, 0.1, 0 , 3000.0,  9000.0, 6000.0,       5.e13,  0,  0,  0,  0.05)

## -----------------------------------------------------------------------------
Fit_res <- phenologyFitter(par.guess=par, 
                           modelfn = PhenoFlex_GDHwrapper,
                           bloomJDays=KA_bloom$pheno[which(KA_bloom$Year %in% c(2003:2008))],
                           SeasonList=SeasonList, lower=lower, upper=upper,
                           control=list(smooth=FALSE, verbose=FALSE, maxit=10,
                                        nb.stop.improvement=5))

## ---- eval=FALSE--------------------------------------------------------------
#  control=list(smooth=FALSE, verbose=TRUE, maxit=1000,
#               nb.stop.improvement=250)

## ---- eval=FALSE--------------------------------------------------------------
#  modelfn=StepChill_Wrapper

## ---- comment=""--------------------------------------------------------------
summary(Fit_res)

## ---- fig.width = 6, fig.height=4---------------------------------------------
plot(Fit_res)

## -----------------------------------------------------------------------------
Fit_res.boot <- bootstrap.phenologyFit(Fit_res, boot.R=10,
                                       control=list(smooth=FALSE, verbose=TRUE, maxit=10,
                                                    nb.stop.improvement=5),
                                       lower=lower, upper=upper, seed=1726354)

## ---- comment=""--------------------------------------------------------------
summary(Fit_res.boot)

## ---- fig.width = 6, fig.height=4---------------------------------------------
plot(Fit_res.boot)

