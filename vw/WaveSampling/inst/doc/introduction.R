## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(WaveSampling)


## -----------------------------------------------------------------------------
# install.packages(sp)
library(sp)
library(sf)
library(sampling)
data("meuse")
data("meuse.riv")
meuse.riv <- meuse.riv[which(meuse.riv[,2] < 334200 & meuse.riv[,2] > 329400),]
meuse_sf <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")

X <- scale(as.matrix(meuse[,1:2]))
pik <- inclusionprobabilities(meuse$copper,30)

## -----------------------------------------------------------------------------
s <- wave(X,pik)
sum(s)

## ----fig.width=5,fig.height=5-------------------------------------------------
library(ggplot2)
p <- ggplot()+
  geom_sf(data = meuse_sf,aes(size=copper),show.legend = 'point',shape = 1,stroke = 0.3)+
  geom_polygon(data = data.frame(x = meuse.riv[,1],y = meuse.riv[,2]),
               aes(x = x,y = y),
               fill = "lightskyblue2",
               colour= "grey50")+
  geom_point(data = meuse,
             aes(x = x,y = y,size = copper),
             shape = 1,
             stroke = 0.3)+
  geom_point(data = meuse[which(s == 1),],
             aes(x = x,y = y,size = copper),
             shape = 16)+
  labs(x = "Longitude",
       y = "Latitude",
       title = NULL,
       size = "Copper",
       caption = NULL)+
  scale_size(range = c(0.5, 3.5))+
  theme_minimal()
p
  

## ----fig.width=7,fig.height=5-------------------------------------------------

library(sp)
library(sampling)
library(ggvoronoi)
data("meuse")
data("meuse.area")

v <- sb_vk(pik,as.matrix(meuse[,1:2]),s)
meuse$v <- v

p <- p + geom_voronoi(data = meuse[which(s == 1),],
               aes(x = x,y = y,fill = v),
               outline =as.data.frame(meuse.area),
               size = 0.1,
               colour = "black")+
  geom_point(data = meuse,
             aes(x = x,y = y,size = copper),
             shape = 1,
             stroke = 0.3)+
  geom_point(data = meuse[which(s == 1),],
             aes(x = x,y = y,size = copper),
             shape = 16)+
  scale_fill_gradient2(midpoint = 1)
p

BalancedSampling::sb(pik,as.matrix(meuse[,1:2]),which(s == 1))



## -----------------------------------------------------------------------------

W <- wpik(X,pik)
W <- W - diag(diag(W))
IB(W,s)


W1 <- wpikInv(X,pik)
IB(W1,s)


