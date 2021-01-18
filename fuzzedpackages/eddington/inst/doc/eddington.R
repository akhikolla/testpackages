## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  fig.align = "center"
)

# Putting this here too so that we don't have pkg startup messages.
library(dplyr)

## ----setup--------------------------------------------------------------------
library(eddington)
head(rides)

## ----xform--------------------------------------------------------------------
library(dplyr)

days <- rides %>%
  group_by(ride_date) %>%
  summarize(n = n(), total = sum(ride_length))

head(days)

## ----summary------------------------------------------------------------------
summary(days)

## ---- echo=FALSE--------------------------------------------------------------
hist(
  as.integer(days$total), 
  breaks = 30, 
  main = "Histogram of Daily Mileages", 
  xlab = "Miles"
)

abline(v = E_num(days$total), col = "darkred")

legend(
  "topright",
  legend = "Eddington Number",
  col = "darkred",
  bty = "n",
  lty = 1L
)

## ----enum---------------------------------------------------------------------
E_num(days$total)

## ----ecum---------------------------------------------------------------------
days$E <- E_cum(days$total)

head(days)

## ----needle, echo=FALSE-------------------------------------------------------
E <- E_num(days$total)

E_contribs <- days[days$total >= E,]


plot(
  y = days$total,
  x = days$ride_date,
  type = "h",
  main = "Ride Mileages in 2009",
  xlab = "Ride Day",
  ylab = "Total Miles",
  bty = "n",
  ylim = c(0, 90)
)


lines(
  y = c(0, days$E),
  x = c(as.Date("2009-01-01"), days$ride_date),
  type = "s",
  col = "darkred"
)

abline(h = E, lty = 2L, col = "darkred")

text(
  E_contribs[,c("ride_date","total")],
  labels = as.integer(E_contribs[["total"]]),
  pos = 3,
  cex = 0.7
  )

legend(
  "topleft",
  title = "Eddington Number",
  legend = c("Cumulative", "Summary"),
  col = "darkred",
  bty = "n",
  lty = c(1L, 2L)
)

## ----enext--------------------------------------------------------------------
E_next(days$total)

## ----ereq---------------------------------------------------------------------
E_req(days$total, 50)

## ----esat---------------------------------------------------------------------
E_sat(days$total, 30)

