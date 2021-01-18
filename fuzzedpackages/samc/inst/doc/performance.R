## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(samc)

## ----fig1, fig.width = 6.5, fig.asp = 1, fig.align = "center", echo = FALSE----
xdat <- seq(0, 1000000, 10000)
ydat <- (xdat * xdat * 8) / (1024^3)

plot(xdat/1000, ydat,
     main = "Transition Matrix Memory Requirements",
     xlab = "Number of Landscape Cells (Thousands)",
     ylab = "RAM (GB)")


## ----fig2, fig.width = 6.5, fig.asp = 1, fig.align = "center", echo = FALSE----
xdat <- seq(0, 1000000, 10000)
ydat <- (xdat * 17 * 8) / (1024^3)

plot(xdat/1000, ydat,
     main = "Transition Matrix Memory Requirements (Optimized)",
     xlab = "Number of Landscape Cells (Thousands)",
     ylab = "RAM (GB)")


## ----fig3, fig.width = 6.5, fig.asp = 1, fig.align = "center", echo = FALSE----
xdat <- seq(0, 50000, 1000)
ydat <- (xdat * xdat * 8) / (1024^3)

plot(xdat/1000, ydat,
     main = "Dense Matrix Memory Requirements",
     xlab = "Number of Landscape Cells (Thousands)",
     ylab = "RAM (GB)")


## ----fig4, out.width = '100%', fig.align = "center", echo = FALSE-------------
knitr::include_graphics("img/bt.png")

## ----fig5, out.width = '100%', fig.align = "center", echo = FALSE-------------
knitr::include_graphics("img/br.png")

## ----fig6, out.width = '100%', fig.align = "center", echo = FALSE-------------
knitr::include_graphics("img/btt.png")

## ----fig7, out.width = '100%', fig.align = "center", echo = FALSE-------------
knitr::include_graphics("img/btr.png")

