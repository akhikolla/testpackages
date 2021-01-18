## ----pkg-load, echo = FALSE, message = FALSE----------------------------------
library (dodgr)

## ----isodists-----------------------------------------------------------------
graph <- weight_streetnet (hampi)
from <- sample (graph$from_id, size = 100)
dlim <- c (1, 2, 5, 10, 20) * 100
d <- dodgr_isodists (graph, from = from, dlim)
dim (d)
knitr::kable (head (d))

## ----isodist-table------------------------------------------------------------
table (d$dlim)

## ----isodist-plot-fakey, eval = FALSE-----------------------------------------
#  from <- sample (graph$from_id, size = 1)
#  dlim <- c (1, 2, 5, 10, 20) * 100
#  d <- dodgr_isodists (graph, from = from, dlim)
#  cols <- terrain.colors (length (dlim))
#  index <- which (d$dlim == max (d$dlim)) # plot max contour first
#  plot (d$x [index], d$y [index], "l", col = cols [1],
#        xlab = "longitude", ylab = "latitude")
#  for (i in seq (dlim) [-1]) {
#      index <- which (d$dlim == rev (dlim) [i])
#      lines (d$x [index], d$y [index], col = cols [i], lwd = i + 1)
#  }

## ----isodist-plot, echo = FALSE-----------------------------------------------
from = "2398957885" # 26 rows
dlim <- c (1, 2, 5, 10, 20) * 100
d <- dodgr_isodists (graph, from = from, dlim)
cols <- terrain.colors (length (dlim))
index <- which (d$dlim == max (d$dlim)) # plot max contour first
plot (d$x [index], d$y [index], "l", col = cols [1],
      xlab = "longitude", ylab = "latitude")
for (i in seq (dlim) [-1]) {
    index <- which (d$dlim == rev (dlim) [i])
    lines (d$x [index], d$y [index], col = cols [i], lwd = i + 1)
}

