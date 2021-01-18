## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = FALSE,
  message = FALSE
)

## ----setup, warning=FALSE-----------------------------------------------------
library(ggip)

## ---- out.width="100%"--------------------------------------------------------
knitr::include_graphics("bits_raw.png")

## ---- out.width="100%"--------------------------------------------------------
knitr::include_graphics("bits_half_reduced.png")

## ---- out.width="100%"--------------------------------------------------------
knitr::include_graphics("bits_reduced.png")

## ----plot_func----------------------------------------------------------------
ordinal_suffix <- function(x) {
    suffix <- c("st", "nd", "rd", rep("th", 17))
    suffix[((x-1) %% 10 + 1) + 10*(((x %% 100) %/% 10) == 1)]
}

plot_curve <- function(curve, curve_order) {
  pixel_prefix <- 32L
  canvas_prefix <- as.integer(pixel_prefix - (2 * curve_order))
  canvas_network <- ip_network(ip_address("0.0.0.0"), canvas_prefix)
  n_pixels <- 2^curve_order
  
  ggplot(data.frame(address = seq(canvas_network))) +
    geom_path(aes(address$x, address$y)) +
    coord_ip(
      canvas_network = canvas_network,
      pixel_prefix = pixel_prefix,
      curve = curve,
      expand = TRUE
    ) +
    theme_ip_light() +
    labs(title = paste0(
      curve_order, ordinal_suffix(curve_order),
      " order (", n_pixels, "x", n_pixels, " grid)"
    ))
}

## ----hilbert, fig.show="hold", out.width="30%"--------------------------------
plot_curve("hilbert", 2)
plot_curve("hilbert", 3)
plot_curve("hilbert", 4)

## ----morton, fig.show="hold", out.width="30%"---------------------------------
plot_curve("morton", 2)
plot_curve("morton", 3)
plot_curve("morton", 4)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  coord_ip(
#    canvas_prefix = ip_network("0.0.0.0/0"),
#    pixel_prefix = 4,
#    curve = "hilbert"
#  )

## ---- fig.align="center", fig.asp=1, fig.width=5------------------------------
curve_order <- 2
pixel_prefix <- 2 * curve_order
vertices <- subnets(ip_network("0.0.0.0/0"), new_prefix = pixel_prefix)

data <- data.frame(ip = network_address(vertices), label = as.character(vertices))
nudge <- c(1, 0, 0, 1, 1, 1, 0, 0, 0, 0, -1, -1, -1, 0, 0, -1)

ggplot(data, aes(ip$x, ip$y)) +
  geom_path() +
  geom_label(aes(label = label), nudge_x = 0.2 * nudge) +
  coord_ip(pixel_prefix = pixel_prefix, expand = TRUE) +
  theme_ip_light() +
  labs(title = paste0("Hilbert curve: ", curve_order, ordinal_suffix(curve_order), " order"))

