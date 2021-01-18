## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.asp = 1,
  fig.width = 5
)

## ----setup, message=FALSE-----------------------------------------------------
library(ggplot2)
library(dplyr)
library(ipaddress)
library(ggip)

## ----before_transform, echo=FALSE---------------------------------------------
tibble(
  label = c("A", "B", "C"),
  address = ip_address(c("0.0.0.0", "192.168.0.1", "255.255.255.255"))
)

## ----after_transform, echo=FALSE----------------------------------------------
tibble(
  label = c("A", "B", "C"),
  address = tibble(
    ip = ip_address(c("0.0.0.0", "192.168.0.1", "255.255.255.255")),
    x = as.integer(c(0, 214, 255)),
    y = as.integer(c(255, 142, 255))
  )
)

## -----------------------------------------------------------------------------
tibble(address = ip_address(c("0.0.0.0", "128.0.0.0", "192.168.0.1"))) %>%
  ggplot(aes(x = address$x, y = address$y, label = address$ip)) +
  geom_point() +
  geom_label(nudge_x = c(10, 0, -10), nudge_y = -10) +
  coord_ip(expand = TRUE) +
  theme_ip_light()

## ---- fig.asp=0.8, fig.width=6.25---------------------------------------------
iana_ipv4 %>%
  ggplot(aes(xmin = network$xmin, ymin = network$ymin, xmax = network$xmax, ymax = network$ymax)) +
  geom_rect(aes(fill = allocation)) +
  scale_fill_brewer(palette = "Accent", name = NULL) +
  coord_ip() +
  theme_ip_dark()

## -----------------------------------------------------------------------------
tibble(address = sample_ipv4(10000)) %>%
  ggplot(aes(ip = address)) +
  stat_summary_address() +
  scale_fill_viridis_c(guide = "none") +
  coord_ip() +
  theme_ip_dark()

