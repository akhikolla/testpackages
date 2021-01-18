## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(data.table)
print(spread::norway_seiiar_noinfected_2017)

## ------------------------------------------------------------------------
vax_measles <- fhidata::norway_childhood_vax[
  year==2016 & 
  stringr::str_detect(location_code,"^municip") & 
  vax=="measles",
  c("location_code","proportion")
  ]

print(vax_measles)

## ------------------------------------------------------------------------
norway_seiiar_measles_noinfected_2017 <- spread::convert_blank_seiiar_with_vax(
  seiiar = spread::norway_seiiar_noinfected_2017, 
  vax = vax_measles
  )

## ------------------------------------------------------------------------
print(norway_seiiar_measles_noinfected_2017)

## ------------------------------------------------------------------------
norway_seiiar_measles_oslo_2017 <- copy(norway_seiiar_measles_noinfected_2017)
norway_seiiar_measles_oslo_2017[location_code == "municip0301", I := 10]
norway_seiiar_measles_oslo_2017[location_code == "municip0301", S := S - I]

print(norway_seiiar_measles_oslo_2017[location_code=="municip0301"])  

