## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, include = TRUE, eval = TRUE)

## ----admission, echo = FALSE, message = FALSE, warning = FALSE----------------
## https://stat.ethz.ch/R-manual/R-patched/library/datasets/html/UCBAdmissions.html
library(dplyr)

adm_ <- as.data.frame(UCBAdmissions) 
adm <- data.frame(Admit = c(), Gender = c(), Dept = c())

for (i in 1:nrow(adm_)) {
  freq_i <- adm_[i, "Freq"]
  for (j in 1:freq_i) {
    adm <- rbind(adm, adm_[i, c("Admit", "Gender", "Dept")])
  }
}

adm <- adm[sample(nrow(adm)), ] %>%
  as_tibble() %>%
  mutate_all(.funs = as.character)

pander::pander(head(adm))

## ----cont_adm, echo = FALSE---------------------------------------------------
pander::pander(ftable(adm))

## ----undirected_graph, echo = FALSE, message = FALSE, warning = FALSE, fig.align = "center", fig.cap = "\\label{fig:DG} $G$: An undirected decomposable graph."----
# library(igraph)
el <- matrix( c("a", "b",
  "b", "c",
  "b", "d",
  "d", "e",
  "d", "c"),
  nc = 2,
  byrow = TRUE
)
G <- igraph::graph_from_edgelist(el, directed = FALSE)
set.seed(7)
plot(G, layout = igraph::layout_nicely(G), 
     edge.arrow.size=1, 
     vertex.label.cex=1.5, 
     # vertex.label.family="Helvetica",
     vertex.label.font=1,
     vertex.shape="circle", 
     vertex.size=40, 
     vertex.label.color="black", 
     edge.width=5)

