## ----echo = TRUE, message = FALSE, warning = FALSE----------------------------
library(dplyr)
library(molic)
y     <- unlist(derma[80, -35]) # a patient with seboreic dermatitis
psor  <- derma %>%
  filter(ES == "psoriasis") %>%
  dplyr::select(-ES)

## ----fit_graph----------------------------------------------------------------
library(ess)
g <- fit_graph(psor, q = 0, trace = FALSE)

## ----color_nodes, echo = TRUE-------------------------------------------------
vs   <- names(adj_lst(g))
vcol <- structure(vector("character", length(vs)), names = vs)
vcol[grepl("c", vs)] <- "tomato"  # clinical attributes
vcol[grepl("h", vs)] <- "#98FB98" # histopathological attributes
vcol["age"]          <- "gray"    # age variable

## ----plot_graph, fig.width=7, fig.height=4.5,fig.show='hold',fig.align='center'----
plot(g, vcol, vertex.size = 10, vertex.label = NA)

## ----comment = NA-------------------------------------------------------------
set.seed(300718)
m <- fit_outlier(psor, g, y)
print(m)

## ----plot_outlier, fig.align='center', fig.width=5, fig.height=3,fig.show='hold'----
plot(m)

## ----echo = TRUE, eval = FALSE, fit_multiple_models, fig.align='center', message = FALSE, warning = FALSE, fig.width=5, fig.height=2.8,fig.show='hold'----
#  set.seed(300718)
#  mm <- fit_multiple_models(derma, y, "ES", q = 0,trace = FALSE)
#  plot(mm)

## ----mult_models, echo = FALSE, out.width = "70%", fig.align = 'center'-------
knitr::include_graphics("multiple_models.png")

