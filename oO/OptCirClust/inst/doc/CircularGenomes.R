## ----  results='hide', message=FALSE, warning=FALSE, echo=FALSE---------------
library(OptCirClust)
library(ape)
library(bazar)
library(knitr)
library(graphics)

opar <- par(mar=c(0,0,2,0))

opts_chunk$set(fig.width=6, fig.height=4) 

Event <- "CG"

K <- 14

# Seq <- read.GenBank("CP019943.1", as.character = TRUE)[[1]]
file <- system.file("extdata", "CP019943.1.fasta", package = "OptCirClust")

Seq <- read.dna(file, format="fasta", as.matrix=FALSE, as.character = TRUE)

Seq <- toupper(concat0(Seq))

V <- gregexpr(Event, Seq)

O <- sort(V[[1]][1:length(V[[1]])])

Circumference <- nchar(Seq)

set.seed(1)

result <- CirClust(O, K, Circumference, method = "FOCC")

plot(result, main = "Optimal circular clustering")

# arrows(.58, - 1.75, 0.48, -1.45, length = 0.125, angle = 30, code = 2, col="orange", lwd=4)
# arrows(0, -10, 0, 0, length = 0.125, angle = 30, code = 2, col="orange", lwd=4)
arrows(0.167, -0.55, 0,-0.145, length = 0.125, angle = 30, code = 1, col="orange", lwd=4)

result_km <- CirClust(O, K, Circumference, method = "HEUC")

plot(result_km, main = "Heuristic circular clustering",)

# arrows(.58, - 1.75, 0.4, -1.5, length = 0.125, angle = 30, code = 2, col="orange", lwd=4)

arrows(0.135, -0.55, 0,-0.145, length = 0.125, angle = 30, code = 1, col="orange", lwd=4)

par(opar)

