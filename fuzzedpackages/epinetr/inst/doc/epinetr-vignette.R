## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

## ----message=FALSE, fig.width=5, fig.asp=1, fig.align='center', fig.cap="An epistatic network generated between 50 QTLs."----
library(epinetr)

# Build a population of size 1000, with 50 QTLs, broad-sense heritability of 0.6,
# narrow-sense heritability of 0.4 and overall trait variance of 50.
pop <- Population(popSize = 1000, map = map100snp, alleleFrequencies = runif(100),
                  QTL = 50, broadH2 = 0.6, narrowh2 = 0.4, traitVar = 50)

# Attach additive effects
pop <- addEffects(pop)

# Attach an epistatic network
pop <- attachEpiNet(pop)

# Plot the network
plot(getEpiNet(pop))

## -----------------------------------------------------------------------------
# Inspect initial phenotypic components
head(getComponents(pop))

## ----message=FALSE, fig.width=7, fig.height=4, fig.align='center', fig.cap="A graphical representation of a simulation run across 250 generations."----
# Run a simulation across 250 generations
pop <- runSim(pop, generations = 250, selection = "ranking")

# Plot the simulation run
plot(pop)

## -----------------------------------------------------------------------------
# Get the allele frequencies
af <- getAlleleFreqRun(pop)

# Get the phased genotypes of the resulting population
geno <- getPhased(pop)

# Get a subset of the resulting population
ID <- getComponents(pop)$ID
ID <- sample(ID, 50)
pop2 <- getSubPop(pop, ID)

## -----------------------------------------------------------------------------
head(map100snp)

## -----------------------------------------------------------------------------
nrow(map100snp)

## -----------------------------------------------------------------------------
length(unique(map100snp[, 2]))

## ----message=FALSE------------------------------------------------------------
pop <- Population(popSize = 200, map = map100snp, QTL = 20,
                  alleleFrequencies = runif(100),
                  broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40)
pop

## ----message=FALSE------------------------------------------------------------
pop <- Population(popSize = 200, map = map100snp, QTL = 20,
                  alleleFrequencies = runif(100))
pop

## ----message=FALSE------------------------------------------------------------
pop <- Population(popSize = 200, map = map100snp,
                  QTL = c(62, 55, 92, 74, 11, 38),
                  alleleFrequencies = runif(100),
                  broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40)
pop

## -----------------------------------------------------------------------------
getQTL(pop)

## -----------------------------------------------------------------------------
dim(geno100snp)

## -----------------------------------------------------------------------------
geno100snp[1, 1:10]

## ----message=FALSE------------------------------------------------------------
pop <- Population(popSize = nrow(geno100snp), map = map100snp, QTL = 20,
                  genotypes = geno100snp,
                  broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40)

## ----message=FALSE------------------------------------------------------------
pop <- Population(popSize = nrow(geno100snp), map = map100snp, QTL = 20,
                  genotypes = geno100snp, literal = FALSE,
                  broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40)

## ----message=FALSE------------------------------------------------------------
pop <- Population(pop, broadH2 = 0.7, traitVar = 30)
pop

## ----message=FALSE------------------------------------------------------------
pop <- Population(pop, popSize = 800)
pop

## ----message=FALSE------------------------------------------------------------
pop <- Population(pop, alleleFrequencies = runif(100))
pop

## ----message=FALSE------------------------------------------------------------
pop <- addEffects(pop)
pop

## ----message=FALSE------------------------------------------------------------
pop <- addEffects(pop, distrib = runif)

## ----message=FALSE------------------------------------------------------------
effects <- c( 1.2,  1.5, -0.3, -1.4,  0.8,
              2.4,  0.2, -0.8, -0.4,  0.8,
             -0.2, -1.4,  1.4,  0.2, -0.9,
              0.4, -0.8,  0.0, -1.1, -1.3)
pop <- addEffects(pop, effects = effects)
getAddCoefs(pop)

## ----message=FALSE------------------------------------------------------------
pop <- Population(pop, narrowh2 = 0.4)
getAddCoefs(pop)

## ----message=FALSE------------------------------------------------------------
pop <- attachEpiNet(pop)
pop

## ----fig.width=5, fig.asp=1, fig.align='center', fig.cap="A random epistatic network generated between 20 QTLs."----
epinet <- getEpiNet(pop)
plot(epinet)

## ----message=FALSE, fig.width=5, fig.asp=1, fig.align='center', fig.cap="An epistatic network generated using the Barabasi-Albert model between 20 QTLs."----
pop <- attachEpiNet(pop, scaleFree = TRUE)
plot(getEpiNet(pop))

## ----message=FALSE, fig.width=5, fig.asp=1, fig.align='center', fig.cap="An epistatic network where 7 QTLs have no epistatic effects."----
pop <- attachEpiNet(pop, scaleFree = TRUE, additive = 7)
plot(getEpiNet(pop))

## ----message=FALSE, fig.width=5, fig.asp=1, fig.align='center', fig.cap="An epistatic network with a minimum of 2 interactions per epistatic QTL."----
pop <- attachEpiNet(pop, scaleFree = TRUE, additive = 7, m = 2)
plot(getEpiNet(pop))

## ----message=FALSE, fig.width=5, fig.asp=1, fig.align='center', fig.cap="An epistatic network featuring 2-way to 7-way interactions"----
pop <- attachEpiNet(pop, scaleFree = TRUE, additive = 7, m = 2, k=2:7)
plot(getEpiNet(pop))

## -----------------------------------------------------------------------------
inc <- getIncMatrix(pop)
dim(inc)

## ----echo=1-------------------------------------------------------------------
inc[, 1:5]
inci <- which(inc[, 1:5] == 1)
incs <- rowSums(inc[, 1:5])
incm <- which(incs == max(incs))

## -----------------------------------------------------------------------------
rincmat100snp

## ----message=FALSE------------------------------------------------------------
pop <- attachEpiNet(pop, incmat = rincmat100snp)

## ---- fig.width=5, fig.asp=1, fig.align='center', fig.cap="An epistatic network derived from a user-defined incidence matrix."----
plot(getEpiNet(pop))

## ----message=FALSE------------------------------------------------------------
# Include the 20th QTL in the first interaction
mm <- rincmat100snp
mm[20, 1] <- 1
pop <- attachEpiNet(pop, incmat = mm)

## ---- fig.width=5, fig.asp=1, fig.align='center', fig.cap="A user-defined epistatic network featuring a single 3-way interaction."----
plot(getEpiNet(pop))

## ----message=FALSE------------------------------------------------------------
# 20 coefficients, 3 of which are 0
coefs <- sample(c(rep(0, 3), rnorm(17)), 20)
pop <- addEffects(pop, effects = coefs)
getAddCoefs(pop)

## ---- echo=FALSE--------------------------------------------------------------
foo <- which(getAddCoefs(pop) == 0)

## ----message=FALSE------------------------------------------------------------
pop <- Population(pop, QTL = 1:15)

## ----message=FALSE------------------------------------------------------------
pop <- attachEpiNet(pop, additive = 1:5)

## -----------------------------------------------------------------------------
getIncMatrix(pop)

## ---- fig.width=5, fig.asp=1, fig.align='center', fig.cap="The epistatic network generated with the first 5 of the 15 QTLs being purely additive."----
plot(getEpiNet(pop))

## ----message=FALSE------------------------------------------------------------
coefs <- rnorm(15)
coefs[6:10] <- 0
pop <- addEffects(pop, effects = coefs)
getAddCoefs(pop)

## ----message=FALSE------------------------------------------------------------
pop <- Population(popSize = 200, map = map100snp, QTL = 20,
                  alleleFrequencies = runif(100),
                  broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40)
pop <- addEffects(pop)
pop <- attachEpiNet(pop, k = 3)

## ---- fig.width=5, fig.asp=1, fig.align='center', fig.cap="An epistatic network consisting of 3-way interactions."----
plot(getEpiNet(pop))

## -----------------------------------------------------------------------------
qtls <- which(getIncMatrix(pop)[, 1] > 0)
qtls

## -----------------------------------------------------------------------------
interaction1 <- getInteraction(pop, 1) # Return first interaction array
interaction1

## -----------------------------------------------------------------------------
interaction1[2, , ]

## -----------------------------------------------------------------------------
interaction1[2, 1, ]

## -----------------------------------------------------------------------------
interaction1[2, 1, 3]

## -----------------------------------------------------------------------------
components <- getComponents(pop)
head(components)

## -----------------------------------------------------------------------------
mean(components$Additive)

## -----------------------------------------------------------------------------
var(components$Additive)

## -----------------------------------------------------------------------------
mean(components$Epistatic)

## -----------------------------------------------------------------------------
var(components$Epistatic)

## -----------------------------------------------------------------------------
mean(components$Environmental)

## -----------------------------------------------------------------------------
var(components$Environmental)

## -----------------------------------------------------------------------------
getAddOffset(pop)

## -----------------------------------------------------------------------------
getEpiOffset(pop)

## -----------------------------------------------------------------------------
mean(components$Phenotype)

## -----------------------------------------------------------------------------
var(components$Phenotype)

## -----------------------------------------------------------------------------
cov(components$Additive, components$Epistatic)
cov(components$Additive, components$Environmental)
cov(components$Environmental, components$Epistatic)

## -----------------------------------------------------------------------------
cor(components$Additive, components$Epistatic)

## -----------------------------------------------------------------------------
geno <- getGeno(pop)

## -----------------------------------------------------------------------------
geno <- geno[, getQTL(pop)$Index]

## -----------------------------------------------------------------------------
additive <- geno %*% getAddCoefs(pop) + getAddOffset(pop)
additive[1:5]

## -----------------------------------------------------------------------------
getComponents(pop)$Additive[1:5]

## -----------------------------------------------------------------------------
head(getEpistasis(pop))

## -----------------------------------------------------------------------------
epistatic <- rowSums(getEpistasis(pop)) + getEpiOffset(pop)
epistatic[1:5]

## -----------------------------------------------------------------------------
getComponents(pop)$Epistatic[1:5]

## ----include=FALSE------------------------------------------------------------
geno2 <- matrix(sample(0:2, 100*5, replace = TRUE, prob = c(0.25, 0.5, 0.25)), 5, 100)

## -----------------------------------------------------------------------------
geno2 <- geno2[, getQTL(pop)$Index]
geno2

## -----------------------------------------------------------------------------
additive2 <- geno2 %*% getAddCoefs(pop) + getAddOffset(pop)
additive2[1:5]

## -----------------------------------------------------------------------------
epistatic2 <- rowSums(getEpistasis(pop, geno = geno2)) + getEpiOffset(pop)
epistatic2

## ----message=FALSE------------------------------------------------------------
popRun <- runSim(pop, generations = 150)

## -----------------------------------------------------------------------------
n <- 10
pmf <- 2 * (n - 1:n + 1) / (n * (n + 1))
pmf

## ----echo=FALSE, fig.width=4, fig.asp=1, fig.align='center', fig.cap="Probability mass function for linear ranking selection on a population of 10 individuals."----
df = data.frame(Individual=as.character(1:n), Probability=pmf)
ggplot2::ggplot(data = df, ggplot2::aes(x=reorder(Individual, -Probability), y=Probability)) + ggplot2::geom_bar(stat="identity") + ggplot2::xlab("Ranked individuals") + ggplot2::ylab("Selection probability")

## ----message=FALSE------------------------------------------------------------
popRunRank <- runSim(pop, generations = 150, selection = "ranking")
popRunBurnIn <- runSim(pop, generations = 150, burnIn = 50,
                       selection = "ranking",
                       truncSire = 0.25, truncDam = 0.25,
                       roundsSire = 5, roundsDam = 5,
                       litterDist = c(0.1, 0.3, 0.4, 0.2),
                       breedSire = 7)

## ----fig.width=7, fig.height=4, fig.align='center', fig.cap="A graphical representation of a simulation run using the default parameters."----
plot(popRun)

## ----fig.width=7, fig.height=4, fig.align='center', fig.cap="A graphical representation of a simulation run using linear ranking selection."----
plot(popRunRank)

## ----fig.width=7, fig.height=4, fig.align='center', fig.cap="A graphical representation of a simulation run using ranking selection and a burn-in period of the first 50 generations."----
plot(popRunBurnIn)

## ----echo=FALSE---------------------------------------------------------------
allGenoFileName <- system.file("extdata", "geno.epi", package = "epinetr")
geno <- loadGeno(allGenoFileName)

## ----eval=FALSE---------------------------------------------------------------
#  popRun <- runSim(pop, generations = 150, allGenoFileName = "geno.epi")
#  geno <- loadGeno("geno.epi")

## -----------------------------------------------------------------------------
geno[1:5, 1:8]

## -----------------------------------------------------------------------------
ped <- getPedigree(popRun)
ped[512:517, ]

## -----------------------------------------------------------------------------
qtl <- getQTL(popRun)$Index
af <- getAlleleFreqRun(popRun)
af[, qtl[1]]

## -----------------------------------------------------------------------------
geno <- getPhased(popRun)
geno[1:6, 1:10]

## -----------------------------------------------------------------------------
geno <- getGeno(popRun)
geno[1:6, 1:5]

## -----------------------------------------------------------------------------
ID <- getComponents(popRun)$ID
ID <- sample(ID, 50)
popRun2 <- getSubPop(popRun, ID)

## ----include=FALSE------------------------------------------------------------
pedData <- getPedigree(popRun)
pedData <- pedData[,1:3]

## -----------------------------------------------------------------------------
pedData[201:210, ]

## ----message=FALSE------------------------------------------------------------
popRunPed <- runSim(pop, pedigree = pedData)

## ----fig.width=7, fig.height=4, fig.align='center', fig.cap="A graphical representation of a simulation run using a pedigree data frame."----
plot(popRunPed)

