context("test-LandsepiParams-Method")

options(warn = -1)

params <- createSimulParams()
params <- setSeed(params, seed = 12345)
params <- setTime(params, Nyears = 3, nTSpY = 120)

land_test <- st_read("test.shp")
land_test$Aire <- NULL
land_test$year_1 <- c(0, 1)
land_test$year_2 <- c(1, 0)
land_test$year_3 <- c(2, 1)
land_test$year_4 <- c(2, 0)

croptypes <- data.frame(croptypeID = c(0, 1, 2), croptypeName = c("crop1", "crop2", "crop3"), Susceptible = c(1.0, 0, 0), Resistant1 = c(0, 1.0, 0), Resistant2 = c(0, 0, 1.0))

gene1 <- loadGene(name = "MG1", type = "majorGene")
gene2 <- loadGene(name = "MG2", type = "majorGene")
genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)

cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
# params <- setCultivarGene(params, "Resistant1", c("MG1"))
# params <- setCultivarGene(params, "Resistant2", c("MG2"))

basic_patho_param <- list(
  name = "rust",
  survival_prob = 1e-4,
  repro_sex_prob = 0, ## probability for an infection that its reproduction is sexual rather than clonal.
  infection_rate = 0.4,
  propagule_prod_rate = 3.125,
  latent_period_exp = 10,
  latent_period_var = 9,
  infectious_period_exp = 24,
  infectious_period_var = 105,
  sigmoid_kappa = 5.333,
  sigmoid_sigma = 3,
  sigmoid_plateau = 1
)

test_that("Default values", {
  expect_equal(normalizePath(dirname(params@OutputDir)), normalizePath(test_path()))
  expect_equal(params@OutputGPKG, "landsepi_landscape.gpkg")
  expect_equal(params@TimeParam$Nyears, 3)
  expect_equal(params@TimeParam$nTSpY, 120)
  expect_equal(params@Seed, 12345)
})

test_that("Landscape", {
  land_test_sp <- as_Spatial(land_test)

  expect_equal(nrow(params@Landscape), 0)
  expect_equal(class(params@Landscape)[1], "sf")

  params <- setLandscape(params, land_test_sp)
  expect_equal(nrow(params@Landscape), 2)
  expect_equal(class(params@Landscape)[1], "sf")

  params <- setLandscape(params, land_test)
  expect_equal(nrow(params@Landscape), 2)
  expect_equal(class(params@Landscape)[1], "sf")
})

test_that("Landscape-Check", {
  params_tmp <- setCroptypes(params, croptypes)
  params_tmp <- setLandscape(params_tmp, land_test)

  expect_true(checkLandscape(params_tmp))

  params_tmp@Landscape$year_5 <- c(2, 0)
  expect_false(res <- checkLandscape(params_tmp))

  params_tmp@Landscape$year_5 <- NULL
  params_tmp@Landscape$year_4 <- NULL

  expect_false(res <- checkLandscape(params_tmp))

  params_tmp@Landscape$year_4 <- c(5, 0)

  expect_false(res <- checkLandscape(params_tmp))
})

test_that("Disp-Check", {
  params_tmp <- params
  params_tmp <- setLandscape(params_tmp, land_test)
  params_tmp@DispHost <- c(0, 0, 0, 0)
  params_tmp@DispPatho <- c(0, 0, 0, 0)

  expect_true(checkDispersalHost(params_tmp))
  expect_true(checkDispersalPathogen(params_tmp))

  params_tmp@DispHost <- c(0, 0, 0)
  params_tmp@DispPatho <- c(0, 0, 0)
  expect_false(checkDispersalHost(params_tmp))
  expect_false(checkDispersalPathogen(params_tmp))

  params_tmp@DispHost <- c(-10, 0, 0, 0)
  params_tmp@DispPatho <- c(-10, 0, 0, 0)
  expect_false(checkDispersalHost(params_tmp))
  expect_false(checkDispersalPathogen(params_tmp))

  params_tmp@DispHost <- c(0, 0, 0, 1.1)
  params_tmp@DispPatho <- c(0, 0, 0, 1.1)
  expect_false(checkDispersalHost(params_tmp))
  expect_false(checkDispersalPathogen(params_tmp))
})

test_that("Genes-Check", {
  params_tmp <- params
  params_tmp@Genes <- genes

  expect_true(checkGenes(params_tmp))

  params_tmp@Genes <- cbind(genes, fakecol <- c(2, 3))
  expect_true(checkGenes(params_tmp))

  params_tmp@Genes <- genes[, -1]
  expect_false(checkGenes(params_tmp))

  params_tmp@Genes <- genes[, -4]
  expect_false(checkGenes(params_tmp))
})

test_that("Croptypes-Check", {
  params_tmp <- params
  params_tmp <- setCroptypes(params_tmp, croptypes)
  params_tmp <- setLandscape(params_tmp, land_test)

  expect_true(checkCroptypes(params_tmp))

  # croptypes not in landscape
  croptypes <- data.frame(croptypeID = c(0, 1, 2, 3)
                          , croptypeName = c("crop1", "crop2", "crop3", "crop4")
                          , Susceptible = c(1.0, 0, 0, 0)
                          , Resistant1 = c(0, 1.0, 0, 0)
                          , Resistant2 = c(0, 0, 1.0, 0))
  params_tmp <- setCroptypes(params_tmp, croptypes)
  expect_false(checkCroptypes(params_tmp))

  # less croptypes than in landscape
  croptypes <- data.frame(croptypeID = c(0, 1)
                          , croptypeName = c("crop1", "crop2")
                          , Susceptible = c(1.0, 0)
                          , Resistant1 = c(0, 1.0)
                          , Resistant2 = c(0, 0))
  params_tmp <- setCroptypes(params_tmp, croptypes)
  expect_true(checkCroptypes(params_tmp))

  # cultivars proportions not equal to 1
  croptypes <- data.frame(croptypeID = c(0, 1, 2)
                          , croptypeName = c("crop1", "crop2", "crop3")
                          , Susceptible = c(1.1, 0, 0)
                          , Resistant1 = c(0, -1.0, 0)
                          , Resistant2 = c(0, 0, 0.5))
  params_tmp <- setCroptypes(params_tmp, croptypes)
  expect_false(checkCroptypes(params_tmp))

  params_tmp <- setCultivars(params_tmp, cultivars)
  croptypes <- data.frame(croptypeID = c(0, 1, 2)
                          , croptypeName = c("crop1", "crop2", "crop3")
                          , Susceptible = c(1.0, 0, 0)
                          , Resistant1 = c(0, 1.0, 0)
                          , Resistant2 = c(0, 0, 1.0)
                          , Resistant3 = c(0, 0, 0))
  params_tmp <- setCroptypes(params_tmp, croptypes)
  expect_false(checkCroptypes(params_tmp))

  croptypes <- data.frame(croptypeID = c(0, 1, 2)
                          , croptypeName = c("crop1", "crop2", "crop3")
                          , Susceptible = c(1.0, 0, 0)
                          , Resistant1 = c(0, 1.0, 0)
                          , Resistant3 = c(0, 0, 1))
  params_tmp <- setCroptypes(params_tmp, croptypes)
  expect_false(checkCroptypes(params_tmp))
})

test_that("Cultivars-Check", {
  params_tmp <- params
  params_tmp <- setCultivars(params_tmp, cultivars)
  params_tmp <- setCroptypes(params_tmp, croptypes)

  expect_true(checkCultivars(params_tmp))

  # Cultivars colnames
  expect_error(params_tmp <- setCultivars(params_tmp, cultivars[, -1]))
  params_tmp@Cultivars <- cultivars[, -1]
  expect_false(checkCultivars(params_tmp))

  params_tmp <- setCultivars(params_tmp, cbind(cultivars, fakecol <- c(1, 2, 3)))
  expect_true(checkCultivars(params_tmp))

  # cultivars nb -1
  params_tmp <- setCultivars(params_tmp, cultivars[-1, ])
  expect_false(checkCultivars(params_tmp))

  # cultivars nb +1
  params_tmp <- setCultivars(params_tmp, rbind(cultivars, loadCultivar(name = "Resistant3", type = "growingHost")))
  expect_true(checkCultivars(params_tmp))
})

test_that("CultivarGene-Check", {
  params_tmp <- params
  params_tmp <- setCultivars(params_tmp, cultivars)
  params_tmp <- setCroptypes(params_tmp, croptypes)
  params_tmp <- setGenes(params_tmp, genes)
  params_tmp <- allocateCultivarGenes(params_tmp, "Resistant1", c("MG1"))
  params_tmp <- allocateCultivarGenes(params_tmp, "Resistant2", c("MG2"))

  expect_true(checkCultivarsGenes(params_tmp))

  params_tmp@CultivarsGenes <- cbind(params_tmp@CultivarsGenes, MG3 <- c(1, 1, 1))
  expect_false(checkCultivarsGenes(params_tmp))

  params_tmp@CultivarsGenes <- params_tmp@CultivarsGenes[-1, -3]
  expect_false(checkCultivarsGenes(params_tmp))
})

test_that("Pathogen-Check", {
  params_tmp <- params

  expect_true(checkPathogen(params_tmp))

  params_tmp@Pathogen <- basic_patho_param
  expect_true(checkPathogen(params_tmp))
})
