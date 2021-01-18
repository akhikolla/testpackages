## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE, message=FALSE-----------------------------------------------
#  library(landsepi)

## -----------------------------------------------------------------------------
myDesign <- data.frame(nTSpY = c(120, 110, 100, 90, 80)
                       , pI0 = c(0, 1E-4, 5E-4, 1E-3, 1E-2)
                       , infection_rate = c(0.5, 0.4, 0.3, 0.2, 0.1)
                       , id_landscape = 1:5
                       , aggreg = c(0.07, 0.07, 0.25, 10, 10)
                       , R_efficiency = c(1.00, 0.90, 0.80, 0.70, 0.60)
                       , growth_rate = c(0.1, 0.2, 0.3, 0.1, 0.1)
                       ## create columns to store outputs
                       , durab_MG1 = NA
                       , durab_MG2 = NA
                       , mean_audpc = NA)

## -----------------------------------------------------------------------------
n <- nrow(myDesign)
myDesign <- cbind(simul = 1:n, seed = sample(n*100, n), myDesign)
myDesign

## ---- eval=FALSE--------------------------------------------------------------
#  simul_params <- createSimulParams(outputDir = getwd())

## ---- eval=FALSE--------------------------------------------------------------
#  for (i in 1:n){
#      print(paste("Running simulation", i, "/", n))
#  
#      ## Set Nyears and nTSpY
#      simul_params <- setTime(simul_params
#                              , Nyears = simul_params@TimeParam$Nyears
#                              , nTSpY = myDesign$nTSpY[i])  ## update nTSpY
#  
#      ## set seed (to run stochastic replicates)
#      simul_params <- setSeed(simul_params, myDesign$seed[i])  ## update seed
#  
#      ## Pathogen parameters
#      basic_patho_param <- loadPathogen("rust")
#      basic_patho_param$infection_rate <- myDesign$infection_rate[i]  ## update inf. rate
#      simul_params <- setPathogen(simul_params, basic_patho_param)
#  
#      ## Initial conditions
#      simul_params <- setInoculum(simul_params, myDesign$pI0[i])  ## update pI0
#  
#      ## Landscape parameters
#      landscape <- loadLandscape(myDesign$id_landscape[i])   ## update landscape
#      simul_params <- setLandscape(simul_params, landscape)
#  
#      ## Dispersal parameters
#      disp_patho <- loadDispersalPathogen(myDesign$id_landscape[i])  ## update dispersal
#      simul_params <- setDispersalPathogen(simul_params, disp_patho)
#  
#      ## Genes
#      gene1 <- loadGene(name = "MG 1", type = "majorGene")
#      gene2 <- loadGene(name = "MG 2", type = "majorGene")
#      genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
#      genes$efficiency <- myDesign$R_efficiency[i]  ## update resistance efficiency
#      simul_params <- setGenes(simul_params, genes)
#  
#      ## Cultivars
#      cultivar1 <- loadCultivar(name = "Susceptible", type = "growingHost")
#      cultivar2 <- loadCultivar(name = "Resistant1", type = "growingHost")
#      cultivar3 <- loadCultivar(name = "Resistant2", type = "growingHost")
#      cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3)
#                              , stringsAsFactors = FALSE)
#      cultivars$growth_rate <- myDesign$growth_rate[i]  ## update growth rate
#      simul_params <- setCultivars(simul_params, cultivars)
#  
#      ## Allocate genes to cultivars
#      simul_params <- allocateCultivarGenes(simul_params, "Resistant1", c("MG 1"))
#      simul_params <- allocateCultivarGenes(simul_params, "Resistant2", c("MG 2"))
#  
#      ## Allocate cultivars to croptypes
#      croptypes <- loadCroptypes(simul_params, names = c("Susceptible crop"
#                                                         , "Resistant crop 1"
#                                                         , "Resistant crop 2"))
#      croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
#      croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
#      croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
#      simul_params <- setCroptypes(simul_params, croptypes)
#  
#      ## Allocate croptypes to landscape
#      rotation_sequence <- croptypes$croptypeID ## No rotation: 1 rotation_sequence element
#      rotation_period <- 0 ## same croptypes every years
#      prop <- c(1/3, 1/3, 1/3) ## croptypes proportions
#      aggreg <- myDesign$aggreg[i]
#      simul_params <- allocateLandscapeCroptypes(simul_params,
#                                                 rotation_period = rotation_period,
#                                                 rotation_sequence = rotation_sequence,
#                                                 rotation_realloc = FALSE,
#                                                 prop = prop,
#                                                 aggreg = aggreg,
#                                                 graphic = FALSE)
#  
#      ## configure outputs
#      outputlist <- loadOutputs(epid_outputs = "audpc", evol_outputs = "durability")
#      simul_params <- setOutputs(simul_params, outputlist)
#  
#      ## Check, (save) and run simulation
#      checkSimulParams(simul_params)
#      # simul_params <- saveDeploymentStrategy(simul_params, overwrite = TRUE)
#      res <- runSimul(simul_params, writeTXT=FALSE, graphic = FALSE)
#  
#      ## Extract outputs
#      myDesign$durab_MG1[i]  <- res$evol_outputs$durability[,"MG 1"]
#      myDesign$durab_MG2[i]  <- res$evol_outputs$durability[,"MG 2"]
#      myDesign$mean_audpc[i] <- mean(res$epid_outputs$audpc$total)
#  
#      ## Create myDesign.txt at first simulation and then append
#      if (i ==1){
#          write.table(myDesign[i,], paste(simul_params@OutputDir, "/myDesign.txt", sep="")
#                      , append=FALSE, row.names = FALSE, col.names = TRUE)
#      } else {
#          write.table(myDesign[i,], paste(simul_params@OutputDir, "/myDesign.txt", sep="")
#                      , append=TRUE, row.names = FALSE, col.names = FALSE)
#      }
#  }
#  
#  myDesign

## ---- eval=FALSE--------------------------------------------------------------
#  ## Disable computation of outputs:
#  outputlist <- loadOutputs(epid_outputs = "", evol_outputs = "")
#  simul_params <- setOutputs(simul_params, outputlist)
#  
#  ## Run simulation and keep raw binary files:
#  runSimul(simul_params, writeTXT=FALSE, graphic = FALSE, keepRawResults = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  ##  retrieve parameters from the object simul_params
#  path <- simul_params@OutputDir
#  Nyears <- simul_params@TimeParam$Nyears
#  nTSpY <- simul_params$nTSpY
#  nTS <- Nyears * nTSpY     ## Total number of time-steps
#  
#  Npoly <- nrow(simul_params@Landscape)
#  Nhost <- nrow(simul_params@Cultivars)
#  Npatho <- Npatho <- prod(simul_params@Genes$Nlevels_aggressiveness)
#  
#  ## Initialise lists
#  H <- as.list(1:nTS)
#  Hjuv <- as.list(1:nTS)
#  P <- as.list(1:nTS)
#  L <- as.list(1:nTS)
#  I <- as.list(1:nTS)
#  R <- as.list(1:nTS)
#  index <- 0
#  
#  ## Read binary files and store values in the lists as matrices or arrays
#  for (year in 1:Nyears) {
#  
#      binfileH <- file(paste(path, sprintf("/H-%02d", year), ".bin", sep = ""), "rb")
#      H.tmp <- readBin(con = binfileH, what = "int", n = Npoly * Nhost * nTSpY, size = 4
#                       , signed = T, endian = "little")
#      close(binfileH)
#  
#      binfileHjuv = file(paste(path, sprintf("/Hjuv-%02d", year), ".bin",sep=""), "rb")
#      Hjuv.tmp <- readBin(con=binfileHjuv, what="int", n=Npoly*Nhost*nTSpY, size = 4
#                          , signed=T,endian="little")
#      close(binfileHjuv)
#  
#      binfileP <- file(paste(path, sprintf("/P-%02d", year), ".bin", sep = ""), "rb")
#      P.tmp <- readBin(con = binfileP, what = "int", n = Npoly * Npatho * nTSpY, size = 4
#                       , signed = T, endian = "little")
#      close(binfileP)
#  
#      binfileL <- file(paste(path, sprintf("/L-%02d", year), ".bin", sep = ""), "rb")
#      L.tmp <- readBin(con = binfileL, what = "int", n = Npoly * Npatho * Nhost * nTSpY
#                       , size = 4 , signed = T, endian = "little")
#      close(binfileL)
#  
#      binfileI <- file(paste(path, sprintf("/I-%02d", year), ".bin", sep = ""), "rb")
#      I.tmp <- readBin(con = binfileI, what = "int", n = Npoly * Npatho * Nhost * nTSpY
#                       , size = 4 , signed = T, endian = "little")
#      close(binfileI)
#  
#      binfileR <- file(paste(path, sprintf("/R-%02d", year), ".bin", sep = ""), "rb")
#      R.tmp <- readBin(con = binfileR, what = "int", n = Npoly * Npatho * Nhost * nTSpY
#                       , size = 4 , signed = T, endian = "little")
#      close(binfileR)
#  
#      ## Convert vectors in matrices or arrays
#      for (t in 1:nTSpY) {
#          H[[t + index]] <- matrix(H.tmp[((Nhost * Npoly) * (t-1)+1):(t * (Nhost * Npoly))]
#                                   , ncol = Nhost, byrow = T)
#          Hjuv[[t + index]] <- matrix(Hjuv.tmp[((Nhost*Npoly)*(t-1)+1):(t*(Nhost*Npoly))]
#                                      , ncol=Nhost,byrow=T)
#          P[[t + index]] <- matrix(P.tmp[((Npatho * Npoly) * (t-1)+1):(t * (Npatho * Npoly))]
#                                   , ncol = Npatho, byrow = T)
#          L[[t + index]] <- array(data = L.tmp[((Npatho * Npoly * Nhost) *
#                                                  (t-1)+1):(t * (Npatho * Npoly * Nhost))]
#                                  , dim = c(Nhost, Npatho, Npoly))
#          I[[t + index]] <- array(data = I.tmp[((Npatho * Npoly * Nhost) *
#                                                  (t-1)+1):(t * (Npatho * Npoly * Nhost))]
#                                  , dim = c(Nhost, Npatho, Npoly))
#          R[[t + index]] <- array(data = R.tmp[((Npatho * Npoly * Nhost) *
#                                                  (t-1)+1):(t * (Npatho * Npoly * Nhost))]
#                                  , dim = c(Nhost, Npatho, Npoly))
#      }
#  
#      index <- index + nTSpY
#  }

