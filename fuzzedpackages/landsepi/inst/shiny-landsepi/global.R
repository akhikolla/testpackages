#library(shinycssloaders)
library(shiny)
library(DT)
library(shinyjs)
#library(slickR)
library(gridExtra)
library(png)
library(grid)
library(future)
library(promises)
library(tools)


library("landsepi")
data(package = "landsepi")

## del all file and directory of a path
cleanDir <- function(path) {
  files <- dir(path, full.names=TRUE, no..=TRUE)
  lapply(files, FUN = function(file){
      if( dir.exists(file) ) cleanDir(file)
      file.remove(file)
    })
}

ROOT_PATH <- getwd()

if(!dir.exists(paste0(ROOT_PATH,"/www/tmp/"))) dir.create(paste0(ROOT_PATH,"/www/tmp/"))
setwd(paste0(ROOT_PATH,"/www/tmp/"))

cleanDir(paste0(ROOT_PATH,"/www/tmp/"))

loadDemoMO <- function(params){
  gene1 <- loadGene(name="MG 1", type="majorGene")
  gene2 <- loadGene(name="MG 2", type="majorGene")
    
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type="growingHost")
  cultivar2 <- loadCultivar(name="Resistant1", type="growingHost")
  cultivar3 <- loadCultivar(name="Resistant2", type="growingHost")
  cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant1", c("MG 1"))
  params <- allocateCultivarGenes(params, "Resistant2", c("MG 2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible crop", "Resistant crop 1", "Resistant crop 2"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

loadDemoMI <- function(params){
  gene1 <- loadGene(name="MG 1", type="majorGene")
  gene2 <- loadGene(name="MG 2", type="majorGene")
  
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type="growingHost")
  cultivar2 <- loadCultivar(name="Resistant1", type="growingHost")
  cultivar3 <- loadCultivar(name="Resistant2", type="growingHost")
  cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant1", c("MG 1"))
  params <- allocateCultivarGenes(params, "Resistant2", c("MG 2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible crop", "Mixture"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Mixture", c("Resistant1","Resistant2"))
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

loadDemoRO <- function(params){
  gene1 <- loadGene(name="MG 1", type="majorGene")
  gene2 <- loadGene(name="MG 2", type="majorGene")
  
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type="growingHost")
  cultivar2 <- loadCultivar(name="Resistant1", type="growingHost")
  cultivar3 <- loadCultivar(name="Resistant2", type="growingHost")
  cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant1", c("MG 1"))
  params <- allocateCultivarGenes(params, "Resistant2", c("MG 2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible crop", "Resistant crop 1", "Resistant crop 2"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 1", "Resistant1")
  croptypes <- allocateCroptypeCultivars(croptypes, "Resistant crop 2", "Resistant2")
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

loadDemoPY <- function(params){
  gene1 <- loadGene(name="MG 1", type="majorGene")
  gene2 <- loadGene(name="MG 2", type="majorGene")
  gene1$mutation_prob <- 1E-4
  gene2$mutation_prob <- 1E-4
  
  genes <- data.frame(rbind(gene1, gene2), stringsAsFactors = FALSE)
  params <- setGenes(params, genes)
  
  cultivar1 <- loadCultivar(name="Susceptible", type="growingHost")
  cultivar2 <- loadCultivar(name="Resistant", type="growingHost")
  cultivars <- data.frame(rbind(cultivar1, cultivar2), stringsAsFactors = FALSE)
  
  params <- setCultivars(params, cultivars)
  
  params <- allocateCultivarGenes(params, "Resistant", c("MG 1", "MG 2"))
  
  croptypes <- loadCroptypes(params, names=c("Susceptible crop", "Pyramid"))
  croptypes <- allocateCroptypeCultivars(croptypes, "Susceptible crop", "Susceptible")
  croptypes <- allocateCroptypeCultivars(croptypes, "Pyramid", "Resistant")
  params <- setCroptypes(params, croptypes)
  
  return(params)
}

# Create a numeric input made for positive integer value
IntegerInput <- function(inputId, label, value, max) {
  shiny::numericInput(
    inputId = inputId,
    label = label,
    value = value,
    min = 0,
    max = max,
    step = 1
  )
}

# Create a numeric input made for percentage value
PercentageInput <- function(inputId, label, value) {
  shiny::numericInput(
    inputId = inputId,
    label = label,
    value = value,
    min = 0.0,
    max = 1.0,
    step = 0.02
  )
}

# Take a table of 3 croptype with 3 cultivars and render it
RenderCroptypes <- function(dt) {
  DT::renderDT(
    dt,
    editable = list(target = "row", disable = list(columns = c(0, 10))),
    rownames = FALSE,
    options = list(
      paging = FALSE,
      searching = FALSE,
      bInfo = FALSE,
      ordering = FALSE,
      columnDefs = list(list(
        className = "dt-center", targets = 0:(length(colnames(dt))-1)
      ))
    ),
    selection = "none",
    colnames = c(colnames(dt)),
    class = "cell-border stripe"
  )
}

# Take a table of 3 cultivars with 8 genes and render it
RenderCultivars <- function(dt) {
  DT::renderDT(
    dt,
    options = list(
      paging = FALSE,
      searching = FALSE,
      bInfo = FALSE,
      ordering = FALSE,
      select = list(info = FALSE),
      columnDefs = list(list(
        className = "dt-center gene", targets = 0:length(colnames(dt))
      ))
    ),
    rownames = TRUE,
    selection = "none",
    colnames = c(colnames(dt)),
    class = "cell-border stripe"
  )
}
