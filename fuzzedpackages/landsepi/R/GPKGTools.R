# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
#                    Julien Papaix <julien.papaix@inrae.fr>
#                    Jean-François Rey <jean-francois.rey@inrae.fr>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,i
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

### GPKG manipulation ###

#' @title createLandscapeGPKG
#' @description Creates a GPKG file from a sf object
#' @details Generates a GPKG file with two layouts.
#' One layout called "croptypeID" containing the croptypes ID for every year.
#' One another layout called "area" containing the area of every polygon.
#' @param landscape a sf object
#' @param outputfile GPKG output file name
#' @return the outputfile name (full path)
createLandscapeGPKG <- function(landscape, outputfile) {
  # create and save landscape in layer croptypeID (year_XX croptype ID)
  st_write(
    landscape,
    outputfile,
    "croptypeID",
    layer_options = "OVERWRITE=yes",
    driver = "GPKG",
    quiet = TRUE
  )

  # Get polygons area
  area <- data.frame(area = st_area(landscape))
  # save polygons area in a layer area
  st_write(
    st_sf(st_geometry(landscape), "area" = area),
    outputfile,
    "area",
    layer_options = "OVERWRITE=yes",
    driver = "GPKG",
    quiet = TRUE
  )

  return(outputfile)
}

#' getGPKGArea
#' @description Gets "area" layer as a vector of a GPKG file
#' @param gpkgfile a GPKG file
#' @return a vector of the area of each polygons
#' @export
getGPKGArea <- function(gpkgfile) {
  area <- st_read(
    gpkgfile,
    "area"
  )
  st_geometry(area) <- NULL
  areadf <- as.data.frame(area)

  return(as.vector(areadf[, 1]))
}


#' getGPKGRotation
#' @description Gets Croptypes ID rotation years as a matrix
#' @param gpkgfile a GPKG file
#' @return a matrix as rows for polygons and cols for years
#' @export
getGPKGRotation <- function(gpkgfile) {
  croptypeID <- st_read(
    gpkgfile,
    "CroptypeID"
  )
  st_geometry(croptypeID) <- NULL
  cdf <- as.data.frame(croptypeID)

  ncol <- length(grep("^year_", colnames(cdf)) %in% colnames(cdf))
  # nrow <- nrow(cdf)

  croptype_matrix <- as.matrix(cdf[, grep("^year_", colnames(cdf))], ncol = ncol)

  return(croptype_matrix)
}


#' @title GPKGAddTables
#' @description Adds non spatial data table definitions (sqlite) into a GPKG file.
#' It adds Cultivar, CultivarList, Gene, GeneList tables
#' @param gpkgfile a GPKF file
#' @return the GPKG file name
GPKGAddTables <- function(gpkgfile) {
  outputDB <-
    DBI::dbConnect(RSQLite::SQLite(), gpkgfile)

  res <- DBI::dbGetQuery(outputDB, "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='Cultivar'")
  if (res[1, 1] != 0) {
    warning(paste0(gpkgfile, " Data tables already built !"))
    return(1)
  }
  ## Define extension for non spatial data for gdal API < 2.1
  DBI::dbExecute(outputDB, "INSERT INTO gpkg_extensions (table_name, column_name, extension_name, definition, scope) VALUES (NULL, NULL, 'gdal_aspatial','http://gdal.org/geopackage_aspatial.html','read-write');")

  DBI::dbExecute(outputDB, "INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('Cultivar','aspatial','Cultivar');")
  DBI::dbExecute(outputDB, "INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('CultivarList','aspatial','CultivarList');")
  DBI::dbExecute(outputDB, "INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('Gene','aspatial','Gene');")
  DBI::dbExecute(outputDB, "INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('GeneList','aspatial','GeneList');")

  # Create a table for Cultivar
  DBI::dbExecute(
    outputDB,
    "CREATE TABLE Cultivar( cultivarID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                                                     cultivarName TEXT NOT NULL,
                                                     initial_density REAL NOT NULL CHECK(initial_density >= 0),
                                                     max_density REAL NOT NULL CHECK(max_density > 0),
                                                     growth_rate REAL NOT NULL CHECK(growth_rate BETWEEN 0 AND 1),
                                                     reproduction_rate REAL NOT NULL CHECK(reproduction_rate BETWEEN 0 AND 1),
                                                     death_rate REAL NOT NULL CHECK(death_rate BETWEEN 0 AND 1),
                                                     yield_H REAL NOT NULL CHECK(yield_H >= 0),
                                                     yield_L REAL NOT NULL CHECK(yield_L >= 0),
                                                     yield_I REAL NOT NULL CHECK(yield_I >= 0),
                                                     yield_R REAL NOT NULL CHECK(yield_R >= 0),
                                                     production_cost REAL NOT NULL CHECK(production_cost >= 0),
                                                     market_value REAL NOT NULL CHECK(market_value >= 0));"
  )

  # Create a table for CultivarList
  DBI::dbExecute(
    outputDB,
    "CREATE TABLE CultivarList(rowid INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
    							croptypeID INTEGER NOT NULL,
                                                          cultivarID INTEGER NOT NULL,
                                                          proportion REAL NOT NULL CHECK(proportion BETWEEN 0 AND 1),
                                                          FOREIGN KEY(cultivarID) REFERENCES Cultivar (cultivarID));"
    # PRIMARY KEY(croptypeID, cultivarID),
  )

  # Create a table for Gene
  DBI::dbExecute(
    outputDB,
    "CREATE TABLE Gene(geneID INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                                                  geneName TEXT NOT NULL,
                                                  efficiency REAL NOT NULL CHECK(efficiency BETWEEN 0 AND 1),
                                                  time_to_activ_exp REAL NOT NULL CHECK(time_to_activ_exp >= 0),
                                                  time_to_activ_var REAL NOT NULL CHECK(time_to_activ_var >= 0),
                                                  mutation_prob REAL NOT NULL CHECK(mutation_prob BETWEEN 0 AND 1),
                                                  Nlevels_aggressiveness INTEGER NOT NULL CHECK(Nlevels_aggressiveness >= 1),
                                                  fitness_cost REAL NOT NULL CHECK(fitness_cost BETWEEN 0 AND 1),
                                                  tradeoff_strength REAL NOT NULL CHECK(tradeoff_strength > 0),
                                                  target_trait TEXT NOT NULL);"
  )

  # Create a table for GeneList
  DBI::dbExecute(
    outputDB,
    "CREATE TABLE GeneList(rowid INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                                                      cultivarID INTEGER NOT NULL,
                                                      geneID INTEGER NOT NULL,
                                                      FOREIGN KEY(cultivarID) REFERENCES Cultivar (cultivarID),
                                                      FOREIGN KEY(geneID) REFERENCES Gene (geneID));"
  )

  # Close the output database
  DBI::dbDisconnect(outputDB)
  
  return(gpkgfile)
}

#' @title GPKGAddInputData
#' @details Adds data 'data' values in the table 'table' using 'data' colnames
#' @param gpkgfile a gpkg filename
#' @param table table name
#' @param data a data to write in BDD, should be the return of a function param2XXXXXBDD
#' @param deleteExistingData if TRUE overwrite data if already present in gpkg file, default FALSE
GPKGAddInputData <- function(gpkgfile, table = "", data = data.frame(), deleteExistingData = FALSE) {
  outputDB <- DBI::dbConnect(RSQLite::SQLite(), gpkgfile)

  # check table exist
  res <- DBI::dbGetQuery(outputDB, paste0("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='", table, "'"))
  if (res[1, 1] != 1) {
    warning(paste0(table, " table does not exist in ", gpkgfile))
    return(1)
  }

  # check table is empty
  res <- DBI::dbGetQuery(outputDB, paste0("SELECT count(*) FROM ", table))
  if (res[1, 1] != 0) {
    warning(paste0(table, " table is not empty in ", gpkgfile))
    if (deleteExistingData == TRUE) {
      warning(paste0("Deleting existing data in ", table, " table"))
      DBI::dbExecute(outputDB, paste0("DELETE FROM ", table))
    }
    else {
      warning(paste0("Will try to add data to table ", table))
      warning("Use deleteExistingData = TRUE parameter to allow data overwritting")
    }
  }

  req <- paste0("INSERT INTO ", table, " VALUES('")
  for (i in 1:nrow(data)) {
    reqsql <- paste0(req, paste(data[i, ], collapse = "','"), "')")
    #print(reqsql)
    DBI::dbExecute(outputDB, reqsql)
  }

  # Close the output database
  DBI::dbDisconnect(outputDB)
}

#' getGPKGCroptypes
#' @description Gets Croptypes and Cultivars proportions from a GPKG file
#' @param inputGPKG a GPKG filename
#' @return a data.frame with croptypeID, CultivarNames, and proportions
getGPKGCroptypes <- function(inputGPKG) {
  inputDB <- DBI::dbConnect(RSQLite::SQLite(), inputGPKG)

  # Select croptypes and cultivars name
  sql <- "SELECT CultivarList.croptypeID, Cultivar.cultivarName, CultivarList.proportion FROM CultivarList, Cultivar WHERE CultivarList.cultivarID = Cultivar.cultivarID;"
  croptypes <- DBI::dbGetQuery(inputDB, sql)

  dfcroptypes <- data.frame(matrix(data = 0.0, nrow = length(unique(croptypes$croptypeID)), ncol = length(unique(croptypes$cultivarName)) + 1))
  colnames(dfcroptypes) <- c("croptypeID", unique(croptypes$cultivarName))
  dfcroptypes$croptypeID <- unique(croptypes$croptypeID)
  rownames(dfcroptypes) <- dfcroptypes$croptypeID
  tmp <- apply(croptypes,
    MARGIN = 1,
    FUN = function(l) {
      dfcroptypes[which(dfcroptypes$croptypeID == l[1]), which(colnames(dfcroptypes) == l[2])] <<- l[3]
    }
  )

  DBI::dbDisconnect(inputDB)

  return(dfcroptypes)
}

#' getGPKGCroptypesRaw
#' @description Gets Croptypes and Cultivars proportions from a GPKG file without refactoring
#' @param inputGPKG a GPKG filename
#' @return a data.frame with croptypeID, CultivarID, and proportions
getGPKGCroptypesRaw <- function(inputGPKG) {
  inputDB <- DBI::dbConnect(RSQLite::SQLite(), inputGPKG)

  # Select croptypes and cultivars name
  sql <- "SELECT croptypeID, cultivarID, proportion FROM CultivarList;"
  croptypes <- DBI::dbGetQuery(inputDB, sql)

  DBI::dbDisconnect(inputDB)

  return(croptypes)
}

#' getGPKGCultivars
#' @description Gets Cultivars from a GPKG file
#' @param inputGPKG a GPKG filename
#' @return a data.frame
getGPKGCultivars <- function(inputGPKG) {
  inputDB <- DBI::dbConnect(RSQLite::SQLite(), inputGPKG)

  # Select croptypes and cultivars name
  sql <- "SELECT * FROM Cultivar;"
  dfcultivars <- DBI::dbGetQuery(inputDB, sql)
  DBI::dbDisconnect(inputDB)

  rownames(dfcultivars) <- dfcultivars$cultivarName

  return(dfcultivars[, -1])
}

#' getGPKGGenes
#' @description Gets Genes from a GPKG file
#' @param inputGPKG a GPKG filename
#' @return a data.frame
getGPKGGenes <- function(inputGPKG) {
  inputDB <- DBI::dbConnect(RSQLite::SQLite(), inputGPKG)

  # Select croptypes and cultivars name
  sql <- "SELECT * FROM Gene;"
  dfgenes <- DBI::dbGetQuery(inputDB, sql)
  DBI::dbDisconnect(inputDB)

  rownames(dfgenes) <- dfgenes$geneName

  return(dfgenes[, -1])
}



#' getGPKGCultivarsGenes
#' @description Gets Cultivars Genes from a GPKG file
#' @param inputGPKG a GPKG filename
#' @return a data.frame
getGPKGCultivarsGenes <- function(inputGPKG) {
  inputDB <- DBI::dbConnect(RSQLite::SQLite(), inputGPKG)

  # Select croptypes and cultivars name

  sql <- "SELECT Cultivar.cultivarName, Gene.geneName FROM GeneList, Cultivar, Gene WHERE GeneList.geneID = Gene.geneID AND GeneList.cultivarID == Cultivar.cultivarID;"
  genes <- DBI::dbGetQuery(inputDB, sql)
  DBI::dbDisconnect(inputDB)

  dfCultivarGene <- data.frame(matrix(data = 0, nrow = length(unique(genes$cultivarName)), ncol = length(unique(genes$geneName))))
  colnames(dfCultivarGene) <- unique(genes$geneName)
  rownames(dfCultivarGene) <- unique(genes$cultivarName)
  tmp <- apply(genes,
    MARGIN = 1,
    FUN = function(l) {
      dfCultivarGene[which(rownames(dfCultivarGene) == l[1]), which(colnames(dfCultivarGene) == l[2])] <<- 1
    }
  )

  return(dfCultivarGene)
}

#' getGPKGGeneIDForCultivar
#' @description Gets GeneID of a cultivar
#' @param inputGPKG a GPKG filename
#' @param cultivarID a cultivarID
#' @return a data.frame of GeneID
getGPKGGeneIDForCultivar <- function(inputGPKG, cultivarID) {
  inputDB <- DBI::dbConnect(RSQLite::SQLite(), inputGPKG)

  geneID <- DBI::dbGetQuery(inputDB, paste("SELECT GeneID FROM GeneList WHERE cultivarID=", cultivarID, sep = ""))
  DBI::dbDisconnect(inputDB)

  return(geneID)
}
