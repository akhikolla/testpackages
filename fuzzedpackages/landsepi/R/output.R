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
#
#' @title Generation of epidemiological and economic model outputs
#' @name epid_output
#' @description Generates epidemiological and economic outputs from model simulations.
#' @param types a character string (or a vector of character strings if several outputs are to be computed)
#' specifying the type of outputs to generate (see details):\itemize{
#' \item "audpc": Area Under Disease Progress Curve
#' \item "gla_abs": Absolute Green Leaf Area
#' \item "gla_rel": Relative Green Leaf Area
#' \item "eco_product": Crop production
#' \item "eco_cost": Crop costs
#' \item "eco_benefit": Crop benefits
#' \item "eco_grossmargin": Gross Margin (benefits - costs)
#' \item "HLIR_dynamics", "H_dynamics", "L_dynamics", "IR_dynamics", "HLI_dynamics", etc.: Epidemic dynamics
#' related to the specified sanitary status (H, L, I or R and all their combinations). 
#' Graphics only, works only if graphic=TRUE.
#' \item "all": compute all these outputs (default).
#' }
#' @param time_param list of simulation parameters:\itemize{
#' \item Nyears = number cropping seasons,
#' \item nTSpY = number of time-steps per cropping season.
#' }
#' @param Npatho number of pathogen genotypes.
#' @param area a vector containing polygon areas (must be in square meters).
#' @param rotation a dataframe containing for each field (rows) and year (columns, named "year_1", "year_2", etc.),
#' the index of the cultivated croptype. Importantly, the matrix must contain 1 more column than the real number
#' of simulated years.
#' @param croptypes a dataframe with three columns named 'croptypeID' for croptype index,
#' 'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar within the croptype.
#' @param cultivars_param list of parameters associated with each host genotype (i.e. cultivars):  
#' \itemize{
#' \item name = vector of cultivar names,
#' \item initial_density = vector of host densities (per square meter) at the beginning of the cropping season,
#' \item max_density = vector of maximum host densities (per square meter) at the end of the cropping season,
#' \item cultivars_genes_list = a list containing, for each host genotype, the indices of carried resistance genes.
#' }
#' @param eco_param a list of economic parameters for each host genotype as if cultivated in pure crop:\itemize{
#' \item yield_perHa = a dataframe of 4 columns for the yield associated with hosts in sanitary status H, L, I and R,
#' and one row per host genotype (yields are expressed in weight or volume units / ha / cropping season),
#' \item production_cost_perHa = a vector of overall production costs (in monetary units / ha / cropping season)
#' including planting costs, amortisation, labour etc.,
#' \item market_value = a vector of market values of the productions (in monetary units / weight or volume unit).
#' }
#' @param GLAnoDis the value of absolute GLA in absence of disease (required to compute economic outputs).
#' @param ylim_param a list of graphical parameters for each required output: bounds for y-axes for
#' audpc, gla_abs, gla_rel, eco_cost, eco_product, eco_benefit, eco_grossmargin.
#' @param writeTXT a logical indicating if the output is written in a text file (TRUE) or not (FALSE).
#' @param graphic a logical indicating if a tiff graphic of the output is generated (only if more than one year is simulated).
#' @param path path of text file (if writeTXT = TRUE) and tiff graphic (if graphic = TRUE) to be generated.
#' @details Outputs are computed every year for every cultivar as well as for the whole landscape. \describe{
#' \item{\strong{Epidemiological outputs.}}{
#' The epidemiological impact of pathogen spread can be evaluated by different measures: \enumerate{
#' \item Area Under Disease Progress Curve (AUDPC): average proportion of diseased hosts (status I + R)
#' relative to the carrying capacity (i.e. disease severity).
#' \item Absolute Green Leaf Area (GLAa): average number of healthy hosts (status H) per time step and per square meter.
#' \item Relative Green Leaf Area (GLAr): average proportion of healthy hosts (status H) relative to the total number
#' of existing hosts (H+L+I+R).
#' }
#'  }
#' \item{\strong{Economic outputs.}}{
#' The economic outcome of a simulation can be evaluated using: \enumerate{
#' \item Crop production: yearly crop production (e.g. grains, fruits, wine) in weight (or volume) units
#' per hectare (depends on the number of productive hosts and associated yield).
#' \item Crop benefits: yearly benefits generated from product sales, in monetary units per hectare
#' (depends on crop production and market value of the product).
#' \item Crop costs: yearly costs associated with crop production (including planting, amortisation, labour, ...)
#'  in monetary units per hectare (depends on initial host density and production cost).
#' \item Crop gross margin, i.e. benefits - costs, in monetary units per hectare.
#' }
#'  }
#'  }
#' @return A list containing, for each required type of output, a matrix summarising the output for each year and cultivar
#' (as well as the whole landscape).
#' Each matrix can be written in a txt file (if writeTXT=TRUE), and illustrated in a graphic (if graphic=TRUE).
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018). Assessing the durability and efficiency of
#' landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
#' @seealso \link{evol_output}
#' @examples
#' \dontrun{
#' demo_landsepi()
#' }
#' @include RcppExports.R Math-Functions.R graphics.R
#' @importFrom sf st_read
#' @importFrom grDevices colorRampPalette dev.off graphics.off gray png tiff
#' @importFrom utils write.table
#' @export
epid_output <- function(types = "all", time_param, Npatho, area, rotation, croptypes, cultivars_param, eco_param,
                        GLAnoDis = cultivars_param$max_density[1] ## true for non-growing host
                        , ylim_param = list(
                          audpc = c(0, 0.38),
                          gla_abs = c(0, 1.48),
                          gla_rel = c(0, 1),
                          eco_cost = c(0, NA),
                          eco_product = c(0, NA),
                          eco_benefit = c(0, NA),
                          eco_grossmargin = c(NA, NA)
                        ),
                        writeTXT = TRUE, graphic = TRUE, path = getwd()) {
  valid_outputs <- c("audpc", "gla_abs", "gla_rel", "eco_product", "eco_benefit", "eco_cost", "eco_grossmargin", "HLIR_dynamics")
  if (types[1] == "all") {
    types <- valid_outputs
  }
  valid_outputs <- c(valid_outputs, paste0(c("H", "L", "I", "R", "HL", "HI", "HR", "LI", "LR", "IR", "HLI", "HLR", "HIR", "LIR", "HLIR"), "_dynamics"))
  if (is.na(sum(match(types, valid_outputs)))) {
    stop(paste("Error: valid types of outputs are", paste(valid_outputs, collapse = ", ")))
  }

  ## Time parameters
  Nyears <- time_param$Nyears
  nTSpY <- time_param$nTSpY
  nTS <- Nyears * nTSpY ## Total number of time-steps

  ## Landscape
  rotation <- data.frame(rotation)
  areaTot <- sum(area)
  Npoly <- length(area)

  ## Host parameters
  cultivar_names <- cultivars_param$name
  initial_density <- cultivars_param$initial_density
  max_density <- cultivars_param$max_density
  cultivars_genes_list <- cultivars_param$cultivars_genes_list
  yield_perIndivPerTS <- eco_param$yield_perHa * 1E-4 / nTSpY / GLAnoDis
  production_cost_perIndiv <- eco_param$production_cost_perHa * 1E-4 / initial_density
  market_value <- eco_param$market_value
  Nhost <- length(initial_density)


  ## Calculation of the carrying capacity and number of exisiting hosts
  K <- array(dim = c(Npoly, Nhost, Nyears)) ## for audpc
  C <- array(dim = c(Npoly, Nhost, Nyears)) ## for cost
  area_host <- matrix(0, nrow = Nyears, ncol = Nhost + 1) ## for cost
  for (y in 1:Nyears) {
    if (ncol(rotation) == 1) {
      index_year <- 1
    } else {
      index_year <- y
    }
    ts_year <- ((y - 1) * nTSpY + 1):(y * nTSpY)

    for (poly in 1:Npoly) {
      indices_croptype <- grep(rotation[poly, index_year], croptypes$croptypeID)
      for (i in indices_croptype) {
        index_host <- croptypes[i, "cultivarID"] + 1 ## +1 because C indices start at 0
        prop <- croptypes[i, "proportion"]
        K[poly, index_host, y] <- floor(area[poly] * max_density[index_host] * prop)
        C[poly, index_host, y] <- area[poly] * initial_density[index_host] * prop
        area_host[y, index_host] <- area_host[y, index_host] + area[poly]
      }
    } ## for poly
  } ## for y
  K_host <- apply(K, c(2, 3), sum, na.rm = TRUE)
  C_host <- apply(C, c(2, 3), sum, na.rm = TRUE)
  C_host[C_host == 0] <- NA
  area_host[, Nhost + 1] <- areaTot

  ## IMPORTATION OF THE SIMULATION OUTPUT (only those required depending on parameter "types")
  requireH <- sum(is.element(substr(types, 1, 3), c("gla", "eco"))) > 0
  requireL <- sum(is.element(types, c("gla_rel", "eco_product", "eco_benefit", "eco_cost", "eco_grossmargin"))) > 0
  requireIR <- sum(is.element(types, c("audpc", "gla_rel", "eco_product", "eco_benefit", "eco_cost", "eco_grossmargin"))) > 0
  if (graphic & sum(substr(types, nchar(types) - 7, nchar(types)) == "dynamics") > 0) {
    requireH <- 1
    requireL <- 1
    requireIR <- 1
  }
  H <- as.list(1:nTS)
  # Hjuv <- as.list(1:nTS)
  # P <- as.list(1:nTS)
  L <- as.list(1:nTS)
  I <- as.list(1:nTS)
  R <- as.list(1:nTS)
  index <- 0

  # print("Reading binary files to compute epidemiological model outputs")
  for (year in 1:Nyears) {
    # print(paste("year", year, "/", Nyears))
    if (requireH) {
      binfileH <- file(paste(path, sprintf("/H-%02d", year), ".bin", sep = ""), "rb")
      H.tmp <- readBin(con = binfileH, what = "int", n = Npoly * Nhost * nTSpY, size = 4, signed = T, endian = "little")
      close(binfileH)
    }
    # binfileHjuv = file(paste(path, sprintf("/Hjuv-%02d", year), ".bin",sep=""), "rb")
    #  Hjuv.tmp <- readBin(con=binfileHjuv, what="int", n=Npoly*Nhost*nTSpY, size = 4, signed=T,endian="little")
    # close(binfileHjuv)
    # binfileP <- file(paste(path, sprintf("/P-%02d", year), ".bin", sep = ""), "rb")
    # P.tmp <- readBin(con = binfileP,
    #                  what = "int",
    #                  n = Npoly * Npatho * nTSpY,
    #                  size = 4,
    #                  signed = T,
    #                  endian = "little")
    # close(binfileP)
    if (requireL) {
      binfileL <- file(paste(path, sprintf("/L-%02d", year), ".bin", sep = ""), "rb")
      L.tmp <- readBin(con = binfileL, what = "int", n = Npoly * Npatho * Nhost * nTSpY, size = 4, signed = T, endian = "little")
      close(binfileL)
    }
    if (requireIR) {
      binfileI <- file(paste(path, sprintf("/I-%02d", year), ".bin", sep = ""), "rb")
      I.tmp <- readBin(con = binfileI, what = "int", n = Npoly * Npatho * Nhost * nTSpY, size = 4, signed = T, endian = "little")
      close(binfileI)

      binfileR <- file(paste(path, sprintf("/R-%02d", year), ".bin", sep = ""), "rb")
      R.tmp <- readBin(con = binfileR, what = "int", n = Npoly * Npatho * Nhost * nTSpY, size = 4, signed = T, endian = "little")
      close(binfileR)
    }

    for (t in 1:nTSpY) {
      if (requireH) {
        H[[t + index]] <- matrix(H.tmp[((Nhost * Npoly) * (t - 1) + 1):(t * (Nhost * Npoly))], ncol = Nhost, byrow = T)
      }
      # Hjuv[[t + index]] <- matrix(Hjuv.tmp[((Nhost*Npoly)*(t-1)+1):(t*(Nhost*Npoly))],ncol=Nhost,byrow=T)
      # P[[t + index]] <- matrix(P.tmp[((Npatho * Npoly) * (t - 1) + 1):(t * (Npatho * Npoly))], ncol = Npatho, byrow = T)
      if (requireL) {
        L[[t + index]] <- array(
          data = L.tmp[((Npatho * Npoly * Nhost) * (t - 1) + 1):(t * (Npatho * Npoly * Nhost))],
          dim = c(Nhost, Npatho, Npoly)
        )
      }
      if (requireIR) {
        I[[t + index]] <- array(
          data = I.tmp[((Npatho * Npoly * Nhost) * (t - 1) + 1):(t * (Npatho * Npoly * Nhost))],
          dim = c(Nhost, Npatho, Npoly)
        )
        R[[t + index]] <- array(
          data = R.tmp[((Npatho * Npoly * Nhost) * (t - 1) + 1):(t * (Npatho * Npoly * Nhost))],
          dim = c(Nhost, Npatho, Npoly)
        )
      }
    } ## for t

    index <- index + nTSpY
  } ## for year

  #### HLIR DYNAMIC
  H_host <- NULL
  L_host <- NULL
  I_host <- NULL
  R_host <- NULL
  if (requireH) {
    for (t in 1:nTS) {
      H_host <- cbind(H_host, apply(H[[t]], 2, sum))
    }
  }
  if (requireL) {
    for (t in 1:nTS) {
      L_host <- cbind(L_host, apply(L[[t]], 1, sum))
    }
  }
  if (requireIR) {
    for (t in 1:nTS) {
      I_host <- cbind(I_host, apply(I[[t]], 1, sum))
      R_host <- cbind(R_host, apply(R[[t]], 1, sum))
    }
  }
  N_host <- H_host + L_host + I_host + R_host


  #####       COMPUTATION OF OUTPUTS
  ## as.numeric to allow memory allocation of long integer
  res <- list()

  ## Colours (+1 to avoid picking white)
  # COL <- c("#FF5555","#4F94CD","darkolivegreen4","#CD950C","black") ## colors: red, blue, green, orange, black
  grad_red <- colorRampPalette(c("#FF5555", "white"))
  RED <- grad_red(Nhost + 1)
  grad_blue <- colorRampPalette(c("#4F94CD", "white"))
  BLUE <- grad_blue(Nhost + 1)
  grad_green <- colorRampPalette(c("darkolivegreen4", "white"))
  GREEN <- grad_green(Nhost + 1)
  grad_grey <- colorRampPalette(c("black", "white")) # "gray30"
  GRAY <- grad_grey(Nhost + 1)

  for (type in types) {

    ## Epidemic dynamics (H, L, I, R)
    if (graphic & substr(type, nchar(type) - 7, nchar(type)) == "dynamics") {
      varToPlot <- strsplit(type, "_")[[1]][1]
      varToLegend <- strsplit(varToPlot, "")[[1]]

      tiff(
        filename = paste(path, "/", type, ".tiff", sep = ""),
        width = 180, height = 110, units = "mm", compression = "lzw", res = 300
      )
      par(xpd = NA, cex = .9, mar = c(4, 4.3, 3, 9))

      plot(0, 0,
        type = "n", xaxt = "n", yaxt = "n", bty = "n", xlim = c(1, nTS), ylim = c(0, 1), main = "Epidemic dynamics",
        ylab = "Proportion of hosts relative to carrying capacity", xlab = ""
      )
      for (host in 1:Nhost) {
        if (grepl("H", varToPlot)) {
          lines(1:nTS, H_host[host, ] / rep(K_host[host, ], each = nTSpY), col = GREEN[host], lty = host)
        }
        if (grepl("L", varToPlot)) {
          lines(1:nTS, L_host[host, ] / rep(K_host[host, ], each = nTSpY), col = BLUE[host], lty = host)
        }
        if (grepl("I", varToPlot)) {
          lines(1:nTS, I_host[host, ] / rep(K_host[host, ], each = nTSpY), col = RED[host], lty = host)
        }
        if (grepl("R", varToPlot)) {
          lines(1:nTS, R_host[host, ] / rep(K_host[host, ], each = nTSpY), col = GRAY[host], lty = host)
        }
      }
      if (Nyears == 1) {
        axis(1, at = round(seq(1, nTS, length.out = 8)))
        title(xlab = "Days")
      } else {
        axis(1, at = seq(1, nTS + 1, nTSpY * ((Nyears - 1) %/% 10 + 1)), labels = seq(0, Nyears, ((Nyears - 1) %/% 10 + 1)))
        title(xlab = "Years")
      }
      axis(2, at = seq(0, 1, .2), las = 1)
      legend(nTS * 1.05, .8,
        title.adj = -.01, bty = "n", title = "Status:", legend = varToLegend, border = NA,
        fill = c(GREEN[1], BLUE[1], RED[1], GRAY[1])[is.element(c("H", "L", "I", "R"), varToLegend)], cex = 0.9
      )
      legend(nTS * 1.05, .3,
        title.adj = -.01, bty = "n", title = "Cultivar:", legend = cultivar_names,
        lty = 1:Nhost, lwd = 2, col = GRAY[1:Nhost], cex = 0.9
      )
      dev.off()
    } else { ## Other epidemiological outputs

      output_matrix <- data.frame(matrix(NA, ncol = Nhost + 1, nrow = Nyears))
      rownames(output_matrix) <- paste("year_", 1:Nyears, sep = "")
      colnames(output_matrix) <- c(cultivar_names, "total")

      for (y in 1:Nyears) {
        ts_year <- ((y - 1) * nTSpY + 1):(y * nTSpY)
        ## metrics for each host
        for (host in 1:Nhost) {
          switch(type,
            "audpc" = {
              numerator <- sum(as.numeric(I_host[host, ts_year]), as.numeric(R_host[host, ts_year]))
              denominator <- sum(as.numeric(K_host[host, y] * length(ts_year)))
            },
            "gla_abs" = {
              numerator <- sum(as.numeric(H_host[host, ts_year]))
              denominator <- nTSpY * areaTot
            },
            "gla_rel" = {
              numerator <- sum(as.numeric(H_host[host, ts_year]))
              denominator <- sum(as.numeric(N_host[host, ts_year]))
            },
            { ## Default instructions (i.e. if substr(type,1,3)=="eco")
              numerator <- sum(
                yield_perIndivPerTS[host, "H"] * as.numeric(H_host[host, ts_year]),
                yield_perIndivPerTS[host, "L"] * as.numeric(L_host[host, ts_year]),
                yield_perIndivPerTS[host, "I"] * as.numeric(I_host[host, ts_year]),
                yield_perIndivPerTS[host, "R"] * as.numeric(R_host[host, ts_year])
              )
              denominator <- 1
            }
          )

          if (denominator > 0 & K_host[host, y] > 0) {
            output_matrix[y, host] <- numerator / denominator
          } ## if >0
        } ## for host

        ## metrics for all landscape
        switch(type,
          "audpc" = {
            numerator_tot <- sum(as.numeric(I_host[, ts_year]), as.numeric(R_host[, ts_year]))
            denominator_tot <- sum(as.numeric(K_host[, y] * length(ts_year)))
            if (denominator_tot > 0) {
              output_matrix[y, "total"] <- numerator_tot / denominator_tot
            }
          },
          "gla_rel" = {
            numerator_tot <- sum(as.numeric(H_host[, ts_year]))
            denominator_tot <- sum(as.numeric(N_host[, ts_year]))
            if (denominator_tot > 0) {
              output_matrix[y, "total"] <- numerator_tot / denominator_tot
            }
          },
          { ## i.e. if gla_abs | substr(type,1,3)=="eco"
            output_matrix[y, "total"] <- sum(output_matrix[y, 1:Nhost], na.rm = TRUE)
          }
        )
      } ## for y

      ## Economic outputs
      if (substr(type, 1, 3) == "eco") {
        eco <- list(eco_product = output_matrix)
        eco[["eco_benefit"]] <- data.frame(rep(market_value, each = Nyears) * eco[["eco_product"]][, 1:Nhost])
        ## data.frame to avoid pb when Nhost=1
        colnames(eco[["eco_benefit"]]) <- cultivar_names ## useful if Nhost=1
        eco[["eco_benefit"]]$total <- apply(eco[["eco_benefit"]], 1, sum, na.rm = TRUE)
        eco[["eco_cost"]] <- data.frame(rep(production_cost_perIndiv, each = Nyears) * t(C_host))
        colnames(eco[["eco_cost"]]) <- cultivar_names ## useful if Nhost=1
        eco[["eco_cost"]]$total <- apply(eco[["eco_cost"]], 1, sum, na.rm = TRUE)
        eco[["eco_grossmargin"]] <- eco[["eco_benefit"]] - eco[["eco_cost"]]

        output_matrix <- eco[[type]] / (area_host * 1E-4) ## normalization by area in Ha
      }
      # product <- t(output_matrix[,1:Nhost])
      # benefit <- market_value %*% product
      # cost <- production_cost_perIndiv %*% C_host
      # grossmargin <- benefit - cost

      if (writeTXT) {
        write.table(output_matrix, file = paste(path, "/", type, ".txt", sep = ""), sep = ",")
      }

      if (graphic & Nyears > 1) {
        ## Graphical parameters
        PCH <- 15:(15 + Nhost - 1)
        PCH.tot <- 0
        LTY.tot <- 1
        COL.tot <- RED[1]

        switch(type, "audpc" = {
          main_output <- "AUDPC"
          ylab_output <- "Proportion of diseased hosts: (I+R)/K"
          round_ylab_output <- 2
        }, "gla_abs" = {
          main_output <- expression("Absolute green leaf area (GLA"[a] * ")") ## expression('title'[2])
          ylab_output <- expression("Equivalent nb of healthy hosts / timestep / m"^2) ## expression('title'^2)
          round_ylab_output <- 2
        }, "gla_rel" = {
          main_output <- expression("Relative green leaf area (GLA"[r] * ")")
          ylab_output <- "Proportion of healthy hosts: H/N"
          round_ylab_output <- 2
        }, "eco_product" = {
          main_output <- "Crop production"
          ylab_output <- "Weight (or volume) units / ha"
          round_ylab_output <- 0
        }, "eco_benefit" = {
          main_output <- "Crop benefits"
          ylab_output <- "Monetary units / ha"
          round_ylab_output <- 0
        }, "eco_cost" = {
          main_output <- "Crop costs"
          ylab_output <- "Monetary units / ha"
          round_ylab_output <- 0
        }, "eco_grossmargin" = {
          main_output <- "Gross margin"
          ylab_output <- "Monetary units / ha"
          round_ylab_output <- 0
        })
        ymin_output <- ylim_param[[type]][1]
        ymax_output <- ylim_param[[type]][2]
        if (is.na(ymin_output)) {
          ymin_output <- floor(min(0, min(output_matrix, na.rm = TRUE)))
        }
        if (is.na(ymax_output)) {
          ymax_output <- ceiling(max(output_matrix, na.rm = TRUE))
        }
        output_mean_host <- apply(output_matrix, 2, mean, na.rm = TRUE) ## Average output per cultivar

        graphics.off()
        tiff(
          filename = paste(path, "/", type, ".tiff", sep = ""),
          width = 180, height = 110, units = "mm", compression = "lzw", res = 300
        )
        # m <- matrix(c(rep(1,5),2),6,1)    # matrix(c(rep(1,5),3,rep(2,5),3),6,2)
        # layout(m)
        # par(xpd=F, cex=.9,mar=c(5,4.1,4.1,2.1))
        par(xpd = NA, cex = .9, mar = c(4, 4.3, 3, 9))
        plot(0, 0,
          type = "n", bty = "n", xaxt = "n", yaxt = "n", xlim = c(1, Nyears), ylim = c(ymin_output, ymax_output),
          main = main_output, xlab = "Years", ylab = ylab_output
        )
        axis(1, at = seq(1, Nyears, by = (Nyears - 1) %/% 10 + 1))
        axis(2, at = round(seq(ymin_output, ymax_output, length.out = 5), round_ylab_output), las = 1)
        for (host in 1:Nhost) {
          lines(1:Nyears, output_matrix[, host], type = "o", lwd = 2, lty = host, pch = PCH[host], col = GRAY[host])
          points(par("usr")[1], output_mean_host[host], pch = PCH[host], col = GRAY[host], xpd = TRUE)
        }
        lines(1:Nyears, output_matrix[, "total"], type = "o", lwd = 2, lty = LTY.tot, pch = PCH.tot, col = COL.tot)
        points(par("usr")[1], output_mean_host["total"], pch = PCH.tot, col = COL.tot, cex = 1, xpd = TRUE)

        if (ymin_output < 0) {
          abline(h = 0, lwd = 1, lty = 3, col = BLUE[1], xpd = FALSE)
        }

        legend(Nyears * 1.05, .5 * ymax_output,
          title.adj = -.01, bty = "n", title = "Cultivar:",
          legend = c(cultivar_names, "Whole landscape"), cex = 0.9, lty = c(1:Nhost, LTY.tot),
          lwd = 2, pch = c(PCH, PCH.tot), pt.cex = c(rep(1, Nhost), 1), col = c(GRAY[1:Nhost], COL.tot), seg.len = 2.5
        )
        dev.off()
      }

      res[[type]] <- output_matrix
    } ## for type
  } ## else if output != HLIR_dynamics

  return(res)
}


#' @title Switch from index of genotype to indices of agressiveness on different components
#' @name switch_patho_to_aggr
#' @description Finds the level of aggressiveness on different components (targeted by different resistance genes)
#' from the index of a given pathogen genotype
#' @param index_patho index of pathogen genotype
#' @param Ngenes number of resistance genes
#' @param Nlevels_aggressiveness vector of the number of adaptation levels related to each resistance gene
#' @return a vector containing the indices of aggressiveness on the different components targeter by the resistance genes
#' @examples
#' switch_patho_to_aggr(5, 3, c(2, 2, 3))
#' @export
switch_patho_to_aggr <- function(index_patho, Ngenes, Nlevels_aggressiveness) {
  ## Index_patho starts at 0
  ## Index aggressiveness starts at 0
  aggr <- numeric()
  remainder <- index_patho
  for (g in 1:Ngenes) {
    prod <- 1
    k <- g
    while (k < Ngenes) {
      k <- k + 1
      prod <- prod * Nlevels_aggressiveness[k]
    }
    aggr[g] <- remainder %/% prod ## Quotient
    remainder <- remainder %% prod
  }
  return(aggr)
}



#' @title Generation of evolutionary model outputs
#' @name evol_output
#' @description Generates evolutionary outputs from model simulations.
#' @param types a character string (or a vector of character strings if several outputs are to be computed) specifying the
#' type of outputs to generate (see details):\itemize{
#' \item "evol_patho": Dynamics of pathogen genotype frequencies
#' \item "evol_aggr": Evolution of pathogen aggressiveness
#' \item "durability": Durability of resistance genes
#' \item "all": compute all these outputs (default)
#' }
#' @param time_param list of simulation parameters:\itemize{
#' \item Nyears = number cropping seasons,
#' \item nTSpY = number of time-steps per cropping season.
#' }
#' @param Npoly number of fields in the landscape.
#' @param cultivars_param list of parameters associated with each host genotype (i.e. cultivars)
#' when cultivated in pure crops:\itemize{
#' \item name = vector of cultivar names,
#' \item cultivars_genes_list = a list containing, for each host genotype, the indices of carried resistance genes.
#' }
#' @param genes_param list of parameters associated with each resistance gene and with the evolution of
#' each corresponding pathogenicity gene:\itemize{
#' \item name = vector of names of resistance genes,
#' \item Nlevels_aggressiveness = vector containing the number of adaptation levels related to each resistance gene (i.e. 1 + number
#' of required mutations for a pathogenicity gene to fully adapt to the corresponding resistance gene),
#' }
#' @param thres_breakdown an integer (or vector of integers) giving the threshold (i.e. number of infections) above which a
#' pathogen genotype is unlikely to go extinct and resistance is considered broken down, used to characterise the time to
#' invasion of resistant hosts (several values are computed if several thresholds are given in a vector).
#' @param writeTXT a logical indicating if the output is written in a text file (TRUE) or not (FALSE).
#' @param graphic a logical indicating if graphics must be generated (TRUE) or not (FALSE).
#' @param path a character string indicating the path of the repository where simulation output files are located and
#' where .txt files and graphics will be generated.
#'
#' @details For each pathogen genotype, several computations are performed: \itemize{
#' \item appearance: time to first appearance (as propagule);
#' \item R_infection: time to first true infection of a resistant host;
#' \item R_invasion: time when the number of infections of resistant hosts reaches a threshold above which
#' the genotype is unlikely to go extinct.}
#' The value Nyears + 1 timestep is used if the genotype never appeared/infected/invaded.
#' @return A list containing, for each required type of output, a matrix summarising the output.
#' Each matrix can be written in a txt file (if writeTXT=TRUE), and illustrated in a graphic (if graphic=TRUE).
#' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018). Assessing the durability and efficiency
#' of landscape-based strategies to deploy plant resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
#' @seealso \link{epid_output}
#' @examples
#' \dontrun{
#' demo_landsepi()
#' }
#' @include RcppExports.R Math-Functions.R graphics.R
#' @importFrom grDevices colorRampPalette dev.off graphics.off gray png tiff
#' @importFrom utils write.table
#' @export
evol_output <- function(types = "all", time_param, Npoly, cultivars_param, genes_param, thres_breakdown = 50000,
                        writeTXT = TRUE, graphic = TRUE, path = getwd()) {

  valid_outputs <- c("evol_patho", "evol_aggr", "durability")
  if (types[1] == "all") {
    types <- valid_outputs
  }
  if (is.na(sum(match(types, valid_outputs)))) {
    stop(paste("Error: valid types of outputs are", paste(valid_outputs, collapse = ", ")))
  }

  ## Time parameters
  Nyears <- time_param$Nyears
  nTSpY <- time_param$nTSpY
  nTS <- Nyears * nTSpY ## Total number of time-steps

  ## Host parameters
  Nhost <- length(cultivars_param$name)
  cultivars_genes_list <- cultivars_param$cultivars_genes_list

  # find the indices of resistant hosts
  Rhosts <- (1:Nhost)[as.logical(lapply(cultivars_genes_list, function(x) {
    length(x) > 0
  }))]

  ## Gene parameters
  Nlevels_aggressiveness <- genes_param$Nlevels_aggressiveness
  Ngenes <- length(Nlevels_aggressiveness)
  Npatho <- prod(Nlevels_aggressiveness)
  

  #### IMPORTATION OF THE SIMULATION OUTPUT
  P <- as.list(1:nTS)
  L <- as.list(1:nTS)
  I <- as.list(1:nTS)
  index <- 0

  # print("Reading binary files to compute evolutionary model outputs")
  for (year in 1:Nyears) {
    binfileP <- file(paste(path, sprintf("/P-%02d", year), ".bin", sep = ""), "rb")
    P.tmp <- readBin(con = binfileP, what = "int", n = Npoly * Npatho * nTSpY, size = 4, signed = T, endian = "little")
    binfileI <- file(paste(path, sprintf("/I-%02d", year), ".bin", sep = ""), "rb")
    I.tmp <- readBin(con = binfileI, what = "int", n = Npoly * Npatho * Nhost * nTSpY, size = 4, signed = T, endian = "little")
    binfileL <- file(paste(path, sprintf("/L-%02d", year), ".bin", sep = ""), "rb")
    L.tmp <- readBin(con = binfileL, what = "int", n = Npoly * Npatho * Nhost * nTSpY, size = 4, signed = T, endian = "little")


    for (t in 1:nTSpY) {
      P[[t + index]] <- matrix(P.tmp[((Npatho * Npoly) * (t - 1) + 1):(t * (Npatho * Npoly))], ncol = Npatho, byrow = T)
      L[[t + index]] <- array(
        data = L.tmp[((Npatho * Npoly * Nhost) * (t - 1) + 1):(t * (Npatho * Npoly * Nhost))],
        dim = c(Nhost, Npatho, Npoly)
      )
      I[[t + index]] <- array(
        data = I.tmp[((Npatho * Npoly * Nhost) * (t - 1) + 1):(t * (Npatho * Npoly * Nhost))],
        dim = c(Nhost, Npatho, Npoly)
      )
    }

    index <- index + nTSpY
    close(binfileP)
    close(binfileL)
    close(binfileI)
  }

  ####
  ####          COMPUTATION OF OUTPUTS
  ####
  res <- list()

  ## Matrix to switch from patho to aggr
  pathoToAggr <- matrix(NA, nrow = Npatho, ncol = Ngenes)
  for (p in 1:Npatho) {
    pathoToAggr[p, ] <- switch_patho_to_aggr(p - 1, Ngenes, Nlevels_aggressiveness) + 1
  } ## -1/+1 because indices start at 0 in function switch

  ####  Summary of P, L, I for each pathogen genotype
  P_patho <- NULL
  L_patho <- NULL
  I_patho <- NULL
  IL_patho_Rhost <- NULL ## sum of L and I on resistant hosts
  for (t in 1:nTS) {
    P_patho <- cbind(P_patho, apply(P[[t]], 2, sum))
    L_patho <- cbind(L_patho, apply(L[[t]], 2, sum))
    I_patho <- cbind(I_patho, apply(I[[t]], 2, sum))
    if (length(Rhosts) > 0 & Npatho > 1) {
      IL_patho_Rhost <- cbind(IL_patho_Rhost, apply(L[[t]][Rhosts, , ] + I[[t]][Rhosts, , ], 1 + (length(Rhosts) > 1), sum))
    } ## length(Rhosts)=2 --> need to use the 2nd dimension
  }

  ## Genotype frequencies
  I_pathoProp <- matrix(NA, nrow = Npatho, ncol = nTS)
  Itot <- apply(matrix(I_patho, ncol = nTS), 2, sum) ## matrix to avoid pb if only 1 row
  for (p in 1:Npatho) {
    for (t in 1:nTS) {
      if (Itot[t] > 0) {
        I_pathoProp[p, t] <- I_patho[p, t] / Itot[t]
      }
    }
  }

  ## Computed variables
  if (length(thres_breakdown) == 1) {
    names(thres_breakdown) <- "R_invasion"
  } else {
    names(thres_breakdown) <- paste("R_invasion_", 1:length(thres_breakdown), sep = "")
  }
  output_variables <- c("appearance", "R_infection", names(thres_breakdown))

  ## "time_to_first" variables for each pathogen (must be computed in any case)
  patho_evol <- matrix(NA, ncol = 3 + length(thres_breakdown), nrow = Npatho)
  colnames(patho_evol) <- c(output_variables, "finalPopSize_LI")
  rownames(patho_evol) <- paste("patho_", 1:Npatho, sep = "")
  colnames(pathoToAggr) <- genes_param$name
  rownames(pathoToAggr) <- paste("patho_", 1:Npatho, sep = "")

  for (patho in 1:Npatho) {
    patho_evol[patho, "appearance"] <- min(which(P_patho[patho, ] > 0), nTS + 1, na.rm = TRUE)
    patho_evol[patho, "R_infection"] <- min(which(IL_patho_Rhost[patho, ] > 0), nTS + 1, na.rm = TRUE)
    patho_evol[patho, "finalPopSize_LI"] <- as.numeric(sum(L_patho[patho, nTS], I_patho[patho, nTS]))
    for (k in 1:length(thres_breakdown)) {
      patho_evol[patho, names(thres_breakdown)[k]] <- min(which(IL_patho_Rhost[patho, ] > thres_breakdown[k]),
        nTS + 1,
        na.rm = TRUE
      )
    }
  }
  if (writeTXT & sum(is.element(types, c("evol_patho"))) > 0) {
    write.table(patho_evol, file = paste(path, "/evol_patho", ".txt", sep = ""), sep = ",")
    write.table(pathoToAggr, file = paste(path, "/evol_patho2aggr", ".txt", sep = ""), sep = ",")
  }
  res[["evol_patho"]] <- patho_evol

  ## "Time_to_first" variables for each level of aggressiveness
  if (sum(is.element(types, c("evol_aggr", "durability"))) > 0) {
    aggr_evol <- list()
    for (g in 1:Ngenes) {
      aggr_evol[[g]] <- matrix(NA, nrow = Nlevels_aggressiveness[g], ncol = length(output_variables))
      colnames(aggr_evol[[g]]) <- output_variables
      rownames(aggr_evol[[g]]) <- paste("level_", 1:Nlevels_aggressiveness[g], sep = "")
      for (a in 1:nrow(aggr_evol[[g]])) {
        for (k in colnames(aggr_evol[[g]])) {
          index_patho <- pathoToAggr[, g] == a
          aggr_evol[[g]][a, k] <- min(patho_evol[index_patho, k])
        } ## for k
      } ## for a
      if (writeTXT & sum(is.element(types, c("evol_aggr"))) > 0) {
        write.table(aggr_evol[[g]], file = paste(path, "/evol_aggr_", genes_param$name[g], ".txt", sep = ""), sep = ",")
      }
    } ## for g
    names(aggr_evol) <- genes_param$name
    res[["aggr_evol"]] <- aggr_evol

    ## Durability relative to a given criterion
    if (sum(is.element(types, c("durability"))) > 0) {
      durability <- matrix(NA, nrow = 1, ncol = Ngenes)
      colnames(durability) <- names(aggr_evol)
      for (g in 1:Ngenes) {
        finalLevel <- nrow(aggr_evol[[g]]) ## last line
        criterion <- ncol(aggr_evol[[g]]) ## last column
        durability[, g] <- aggr_evol[[g]][finalLevel, criterion]
      }
      if (writeTXT & sum(is.element(types, c("durability"))) > 0) {
        write.table(durability, file = paste(path, "/durability", ".txt", sep = ""), sep = ",")
      }
      res[["durability"]] <- durability
    } ## if durability
  } ## if aggr_evol


  ####
  ####               GRAPHICS
  ####
  if (graphic) {
    graphics.off()

    ## Dynamics of pathotype frequencies (gray scales)
    I_aggrProp <- list()
    for (g in 1:Ngenes) {
      I_aggrProp[[g]] <- matrix(0, nrow = Nlevels_aggressiveness[g], ncol = nTS)
      for (t in 1:nTS) {
        for (p in 1:Npatho) {
          aggr <- pathoToAggr[p, g]
          I_aggrProp[[g]][aggr, t] <- I_aggrProp[[g]][aggr, t] + I_pathoProp[p, t]
        }
      }

      if (Nlevels_aggressiveness[g] > 1) {
        tiff(
          filename = paste(path, "/freqPatho_gene", g, ".tiff", sep = ""), width = 100, height = 100,
          units = "mm", compression = "lzw", res = 300
        )
        par(xpd = F, mar = c(4, 4, 2, 2))
        plot_freqPatho(genes_param$name[g], Nlevels_aggressiveness[g], I_aggrProp[[g]], nTS, Nyears, nTSpY)
        dev.off()
      }
    }

    ## Dynamics of genotypes frequencies (curves)
    # COL <- c("#FF5555","#4F94CD","darkolivegreen4","#CD950C","black") ## colors: red, blue, green, orange, black
    tiff(
      filename = paste(path, "/freqPathoGenotypes.tiff", sep = ""), width = 180, height = 110,
      units = "mm", compression = "lzw", res = 300
    )
    par(xpd = FALSE, mar = c(4, 4, 0, 9))
    plot(0, 0, type = "n", xlim = c(1, nTS), bty = "n", las = 1, ylim = c(0, 1), xaxt = "n", xlab = "", ylab = "Frequency")
    if (Nyears == 1) {
      axis(1, at = round(seq(1, nTS, length.out = 8)), las = 1)
      title(xlab = "Evolutionnary time (days)")
    } else {
      axis(1,
        las = 1, at = seq(1, nTS + 1, nTSpY * ((Nyears - 1) %/% 10 + 1)),
        labels = seq(0, Nyears, ((Nyears - 1) %/% 10 + 1))
      )
      title(xlab = "Evolutionnary time (years)")
    }

    for (patho in 1:Npatho) {
      lines(I_pathoProp[patho, ], col = patho, lty = patho, lwd = 1.5)
    }

    par(xpd = TRUE)
    legend(nTS * 1.05, .5,
      title.adj = -.01, bty = "n", title = "Pathogen genotype:",
      legend = 1:Npatho, lty = 1:Npatho, lwd = 2, col = 1:Npatho
    )
    dev.off()
  } ## if graphic

  return(res[match(types, valid_outputs)])
}


#' @title Generation of a video
#' @name video
#' @description Generates a video showing the epidemic dynamics on a map representing the cropping landscape.
#' (requires ffmpeg library).
#' @param audpc A dataframe containing audpc outputs (generated through epid_output). 1 year per line and
#' 1 column per cultivar, with an additional column for the average audpc in the landscape.
#' @param time_param list of simulation parameters:\itemize{
#' \item Nyears = number cropping seasons,
#' \item nTSpY = number of time-steps per cropping season.
#' }
#' @param Npatho number of pathogen genotypes.
#' @param landscape a sp object containing the landscape.
#' @param area a vector containing polygon areas (must be in square meters).
#' @param rotation a dataframe containing for each field (rows) and year (columns, named "year_1", "year_2", etc.),
#' the index of the cultivated croptype. Importantly, the matrix must contain 1 more column than the real number
#' of simulated years.
#' @param croptypes a dataframe with three columns named 'croptypeID' for croptype index,
#' 'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar within the croptype.
#' @param croptype_names a vector of croptype names (for legend).
#' @param cultivars_param a list of parameters associated with each host genotype (i.e. cultivars)
#' when cultivated in pure crops:\itemize{
#' \item name = vector of cultivar names,
#' \item max_density = vector of maximum host densities (per square meter) at the end of the cropping season,
#' \item cultivars_genes_list = a list containing, for each host genotype, the indices of carried resistance genes.
#' }
#' @param audpc100S AUDPC of a 100% susceptible landscape, used as a reference.
#' @param keyDates a vector of times (in timesteps) where to draw vertical lines in the AUDPC graphic. Usually
#' used to delimit durabilities of the resistance genes. No line is drawn if keyDates=NULL (default).
#' @param nMapPY an integer specifying the number of epidemic maps per year to generate.
#' @param path path where binary files are located and where the video will be generated.
#' @details The left panel shows the year-after-year dynamics of AUDPC, relative to a fully susceptible landscape,
#' for each cultivar as well as the global average. The right panel illustrates the landscape,
#' where fields are hatched depending on the cultivated croptype, and coloured depending on the prevalence of the disease.
#' Note that up to 9 different croptypes can be represented properly in the right panel.
#' @return A video file of format mp4.
#' @examples
#' \dontrun{
#' demo_landsepi()
#' }
#' @include RcppExports.R graphics.R
#' @importFrom sf st_read
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach "%dopar%"
#' @importFrom foreach foreach
#' @importFrom parallel detectCores
#' @importFrom grDevices colorRampPalette dev.off graphics.off gray png tiff
#' @export
video <- function(audpc, time_param, Npatho, landscape, area, rotation, croptypes, croptype_names = c(), cultivars_param, audpc100S, keyDates = NULL,
                  nMapPY = 5, path = getwd()) {

  if (system("ffmpeg -version", ignore.stdout = TRUE) == 127) {
    stop("You need to install ffmpeg before generating videos. Go to https://ffmpeg.org/download.html")
  }

  ## Time & graphic parameters
  Nyears <- time_param$Nyears
  nTSpY <- time_param$nTSpY
  nTS <- Nyears * nTSpY ## Total number of time-steps

  ## Landscape
  rotation <- data.frame(rotation)
  Npoly <- length(area)

  ## croptype parameters (for legend)
  if (length(croptype_names) == 0) {
    croptype_names <- paste("Croptype", unique(croptypes$croptypeID))
  }
  Ncroptypes <- length(croptype_names)

  ## Host parameters
  cultivar_names <- cultivars_param$name
  cultivars_genes_list <- cultivars_param$cultivars_genes_list
  max_density <- cultivars_param$max_density
  Nhost <- length(max_density)

  ## Calculation of the carrying capacity
  K <- array(dim = c(Npoly, Nhost, Nyears))
  for (y in 1:Nyears) {
    if (ncol(rotation) == 1) {
      index_year <- 1
    } else {
      index_year <- y
    }
    ts_year <- ((y - 1) * nTSpY + 1):(y * nTSpY)

    for (poly in 1:Npoly) {
      indices_croptype <- grep(rotation[poly, index_year], croptypes$croptypeID)
      for (i in indices_croptype) {
        index_host <- croptypes[i, "cultivarID"] + 1 ## +1 because C indices start at 0
        prop <- croptypes[i, "proportion"]
        K[poly, index_host, y] <- floor(area[poly] * max_density[index_host] * prop)
      }
    } ## for poly
  } ## for y

  #### IMPORTATION OF THE SIMULATION OUTPUT
  I <- as.list(1:nTS)
  R <- as.list(1:nTS)
  index <- 0

  for (year in 1:Nyears) {
    binfileI <- file(paste(path, sprintf("/I-%02d", year), ".bin", sep = ""), "rb")
    I.tmp <- readBin(con = binfileI, what = "int", n = Npoly * Npatho * Nhost * nTSpY, size = 4, signed = T, endian = "little")
    binfileR <- file(paste(path, sprintf("/R-%02d", year), ".bin", sep = ""), "rb")
    R.tmp <- readBin(con = binfileR, what = "int", n = Npoly * Npatho * Nhost * nTSpY, size = 4, signed = T, endian = "little")

    for (t in 1:nTSpY) {
      I[[t + index]] <- array(
        data = I.tmp[((Npatho * Npoly * Nhost) * (t - 1) + 1):(t * (Npatho * Npoly * Nhost))],
        dim = c(Nhost, Npatho, Npoly)
      )
      R[[t + index]] <- array(
        data = R.tmp[((Npatho * Npoly * Nhost) * (t - 1) + 1):(t * (Npatho * Npoly * Nhost))],
        dim = c(Nhost, Npatho, Npoly)
      )
    }

    index <- index + nTSpY
    close(binfileI)
    close(binfileR)
  }


  print("Generate video (and create directory maps)")
  graphics.off()
  dir.create("maps")

  ## Standard colours
  RED <- "#FF5555"
  BLUE <- "#4F94CD"
  grad_grey <- colorRampPalette(c("black", "white"))
  GRAY <- grad_grey(Nhost + 1)

  ## Colour palette for prevalence
  nCol <- 21
  grad_whiteYellowRed <- colorRampPalette(c("white", "#FFFF99", "#990000"))
  col_prev <- grad_whiteYellowRed(nCol)

  PCH <- 15:(15 + Nhost - 1)
  PCH.tot <- 0
  LTY.tot <- 1
  COL.tot <- RED

  title <- "Proportion of diseased hosts"
  maxColorScale <- 1
  intvls <- maxColorScale * sort(seq(0, 1, l = nCol), decreasing = FALSE)
  density_hatch <- c(0, 8, 16, 4, 12, 2, 6, 10, 14)[1:Ncroptypes]
  angle_hatch <- c(0, 30, 120, 150, 60, 90, 180, 45, 135)[1:Ncroptypes]
  # legend_prev <- sprintf("%.3f", intvls)
  legend_prev <- paste(intvls * 100, "%", sep = "")

  registerDoParallel(parallel::detectCores())
  foreach(y = 1:Nyears) %dopar% {
    print(paste("Year", y, "/", Nyears))
    if (ncol(rotation) == 1) {
      index_year <- 1
    } else {
      index_year <- y
    }

    K_poly <- apply(as.data.frame(K[, , index_year]), 1, sum, na.rm = TRUE) ## as.data.frame in case Nhost==1 & Npatho==1
    for (d in round(seq(1, nTSpY, length.out = nMapPY))) {
      subtitle <- paste("Year =", y, "   Day =", d)
      ## Proportion of each type of host relative to carrying capacity
      ts <- (y - 1) * nTSpY + d

      propI <- apply(I[[ts]], 3, sum) / K_poly
      propR <- apply(R[[ts]], 3, sum) / K_poly
      propIR <- (propI + propR)
      intvlsIR <- findInterval(propIR, intvls)

      png(
        filename = paste(path, "/maps/HLIR_", sprintf("%02d", y), "-", sprintf("%03d", d), ".png", sep = ""),
        width = 700, height = 350, units = "mm", res = 72, type = "cairo", bg = "white", antialias = "default"
      )
      par(mfrow = c(1, 2), cex = 2, xpd = NA, mar = c(9.5, 5, 4, 2))

      ## moving AUDPC
      plot(0, 0,
        type = "n", bty = "n", xlim = c(1, Nyears), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
        main = "Disease severity relative to a fully susceptible landscape"
      )
      axis(1, at = seq(1, Nyears, by = (Nyears - 1) %/% 10 + 1))
      ticksMarks <- round(seq(0, 1, length.out = 5), 2)
      axis(2, at = ticksMarks, labels = paste(100 * ticksMarks, "%"), las = 1)
      title(xlab = "Years", mgp = c(2, 1, 0))
      for (host in 1:Nhost) {
        lines(1:y, audpc[1:y, host] / audpc100S, type = "o", lwd = 2, lty = host, pch = PCH[host], col = GRAY[host])
      }
      lines(1:y, audpc[1:y, "total"] / audpc100S, type = "o", lwd = 2, lty = LTY.tot, pch = PCH.tot, col = COL.tot)
      ## Add lines for durability
      if (!is.null(keyDates)) {
        for (k in 1:length(keyDates)) {
          date <- ceiling(keyDates[k] / nTSpY)
          if (!is.na(date) & date <= y) {
            abline(v = date, col = BLUE, lty = 2, lwd = 4, xpd = FALSE)
          }
        }
      }

      legend(Nyears / 3, -1 / 5,
        bty = "n", legend = c(cultivar_names, "Whole landscape"),
        cex = 1, lty = c(1:Nhost, LTY.tot), lwd = 2, pch = c(PCH, PCH.tot), pt.cex = c(rep(1, Nhost), 1),
        col = c(GRAY[1:Nhost], COL.tot), seg.len = 2.5
      )

      ## Map dynamic of diseased hosts
      plotland(
        landscape, col_prev[intvlsIR], density_hatch[rotation[, index_year] + 1], angle_hatch[rotation[, index_year] + 1],
        col_prev, density_hatch, angle_hatch, title, subtitle,
        croptype_names,
        legend_prev, " "
      )
      dev.off()
    } ## for d
  } ## for y

  ## Convert png files to mp4 video
  # cat *.png | ffmpeg -f image2pipe -framerate 3 -i - -vcodec libx264 -vb 1024k  video.mp4
  # fast
  # ffmpeg -f image2 -framerate 2.5 -i HLIR_%04d.png  -vcodec libx264 -vb 1024k -crf 25 -pix_fmt yuv420p video.mp4
  # slow
  # ffmpeg -f image2 -framerate 2.5 -i HLIR_%04d.png -vcodec vp9 -vb 1024k -crf 25 -pix_fmt yuv420p video.webm
  setwd("maps")
  # rename png file to HLIR_XXXX.png
  lfiles <- list.files(pattern = "*.png")
  file.rename(lfiles, paste0("HLIR_", sprintf("%04d", 1:length(lfiles)), ".png"))
  title_video <- paste("video.mp4", sep = "")
  command_line <- paste("ffmpeg -f image2 -framerate ", nMapPY / 2,
    " -i HLIR_%04d.png -vcodec libx264 -vb 1024k -crf 25 -pix_fmt yuv420p 2>&1 > /dev/null ",
    title_video,
    sep = ""
  )
  system(command_line, ignore.stdout = TRUE, ignore.stderr = TRUE)

  system(paste("mv ", title_video, " ", "../."))
  setwd("../")
  print("remove directory maps")
  system("rm -rf maps")
}
