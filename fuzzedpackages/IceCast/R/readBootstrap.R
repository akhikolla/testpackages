#' Read in individual binary files of monthly observation data. The observations
#' are from the monthly sea ice concentration obtained from the National
#' Aeronautics and Space Administration (NASA) satellites Nimbus-7
#' SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm.
#' The results are distributed by the National Snow and Ice Data Center (NSIDC)
#' (Comiso 2017). Functions assume file name conventions are the
#' same as used by NSIDC.
#' @title Read individual bootstrap binary file
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center.
#'             doi: \url{https://doi.org/10.5067/7Q8HCCWS4I0R}
#'
#' @param file_name file name for binary bootstrap data
#' @param nX dimension in the x (defaults to value for Northern Polar
#'           stereographic grid: 304)
#' @param nY dimension in the y (defaults to value for Northern Polar
#'           stereographic grid: 448)
#' @return numeric vector of concentrations
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' #fileName should be the binary file
#' rawData <- read_bootstrap(file_name)
#'}
read_bootstrap <- function(file_name, nX = 304, nY = 448) {
  to.read <- file(file_name, "rb")
  dat <- readBin(to.read, integer(), n = nX*nY, size = 2, endian = "little")/10
  close(to.read)
  return(dat)
}

#' Function to process monthly bootstrap data over multiple years. The
#' observations are from the monthly sea ice concentration obtained from the
#' National Aeronautics and Space Administration (NASA) satellites Nimbus-7
#' SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm. The
#' resultsare distributed by the National Snow and Ice Data Center (NSIDC)
#' (Comiso 2017). Functions assume file name conventions are the same as used
#'by NSIDC.
#' @title Read in a set of bootstrap observations over a set of year
#' @references Bootstrap sea ice concentration:
#'             Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             {Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center}
#' @param start_year first year to read in
#' @param end_year last year to read in
#' @param file_folder folder in which binary files are stored
#' @param version either 2 or 3 indicating which version of the bootstrap data
#'                you are using
#' @param nX longitude dimension
#' @param nY latitude dimension
#' @details raw binary files for 2012-2013 are included in the package as an
#'          example
#' @export
#' @return bootstrap observations sorted into array of dimension: year x month x
#'         lon x lat
#' @examples
#' \dontrun{
#' #my_file_path should be a file path where the 1983 binary files are stored
#' observed_demo <- read_monthly_BS(start_year = 1983, end_year = 1983,
#'                              file_folder = my_file_path)
#' }
read_monthly_BS <- function(start_year, end_year, file_folder, version, nX = 304,
                          nY = 448) {
  years <- start_year:end_year; nYears <- length(years)
  obs <- array(dim = c(nYears, 12, nX, nY))
  stopifnot(version == 2 || version == 3 || version == 3.1 )
  for (i in 1:nYears) {
    for (j in 1:12) {
      if (version == 2) { #no missing data in V2
        file_name <- Sys.glob(paste(file_folder, sprintf('bt_%i%02d_*_v02_n.bin',
                                                       years[i], j), sep = ""))
        obs[i, j, ,nY:1] <- read_bootstrap(file_name)
      } else if (version == 3 & !(j == 12 & years[i] == 1987)
                 & !(j == 1 & years[i] == 1988)) {
        #missing Dec 1987 and Jan 1988 in V3 because "major data gap in the
        #SSM/I data occurs from 03 December 1987 to 13 January 1988"
        file_name <- Sys.glob(paste(file_folder,
                                   sprintf('bt_%i%02d_*_v03_n.bin', years[i], j),
                                   sep = ""))
        obs[i, j, ,nY:1] <- read_bootstrap(file_name)
      } else if (version == 3.1 & !(j == 12 & years[i] == 1987)
                 & !(j == 1 & years[i] == 1988)) {
        #missing Dec 1987 and Jan 1988 in V3 because of "major data gap in the
        #SSM/I data occurs from 03 December 1987 to 13 January 1988"
        file_name <- Sys.glob(paste(file_folder,
                                   sprintf('bt_%i%02d_*_v3.1_n.bin', years[i], j),
                                   sep = ""))
        obs[i, j, ,nY:1] <- read_bootstrap(file_name)
      } else {
        stopifnot(version == 3 || version == 3.1)
        obs[i, j, ,nY:1] <- NA
      }
    }
  }
  return(obs)
}
