#' Example library
#'
#' @description Example spectral library to be used in microplastic evaluations. An example on how to format a \code{\link[=SpectralReference-class]{SpectralReference}} object.
#' 
#' It was taken form:
#' 
#' Primpke, S., Wirth, M., Lorenz, C., Gerdts, G. 2018. Reference database design for the automated analysis of microplastic samples based on Fourier transform infrared (FTIR) spectroscopy. Analytical and Bioanalytical Chemistry 410: 5131-5141. You might access to the article here \url{https://doi.org/10.1007/s00216-018-1156-x}
#'
#' @format an S4 \code{\link[=SpectralReference-class]{SpectralReference}} object with:
#' \describe{
#'   \item{substances}{character vector with the 270 names of the polymers included in Spectra}
#'   \item{clusterlist}{integer vector with the corresponding cluster of each polymer described}
#'   \item{clusternames}{the name of the cluster, sorted in ascending order according to clusterlist}
#'   \item{Spectra}{matrix with the absorbance readings of each polymer (270 x 1176)}
#'   \item{wavenumbers}{read wavenumbers. Correspond to Spectra ncol}
#' }
#' @source \url{https://doi.org/10.1007/s00216-018-1156-x}
"primpke"
