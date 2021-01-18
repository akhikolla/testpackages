
#' Machine temperature data.
#'
#' @description Temperature sensor data of an internal component of a large, industrial machine. The data contains three known anomalies.
#' The first anomaly is a planned shutdown of the machine. The second anomaly is difficult to detect and directly led to the third anomaly, a catastrophic failure of the machine.
#' The data consists of 22695 observations of machine temperature recorded at 5 minute intervals along with the date and time of the measurement. The data was obtained from the Numenta
#' Anomaly Benchmark \insertCite{AHMAD2017134}{anomaly}, which can be found at https://github.com/numenta/NAB.
#' 
#' @name machinetemp
#' 
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(machinetemp)
#'
#' @rdname machine-temperature-data
#'
#' @format A dataframe with 22695 rows and 2 columns. The first column contains the date and time of the temperature measurement. The second column contains the
#' machine temperature.
#'
#' @references \insertRef{AHMAD2017134}{anomaly}
#' 
NULL

#' Kepler Lightcurve  data.
#' @description One of the most successful approaches for the detection of exoplanets is the so called transit method: A star's brightness is continuously measured over time by a
#' powerful telescope. If one or multiple planets orbit this star the recorded luminosity of the star will exhibit periodically recurring dips due to the transits of the planet in
#' front of the telescope's lens -- an effect comparable to that of an eclipse. Given how small planets are compared to stars the transit signals are known to be very weak. 
#'
#' The stars included in this file all have known exoplanets with the following periods: 
#'
#' Kepler 1871056: 2 planets with orbital periods of 40.8 and 140.1 days \cr
#' Kepler 2307415: 2 planets with orbital periods of 4.61 and 12.12 days \cr
#' Kepler 3102384: 2 planets with orbital periods of 10.57 and 523.9 days \cr
#' Kepler 3231341: 4 planets with orbital periods of 4.24, 8.15, 12.33, and 19.00 days \cr
#' Kepler 3447722: 3 planets with orbital periods of 10.30, 16.09, and 35.68 days \cr
#' Kepler 4139816: 4 planets with orbital periods of 3.34, 7.82, 20.06, and 46.18 days \cr
#' Kepler 10965588: 1 planet with orbital period of 62.89 days \cr
#' 
#' More information about the exoplanets above and more data can be found at \cr https://exoplanetarchive.ipac.caltech.edu/index.html  
#' 
#' This research has made use of the NASA Exoplanet Archive, which is operated by the California Institute of Technology, under contract with the National Aeronautics and
#' Space Administration under the Exoplanet Exploration Program.
#'
#' @name Lightcurves
#' 
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(Lightcurves)
#'
#' @rdname lightcurves-data
#'
#' @format A list of seven dataframes named "Kepler1871056", "Kepler2307415", "Kepler3102384", "Kepler3231341", "Kepler3447722", "Kepler4139816", and "Kepler10965588". Each dataframe
#' consists of two columns called "Brightness" and "Day", containing measurements of a star's brightness and the measurement's timestamp respectively.
#'
#' 
NULL


