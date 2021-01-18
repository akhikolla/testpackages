#' Computed mappings for predictions for September 1993-2008
#'
#' Output of the \code{create_mapping} function for September 1993-2008 using
#' predictions from the European Center for Medium-Range Weather Forecasts
#' (ECMWF) ensemble converted to a Polar stereographic grid.
#'
#' @docType data
#' @keywords datasets
#' @usage data(pred_maps)
#' @references Copernicus Climate Change Service (2019). Description of the c3s seasonal multi-system.\url{https://confluence.ecmwf.int/display/COPSRV/Description+of+the+C3S+seasonal+multi-system}
#'
#'             Sea Ice Prediction Network (2019). Sea ice prediction network
#'             predictability portal. \url{https://atmos.uw.edu/sipn/.}
#' @examples
#' data(pred_maps)
"pred_maps"

#' Discrepancy maps for September 1993-2007 (lead time 2.5-months)
#'
#' The object \code{discrep} is obtained from running the \code{createMapping}
#' function for September 1993-2007. The predictions used are from European Center
#' for Medium-Range Weather Forecasts (ECMWF) at a 2.5-month lead time and are
#' converted to a Polar Stereographic grid. Model output is available from
#' the Sea Ice Prediction Network Predicatability Portal or the Copernicus
#' Climate Change Service data store. The observations used are from the monthly
#' sea ice concentration obtained fromthe National Aeronautics and Space
#' Administration (NASA) satellites Nimbus-7 SMMR and DMSP SSM/I-SSMIS and
#' processed by the bootstrap algorithm. The results are distributed by the
#' National Snow and Ice Data Center (Comiso 2017).
#' @docType data
#' @format Object obtained from the \code{createMapping} function (see details)
#'
#' @details The object \code{discrep} is obtained from running the
#'          \code{createMapping} function. It is a list of four objects where
#'          \code{startYear} and \code{endYear} give the
#'          first year and last year that were mapped. The variables
#'          \code{obsList} and \code{predList} are lists of arrays with one
#'          3-dimensional array for each region. The first dimension is for the
#'          year. The other two dimensions are for the fixed points'
#'          y-coordinates, the mapped points' x-coordinates, the mapped points'
#'          y-coordinates, the length of the mapping vectors in the x-direction,
#'          the length of the vectors in the y-direction, and the angles of the
#'          mapping vectors.
#'
#' @keywords datasets
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center. doi: \url{https://doi.org/10.5067/7Q8HCCWS4I0R}
#'
#'             Copernicus Climate Change Service (2019). Description of the c3s
#'             seasonal multi-system.\url{https://confluence.ecmwf.int/display/COPSRV/Description+of+the+C3S+seasonal+multi-system}
#'
#'             Sea Ice Prediction Network (2019). Sea ice prediction network
#'             predictability portal. \url{https://atmos.uw.edu/sipn/.}
#' @examples
#' data(discrep)
#' names(discrep)
"discrep"


#' Observed sea ice September 2006-2007
#'
#' The object \code{observed} is an array obtained from the function
#' \code{readMonthlyBS}for startYear = 2006 and endYear = 2007. It gives the
#' observed sea ice concentrations arranged in an array of dimension of year x
#' month x lon x lat. The observations are from the monthly sea ice
#' concentration obtained from the National Aeronautics and Space Administration
#' (NASA) satellites Nimbus-7 SMMR and DMSP SSM/I-SSMIS and processed by the
#' bootstrap algorithm. The results are distributed by the National Snow and Ice
#'Data Center (Comiso 2017).
#' @docType data
#' @format array of dimension of 2 x 12 x 304 x 448 (year x month x longitude
#' x latitude)
#' @keywords datasets
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             {Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center}
#' @examples
#' data(obsSep2006_2007)
#' dim(obsSep2006_2007)
"obsSep2006_2007"

#' Observed sea ice September 2008
#'
#' The object \code{observed} is an binary matrix of dimension lon x lat that
#' indicates whether sea ice concentration was at least 15\%. The observations
#' are from the monthly sea ice concentration obtained from the National
#' Aeronautics and Space Administration (NASA) satellites Nimbus-7 SMMR and
#' DMSP SSM/I-SSMIS and processed by the bootstrap algorithm. The results are
#' distributed by the National Snow and Ice Data Center (Comiso 2017).
#' @docType data
#' @format array of dimension of 2 years x 12 months x longitude x latitude
#' @keywords datasets
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center.
#'             doi: \url{https://doi.org/10.5067/7Q8HCCWS4I0R}
#' @examples
#' data(obsSep2008)
#' dim(obsSep2008)
"obsSep2008"

#' Ensemble estimated sea ice probability September 2008 (lead time 2.5 months)
#'
#' The object \code{sip2006_2007} is an array of the sea ice probability
#' predicted  from the European Center for Medium-Range Weather Forecasts
#' (ECMWF) ensemble converted to a Polar stereographic grid.
#'
#' @docType data
#' @format matrix of dimension 304 x 448 (longitude x latitude)
#' @keywords datasets
#' @references
#'             Copernicus Climate Change Service (2019). Description of the c3s
#'             seasonal multi-system. \url{https://confluence.ecmwf.int/display/COPSRV/Description+of+the+C3S+seasonal+multi-system}
#'
#'             Sea Ice Prediction Network (2019). Sea ice prediction network
#'             predictability portal. \url{https://atmos.uw.edu/sipn/.}#' @examples
#' data(sipSep2008)
#' dim(sipSep2008)
"sipSep2008"



#' Ensemble sea ice probability for September 2006-2007 (lead time 2.5 months)
#'
#' The object \code{sipSep2006_2007} is an array of the proportion of ensemble
#' members that have sea ice concentrations of at least 15\%. The predictions are
#' are from the European Center for Medium-Range Weather Forecasts (ECMWF)
#' ensemble and converted to a Polar stereographic grid.
#' @docType data
#' @format array of dimension of 2 x 304 x 448 (corresponding to year x
#'         longitude x latitude)
#' @keywords datasets
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center.
#'             doi: \url{https://doi.org/10.5067/7Q8HCCWS4I0R}
#'
#'             Copernicus Climate Change Service (2019). Description of the c3s
#'             seasonal multi-system. \url{https://confluence.ecmwf.int/display/COPSRV/Description+of+the+C3S+seasonal+multi-system}
#'
#'             Sea Ice Prediction Network (2019). Sea ice prediction network
#'             predictability portal. \url{https://atmos.uw.edu/sipn/.}
#' @examples
#' data(sipSep2006_2007)
#' dim(sipSep2006_2007)
"sipSep2006_2007"


#' Ensemble sea ice probability for September 2008 (lead time 2.5 months)
#'
#' The object \code{sipSep2008} is an array of the proportion of ensemble
#' members that have sea ice concentrations of at least 15\%. The predictions are
#' from the European Centerfor Medium-Range Weather Forecasts (ECMWF) ensemble
#' converted to the Polar stereographic grid.
#' @docType data
#' @format matrix of dimension 304 x 448 (longitude x latitude)
#' @references
#'             Copernicus Climate Change Service (2019). Description of the c3s
#'             seasonal multi-system.\url{https://confluence.ecmwf.int/display/COPSRV/Description+of+the+C3S+seasonal+multi-system}
#'
#'             Sea Ice Prediction Network (2019). Sea ice prediction network
#'             predictability portal. \url{https://atmos.uw.edu/sipn/.}
#' @keywords datasets
#' data(sipSep2008)
#' dim(sipSep2008)
"sipSep2008"

#' Example of mapped points
#'
#' Example of a set of mapped points organized as an n x 2 matrix of coordinates.
#' This is used to demonstrate the \code{makePolygons} function.
#' @docType data
#' @format matrix of 1027 x 2
#' @keywords datasets
#' @examples
#' data(mappedPoints)
#' head(mappedPoints)
#' plot(mappedPoints, type = "l")
"mappedPoints"

#' Example of a line that contains self-intersections
#'
#' Example of a line that contains self-intersections. We will use it to
#' demonstrate the functions that address these intersections.
#' @docType data
#' @format n x 2 matrix of coordinates
#' @keywords datasets
#' @examples
#' data(interEx)
#' plot(interEx)
"interEx"

#' Coordinates of an observed line segment
#'
#' Example of the coordinates for an observed line segment. We will use it to
#' demonstrate the \code{intLine} function.
#' @docType data
#' @format n x 2 matrix of coordinates
#' @keywords datasets
#' @examples
#' data(obsLEx)
#' head(obsLEx)
"obsLEx"

#' Coordinates of a predicted line segment
#'
#' Example of the coordinates for a predicted line segment. We will use it to
#' demonstrate the \code{intLine} function.
#' @docType data
#' @format n x 2 matrix of coordinates
#' @keywords datasets
#' @examples
#' data(predLEx)
#' head(predLEx)
"predLEx"

#' Coordinates of a line segment with self-intersections
#'
#' Example of a line segment with self-intersections. We will use it
#' to demonstrate the \code{untwistSec} function.
#' @docType data
#' @format n x 2 matrix of coordiantes
#' @keywords datasets
#' @examples
#' data(currSecEx)
#' head(currSecEx)
"currSecEx"

#' Binary matrix indicating where there is land
#'
#' Binary matrix of dimension 304 x 448 with value for 1 for land grid boxes and
#' 0 otherwise. Data are on a north Polar Stereographic grid with the land mask
#' simplified to match model output from the CM2.5 Forecast-oriented Low-Ocean
#' Resolution (FLOR) model produced by the National Oceanic and Atmospheric
#' Administration’s Geophysical Fluid Dynamics Laboratory converted to a Polar
#' Stereographic grid (Vecchi et al. 2014; Msadek et al. 2014).
#' Weights for converting to a polar stereograhic grid were obtained
#' from the spherical coordinate remapping and interpolation package
#'(SCRIP) (Jones 1997).
#' @docType data
#' @format 304 x 448 matix
#' @keywords datasets
#' @references Vecchi, Gabriel A., et al.
#'             \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical} cyclone activity."
#'             Journal of Climate 27.21 (2014): 7994-8016.
#'
#'             Msadek, R., et al.
#'             \href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditionsin seasonal predictions of Arctic sea ice extent."}
#'             Geophysical Research Letters 41.14 (2014): 5208-5215.
#'
#'             National Center for Atmospheric Research, 2017: Earth system grid
#'             at NCAR. \url{https://www.earthsystemgrid.org/home.html}.

#' @examples
#' data(land_mat)
#' image(land_mat, xaxt = "n", yaxt = "n")
"land_mat"

#' Binary predictions from ECMWF ensemble, September 1993-2018
#'
#' The object \code{ensemble_bin} is a binary array indicating if at least half of
#' the ensemble members have sea ice concentrations of at least 15\% from
#' September 1993-2018. The predictions are from the European Center
#' for Medium-Range Weather Forecasts (ECMWF) ensemble at a 2.5-month
#' lead time. They have been converted to a Polar stereographic grid.
#' @docType data
#' @format array of dimension of 26 x 304 x 448 (corresponding to year x
#'         longitude x 448 latitude)
#' @keywords datasets
#' @references Copernicus Climate Change Service (2019). Description of the c3s
#'             seasonal multi-system. \url{https://confluence.ecmwf.int/display/COPSRV/Description+of+the+C3S+seasonal+multi-system}
#'
#'             Sea Ice Prediction Network (2019). Sea ice prediction network
#'             predictability portal. \url{https://atmos.uw.edu/sipn/.}
#' @examples
#' data(ecmwf_bin)
"ecmwf_bin"


#' Sample parameter information for generating a contour
#'
#' Example list with two elements, \code{mu_est} and \code{sigma_est}, which
#' give the mean and covariance from which an example contour can be generated
#'
#' @docType data
#' @keywords datasets
#' @usage data(pars_1)
#' @examples
#' data(pars_1)
"pars_1"

#' Sample list of contours
#'
#' Example list of ten contours in the form of \code{SpatialPolygons} objects
#'
#' @docType data
#' @keywords datasets
#' @usage data(merged)
#' @examples
#' data(merged)
"merged"



#' Output from \code{create_mapping} function
#'
#' Sample output from the \code{create_mapping} function for observations
#' from September 1993-2007. It is a list of four objects. The first two items
#' in the list, \code{start_year} and \code{end_year}, give the first and last
#' year that were mapped. The second two items, \code{obs_list} and
#' \code{pred_list}, are lists of arrays with one 3-dimensional array for each
#' region. The first dimension is for the  year. The other two dimensions are
#' for the fixed points' y-coordinates, the mapped points' x-coordinates, the
#' mapped points' y-coordinates, the length of the mapping vectors in the
#' x-direction, the length of the vectors in the y-direction, and the angles of
#' the mapping vectors.
#'
#'
#' @docType data
#' @keywords datasets
#' @usage data(obs_maps)
#' @examples
#' data(obs_maps)
"obs_maps"


#' Sample collection of completely ice-filled regions
#'
#' Example of a collection of several completely ice-filled regions stored
#' as a single \code{SpatialPolygons} object
#' @docType data
#' @keywords datasets
#' @usage data(full)
#' @examples
#' data(full)
"full"

#' Sample output from \code{get_map}
#'
#' Example output from the \code{get_map} function for the Central Arctic region.
#'
#' @docType data
#' @keywords datasets
#' @usage data(map_curr_1)
#' @examples
#' data(map_curr_1)
"map_curr_1"


#' Observed sea ice for September 2008
#'
#' Array of dimension longitude by latitude. Binary indicate of whether
#' sea ice concentration of at least 15\% was observed. Computed from NASA
#' Boostrap sea ice concentration product (Comiso 2017).
#'
#' @docType data
#' @format array
#' @keywords datasets
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center.
#'             doi: \url{https://doi.org/10.5067/7Q8HCCWS4I0R}
#' @usage data(obs_9_2008)
"obs_9_2008"

#' Observed sea ice for September 2005-2007
#'
#' Array of dimension year x longitude by latitude. Binary indicate of whether
#' sea ice concentration of at least 15\% was observed. Computed from NASA
#' Boostrap sea ice concentration product (Comiso 2017).
#'
#' @docType data
#' @format array
#' @keywords datasets
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center.
#'             doi: \url{https://doi.org/10.5067/7Q8HCCWS4I0R}
#' @usage data(obs_9_2005_2007)
"obs_9_2005_2007"

#' Post-procesed ensemble forecast for 2005-2007
#'
#' Array of dimension year by longitude by latitude that gives example
#' forecasts post-processed with a contour model. The initial forecasts
#' are from the European Center for Medium-Range Weather Forecasts (ECMWF)
#' ensemble at a 2.5-month lead time. They have been converted to a Polar
#' stereographic grid. Model output is available from
#' the Sea Ice Prediction Network Predicatability Portal or the Copernicus
#' Climate Change Service data store.
#'
#' @docType data
#' @format array
#' @references Copernicus Climate Change Service (2019). Description of the c3s
#'             seasonal multi-system. \url{https://confluence.ecmwf.int/display/COPSRV/Description+of+the+C3S+seasonal+multi-system}
#'
#'             Sea Ice Prediction Network (2019). Sea ice prediction network
#'             predictability portal. \url{https://atmos.uw.edu/sipn/.}
#' @keywords datasets
#' @usage data(ppe_9_2005_2007)
"ppe_9_2005_2007"


#' Post-procesed ensemble forecast for 2008
#'
#' Array of dimension year by longitude by latitude.The initial forecasts
#' are from the European Center for Medium-Range Weather Forecasts (ECMWF)
#' ensemble at a 2.5-month lead time. They have been converted to a Polar
#' stereographic grid. Model output is available from
#' the Sea Ice Prediction Network Predicatability Portal or the Copernicus
#' Climate Change Service data store.
#'
#' @docType data
#' @format array
#' @keywords datasets
#' @references Copernicus Climate Change Service (2019). Description of the c3s
#'             seasonal multi-system. \url{https://confluence.ecmwf.int/display/COPSRV/Description+of+the+C3S+seasonal+multi-system}
#'
#'             Sea Ice Prediction Network (2019). Sea ice prediction network
#'             predictability portal. \url{https://atmos.uw.edu/sipn/.}
#' @usage data(ppe_9_2008)
"ppe_9_2008"


#' Climatology forecast for 2005-2007.
#'
#' Proportion of times in the preceding ten years that sea ice concentration of
#' at least 15\%  was observed in each  grid box. Array of dimension year by
#' longitude by latitude. Computed from NASA Boostrap sea ice concentration
#' product.
#'
#' @docType data
#' @format array
#' @keywords datasets
#' @usage data(clim_9_2005_2007)
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center.
#'             doi: \url{https://doi.org/10.5067/7Q8HCCWS4I0R}
"clim_9_2005_2007"

#' Climatology forecast for 2008.
#'
#' Proportion of times in the preceding ten years that sea ice concentration of
#' at least 15\%  was observed in each  grid box. Array of dimension year by
#' longitude by latitude. Computed from NASA Boostrap sea ice concentration
#' product.
#'
#' @docType data
#' @format array
#' @keywords datasets
#' @usage data(clim_9_2008)
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center.
#'             doi: \url{https://doi.org/10.5067/7Q8HCCWS4I0R}
"clim_9_2008"


#' Proportion of total area by grid box
#'
#' Matrix of dimension longitude by latitude. Elements give the proportion of
#' the total area (within the seas of the Arctic) in each grid box. The sum
#' of all elements is 1.
#'
#' @docType data
#' @format array
#' @keywords datasets
#' @usage data(prop_area)
"prop_area"







