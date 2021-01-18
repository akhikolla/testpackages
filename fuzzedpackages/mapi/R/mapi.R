#' @title MAPI, general presentation
#' @name mapi
#' @aliases mapi-package
#' @import data.table
#' @import sf
#' @import stats
#' @import parallel
#' @import pbapply
## usethis namespace: start
#' @useDynLib mapi
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
#' @docType package
#' 
#' @description 
#' MAPI is an exploratory method providing graphical representations of the spatial variation of 
#' pairwise metrics (eg. distance, similarity coefficient, ...) computed between georeferenced samples.
#' 
#' \subsection{Principle}{
#' As schematically illustrated Figure 1, MAPI relies on spatial joins between a hexagonal grid and a network
#' of georeferenced samples connected by ellipses, i.e. polygons with 32 segments approaching an elliptical shape.
#'
#' The shape of the ellipses can be controlled through the eccentricity value and the sample locations can be 
#' "blurred" by applying an error circle of a given radius on the geographic coordinates.
#' Each elliptical polygon is associated to 1) the value of the pairwise metric computed between the samples 
#' it connects and 2) a weight corresponding to the inverse of its area (i.e. larger ellipses have lower weights).
#' 
#' Each cell of the grid receives the weighted mean of the pairwise metric values associated to the ellipses intersecting the cell.
#' }
#' \figure{fig3D.png}\cr
#' \emph{Figure 1: Schematic principle of the MAPI method from Piry et al. 2016.}
#' 
#' 
#' \subsection{Input data}{
#'   The analysis requires two tables (data.frame or data.table):
#' 
#'   1) Information on samples: table with three mandatory columns and column names: 'ind' (sample name), 'x' and 'y' (projected coordinates). 
#'   An optional column 'errRad' (radius of error circle on sample coordinates) can be provided.
#' 
#'   MAPI requires cartesian coordinates (ie. projected, such as UTM or Lambert) NOT (yet?) angular coordinates (eg. latitude/longitude).
#'   The package \pkg{sf} provides the \code{st_transform} function for coordinates transformation and projection.
#'   GIS software such as QGis can also help with datum transformation.
#'   
#'   Example of 'samples' data:
#'   \tabular{llrrr}{
#'      \tab ind   \tab x      \tab y      \tab errRad \cr
#'   1  \tab 2_42  \tab 12000  \tab 5000   \tab     10 \cr
#'   2  \tab 2_47  \tab 17000  \tab 5000   \tab     10 \cr
#'   3  \tab 1_82  \tab  2000  \tab 9000   \tab     10 \cr
#'   4  \tab 2_100 \tab 20000  \tab 10000  \tab     10 \cr
#'   5  \tab 2_87  \tab 17000  \tab 9000   \tab     10 \cr
#'   6  \tab 1_11  \tab 1000   \tab 2000   \tab     10 \cr
#'   ...  \tab ...  \tab ...   \tab ...    \tab   ...  \cr
#'   }
#'  
#'   2) Values of the pairwise metric computed between samples provided, either, as a complete matrix with the same number of columns and rows 
#'   (column and row names must match the sample names provided in the 'samples' data) or as a table with three mandatory columns and 
#'   column names: 'ind1', 'ind2' (sample names) and 'value' (pairwise metric values).
#'   
#'   Example of 'metric' data:
#'   \tabular{lllr}{
#'       \tab ind1  \tab ind2 \tab   value \cr
#'   1   \tab  1_1  \tab 1_2 \tab 0.055556 \cr
#'   2   \tab  1_1  \tab 1_3 \tab 0.020833 \cr
#'   3   \tab  1_1  \tab 1_4 \tab 0.125000 \cr
#'   4   \tab  1_1  \tab 1_5 \tab 0.125000 \cr
#'   5   \tab  1_1  \tab 1_6 \tab 0.020833 \cr
#'   6   \tab  1_1  \tab 1_7 \tab 0.090278 \cr
#'   ... \tab ...   \tab ... \tab ...      \cr
#'   }
#' 
#' 
#' }
#' 
#' \subsection{Try it}{
#'  Using the test dataset ('samples' and 'metric') included in the package, let's run an (almost) automatic MAPI analysis
#' 
#'  Test data result from population genetic simulations in which two panmictic populations are
#'  separated by a barrier to dispersal. As we use dummy coordinates, there is no appropriate crs, so we just use 'crs=3857' 
#'  (a pseudo-mercator projection). Of course, in this case, sample position on earth is meaningless. 
#'  For a real dataset, 'crs' must be the EPSG code of the projection of your cartesian coordinates.
#'  
#'  \preformatted{
#'  # Load the package
#'  library(mapi)
#'  
#'  # Load 'samples' data
#'  data("samples")
#'  
#'  # Load 'metric' data. For our simulated data set the parwise metric 
#'  # computed between samples is the individual genetic distance â of Rousset (2000).
#'  data("metric")
#'  
#'  # Run MAPI the lazy way (automatic) with 1000 permutations 
#'  # for detection of significant (dis)continuous areas.
#'  # As crs must be set, we go with crs=3857 even if we use dummy coordinates.
#'  # Of course, this have no geographical meaning.
#'  # As we have a regular sampling, we use beta=0.5
#'  my.results <- MAPI_RunAuto(samples, metric, crs=3857, beta=0.5, nbPermuts=1000)
#'  
#'  # Get significant areas with a FDR control at alpha=0.05 (5\%, by default)
#'  my.tails <- MAPI_Tails(my.results, alpha=0.05)
#'  
#'  # Look at the result Figure 2.
#'  MAPI_Plot2(my.results, tails=my.tails)
#'  }
#'  Spatial variation of the genetic distance is represented with a color scale from dark brown (lowest values) to dark blue (higher value). The central blue area identified as a significant area of discontinuity corresponds
#'   to the position of the simulated barrier. Note that due to the permutation procedure, delineation of the significant areas may vary slightly among runs.
#'  
#' }
#' \figure{mapiPlotOutput.png}\cr
#' \emph{Figure 2: MAPI graphical Output produced using the MAPI_Plot2 function.}
#' 
#' 
#' 
#' \subsection{To go deeper}{
#'  \code{\link{MAPI_RunAuto}} is a wrapper which calls three other functions: \code{\link{MAPI_CheckData}}, \code{\link{MAPI_GridAuto}} and \code{\link{MAPI_RunOnGrid}}.
#' 
#'  \code{\link{MAPI_GridAuto}} is itself another wrapper around \code{\link{MAPI_EstimateHalfwidth}} and \code{\link{MAPI_GridHexagonal}}.
#' 
#'  Typically, a "manual" MAPI analysis will involve the following ordered steps:
#'  \enumerate{
#'    \item{ \code{\link{MAPI_CheckData}} }
#'    \item{ \code{\link{MAPI_EstimateHalfwidth}} }
#'    \item{ \code{\link{MAPI_GridHexagonal}} }
#'    \item{ \code{\link{MAPI_RunOnGrid}} }
#'    \item{ \code{\link{MAPI_Tails}} }
#'    \item{ \code{\link{MAPI_Plot2}} }
#'  }
#'  
#'  Within this general framework, you may, for example:
#'  \itemize{
#'  \item{set your own value for 'halfwidth' (ignore step 2)}
#'  \item{use your own grid, or reuse one from another run (ignore steps 2 & 3)}
#'  \item{tweak some MAPI parameters (such as dMin or dMax for filtering on geographic distances between samples)}
#'  \item{discard poorly supported cells prior detecting significant areas of (dis)continuity  (parameter \code{minQ}) and/or change significance level (parameter \code{alpha} in \code{\link{MAPI_Tails}})}
#'  \item{build your MAPI maps with a GIS software (ignore step 6). See 'Export results' section below}
#'  }
#' 
#' }
#' 
#' 
#' 
#' \subsection{Export results}{
#' Output tables (weighted mean of the pairwise metric within cell and polygons delineating significant areas of (dis)continuity) are spatial objects built using the package \pkg{sf}. 
#' Refer to \pkg{sf} documentation to export MAPI results in various format.
#' Below is an example of how MAPI results can be exported as ESRI Shapefiles:\cr
#' \preformatted{
#' library(sf)
#' # Export results for our test dataset
#' st_write(my.results, dsn=".", layer="myFirstMapiResult", 
#'    driver="ESRI Shapefile", update=TRUE, delete_layer=TRUE)
#' st_write(my.tails, dsn=".", layer="myFirstMapiResultTails", 
#'    driver="ESRI Shapefile", update=TRUE, delete_layer=TRUE)
#' }
#' You may now open these files \file{myFirstMapiResult.shp} and \file{myFirstMapiResultTails.shp} in a GIS software such as QGis and customize the layout.
#' 
#' Overlaying MAPI results with landscape layouts can help in analyzing the relationship between environmental features and spatial genetic patterns (eg. Piry & al., 2016; Piry & al., 2018).
#' }
#' 
#' 
#' @references 
#' \subsection{Description of the MAPI method}{\cr
#'   Piry, S., Chapuis, M.-P., Gauffre, B., Papaïx, J., Cruaud, A., and Berthier, K. \bold{(2016)}. 
#'   Mapping Averaged Pairwise Information (MAPI): a new exploratory tool to uncover spatial structure. 
#'   \emph{Methods in Ecology and Evolution} \bold{7}:(12), 1463–1475. 
#'   \href{https://doi.org/10.1111/2041-210X.12616}{doi: 10.1111/2041-210X.12616}
#' }
#' \subsection{Applications of MAPI in Landscape Genetics}{\itemize{
#'   \item{Piry, S., Berthier, K., Streiff, R., Cros-Arteil, S., Tatin, L., Foucart, A., Bröder, L., Hochkirch, A., and Chapuis, M.-P. \bold{(2018)}. 
#'   Fine-scale interactions between habitat quality and genetic variation suggest an impact of grazing on the critically endangered Crau Plain grasshopper (Pamphagidae: \emph{Prionotropis rhodanica}). 
#'   \emph{Journal of Orthoptera Research} \bold{27}, 61–73.
#'   \href{https://doi.org/10.3897/jor.27.15036}{doi: 10.3897/jor.27.15036}}
#'   \item{Dellicour S, Prunier JG, Piry S, \& al. \bold{(2019)}
#'   Landscape genetic analyses of \emph{Cervus elaphus} and \emph{Sus scrofa}: comparative study and analytical developments. 
#'   \emph{Heredity}.
#'   \href{https://doi.org/10.1038/s41437-019-0183-5}{doi: 10.1038/s41437-019-0183-5}}
#' }}
#'
NULL


