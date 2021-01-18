# General Clases ----

#' Spectral info
#' 
#' @description General class that holds basic information about an image. Usually used when working with mosaics, as in this case the image is not read/loaded to R until is processed (e.g. by running a mosaic_sam process). Currently, the only function that returns an object of this class is he \code{\link{mosaic_info}} function. The object can by used later to read/load to R a single mosaic chuck using the \code{\link{mosaic_chunk}} function, or to process the whole mosaic image by calling the \code{\link{mosaic_sam}} function.
#' 
#' The class \code{\link[=Tile-class]{Tile}} contains the "Spectral_info" class. By this means, an object of class \code{\link[=Tile-class]{Tile}} holds the same information about the reading, and an extra slot called "Spectra" that holds the actual readings. 
#'
#' @slot file character.  The full path to the .dmt or .bsp file.
#' @slot date POSIXct.  The date in which the image (mosaic or tile) was taken.
#' @slot fpasize numeric. The size of the focal plane array used when measuring.
#' @slot wavenumbers numeric. The wavenumbers at which the tarject image (mosaic or tile) was taken.
#' @slot path character.  The path to the folder that holds all the uFTIR Microscope output files.
#'
#' @return
#' An S4 object of class "Spectral_info"
#' @seealso 
#' \code{\link[=Tile-class]{Tile}}
#' 
#' @export
#' @examples
#' NULL
setClass("SpectralInfo",
         slots = c(file = "character",
                   date = "POSIXct",
                   fpasize = "numeric",
                   wavenumbers = "numeric",
                   path = "character"
                   )
)

# Single tile Clases ----

#' Tile
#' 
#' @description Class to hold the reading of a single tile. It can also store a single tile of a mosaic, to which we will refer to as chunks. It is an extended class for \code{\link[=SpectralInfo-class]{SpectralInfo}}, adding an extra slot to hold the uFTIR Microscope readings -called Spectra. The functions \code{\link{tile_read}} and \code{\link{mosaic_chunk}} return an object of this class. The class has a defined method to plot, which can be called either by the generic \code{\link[graphics]{plot}} or by \code{\link{plot_tile}}. You can use an object of this class in the following functions:
#' \describe{
#'     \item{\code{\link{plot_tile}}}{To plot a heat map of the reading. Available too through the generic \code{\link[graphics]{plot}}.}
#'     \item{\code{\link{tile_base_corr}}}{A preprocessing function that replace negative values with zeros.}
#'     \item{\code{\link{wavealign}}}{A resampling method that allows to match the wavenumbers of a given tile/chunk with the wavenumbers of a \code{\link[=SpectralReference-class]{SpectralReference}}.}
#'     \item{\code{\link{tile_sam}}}{When a Tile object is linked to a \code{\link[=SpectralReference-class]{SpectralReference}} by a former call to \code{\link{wavealign}}, this function performs a spectral angle mapper match between the Tile and the SpectralReference.}
#'     \item{\code{\link{get_profile_tile}}}{To get the pixels spectra that matched a given polymer/cluster.}
#' }
#'
#' @slot Spectra array. The aquired spectra (.dat or dmd file)
#'
#' @return
#' An S4 object of class "Tile".
#' 
#' @export
#' @examples
#' NULL
setClass("Tile",
         slots = c(Spectra = "array"),
         contains = "SpectralInfo"
         )

# Here we check that the wavenumbers match with the third dimension of the array.
setValidity("Tile",
            function(object){
              isTRUE(all.equal(length(object@wavenumbers),
                               dim(object@Spectra)[length(dim(object@Spectra))])
                     )
            }
)

# Other Clasess ----

#' Spectral Reference Class
#'
#' @description I made this class to hold spectral libraries. This libraries can be used to match a know spectrum with a pixel reading. \code{\link[=SpectralReference-class]{SpectralReference}} objects are used by several functions in the package, at all analysis stages (i.e. preprocessing, processing, and summarizing/plotting). You can see an example in the data \code{\link{primpke}} (you can load it by calling data("primpke")). The format includes character/integer vectors to hold the names of the polymers, cluster IDs, and cluster names, all of which should correspond to one another. It also hold a matrix with the spectra of each polymer, with rows matching the vectors provided and a fourth vector holding the wavelengts (that should be equal to the matrix cols).
#' 
#' The example provided in the package was taken form Primpke, S., Wirth, M., Lorenz, C., Gerdts, G. 2018. Reference database design for the automated analysis of microplastic samples based on Fourier transform infrared (FTIR) spectroscopy. Analytical and Bioanalytical Chemistry 410: 5131-5141. You might access to the article here \url{https://doi.org/10.1007/s00216-018-1156-x}.
#'  
#' If you have reference data that was taken considering a different range of wavenumbers, you should resample it first for the wavenumbers to match. 
#'
#' @slot substances A character vector identifying by name the polymers whose spectrum is provided in Spectra. The order of the elements should be consistent with clusterlist, clusternames and Spectral rows. Along the same lines, the object expects length(substances) == nrow(Spectra).
#' @slot clusterlist If the polymers are aggregated in clusters this slot should hold an integer vector with the correspondent cluster for each polymer individuated in substances. If you don't want to use clusters or you don't have clusters for your polymers/library you can place a sequence of integers from 1:nrow(Spectra) (e.g. by calling seq_along(substances)) in this slot.
#' @slot clusternames Character vector holding reference names for each cluster, sorted according to clusterlist. In other words, if you want to name cluster 1L "polystyrene" the first position in clusternames should be equal to "polystyrene". If you don't want to use clusters or you don't have clusters for your polymers/library you can duplicate substances. The program expects that the length of the clusternames vector equals the length of unique elements in clusterlist.
#' @slot Spectra A matrix holding row-wise the spectrum of each substance. Each row should correspond to the spectrum of one substance. In the same lines, columns should hold recorded measures for each wavenumber. The program expects that the length of the substance vector equals the number of rows of Spectra and that the number of colums of Spectra equals the length of the wavenumbers vector.
#' @slot wavenumbers A numeric vector with the wavenumbers.
#'
#' @return 
#' A S4 object of class "SpectralReference".
#' @export
#' @examples
#' data("primpke")
setClass("SpectralReference",
         slots = c(substances = "character",
                   clusterlist = "integer",
                   clusternames = "character",
                   Spectra = "matrix",
                   wavenumbers = "numeric")
         )

setValidity("SpectralReference",
            function(object){
              var <- FALSE
              errors <- c()
              if(length(dim(object@Spectra)) != 2){
                errors <- c(errors, "Spectra should be a 2x2 matrix\n")
              }
              if(!isTRUE(all.equal(length(object@wavenumbers), ncol(object@Spectra)))){
                errors <- c(errors, "Column length should be equal to the length of wavenumbers\n")
              }
              if(!isTRUE(all.equal(length(object@clusterlist), length(object@substances)))){
                errors <- c(errors, "The cluster list should correspond with the substances list\n")
              }
              if(length(errors)==0) var <- TRUE
              var
            }
)

#' Spectral Pack
#'
#' @description A call to the function \code{\link{wavealign}} (specifically a call to \code{\link{TileRead.wavealign}}) returns an object of class \code{\link[=SpectralPack-class]{SpectralPack}}. This class reunites the uFTIR Microscope readings loaded to R by a call to \code{\link{tile_read}} with a tarjet spectral library of class \code{\link[=SpectralReference-class]{SpectralReference}}. Since this object is an output of a call to \code{\link{wavealign}} both the \code{\link[=Tile-class]{Tile}} object and \code{\link[=SpectralReference-class]{SpectralReference}} have congruent wavenumbers. See \code{\link{wavealign}} how the resampling is done when calling different methods. A \code{\link[=SpectralPack-class]{SpectralPack}} object is used as imput when calling \code{\link{tile_sam}} to run the Spectral Angle Mapping algorithm.
#'
#' @slot Readings an object of class \code{\link[=Tile-class]{Tile}}.
#' @slot Reference an object of class \code{\link[=SpectralReference-class]{SpectralReference}}.
#'
#' @return
#' An S4 object of class "SpectralPack" ready for spectral angle mapping.
#'
#' @export
#' @examples
#' NULL
setClass("SpectralPack",
         slots = c(Readings = "Tile",
                   Reference = "SpectralReference")
)

#' Sam Prediction
#'
#' @description The class holds the resuls of the Spectral Angle Mapping algorithm and it originates by a call to \code{\link{tile_sam}} or to \code{\link{mosaic_compose}}. It has three slots with congruent cubes. Each cube has its first two dimensions equal to the FPA size of the uFTIR Microscope image processed and the third dimension equal to the number of substances in the \code{\link[=SpectralReference-class]{SpectralReference}} passed when calling \code{\link{tile_sam}} or \code{\link{mosaic_sam}}. This class have methods to plot, summarize, and extract the readings of pixel matching a given SpectralReference substance. You can use the generics \code{\link[graphics]{plot}}, \code{\link[base]{summary}}, and \code{\link{get_profile}} but beware, you might have to provide extra arguments! (so, please don't see only the generic but the method for this class, see the section see also)
#' 
#' The \code{\link[=SAM-class]{sam}} object can be postprocessed by calling \code{\link{smooth_sam}} and/or clipped to a given extent by calling \code{\link{clipper}} (currently only circular extents are supported).
#'
#' @slot raw_sam a cube -array-. The numeric score for each \code{\link[=SpectralReference-class]{SpectralReference}} substance passed by the call (slices), for each pixel of the FPA array of the \code{\link[=Tile-class]{Tile}} (rows and cols).
#' @slot substances a cube -array-. The slices of the cube hold the substances (as number) of the \code{\link[=SpectralReference-class]{SpectralReference}} object passed when calling \code{\link{tile_sam}} or \code{\link{mosaic_sam}}. The substance number are sorted in a decreacing order. In this way, the first slide is the best match, the second the second-best, and so on according to the Spectral Angle Mapping algorithm.
#' @slot clusters a cube -array-. Equivalent to the cube in the substances slot but instead of a substance number it has per ech slice the cluster number that corresponds to the clusterlist slot of the \code{\link[=SpectralReference-class]{SpectralReference}} object.
#'
#' @return
#' An S4 object of class SAM
#' @export
#' 
#' @seealso 
#' \code{\link{plot_tile}} \code{\link{summary_sam}} \code{\link{get_profile_tile}} \code{\link{get_profile_sinfo}}
#' @examples
#' NULL
setClass("SAM",
         slots = c(raw_sam = "array",
                   substances = "array",
                   clusters = "array"
         )
)

#' Smooth
#' 
#' @description An object of this class is returned by a call to \code{\link{smooth_sam}}. It holds a smoother version of a SAM object, keeping only the clusters slot and slicing it to a user-given number of slices. This class have methods to plot, summarize, and extract the readings of pixel matching a given SpectralReference substance. You can use the generics \code{\link[graphics]{plot}}, \code{\link[base]{summary}}, and \code{\link{get_profile}} but beware, you might have to provide extra arguments! (so, please don't see only the generic but the method for this class, see the section see also).
#'
#' @slot smooth a cube -array- (even if only has 1 slice) that holds a smooth version of a SAM clusters slot.
#'
#' @return
#' An S4 object of class "Smooth.
#' @export
#' 
#' @seealso 
#' \code{\link{plot_tile}} \code{\link{summary_sam}} \code{\link{get_profile_tile}} \code{\link{get_profile_sinfo}} 
#' @examples
#' NULL
setClass("Smooth",
         slots = c(smooth = "array")
)

#' Class holding info for the clipper
#'
#' @description The filters we used to support our samples under the uFTIR microscope are round. This means that a cube might not be the best representation of our data, as they have noisy borders that might be missinterpreted by the program. Then, I implemented a \code{\link{clipper}} function to clip objects of class \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}} (even an S3 matrix) before summarizing them. Although this can be done step by step, you can also create a clip-mask for the \code{\link{summary_sam}} function to clip them for you. To crate this clip-mask, you can call the function \code{\link{toClip}} which returns an object of this class.
#' 
#' The class has yet another use. Since the samples are placed under the microscope by hand, the cropping area is not always the same and (usually) it has to be adjusted. To have a visual aid, you can use the function \code{\link{toClip}}, to then plot a \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}} and overlay the croppling circle by calling \code{\link[graphics]{polygon}} and the xycoords slot of an object of this class (returned by \code{\link{toClip}}). Using that process you can test manually different center points and radius for the clipping circle. 
#' 
#' @slot xycoords polygon (circle) coodinates to plot. matrix with two columns (x,y).
#' @slot rad cicle radius. integer.
#' @slot centre circle centre. Numeric vector of two elements (x,y). 
#'
#' @return
#' An S4 object of class "clipmask" to be used (internally) by the \code{\link{clipper}} function.
#' @export
#' 
#' @seealso 
#' \code{\link{toClip}} \code{\link{clipper}}
#' @examples
#' NULL
setClass("clipmask", 
         slots = c(
           xycoords = "matrix",
           rad = "integer",
           centre = "numeric")
)

# a toy class to allow ploting the out of the clipper function
setOldClass("clipper")
