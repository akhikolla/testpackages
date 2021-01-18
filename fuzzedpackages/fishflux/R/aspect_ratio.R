#' A function to find aspect ratio
#'
#' A function to find aspect ratio of a species on either species or genus
#' level using rfishbase. It returns a data frame containing the aspect ratio
#' and the level at which the aspect ratio was found (species or genus).
#'
#' @param sp A character value containing the species name
#' 
#' @keywords fish aspect-ratio fishbase
#' 
#' @importFrom rfishbase morphometrics species_list
#' @importFrom dplyr select summarise group_by select
#'
#' @examples
#' \donttest{
#' library(fishflux)
#' library(plyr)
#' aspect_ratio("Lutjanus griseus")
#' ldply(lapply(c("Chlorurus spilurus","Zebrasoma scopas"), aspect_ratio))
#' }
#' 
#' @export
aspect_ratio <- function(sp) {
    check_name_fishbase(sp)
    ma <- morphometrics(sp)
    if (length(ma) == 0) {
        genus <- strsplit(sp, " ")[[1]][1]
        gn    <- species_list(Genus = genus)
        ma    <- morphometrics(gn)
        asp   <- select(ma, Species, AspectRatio, SL, TL)
        asp   <- summarise(group_by(asp, Species), AspectRatio = mean(AspectRatio))
        asp   <- mean(asp$AspectRatio, na.rm = TRUE)
        level <- "genus"
    } else {
        asp   <- select(ma, Species, AspectRatio, SL, TL)
        asp   <- mean(asp$AspectRatio, na.rm = TRUE)
        level <- "species"
    }
    data.frame(species = sp,
               aspect_ratio = asp,
               level = level)
}
