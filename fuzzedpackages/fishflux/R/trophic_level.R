#' A function to find trophic level
#'
#' A function to find trophic level of a species on either species or genus
#' level using rfishbase. It returns a data frame containing the trophic level
#' and the level at which the trophic level was found (species or genus).
#'
#' @param sp A character value containing the species name
#' 
#' @keywords fish trophic level fishbase
#' 
#' @importFrom rfishbase ecology species species_list
#' @importFrom dplyr select
#' 
#' @examples
#' \donttest{
#' library(fishflux)
#' library(plyr)
#' trophic_level("Lutjanus griseus")
#' ldply(lapply(c("Chlorurus spilurus","Zebrasoma scopas"), trophic_level))}
#' 
#' @export
trophic_level <- function(sp) {
  check_name_fishbase(sp)
  ecogn <- ecology(sp)
  if (length(ecogn) == 0) {
    genus <- species(sp)$Genus
    gn    <- species_list(Genus = genus)
    ecogn <- ecology(gn)
    level <- "genus"
  } else {
    level <- "species"
  }
  troph <- select(ecogn, Species, DietTroph, FoodTroph)
  troph$troph <- rowMeans(troph[, 2:3], na.rm = TRUE)
  troph <- mean(troph$troph, na.rm = TRUE)
  data.frame(species = sp,
             trophic_level = troph,
             level = level)
}
