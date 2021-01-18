#' A function to find growth parameters on fishbase
#'
#' A function to find growth parameters of a species using rfishbase.
#' It returns a data frame containing K, t0 and Linf, the source.
#' This function is useful to see what is available on fishbase.
#' Nevertheless, we strongly recommend to check the source and only use otolith based studies.
#'
#' @param sp A charachter value containing the species name
#' @param otolith A logical value. If TRUE, only results from otolith analysis are returned. If false, all growth studies will be returned.
#' 
#' @keywords fish trophic level fishbase
#' 
#' @importFrom rfishbase popgrowth
#' @importFrom dplyr select filter
#' 
#' @examples
#' \donttest{library(fishflux)
#' growth_params("Lutjanus griseus")}
#' 
#' @export
growth_params <- function(sp, otolith = TRUE) {
  if (length(suppressMessages(name_errors(sp))) > 0) {
    stop("Species name is incorrect")
  }
 pop <- popgrowth(sp)
 if (length(pop) == 0) {
   growth <- data.frame(species = NA, Locality = NA, k = NA,
                        Linf = NA, t0 = NA, method = NA, comments = NA)
 } else {
   growth <- select(pop, species = Species, Locality,
                    k = K, Linf = Loo, t0 = to,
                    method = Data,
                    comments = Comment)
    if (otolith) {
      growth <- filter(growth, method %in%
                  c("annuli on otoliths", "annuli on many otoliths"))
    }
    if (nrow(growth) < 1) {
      stop("No otolith based parameters available in fishbase.
           Try otolith = FALSE or search literature.")
    }
  }
 growth
}
