#' Simulate genotypes for pedigree for multiple chromosomes
#'
#' Simulate genotypes along all chromosomes for a pedigree. This is a
#' wrap up of sim_from_pedigree.
#'
#' @inheritParams sim_from_pedigree
#' @param map marker locations, a list with elements for each
#' chromosome
#'
#' @return A list with each component being the result from
#' `sim_from_pedigree`, of length same as `map`.
#'
#' @export
#' @keywords datagen
#' @seealso [check_pedigree()],
#' [sim_ril_pedigree()], [sim_ail_pedigree()]
#' [sim_from_pedigree()]
#'
#' @examples
#' library(qtl)
#' # marker map
#' map <- sim.map(len=rep(100, 19), n.mar=10, include.x=FALSE)
#' # simulate AIL pedigree
#' tab <- sim_ail_pedigree(12, 30)
#' # simulate data from that pedigree
#' dat <- sim_from_pedigree_allchr(tab, map)
#'
sim_from_pedigree_allchr <- function(pedigree, map, m=10, p=0, obligate_chiasma=FALSE){
  chr.len <- sapply(map, max)
  is.xchr <- sapply(map, class) == "X"
  dat <- NULL
  for(chr in 1:length(map))
      dat[[chr]] <- sim_from_pedigree(pedigree, L=chr.len[chr], xchr=is.xchr[chr],
                                      m=m, p=p, obligate_chiasma=obligate_chiasma)
  names(dat) <- names(map)
  return(dat)
}
