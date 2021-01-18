# sim_4way_pedigree
#
#' Simulate pedigree for 4-way intercross
#'
#' Simulate a 4-way cross, among four inbred lines (a table of
#' individual, mom, dad, sex)
#'
#' @param ngen Number of intercross generations (1 or 2)
#' @param nsibs Vector with number of siblings in the sibships in the
#' last generation.
#'
#' @return A data frame with five columns: individual ID, mother ID,
#' father ID, sex, and generation.  Founders have `0` for mother
#' and father ID. Sex is coded 0 for female and 1 for male.
#'
#' @details We start with a set of 4 individuals (representing four
#' inbred lines), and make a pair of crosses to generate a pair of
#' heterozygous individuals.  These are then crosses to generate a set
#' of F1 individuals. If `ngen==1`, we stop there, with
#' `sum(nsibs)` individuals in this last generation.  If
#' `gen==2`, we generate `length(nsibs)` male/female pairs
#' of F1 offspring; these are intercrossed to generate a set of
#' sibships, with lengths defined by the values in `nsibs`.
#' Individuals in the last generation are alternating female/male.
#'
#' @export
#' @keywords datagen
#' @seealso [sim_from_pedigree()],
#' [sim_ril_pedigree()], [sim_do_pedigree()],
#' [sim_ail_pedigree()]
#'
#' @examples
#' # 100 F1s between heterozygous parents
#' tab <- sim_4way_pedigree(1, 100)
#' # could also do this
#' tab2 <- sim_4way_pedigree(1, rep(10, 10))
#'
#' # 120 F2s in 10 sibships each of size 12
#' tab3 <- sim_4way_pedigree(ngen=2, rep(12, 10))
sim_4way_pedigree <-
    function(ngen=1, nsibs=100)
{
    if(ngen != 1 && ngen != 2)
        stop("ngen must be 1 or 2")
    if(any(is.na(nsibs) | nsibs <= 0))
        stop("nsibs can't have missing or non-positive values")

    # initial generations
    id <-  c(1,2,3,4, 5,6)
    mom <- c(0,0,0,0, 1,3)
    dad <- c(0,0,0,0, 2,4)
    sex <- c(0,1,0,1, 0,1)
    gen <- c(0,0,0,0, 1,1)

    if(ngen==1) nf1 <- sum(nsibs)
    if(ngen==2) nf1 <- 2*length(nsibs)

    id <- c(id, (1:nf1)+6)
    mom <- c(mom, rep(5, nf1))
    dad <- c(dad, rep(6, nf1))
    sex <- c(sex, rep(c(0,1), ceiling(nf1/2))[1:nf1])
    gen <- c(gen, rep(2, nf1))

    if(ngen==2) {
        totsibs <- sum(nsibs)

        prev <- length(id)
        for(i in seq(along=nsibs)) {
            id <- c(id, prev+(1:nsibs[i]))
            mom <- c(mom, rep(i*2+5, nsibs[i]))
            dad <- c(dad, rep(i*2+6, nsibs[i]))
            prev <- length(id)
        }
        sex <- c(sex, rep(c(0,1), ceiling(totsibs/2))[1:totsibs])
        gen <- c(gen, rep(3, totsibs))
    }

    data.frame(id=id, mom=mom, dad=dad, sex=sex, gen=gen)
}
