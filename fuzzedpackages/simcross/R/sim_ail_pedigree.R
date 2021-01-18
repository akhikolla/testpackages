# sim_ail_pedigree
#
#' Simulate AIL pedigree
#'
#' Simulate a pedigree for advanced intercross lines (a table of
#' individual, mom, dad, sex)
#'
#' @param ngen Number of generations of outbreeding
#' @param npairs Number of breeding pairs at each generation
#' @param nkids_per Number of offspring per pair for the last
#' generation
#' @param design How to choose crosses: either random but avoiding
#' siblings, or completely at random
#'
#' @return A data frame with five columns: individual ID, mother ID,
#' father ID, sex, and generation.  Founders have `0` for mother
#' and father ID. Sex is coded 0 for female and 1 for male.
#'
#' @details Advanced intercross lines (AIL) are generated from a pair of inbred lines.
#' We cross them and then cross the F1 to generate `npair` breeding pairs.
#' The subsequent `ngen` outbreeding generations then proceed by
#' crossing a male and female from the preceding generation (mated
#' completely at random, with `design="random"`, or avoiding
#' siblings, with `design="nosib"`). Each breeding pair gives a
#' single female and a single male to the next generation, except at
#' the last generation `nkids_per` offspring are mated, in equal
#' numbers male and female. (If `nkids_per` is an odd number, the
#' number of males and females in each sibship will differ by one,
#' alternating between sibships, with one additional female and then
#' one additional male.
#'
#' @export
#' @keywords datagen
#' @seealso [sim_from_pedigree()],
#' [sim_ril_pedigree()], [sim_do_pedigree()],
#' [sim_4way_pedigree()]
#'
#' @examples
#' tab <- sim_ail_pedigree(12, 30)
sim_ail_pedigree <-
    function(ngen=12, npairs=30, nkids_per=5, design=c("nosib", "random"))
{
    design <- match.arg(design)

    # initial generations
    id <- c(1,2,3,4, (1:(npairs*2)+4))
    mom <- c(0,0,1,1,rep(3, npairs*2))
    dad <- c(0,0,2,2,rep(4, npairs*2))
    sex <- c(0,1,0,1,rep(c(0,1), npairs))
    gen <- c(0,0,1,1,rep(2, npairs*2))

    if(ngen <= 2)
        stop("No. generations should be > 2")

    cur_nind <- length(id)
    moms <- seq(5, length(id), by=2)
    dads <- seq(6, length(id), by=2)
    for(i in 3:ngen) {
        if(i > 3) {
            dads <- sample(dads)
            while(design=="nosib" && any(dads - moms == 1)) { # sample until no sibs
                dads <- sample(dads)
            }
        }

        if(i < ngen) {
            kids <- 1:(npairs*2)+max(id)

            id <- c(id, kids)
            mom <- c(mom, rep(moms, each=2))
            dad <- c(dad, rep(dads, each=2))
            sex <- c(sex, rep(c(0,1), npairs))
            gen <- c(gen, rep(i, npairs*2))

            moms <- kids[seq(1, length(kids), by=2)]
            dads <- kids[seq(2, length(kids), by=2)]
        }
        else { # last generation: expand
            kids <- 1:(npairs*nkids_per)+max(id)

            id <- c(id, kids)
            mom <- c(mom, rep(moms, each=nkids_per))
            dad <- c(dad, rep(dads, each=nkids_per))
            sex <- c(sex, rep(c(0,1), ceiling(npairs*nkids_per/2))[1:(npairs*nkids_per)])
            gen <- c(gen, rep(i, npairs*nkids_per))
        }
    }

    data.frame(id=id, mom=mom, dad=dad, sex=sex, gen=gen)
}
