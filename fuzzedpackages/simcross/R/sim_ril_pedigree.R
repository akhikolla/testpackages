# sim_ril_pedigree
#
#' Generate a ril pedigree
#'
#' Generate a pedigree for multi-way recombinant inbred lines (a table
#' of individual, mom, dad, sex)
#'
#' @param ngen Number of generations of inbreeding
#' @param selfing If TRUE, use selfing
#' @param parents Vector of the parents' IDs. Should be integers, and
#'     length must be a power of 2 (i.e., 2, 4, 8, ...)
#' @param firstind Positive integer to assign to the first child. Must
#'     be greater than `max(parents)`.
#'
#'
#' @return A data frame with five columns: individual ID, mother ID,
#' father ID, sex, and generation.  Founders have `0` for mother
#' and father ID. Sex is coded 0 for female and 1 for male.
#'
#' @export
#' @seealso [sim_from_pedigree()],
#' [sim_ail_pedigree()], [sim_do_pedigree()],
#' [sim_4way_pedigree()]
#' @keywords datagen
#'
#' @examples
#' tab <- sim_ril_pedigree(7)
sim_ril_pedigree <-
    function(ngen=20, selfing=FALSE, parents=1:2, firstind=max(parents)+1)
{
    if(!is.numeric(parents))
        stop("parents should be vector of integers")

    nparents <- length(parents)
    if(nparents %% 2 != 0 || nparents<2)
        stop("number of parents must = 2^k for some positive k")

    if(firstind <= max(parents))
        stop("firstind should be > max(parents); otherwise there will be ID conflicts")

    # first generation of pedigree
    ped <- data.frame(id=parents,
                      mom=0,
                      dad=0,
                      sex=c(0,1),
                      gen=0)
    nextind <- firstind
    generation <- 0

    # parents -> matrix of couples
    parents <- matrix(parents, ncol=2, byrow=TRUE)

    # pedigree up to the last outcross
    while(length(parents) > 2) {
        generation <- generation + 1
        offspring <- seq(nextind, by=1, length=nrow(parents))
        nextind <- max(offspring) + 1
        ped <- rbind(ped,
                     data.frame(id=offspring,
                                mom=parents[,1],
                                dad=parents[,2],
                                sex=c(0,1),
                                gen=generation) )

        parents <- matrix(offspring, ncol=2, byrow=TRUE)
    }

    # parents now has length 2
    generation <- generation + 1

    if(selfing) {
        # one last cross
        ped <- rbind(ped,
                     data.frame(id=nextind,
                                mom=parents[,1],
                                dad=parents[,2],
                                sex=0,
                                gen=generation))
        # then start selfing
        nextind <- nextind + 1
        for(g in generation + 1:ngen) {
            ped <- rbind(ped,
                         data.frame(id=nextind,
                                    mom=nextind-1,
                                    dad=nextind-1,
                                    sex=0,
                                    gen=g))
            nextind <- nextind + 1
        }
    }
    else {
        # one last cross
        ped <- rbind(ped,
                     data.frame(id=c(nextind,nextind+1),
                                mom=parents[,1],
                                dad=parents[,2],
                                sex=c(0,1),
                                gen=generation))
        # then start selfing
        nextind <- nextind + 2
        for(g in generation + 1:ngen) {
            ped <- rbind(ped,
                         data.frame(id=c(nextind,nextind+1),
                                    mom=nextind-2,
                                    dad=nextind-1,
                                    sex=c(0,1),
                                    gen=g))
            nextind <- nextind + 2
        }

    }
    ped
}
