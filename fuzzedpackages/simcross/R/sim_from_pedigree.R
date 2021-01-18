# sim_from_pedigree
#
#' Simulate genotypes for pedigree
#'
#' Simulate genotypes along one chromosome for a pedigree
#'
#' @param pedigree Matrix or data frame describing a pedigree, with first four
#' columns being individual ID, mom ID, dad ID, and sex (female as
#' `0`, male as `1`).
#' @param L Length of chromosome in cM (or a vector of chromosome lengths)
#' @param xchr If TRUE, simulate X chromosome. (If `L` is a
#' vector, this should be a vector of TRUE/FALSE values, of the same
#' length as `L`, or a character string with the name of the X
#' chromosome, in `L`.)
#' @param m Crossover interference parameter, for chi-square model
#' (m=0 corresponds to no interference).
#' @param p proportion of crossovers coming from no-interference process
#' @param obligate_chiasma If TRUE, require an obligate chiasma on the
#' 4-strand bundle at meiosis.
#'
#' @return A list with each component being the data for one
#' individual, as produced by the [cross()] function.  Those
#' results are a list with two components, corresponding to the
#' maternal and paternal chromosomes. The chromosomes are represented
#' as lists with two components: an integer vector of alleles in
#' chromosome intervals, and a numeric vector of locations of the
#' right-endpoints of those intervals; these two vectors should have
#' the same length.
#'
#' If the input `L` is a vector, in order to simulate multiple
#' chromosomes at once, then the output will be a list with length
#' `length(L)`, each component being a chromosome and having the
#' form described above.
#'
#' @export
#' @keywords datagen
#' @seealso [check_pedigree()],
#' [sim_ril_pedigree()], [sim_ail_pedigree()],
#' [sim_ril_pedigree()]
#'
#' @examples
#' # simulate AIL pedigree
#' tab <- sim_ail_pedigree(12, 30)
#' # simulate data from that pedigree
#' dat <- sim_from_pedigree(tab)
#' # simulate multiple chromosomes
#' dat <- sim_from_pedigree(tab, c("1"=100, "2"=75, "X"=100), xchr="X")
sim_from_pedigree <-
    function(pedigree, L=100, xchr=FALSE, m=10, p=0, obligate_chiasma=FALSE)
{
    if(length(L) > 1) { # multiple chromosomes
        result <- vector("list", length(L))
        if(is.null(names(L))) names(L) <- seq(along=L)
        names(result) <- names(L)

        if(is.character(xchr)) # xchr is chromosome names
            xchr <- names(L) %in% xchr

        if(is.null(xchr))
            xchr <- rep(FALSE, length(L))

        if(length(xchr) == 1) # if single value, apply to all chromosomes
            xchr <- rep(xchr, length(L))

        if(length(xchr) != length(L))
            stop("length(xchr) != length(L)")

        for(i in seq(along=L))
            result[[i]] <- sim_from_pedigree(pedigree, L[i], xchr[i],
                                             m, p, obligate_chiasma)
        return(result)
    }

    if(length(unique(pedigree[,1])) != nrow(pedigree))
        stop("IDs must be unique")
    rownames(pedigree) <- pedigree[,1]

    result <- vector("list", nrow(pedigree))
    names(result) <- as.character(pedigree[,1])

    if(obligate_chiasma) Lstar <- calc_Lstar(L, m, p)
    else Lstar <- L

    for(i in 1:nrow(pedigree)) {
        if(pedigree[i,2]==0 || pedigree[i,3]==0) # founder
            result[[i]] <- create_parent(L, allele=pedigree[i,1])
        else {
            mom <- which(pedigree[,1]==pedigree[i,2])
            dad <- which(pedigree[,1]==pedigree[i,3])

            if(length(mom) < 1 || length(dad) < 1) {
                print(pedigree[i,,drop=FALSE])
                stop("parents not found")
            }
            if(mom >= i || dad >= i) {
                print(pedigree[i,,drop=FALSE])
                stop("Pedigree problem: parents follow individual")
            }

            result[[i]] <- cross(result[[mom]], result[[dad]],
                                 m=m, p=p,
                                 xchr=xchr, male=pedigree[i,4]==1,
                                 obligate_chiasma=obligate_chiasma,
                                 Lstar=Lstar)
        }
    }
    result
}
