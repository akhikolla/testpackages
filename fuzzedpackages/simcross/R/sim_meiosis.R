## sim_meiosis.R

# create_parent
#
#' Create a parent object
#'
#' Create a parent object
#'
#' @param L chromosome length in cM
#' @param allele vector of integers for alleles, of length 1 or 2
#'
#' @return A list with two components, for the individual's two
#' chromosomes.  Each is a list with alleles in chromosome intervals
#' (as integers) and locations of the right endpoints of those
#' intervals.
#'
#' @keywords datagen
#' @export
#' @seealso [cross()], [sim_meiosis()]
#'
#' @examples
#' create_parent(100, 1)
#' create_parent(100, 1:2)
create_parent <-
    function(L, allele=1)
{
    if(length(allele) == 1) allele <- rep(allele,2)
    if(length(allele) != 2)
        stop("allele should be of length 1 or 2")
    if(!is.integer(allele)) {
        if(!is.numeric(allele))
            stop("allele should be a vector with 1 or 2 integers")
        storage.mode(allele) <- "integer"
    }

    list(mat=list(alleles=allele[1], locations=L),
         pat=list(alleles=allele[2], locations=L))
}


# check that data for an individual conforms to expected format
check_individual <-
    function(ind, tol=1e-12)
{
    # list with two components, named "mat" and "pat"
    if(!is.list(ind) || length(ind)!=2 ||
       !all(names(ind) == c("mat", "pat")))
        stop('ind should be list with "mat" and "pat"')

    # check each chromosome
    for(i in 1:2) {
        # list with "alleles" and "locations" components
        if(!is.list(ind[[i]]) || length(ind[[i]])!=2 ||
           !all(names(ind[[i]]) == c("alleles", "locations")))
            stop('chromosome should be list with "alleles" and "locations"')

        alleles <- ind[[i]]$alleles
        locations <- ind[[i]]$locations

        # locations numeric
        if(!is.numeric(locations))
            stop("locations should be numeric")

        # alleles integer
        if(!is.integer(alleles))
            stop("alleles should be integers")

        # locations non-decreasing
        if( length(locations) > 1 && min(diff(locations)) < 0 )
            stop("locations should be non-decreasing")
    }

    # start and end positions are the same
    starts <- vapply(ind, function(a) a$location[1], 0.0)
    ends <- vapply(ind, function(a) a$location[length(a$location)], 0.0)
    stopifnot(abs(diff(starts)) < tol, abs(diff(ends)) < tol)

    TRUE
}

#' Simulate crossover locations using the Stahl model
#'
#'
#' Simulate crossover locations on a single meiotic product using the
#' Stahl model.
#'
#' @details Chiasma locations are a superposition of two
#' processes: a proportion p exhibiting no interference, and a
#' proportion (1-p) following the chi-square model with interference
#' parameter m.  Crossover locations are derived by thinning the
#' chiasma locations with probability 1/2.
#'
#' @param L length of chr in cM
#' @param m Interference parameter (`m=0` is no interference)
#' @param p Proportion of chiasmata from no-interference mechanism
#' (`p=0` gives pure chi-square model)
#' @param obligate_chiasma If TRUE, require an obligate chiasma on the
#' 4-strand bundle at meiosis.
#' @param Lstar Adjusted chromosome length, if
#' `obligate_chiasma=TRUE`. Calculated if not provided.
#'
#' @return Numeric vector of crossover locations, in cM
#'
#' @keywords datagen
#'
#' @details Simulations are under the Stahl model with the
#' interference parameter being an integer. This is an extension of
#' the chi-square model, but with chiasmata being the superposition of
#' two processes, one following the chi-square model and the other
#' exhibiting no interference.
#'
#' @examples
#' x <- sim_crossovers(200, 10, 0)
#' x <- sim_crossovers(200, 10, 0.04)
#' x <- sim_crossovers(100, 0, 0, obligate_chiasma=TRUE)
#'
#' @references
#' Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
#' interference in arabidopsis.  \emph{Genetics} \bold{160}, 1631--1639.
#'
#' Foss, E., Lande, R., Stahl, F. W. and Steinberg, C. M. (1993) Chiasma
#' interference as a function of genetic distance. \emph{Genetics}
#' \bold{133}, 681--691.
#'
#' Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis
#' of crossover interference using the chi-square model.  \emph{Genetics}
#' \bold{139}, 1045--1056.
#'
#' @export
#'
sim_crossovers <-
    function(L, m=10, p=0.0, obligate_chiasma=FALSE, Lstar=NULL)
{
    if(is.null(Lstar)) {
        if(obligate_chiasma)
            Lstar <- calc_Lstar(L, m, p)
        else Lstar <- L
    }

    .sim_crossovers(L, m, p, obligate_chiasma, Lstar)
}


#' Simulate meiosis
#'
#' Output a random meiotic product from an input individual.
#'
#' @param parent An individual object, as output by
#' [create_parent()] or [cross()]
#' @param m interference parameter for chi-square model
#' @param p Proportion of chiasmata coming from no-interference process.
#' @param obligate_chiasma If TRUE, require an obligate chiasma on the
#' 4-strand bundle at meiosis.
#' @param Lstar Adjusted chromosome length, if
#' `obligate_chiasma=TRUE`. Calculated if not provided.
#'
#' @return A list with alleles in chromosome intervals (as integers)
#' and locations of the right endpoints of those intervals.
#'
#' @details Simulations are under the Stahl model with the
#' interference parameter being an integer. This is an extension of
#' the chi-square model, but with chiasmata being the superposition of
#' two processes, one following the chi-square model and the other
#' exhibiting no interference.
#'
#' @references
#' Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
#' interference in arabidopsis.  \emph{Genetics} \bold{160}, 1631--1639.
#'
#' Foss, E., Lande, R., Stahl, F. W. and Steinberg, C. M. (1993) Chiasma
#' interference as a function of genetic distance. \emph{Genetics}
#' \bold{133}, 681--691.
#'
#' Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis
#' of crossover interference using the chi-square model.  \emph{Genetics}
#' \bold{139}, 1045--1056.
#'
#' @keywords datagen
#' @export
#' @seealso [create_parent()], [cross()],
#' [sim_crossovers()], [calc_Lstar()]
#'
#' @examples
#' ind <- create_parent(100, 1:2)
#' prod <- sim_meiosis(ind)
sim_meiosis <-
    function(parent, m=10, p=0.0, obligate_chiasma=FALSE, Lstar=NULL)
{
    if(is.null(Lstar)) {
        L <- max(parent$mat$locations)
        if(obligate_chiasma)
            Lstar <- calc_Lstar(L, m, p)
        else Lstar <- L
    }

    .sim_meiosis(parent, m, p, obligate_chiasma, Lstar)
}


# cross
#
#' Cross two individuals
#'
#' Simulate the cross of two individuals to create a
#' single progeny
#'
#' @param mom An individual object, as produced by
#' [create_parent()] or this function.
#' @param dad An individual object, as produced by
#' [create_parent()] or this function.
#' @param m interference parameter for chi-square model
#' @param p proportion of crossovers coming from no-interference process
#' @param xchr If TRUE, simulate X chromosome
#' @param male If TRUE, simulate a male (matters only if
#' `xchr=TRUE`)
#' @param obligate_chiasma If TRUE, require an obligate chiasma on the
#' 4-strand bundle at meiosis.
#' @param Lstar Adjusted chromosome length, if
#' `obligate_chiasma=TRUE`. Calculated if not provided.
#'
#' @details Simulations are under the Stahl model with the
#' interference parameter being an integer. This is an extension of
#' the chi-square model, but with chiasmata being the superposition of
#' two processes, one following the chi-square model and the other
#' exhibiting no interference.
#'
#' @return A list with two components, for the individual's two
#' chromosomes.  Each is a list with alleles in chromosome intervals
#' (as integers) and locations of the right endpoints of those
#' intervals.
#'
#' @keywords datagen
#' @export
#' @seealso [create_parent()], [sim_meiosis()],
#' [sim_crossovers()], [calc_Lstar()]
#'
#' @examples
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 1:2)
#' child <- cross(mom, dad)
cross <-
    function(mom, dad, m=10, p=0, xchr=FALSE, male=FALSE,
             obligate_chiasma=FALSE, Lstar=NULL)
{
    if(is.null(Lstar)) {
        L <- max(mom$mat$locations)
        if(obligate_chiasma)
            Lstar <- calc_Lstar(L, m, p)
        else Lstar <- L
    }

    if(!xchr) {
        return(list(mat=sim_meiosis(mom, m, p, obligate_chiasma, Lstar),
                    pat=sim_meiosis(dad, m, p, obligate_chiasma, Lstar)))
    }
    else {
        if(male)
            return(list(mat=sim_meiosis(mom, m, p, obligate_chiasma, Lstar),
                        pat=dad$pat))
        else
            return(list(mat=sim_meiosis(mom, m, p, obligate_chiasma, Lstar),
                        pat=dad$mat))
    }
}

#' Calculate adjusted chromosome length for obligate chiasma
#'
#' Calculate the reduced chromosome length that will give the target
#' expected number of chiasmata when conditioning on there being at
#' least one chiasma on the four-strand bundle.
#'
#' @param L Length of chromosome (in cM); must be > 50
#' @param m Interference parameter for chi-square model
#' @param p Proportion of chiasmata coming from no-interference
#' process
#'
#' @return Adjusted length of chromosome
#'
#' @keywords utilities
#' @importFrom stats uniroot dpois
#' @export
#' @seealso [cross()], [sim_meiosis()],
#' [sim_crossovers()]
#'
#' @examples
#' calc_Lstar(100, 0, 0)
#' calc_Lstar(60, 10, 0.1)
calc_Lstar <-
    function(L, m=0, p=0)
{
    if(L <= 50) stop("Must have L > 50")
    if(m < 0) stop("Must have m==0")
    if(!is.integer(m) && m %% 1 > 1e-8) {
        warning("m must be an non-negative integer; rounding")
        m <- round(m)
    }
    if(p < 0 || p > 1)
        stop("p must be in [0, 1]")
    if(p==1) { # if p == 1, might as well take m=0, p=0
        m <- 0
        p <- 0
    }

    func_to_zero <- function(Lstar, L, m=0, p=0) {
        if(m==0)
            denom <- 1 - exp(-Lstar/50)
        else {
            lambda1 <- Lstar/50 * (m+1) * (1-p)
            lambda2 <- Lstar/50 * p
            denom <- 1 - sum(dpois(0:m, lambda1) * (m+1 - (0:m))/(m+1)) * exp(-lambda2)
        }

        2*L - 2*Lstar / denom
    }

    uniroot(func_to_zero, c(1e-8, L), L=L, m=m, p=p,
            tol=sqrt(.Machine$double.eps))$root
}
