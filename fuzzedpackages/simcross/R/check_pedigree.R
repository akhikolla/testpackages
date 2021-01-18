# check_pedigree
#
#' Check a pedigree for errors
#'
#' Perform a series of checks on the tabular data for a pedigree,
#' checking for problems
#'
#' @param pedigree Numeric matrix or data frame with four columns: ID,
#' mom ID, dad ID, sex. Sex is coded as `0`=female,
#' `1`=male. There can be additional columns, but they'll be
#' ignored.
#' @param ignore_sex If TRUE, ignore the sex values completely
#' (appropriate for hermaphroditic species.)
#'
#' @return TRUE (invisibly) if everything is okay; otherwise gives an
#' error.
#'
#' @details The parents should be listed before any of their
#' offspring. Founders should have 0's for mother and father; all
#' others should have non-zero values for the parents, and the parents
#' should appear in the pedigree. Father should be male and mothers
#' should be female (unless `ignore_sex=TRUE`). Individual
#' identifiers should be unique and non-zero. There should be no
#' missing values anywhere. (`NA`s are allowed in the sex column
#' if `ignore_sex=TRUE`.)
#'
#' @export
#' @keywords utilities
#' @seealso [sim_from_pedigree()],
#' [sim_ril_pedigree()]
#'
#' @examples
#' tab <- sim_ril_pedigree(7)
#' check_pedigree(tab)
#'
check_pedigree <-
    function(pedigree, ignore_sex=FALSE)
{
    if(ncol(pedigree) < 4)
        stop("pedigree must have at least 4 columns")
    id <- pedigree[,1]
    mom <- pedigree[,2]
    dad <- pedigree[,3]
    sex <- pedigree[,4]

    # check for missing values
    if(any(is.na(id)))
        stop("missing values not allowed in id column")
    if(any(is.na(mom)))
        stop("missing values not allowed in mom column")
    if(any(is.na(dad)))
        stop("missing values not allowed in dad column")
    if(!ignore_sex && any(is.na(sex)))
        stop("missing values not allowed in sex column")

    # check numeric IDs: non-zero, unique
    if(any(id==0))
        stop("0 is not a valid ID")
    if(length(id) != length(unique(id)))
        stop("IDs must be unique")

    # mom and dad are both 0 or neither 0
    if(any((mom==0 & dad != 0) | (dad==0 & mom != 0)))
        stop("mom and dad must be both 0 or neither 0")

    # parents in pedigree
    match_mom <- match(mom, id)
    match_dad <- match(dad, id)
    if(any(mom!=0 & is.na(match_mom))) {
        not_found <- mom[mom!=0 & is.na(match_mom)]
        stop("Some moms not found: ", paste(not_found, collapse=" "))
    }
    if(any(dad!=0 & is.na(match_dad))) {
        not_found <- dad[dad!=0 & is.na(match_dad)]
        stop("Some dads not found: ", paste(not_found, collapse=" "))
    }
    stopifnot(dad==0 | !is.na(match(dad, id)))

    # parents precede offspring
    match_mom <- match_mom[mom != 0]
    momsub <- which(mom != 0)
    if(any(match_mom >= momsub)) {
        problems <- momsub[match_mom >= momsub]
        stop("Some moms appear after their offspring: ",
             paste(problems, collapse=" "))
    }
    match_dad <- match_dad[dad != 0]
    dadsub <- which(dad != 0)
    if(any(match_dad >= dadsub)) {
        problems <- dadsub[match_dad >= dadsub]
        stop("Some dads appear after their offspring: ",
             paste(problems, collapse=" "))
    }

    if(ignore_sex) return(TRUE) # skip the rest

    # check sex
    if(any(sex[match_mom] != 0)) {
        male_moms <- unique(pedigree$mom[momsub[sex[match_mom] != 0]])
        stop("Some moms are male: ",
             paste(male_moms, collapse=" "))
    }
    if(any(sex[match_dad] != 1)) {
        female_dads <- unique(pedigree$dad[dadsub[sex[match_dad] != 1]])
        stop("Some dads are female: ",
             paste(female_dads, collapse=" "))
    }

    TRUE
}
