#' Simulate AIL pedigree with fixed n
#'
#' Simulate a pedigree for advanced intercross lines (a table of
#' individual, mom, dad, sex) so that the last generation reaches a
#' desired sample size n
#'
#' @param ngen Number of generations of outbreeding
#' @param npairs Number of breeding pairs at each generation. If
#' missing, we use 30 when `method="last2"` and 300 when
#' `method="sub2"`.
#' @param nkids_per Number of offspring per pair for the last
#' generation
#' @param design How to choose crosses: either random but avoiding
#' siblings, or completely at random
#' @param method Method used to generate pedigree: either expand at the last
#' two generations or generate a pedigree with a large number of pairs and
#' select a subset to have the desired sample size.
#' @param nsample_ngen Number of individuals desired at the last
#' generation
#'
#' @details The default value for `npairs` depends on the choice of `method`.
#' For `method="last2"`, we use a default of `npairs=30`; for
#' `method="sub2"`, we use a default of `npairs=300`.
#'
#' @return A data frame with five columns: individual ID, mother ID,
#' father ID, sex, and generation.  Founders have `0` for mother
#' and father ID. Sex is coded 0 for female and 1 for male.
#'
#'
#' @export
#' @keywords datagen
#' @seealso [sim_from_pedigree()],
#' [sim_ril_pedigree()], [sim_ail_pedigree()],
#' [sim_do_pedigree()], [sim_4way_pedigree()],
#' [sim_do_pedigree_fix_n()]
#'
#' @examples
#' tab <- sim_ail_pedigree_fix_n(12)
sim_ail_pedigree_fix_n <- function(ngen=12, nkids_per=5,
                                  nsample_ngen=150, npairs=NULL,
                                  method=c("last2", "sub2"),
                                  design=c("nosib", "random")){
  method <- match.arg(method)
  design <- match.arg(design)

  if(method =="last2"){
    if(is.null(npairs))
      npairs <- 30
    ped <- sim_ail_pedigree_last2(ngen = ngen, npairs = npairs, nkids_per = nkids_per,
                                  design = design, nsample_ngen = nsample_ngen)
    id <- which(ped[, "gen"] == ngen)
  } else if(method =="sub2"){
    if(is.null(npairs))
      npairs <- 300
    npairs.selc <- nsample_ngen / nkids_per
    ped <- sim_ail_pedigree(ngen = ngen, npairs = npairs,
                            nkids_per = nkids_per, design = design)
    selc.dad <- sample(which(ped[, "gen"] == ngen-1 & ped[, "sex"] == 1),
                       size=npairs.selc)
    id <- which(ped[, "gen"] == ngen & ped[,"dad"] %in% selc.dad)
  }
  attr(ped, "last.gen.id") <- id
  ped
}

sim_ail_pedigree_last2 <- function(ngen=12, npairs=30, nkids_per=5, design=c("nosib", "random"),
                                   nsample_ngen=150)
{
  design <- match.arg(design)
  npairs_la2 <- ceiling(nsample_ngen/nkids_per)
  nkids_la2 <- ceiling(npairs_la2*2/npairs)

  ## second-to-last generation
  ped <- sim_ail_pedigree(ngen=ngen-1, npairs=npairs, nkids_per=nkids_la2,
                          design=design)
  id <- ped[,"id"]
  mom <- ped[,"mom"]
  dad <- ped[,"dad"]
  sex <- ped[,"sex"]
  gen <- ped[,"gen"]

  ## last generation
  n.last <- npairs_la2*nkids_per
  kids <- 1:(n.last)+max(id)

  wh <- which(ped[,"gen"] == ngen-1 & ped[, "sex"] == 0)
  moms <- ped[wh, "id"]
  wh <- which(ped[,"gen"] == ngen-1 & ped[, "sex"] == 1)
  dads <- ped[wh, "id"]
  rownames(ped) <- ped[, "id"]

  while(design=="nosib") { # sample until no sibs
    dads <- sample(dads)
    if(all(ped[dads, "dad"] != ped[moms, "dad"]))
        break
  }

  wh <- sort(sample(n.last, size=nsample_ngen))
  id <- c(id, kids[wh])
  mom <- c(mom, rep(moms, each=nkids_per)[wh])
  dad <- c(dad, rep(dads, each=nkids_per)[wh])
  sex <- c(sex, rep_len(c(0,1), length.out=n.last)[wh])
  gen <- c(gen, rep(ngen, n.last)[wh])

  data.frame(id=id, mom=mom, dad=dad, sex=sex, gen=gen)
}
