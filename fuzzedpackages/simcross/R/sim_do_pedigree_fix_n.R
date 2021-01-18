#' Simulate a pedigree for Diversity Outbreds for a target sample size
#'
#' Simulate a pedigree for Diversity Outbred (DO) mice (a table of
#' individual, mom, dad, sex) so that the last generation reaches a
#' desired sample size.
#'
#' @param ngen Number of generations of outbreeding
#' @param npairs Number of breeding pairs at each generation. If
#' missing, we use 30 when `method="last2"` and 300 when
#' `method="sub2"`.
#' @param nkids_per Number of offspring per pair for the last
#' generation
#' @param design How to choose crosses: either random but avoiding
#' siblings, or completely at random
#' @param method Method used to generate the pedigree: either expand
#' at the last two generations or generate a pedigree with a large
#' number of pairs and then select a subset to have the desired sample
#' size. With `method="fixcc"`, we use the pre-CC generations as performed at the
#' Jackson Lab.
#' @param selc.method Method used to select the individuals from last
#' generation.
#' @param nsample_ngen Number of individuals desired at the last
#' generation
#' @param nccgen The number of generations for each CC line, only used
#' when `method` is not `"fixcc"`.
#'
#' @details The default number of breeding pairs depends on the chosen
#' `method`. With `method="last2"`, the default is `npairs=30`;
#' with `method="sub2"`, the default is `npairs=300`;
#' with `method="fixcc"`, `npairs` is ignored and is fixed at 144.
#'
#' @return A data frame with six columns: individual ID, mother ID, father
#' ID, sex, generation, and TRUE/FALSE indicator for whether DO or pre-DO.
#' Founders have `0` for mother and father ID. Sex is coded 0 for
#' female and 1 for male.
#'
#' @export
#' @seealso [sim_from_pedigree()],
#' [sim_ril_pedigree()], [sim_ail_pedigree()],
#' [sim_do_pedigree()], [sim_4way_pedigree()],
#' [sim_ail_pedigree_fix_n()]
#'
#' @examples
#' tab <- sim_do_pedigree_fix_n(8)
sim_do_pedigree_fix_n <- function(ngen=12, nkids_per=5, nccgen=15,
                                  nsample_ngen=150, npairs=NULL,
                                  method=c("last2", "sub2", "fixcc"),
                                  design=c("nosib", "random"),
                                  selc.method=c("byfamily","byindiv")){
  method <- match.arg(method)
  design <- match.arg(design)
  selc.method <- match.arg(selc.method)

  if(method =="last2"){
    if(is.null(npairs))
      npairs <- 30
    ped <- sim_do_pedigree_last2(ngen = ngen, npairs = npairs, nkids_per = nkids_per,
                                 nccgen = nccgen, design = design,
                                 nsample_ngen = nsample_ngen)
  } else if(method =="sub2"){
    if(is.null(npairs))
      npairs <- 300
    ped <- sim_do_pedigree(ngen = ngen, npairs = npairs,
                           ccgen = rep(nccgen, npairs),
                           nkids_per = nkids_per, design = design)
  } else if(method=="fixcc"){
    ## use fixed ccgen value
    alpha <- c(21, 64, 24, 10, 5, 9, 5, 3, 3)
    names(alpha) <- 4:12
    npairs <- sum(alpha)
    ccgen <- rep(as.numeric(names(alpha)), alpha)
    ped <- sim_do_pedigree(ngen = ngen, npairs = npairs, ccgen=ccgen,
                           nkids_per = nkids_per, design = design)
  }

  ## selecting samples...
  if(method == "last2"){
    id <- which(ped[, "gen"] == ngen & ped[,"do"]==1)
  } else {
    if(selc.method == "byfamily"){
      npairs.selc <- nsample_ngen / nkids_per
      selc.dad <- sample(which(ped[, "gen"] == ngen-1 &
                               ped[, "sex"] == 1 & ped[,"do"]==1),
                         size=npairs.selc)
      id <- which(ped[,"dad"] %in% selc.dad &
                  ped[, "gen"] == ngen & ped[,"do"]==1)
    }else{
      id <- ped[ped[, "gen"] == ngen & ped[, "do"] == 1, "id"]
      id <- sample(id, nsample_ngen)
    }
  }

  attr(ped, "last.gen.id") <- id
  ped
}


sim_do_pedigree_last2 <- function(ngen=12, npairs=30, nkids_per=5, nccgen=15,
                                  design=c("nosib", "random"),
                                  nsample_ngen=150)
{
  design <- match.arg(design)
  npairs_la2 <- ceiling(nsample_ngen/nkids_per)
  nkids_la2 <- ceiling(npairs_la2*2/npairs)

  ## second-to-last generation
  ped <- sim_do_pedigree(ngen = ngen-1, npairs = npairs,
                         ccgen = rep(nccgen, npairs),
                         nkids_per=nkids_la2, design=design)
  id <- ped[,"id"]
  mom <- ped[,"mom"]
  dad <- ped[,"dad"]
  sex <- ped[,"sex"]
  gen <- ped[,"gen"]
  do <- ped[,"do"]

  ## last generation
  n.last <- npairs_la2*nkids_per
  kids <- 1:(n.last)+max(id)

  wh <- which(ped[,"gen"] == ngen-1 & ped[, "sex"] == 0 & ped[, "do"] == 1)
  moms <- ped[wh, "id"]
  wh <- which(ped[,"gen"] == ngen-1 & ped[, "sex"] == 1 & ped[, "do"] == 1)
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
  do <- c(do, rep(1, length(wh)))

  data.frame(id=id, mom=mom, dad=dad, sex=sex, gen=gen, do=do)
}
