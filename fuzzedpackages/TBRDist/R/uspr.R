#' Calculate SPR, TBR and Replug distances on unrooted trees
#'
#' Calculate SPR, TBR and Replug distances on unrooted trees, and the
#' information content of the maximum agreement forest.
#'
#' Note that these distances are NP-hard to compute, so the running time of the
#' algorithms used in this software scale exponentially with the distance
#' computed.
#' The version of 'uspr' linked in this package is aimed at trees with up to
#' 50 leaves and uSPR distances up to 14.
#'
#' If you are interested in comparing rooted trees in terms of SPR operations,
#' you should use '[rspr](https://github.com/cwhidden/rspr/)' instead.
#' 'rspr' is also much more efficient and can easily handle pairs of binary
#' rooted trees with 200+ leaves and distances > 50.
#' rspr is not yet incorporated in this R package; please
#' [contact the maintainer](https://github.com/ms609/TBRDist/issues/2/)
#' if this would be useful to you.
#'
#'
#' @param tree1,tree2 Trees of class `phylo`, or lists thereof.
#' @param checks Logical specifying whether to check that trees are properly
#' formatted and labelled.  Specify `FALSE` at your peril, as improper
#' input is likely to crash R.
#' @param allPairs Logical; if `TRUE`, compare each tree in `tree1` with each
#' tree in `tree2`; if `FALSE`, compare each tree in `tree1` only with the
#' tree at the corresponding index in `tree2`.  If `tree2` is not specified,
#' each tree in `tree1` will be compared with each other tree in `tree1`.
#'
#' @param useTbrApproxEstimate,useTbrEstimate,useReplugEstimate Logical
#' specifying whether to use approximate TBR distance, TBR distance or Replug
#' distance to help estimate the SPR distance.
#'
#' @return `USPRDist()` returns a vector of SPR distances between each pair of
#' unrooted trees.
#'
#' @examples
#' tree1 <- TreeTools::BalancedTree(6)
#' tree2 <- TreeTools::PectinateTree(6)
#'
#' # SPR distance
#' USPRDist(tree1, tree2)
#'
#' @name TreeRearrangementDistances
#' @author
#' Algorithms implemented by Chris Whidden (<cwhidden@fredhutch.org>)
#'
#' R wrappers by Martin R. Smith (<martin.smith@durham.ac.uk>)
#'
#' @references
#' If you use these functions in your research, please cite:
#'
#' * Chris Whidden and Frederick A. Matsen IV. Calculating the Unrooted
#' Subtree-Prune-and-Regraft Distance.
#' arXiv:[1511.07529](https://arxiv.org/abs/1511.07529).
#'
#' @export
USPRDist <- function (tree1, tree2 = NULL, allPairs = is.null(tree2),
                      checks = TRUE,
                      useTbrApproxEstimate = TRUE,
                      useTbrEstimate = TRUE,
                      useReplugEstimate = TRUE) {
  treeLists <- .PrepareTrees(tree1, tree2, allPairs, checks, unroot = TRUE)
  ret <- uspr_dist(treeLists[[1]], treeLists[[2]],
                   useTbrApproxEstimate = useTbrApproxEstimate,
                   useTbrEstimate = useTbrEstimate,
                   useReplugEstimate = useReplugEstimate)
  .DistReturn(ret, tree1, tree2, allPairs)
}

#' @rdname TreeRearrangementDistances
#' @param maf Logical specifying whether to report a maximum agreement forest
#' corresponding to the optimal score.
#'
#' @return `ReplugDist()` returns a vector of Replug distances between each pair
#' of trees, or (if `maf = TRUE`) a named list whose second and third elements
#' list a vector of maximum agreement forests for each pair of trees.
#'
#' @examples
#' # Replug distance
#' ReplugDist(tree1, tree2)
#' ReplugDist(tree1, tree2, maf = TRUE)
#'
#' @export
ReplugDist <- function (tree1, tree2 = NULL, allPairs = is.null(tree2),
                        checks = TRUE, maf = FALSE) {
  treeLists <- .PrepareTrees(tree1, tree2, allPairs, checks, unroot = TRUE)
  ret <- replug_dist(treeLists[[1]], treeLists[[2]])
  if (maf) {
    names(ret) <- c('replug_dist', 'maf_1', 'maf_2')
    .DistReturn(ret, tree1, tree2, allPairs)
  } else {
    .DistReturn(ret[[1]], tree1, tree2, allPairs)
  }
}

#' @rdname TreeRearrangementDistances
#' @param exact Logical specifying whether to calculate the exact TBR distance.
#' @param approximate Logical specifying whether to calculate the approximate
#' TBR distance.  By default, is set to the opposite of `exact`; either
#' `approximate` or `exact` should usually be set to `TRUE` if a distance is
#' required.
#' @param countMafs Logical specifying whether to count the number of maximum
#' agreement forests found.
#' @param printMafs Logical specifying whether to print maximum agreement
#' forests to stdout whilst counting.
#' Use [`capture.output`]`(TBRDist(tree1, tree2, printMafs = TRUE))` to access
#' these in R.
#' @param optimize Logical specifying whether to use the default optimizations.
#' @param protectB Logical specifying whether to use the 'PROTECT_B'
#' optimization.
#' Overrides value inherited from `optimize`.
#' @return `TBRDist()` returns a named list, each element of which bears a
#' vector corresponding to the requested value for each tree pair.  If only the
#' exact value is requested (`exact = TRUE`), an unnamed vector of distances is
#' returned.
#'
#' @examples
#' # TBR distance between two trees
#' TBRDist(tree1, tree2, exact = TRUE)
#'
#' # Compare a list against one tree, using approximate distances
#' TBRDist(list(tree1, tree2), tree2, exact = FALSE)
#'
#' # Compare all pairs in two lists
#' TBRDist(list(tree1, tree2), list(tree1, tree2, tree2), allPairs = TRUE,
#'         exact = FALSE)
#'
#' # Compare each tree in a list against each other
#' TBRDist(list(one = tree1, two = tree2, twoAgain = tree2))
#'
#' # Compare each pair in two lists
#' TBRDist(list(tree1, tree2, tree2),
#'         list(tree2, tree1, tree2),
#'         exact = TRUE, approximate = TRUE, countMafs = TRUE)
#'
#' # Capture maximum agreement forests
#' mafs <- capture.output(TBRDist(tree1, tree2, approximate = FALSE,
#'                         printMafs = TRUE))
#' head(mafs)
#'
#' MAFInfo(tree1, tree2)
#' MAFInfo(list(tree2, tree1), list(tree1, tree2))
#' @export
TBRDist <- function (tree1, tree2 = NULL, allPairs = is.null(tree2),
                     checks = TRUE,
                     maf = FALSE, countMafs = FALSE, printMafs = FALSE,
                     exact = maf, approximate = !exact,
                     optimize = TRUE, protectB = TRUE) {
  if (!exact && !approximate && !countMafs && !printMafs) {
    warning("Nothing to do in TBRDist.")
  }
  if (maf && !exact) {
    warning("Maximum agreeement forest requires exact = TRUE")
  }
  treeLists <- .PrepareTrees(tree1, tree2, allPairs, checks, unroot = TRUE,
                             keepLabels = maf)

  whichRets <- c(exact, rep(approximate, 2L), countMafs,
                 rep(exact && maf, 2L))

  ret <- tbr_dist(treeLists[[1]], treeLists[[2]],
                  printMafs = printMafs, countMafs = countMafs,
                  optimize = optimize, protectB = protectB,
                  exact = exact, approximate = approximate)[whichRets]
  if (!any(whichRets)) {
    invisible()
  } else if (exact && sum(whichRets) == 1) {
    .DistReturn(ret[[1]], tree1, tree2, allPairs)
  } else {
    names(ret) <- c('tbr_exact', 'tbr_min', 'tbr_max', 'n_maf',
                    'maf_1', 'maf_2')[whichRets]
    .DistReturn(ret, tree1, tree2, allPairs)
  }
}

.Entropy <- function (p) -sum(p[p > 0] * log2(p[p > 0]))

#' @rdname TreeRearrangementDistances
#' @return `MAFInfo()` returns the information content of the maximum agreement
#' forest, in bits.  This is defined as the sum of the phylogenetic information
#' content of each constituent subtree, plus the entropy of the clusters implied
#' by the division of the tree into subtrees.  Note that as there is no
#' guarantee that the most informative MAF will be encountered,
#' this measure is approximate only.  `exact` will only serve to guarantee
#' that a MAF corresponding to the exact TBR distance is among those sampled.
#' @importFrom TreeDist .TreeDistance
#' @export
MAFInfo <- function(tree1, tree2 = tree1, exact = FALSE) {
  .TreeDistance(.MAFInfo, tree1, tree2, exact = exact)
}

#' @importFrom utils capture.output
#' @importFrom TreeTools LnUnrooted
.MAFInfo <- function (tree1, tree2 = NULL, exact = FALSE, ...) {
  mafs <- capture.output(TBRDist(tree1, tree2, exact = exact,
                                 approximate = FALSE, printMafs = TRUE))
  # c(F,T) ensures that output of exact = TRUE is disregarded.
  mafs <- unique(mafs[rep(c(FALSE, TRUE), length.out = length(mafs))])
  sizes <- lapply(strsplit(mafs, ';'), vapply, function (tr)
    lengths(regmatches(tr, gregexpr(",", tr))) + 1L, integer(1))

  phylogeneticInfo <- vapply(sizes, function (x) sum(LnUnrooted(x)), 0) / log(2)
  clusteringInfo <- vapply(sizes, function (x) .Entropy(x / sum(x)), 0)
  totalInfo <- phylogeneticInfo + clusteringInfo

  # Return:
  max(totalInfo)
}

#' Prepare trees for passing to uspr
#'
#' Converts trees to Newick strings, as expected by the 'uspr' C++ library.
#'
#' @param keepLabels Logical specifying whether to pass text labels to distance
#' calculator.  This is slower, but makes maximum agreement forests easier
#' to read.
#' @return `.PrepareTrees()` returns a two-element list, each entry of which is
#' a character vector listing trees in Newick format.
#'
#' @template MRS
#' @importFrom TreeTools as.Newick RenumberTips UnrootTree
#' @importFrom ape write.tree
#' @keywords internal
#' @export
.PrepareTrees <- function (tree1, tree2, allPairs = FALSE, checks = TRUE,
                           unroot = FALSE, keepLabels = FALSE) {
  if (unroot) {
    tree1 <- UnrootTree(tree1)
    tree2 <- UnrootTree(tree2)
  }

  if (checks) {

    if (inherits(tree1, 'phylo')) tree1 <- list(tree1)
    if (inherits(tree2, 'phylo')) tree2 <- list(tree2)

    if (allPairs) {
      if (is.null(tree2)) {
        nTree <- length(tree1)
        if (length(tree1) < 2L) {
          warning("Comparing a tree with itself only.")
          tree2 <- tree1
        } else {
          selector <- matrix(seq_len(nTree), nTree, nTree)
          tree2 <- tree1[t(selector)[lower.tri(selector)]]
          tree1 <- tree1[selector[lower.tri(selector)]]
        }
      } else {
        nTree1 <- length(tree1)
        tree1 <- rep(tree1, each = length(tree2))
        tree2 <- rep(tree2, nTree1)
      }
    } else {
      if (length(tree1) != length(tree2)) {
        if (length(tree1) == 1L) {
          tree1 <- rep(tree1, length(tree2))
        } else if (length(tree2) == 1L) {
          tree2 <- rep(tree2, length(tree1))
        } else {
          stop("if allPairs = FALSE, tree1 and tree2 must contain the same ",
               "number of trees, or a single tree.")
        }
      }
    }

    vapply(seq_along(tree1),
           function (i) .CatchBadPair(i, tree1[[i]], tree2[[i]]),
           integer(0))
  }
  if (keepLabels) {
    class(tree1) <- class(tree2) <- 'multiPhylo'
    list(write.tree(tree1), write.tree(tree2))
  } else {
    tree1[-1] <- lapply(tree1[-1], RenumberTips, tree1[[1]]$tip.label)
    tree2 <- lapply(tree2, RenumberTips, tree1[[1]]$tip.label)
    list(as.Newick(tree1), as.Newick(tree2))
  }
}


#' Internal functions
#'
#' These helper functions are unlikely to be of day-to-day use, but are
#' exported in case they are valuable to package developers.
#'
#' @template MRS
#' @name internals

#' @rdname internals
#' @param i Integer iterating tree pair.
#' @param tree1,tree2 Trees of class `phylo`.
#' @return `.CatchBadPair()` returns `integer(0)`, but will `stop` R with an
#' error message if the labels of `tree1` and `tree2` differ in anything but
#' order.
#' @keywords internal
#' @export
.CatchBadPair <- function (i, tree1, tree2) {
  lab1 <- tree1$tip.label
  lab2 <- tree2$tip.label
  if (length(lab1) != length(lab2)) {
    stop("Problem with tree pair " , i, ": Trees must be the same size.")
  }
  if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0){
    stop("Problem with tree pair " , i, ": Trees must bear identical labels.")
  }
  integer(0)
}

#' @rdname internals
#' @param ret A list containing the results of a tree comparison, probably
#' performed in C.  Each element of the list will correspond to an aspect
#' of tree distance (e.g. TBR distance, MAF composition).
#' @param allPairs Logical specifying whether all trees were compared with all
#' other trees, in which case a structure of class `dist` will be returned for
#' numeric entries, and a symmetric matrix (with diagonal marked `NA`) returned
#' for non-numeric entries.
#' @return `.DistReturn()` returns a structure of class `numeric`, `matrix` or
#' `dist` containing the distances between each pair of trees.
#' @keywords internal
#' @export
.DistReturn <- function (ret, tree1, tree2, allPairs) {
  .ReturnMatrix <- function (dat) {
    matrix(dat, length(tree1), length(tree2), byrow = TRUE,
           dimnames = list(names(tree1), names(tree2)))
  }
  .ReturnDist <- function (dat, nTree) {
    if (is.numeric(dat)) {
      ret <- structure(dat, Size = nTree, Diag = FALSE, Upper = FALSE,
                       Labels = names(tree1), class = 'dist')
    } else {
      ret <- matrix(NA, nTree, nTree)
      ret[lower.tri(ret)] <- dat
      ret[upper.tri(ret)] <- t(dat)
    }

    # Return:
    ret
  }

  if (inherits(tree1, 'phylo')) tree1 <- list(tree1)

  if (allPairs) {
    if (is.null(tree2)) {
      nTree1 <- length(tree1)
      if (nTree1 < 2L) {
        ret
      } else {
        if (is.list(ret)) {
          lapply(ret, .ReturnDist, nTree1)
        } else {
          .ReturnDist(ret, length(tree1))
        }
      }
    } else {
      names1 <- names(tree1)
      names2 <- names(tree2)
      if (is.list(ret)) {
        lapply(ret, .ReturnMatrix)
      } else {
        .ReturnMatrix(ret)
      }
    }
  } else {
    ret
  }
}
