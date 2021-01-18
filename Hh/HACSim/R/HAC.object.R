## Master class, both functions return an object of this type
HACClass <- function(input.seqs = NULL,
                     subset.seqs = NULL,
                     prop.seqs = NULL,
                     prop.haps = NULL,
                     subset.haps = NULL,
                     N = NA,
                     Hstar = NA,
                     probs = NA,
                     p = NA,
                     perms = NA,
                     conf.level = NA,
                     num.iters = NA,
                     progress = NA,
                     filename = NULL) {
  HACObject <- list(
    input.seqs = input.seqs,
    subset.seqs = subset.seqs,
    prop.seqs = prop.seqs,
    prop.haps = prop.haps,
    subset.haps = subset.haps,
    N = N,
    Hstar = Hstar,
    probs = probs,
    p = p,
    perms = perms,
    conf.level = conf.level,
    num.iters = num.iters,
    progress = progress,
    filename = filename
  )

  ## Set the name for the class
  class(HACObject) <- append(class(HACObject), "HAC")
  return(HACObject)
}

## This uses DNA sequences from a FASTA file
HACReal <- function(perms = 10000,
                    p = 0.95,
                    conf.level = 0.95,
                    subsample = FALSE,
                    prop = NULL,
                    progress = TRUE,
                    num.iters = NULL,
                    filename = NULL) {
  # input.seqs <- TRUE # analyze DNA sequence file?
  # subset.haps <- NULL # subset haplotypes?
  # prop.haps <- NULL # proportion of haplotypes to subsample
  # subset.seqs <- TRUE # subset DNA sequences?

  if (subsample == TRUE) {
    subset.seqs <- TRUE
    prop.seqs <- prop # proportion of DNA sequences to subsample
  } else {
    subset.seqs <- FALSE
    prop.seqs <- NA # proportion of haplotypes to subsample
  }

  objectHAC <- HACClass(
    input.seqs = TRUE,
    subset.seqs = subset.seqs,
    prop.seqs = prop.seqs,
    perms = perms,
    p = p,
    conf.level = conf.level,
    num.iters = num.iters,
    progress = progress,
    filename = filename
  )

  return(objectHAC)
}

## For simulating custom data distributions
HACHypothetical <- function(N,
                            Hstar,
                            probs,
                            perms = 10000,
                            p = 0.95,
                            conf.level = 0.95,
                            subsample = FALSE,
                            prop = NULL,
                            progress = TRUE,
                            num.iters = NULL,
                            filename = NULL) {
  if (missing(N)) {
    stop("Please provide a value for N")
  }
  if (missing(Hstar)) {
    stop("Please provide a value for Hstar")
  }
  if (missing(probs)) {
    stop("Please provide a value for probs")
  }

  # Representation of what the parameters will be

  # input.seqs <- FALSE # subset DNA sequences?
  # subset.seqs <- FALSE # subset DNA sequences?
  # prop.seqs <- NULL # proportion of DNA sequences to subsample
  # subset.haps <- NULL # subset haplotypes?

  if (subsample == TRUE) {
    prop.haps <- prop # proportion of haplotypes to subsample
    subset.haps <- sort(sample(Hstar, size = ceiling(prop.haps * Hstar), replace = FALSE))
  } else {
    prop.haps <- NULL # proportion of haplotypes to subsample
    subset.haps <- NULL
  }
  objectHAC <- HACClass(
    N = N,
    Hstar = Hstar,
    probs = probs,
    input.seqs = FALSE,
    subset.seqs = FALSE,
    prop.haps = prop.haps,
    subset.haps = subset.haps,
    perms = perms,
    p = p,
    conf.level = conf.level,
    num.iters = num.iters,
    progress = progress,
    filename = filename
  )
  return(objectHAC)
}
