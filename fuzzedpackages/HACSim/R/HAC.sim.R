### HACSim: Haplotype Accumulation Curve Simulator ###

##########

# Author: Jarrett D. Phillips
# Last modified: December 21, 2019

##########

## Best run in RStudio ##
## DO NOT change order of code (can throw errors)! ##

#####

## Input parameters ###

# Required #

# N = Number of specimens (DNA sequences)
# Hstar = Number of observed unique haplotypes
# probs = Probability frequency distribution of haplotypes

# Optional #

# p = Proportion of unique haplotypes to recover
# perms = Number of permutations (replications)
# input.seqs = Analyze inputted aligned/trimmed FASTA DNA sequence file (TRUE / FALSE)?
# subset.seqs = Subset of DNA sequences to sample
# prop.seqs = Proportion of DNA sequences to sample
# prop.haps = Proportion of haplotypes to sample
# subset.haps = Subset of haplotypes to sample
# conf.level = Confidence level for accumulation curve and confidence intervals
# num.iters = Number of iterations to perform (1 or NULL (i.e., all))
# progress = Print results to console (TRUE / FALSE)

#####

HAC.sim <- function(N,
                    Hstar,
                    probs,
                    perms = 10000,
                    K = 1, # DO NOT CHANGE
                    p = 0.95,
                    subset.haps = NULL,
                    prop.haps = NULL,
                    input.seqs = FALSE,
                    subset.seqs = FALSE,
                    prop.seqs = NULL,
                    conf.level = 0.95,
                    df = NULL, # dataframe
                    num.iters = NULL,
                    progress = TRUE) {
  if ((is.null(num.iters)) || (num.iters == 1)) {
    cat("\n \n")

    ## Display progress bar ##

    if (progress == TRUE) {
      pb <- utils::txtProgressBar(min = 0, max = 1, style = 3)
    }

    ## Load DNA sequence data and set N, Hstar and probs ##

    if (input.seqs == TRUE) {
      seqs <- read.dna(file = file.choose(), format = "fasta")

      bf <- base.freq(seqs, all = TRUE)[5:17]

      if (any(bf > 0)) {
        warning("Inputted DNA sequences contain missing and/or ambiguous 
	    nucleotides, which may lead to overestimation of the number of 
	    observed unique haplotypes. Consider excluding sequences or alignment 
	    sites containing these data. If missing and/or ambiguous bases occur 
	    at the ends of sequences, further alignment trimming is an option.")
      }

      assign("ptm", proc.time(), envir = envr)

      if (subset.seqs == TRUE) { # take random subset of sequences (e.g., prop.seqs = 0.10 (10%))
        seqs <- seqs[sample(nrow(seqs), size = ceiling(prop.seqs * nrow(seqs)), replace = FALSE), ]
        seqsfile <- tempfile(fileext = ".fas")
        write.dna(seqs, file = seqsfile, format = "fasta")
      }

      assign("N", dim(seqs)[[1]], envir = envr)
      h <- sort(haplotype(seqs), decreasing = TRUE, what = "frequencies")
      rownames(h) <- 1:nrow(h)
      assign("Hstar", dim(h)[[1]], envir = envr)
      assign("probs", lengths(attr(h, "index")) / N, envir = envr)
    }

    ## Error messages ##

    if (N < Hstar) {
      stop("N must be greater than or equal to Hstar")
    }

    if (N == 1) {
      stop("N must be greater than 1")
    }

    if (Hstar == 1) {
      stop("H* must be greater than 1")
    }

    if (!isTRUE(all.equal(1, sum(probs), tolerance = .Machine$double.eps^0.25))) {
      stop("probs must sum to 1")
    }

    if (perms == 1) {
      stop("perms must be greater than 1")
    }

    if ((p <= 0) || (p > 1)) {
      stop("p must be greater than 0 and less than or equal to 1")
    }

    ## Set up container to hold the identity of each individual from each permutation ##

    num.specs <- N

    ## Create an ID for each haplotype ##

    if (is.null(subset.haps)) {
      haps <- 1:Hstar
    } else {
      subset.haps <- subset.haps
    }

    ## Assign individuals (N) ##

    specs <- 1:num.specs

    ## Generate permutations. Assume each permutation has N individuals, and sample those
    ## individuals' haplotypes from the probabilities ##

    gen.perms <- function() {
      if (is.null(subset.haps)) {
        sample(haps, size = num.specs, replace = TRUE, prob = probs)
      } else {
        resample <- function(x, ...) x[sample.int(length(x), ...)]
        resample(subset.haps, size = num.specs, replace = TRUE, prob = probs[subset.haps])
      }
    }

    pop <- array(dim = c(perms, num.specs, K))

    for (i in 1:K) {
      pop[, , i] <- replicate(perms, gen.perms())
    }

    ## Perform haplotype accumulation ##

    HAC.mat <- accumulate(pop, specs, perms, K)
    # HAC.mat <- drop(HAC.mat)

    ## Update progress bar ##

    if (progress == TRUE) {
      utils::setTxtProgressBar(pb, i)
    }

    ## Calculate the mean and CI for number of haplotypes recovered over all permutations

    means <- apply(HAC.mat, MARGIN = 2, mean)
    # means <- colMeans2(HAC.mat)
    # means <- colmeans(HAC.mat)
    sds <- apply(HAC.mat, MARGIN = 2, sd)
    # sds <- colSds(HAC.mat)
    # sds <- colVars(HAC.mat, std = TRUE)

    lower <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, (1 - conf.level) / 2))
    # lower <- colQuantiles(HAC.mat, probs = (1 - conf.level) / 2)
    # lower <- colQuantile(HAC.mat, probs = (1 - conf.level) / 2)
    upper <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, (1 + conf.level) / 2))
    # upper <- colQuantiles(HAC.mat, probs = (1 + conf.level) / 2)
    # upper <- colQuantile(HAC.mat, probs = (1 + conf.level) / 2)

    ## Make data accessible to user ##

    assign("d", data.frame(specs, means, sds), envir = envr)

    ## Compute simple summary statistics and display output ##
    ## tail() is used here instead of max() because curves will not be monotonic if perms is not set high enough. When perms is large (say 10000), tail() is sufficiently close to max()

    P <- tail(means, n = 1)
    num <- N * Hstar
    num.haps <- length(subset.haps)

    if (is.null(subset.haps)) {
      Q <- Hstar - P
      assign("R", P / Hstar, envir = envr)
      S <- Q / Hstar
      assign("Nstar", num / P, envir = envr)
      assign("X", (num / P) - N, envir = envr)
    } else {
      Q <- num.haps - P
      assign("R", P / num.haps, envir = envr)
      S <- Q / num.haps
      assign("Nstar", (N * num.haps) / P, envir = envr)
      assign("X", ((N * num.haps) / P) - N, envir = envr)
    }

    if (envr$X < 0) {
      envr$X <- 0 # to ensure non-negative result
    }

    moe <- (qnorm((1 + conf.level) / 2) * (tail(envr$d$sds, n = 1) / tail(envr$d$means, n = 1)) * sqrt(N))

    assign("low", signif(N - moe), envir = envr)
    assign("high", signif(N + moe), envir = envr)

    ## Output results to R console and CSV file ##
    if (progress == TRUE) {
      cat(
        "\n \n --- Measures of Sampling Closeness --- \n \n",
        "Mean number of haplotypes sampled: ", P,
        "\n Mean number of haplotypes not sampled: ", Q,
        "\n Proportion of haplotypes sampled: ", envr$R,
        "\n Proportion of haplotypes not sampled: ", S,
        "\n \n Mean value of N*: ", envr$Nstar,
        "\n Mean number of specimens not sampled: ", envr$X
      )

      ## Plot the mean haplotype accumulation curve (averaged over perms number of curves) and haplotype frequency barplot ##
      par(mfrow = c(1, 2))
      
      if (is.null(subset.haps)) {
        plot(specs, means, type = "n", xlab = "Specimens sampled", ylab = "Unique haplotypes", ylim = c(1, Hstar), main = "Haplotype accumulation curve")
      } else {
        plot(specs, means, type = "n", xlab = "Specimens sampled", ylab = "Unique haplotypes", ylim = c(1, length(subset.haps)), main = "Haplotype accumulation curve")
      }
      
      polygon(x = c(specs, rev(specs)), y = c(lower, rev(upper)), col = "gray")
      lines(specs, means, lwd = 2)
      
      if (is.null(subset.haps)) {
        abline(h = envr$R * Hstar, v = max(envr$d$specs), lty = 2) # dashed line
        abline(h = p * Hstar, lty = 3) # dotted line
        barplot(num.specs * probs, xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = haps, main = "Haplotype frequency distribution")
      } else {
        abline(h = envr$R * length(subset.haps), v = max(envr$d$specs), lty = 2) # dashed line
        abline(h = p * length(subset.haps), lty = 3) # dotted line
        barplot(num.specs * (probs[subset.haps] / sum(probs[subset.haps])), xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = subset.haps, main = "Haplotype frequency distribution")
      }
    }

    df[nrow(df) + 1, ] <- c(P, Q, envr$R, S, envr$Nstar, envr$X)
    df
  }
} # end HAC.sim
