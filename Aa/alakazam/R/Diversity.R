# Clonal diversity analysis

#' @include Classes.R
NULL

#### Coverage functions ####

#' Calculate sample coverage
#' 
#' \code{calcCoverage} calculates the sample coverage estimate, a measure of sample 
#' completeness, for varying orders using the method of Chao et al, 2015, falling back 
#' to the Chao1 method in the first order case.
#'
#' @param    x  numeric vector of abundance counts.
#' @param    r  coverage order to calculate.
#' 
#' @return   The sample coverage of the given order \code{r}.
#' 
#' @references
#' \enumerate{
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by 
#'            generalizing the Good-Turing sample coverage theory. 
#'            Ecology. 2015 96, 11891201.
#' }
#' 
#' @seealso  
#' Used by \link{alphaDiversity}.
#'           
#' @examples
#' # Calculate clone sizes
#' clones <- countClones(ExampleDb, groups="sample_id")
#' 
#' # Calculate 1first order coverage for a single sample
#' calcCoverage(clones$seq_count[clones$sample_id == "+7d"])
#'
#' @export
calcCoverage <- function(x, r=1) {
    # Use traditional calculation for 1st order coverage
    if (r == 1) { return(calcChao1Coverage(x)) }
    
    # Use general form for 2nd order and higher coverage
    x <- x[x >= 1]
    n <- sum(x)
    fr <- sum(x == r)
    fs <- sum(x == r + 1)
    
    if (fr == 0) {
        stop("Cannot calculate coverage of order ", r, ". No abundance data with count=", r, ".")
    }
    if (fs == 0) {
        stop("Cannot calculate coverage of order ", r, ". No abundance data with count=", r + 1, ".")
    }
    
    a <- factorial(r)*fr / sum(x[x >= r]^r)
    b <- ((n - r)*fr / ((n - r)*fr + (r + 1)*fs))^r
    rC <- 1 - a*b
    
    return(rC)
}


# Calculate first order coverage
#
# @param    x  a numeric vector of species abundance as counts
#
# @returns  Coverage estimate.
calcChao1Coverage <- function(x) {
    x <- x[x >= 1]
    n <- sum(x)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    
    if (f2 > 0) {
        rC1 <- 1 - (f1 / n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2))
    } else {
        rC1 <- 1 - (f1 / n) * (((n - 1) * (f1 - 1)) / ((n - 1) * (f1 - 1) + 2))
    }
    
    return(rC1)
}


# Calculates diversity under rarefaction
# 
# Calculates Hill numbers under rarefaction
#
# @param    x  vector of observed abundance counts.
# @param    m  the sequence count to rarefy to.
#
# @return   The first order coverage estimate
inferRarefiedCoverage <- function(x, m) {
    x <- x[x >= 1]
    n <- sum(x)
    if (m > n) {
        stop("m must be <= the total count of observed sequences.")
    }
    
    # Unrarefied case
    if (m == n) {
        return(calcCoverage(x, r=1))
    }
    
    # Calculate rarefied coverage
    # TODO: Read up on this and fix
    #rC1 <- iNEXT:::Chat.Ind(x, m)
    y <- x[(n - x) >= m]
    rC1 <- 1 - sum(y/n * exp(lgamma(n - y + 1) - lgamma(n - m - y + 1) - lgamma(n) + lgamma(n - m)))
    
    return(rC1)
}


#### Abundance functions ####

# Calculate undetected species
# 
# Calculates the lower bound of undetected species counts using the Chao1 estimator.
#
# @param    x  vector of observed abundance counts.
# 
# @return   The count of undetected species.
inferUnseenCount <- function(x) {
    x <- x[x >= 1]
    n <- sum(x)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    
    if (f2 > 0) {
        f0 <- ceiling(((n - 1) * f1^2) / (n * 2 * f2))
    } else {
        f0 <- ceiling(((n - 1) * f1 * (f1 - 1)) / (n * 2))
    }
    
    return(f0)
}


# Define undetected species relative abundances
#
# @param    x  vector of detected species abundance counts.
# 
# @return   An adjusted detected species relative abundance distribution.
inferUnseenAbundance <- function(x) {
    x <- x[x >= 1]
    
    # Coverage
    rC1 <- calcCoverage(x, r=1)
    
    # Unseen count
    f0 <- inferUnseenCount(x)
    
    # Assign unseen relative abundance
    p <- rep((1 - rC1) / f0, f0)
    
    return(p)
}


# Adjustement to observed relative abundances
#
# @param    x  vector of observed abundance counts
#
# @return   An adjusted observed species relative abundance distribution.
adjustObservedAbundance <- function(x) {
    x <- x[x >= 1]
    n <- sum(x)
    
    # Coverage
    rC1 <- calcCoverage(x, r=1)
    
    # Calculate tuning parameter
    lambda <- (1 - rC1) / sum(x/n * exp(-x))
    
    # Define adjusted relative abundance
    p <- x/n * (1 -  lambda * exp(-x))
    
    return(p)
}


# Combined unseen inferrence and observed abundance adjustment
#
# @param    x  named vector of observed abundance counts by clone.
#
# @return   A vector containing the complete inferred abundance distribution.
#           Unseen species will be denote by a clone name starting with "U".
inferCompleteAbundance <- function(x) {
    # Infer complete abundance distribution
    p1 <- adjustObservedAbundance(x)
    p2 <- inferUnseenAbundance(x)
    names(p2) <- if (length(p2) > 0) { paste0("U", 1:length(p2)) } else { NULL }
    
    return(c(p1, p2))
}

#' Tabulates clones sizes
#' 
#' \code{countClones} determines the number of sequences and total copy number of 
#' clonal groups.
#'
#' @param    data    data.frame with Change-O style columns containing clonal assignments.
#' @param    groups  character vector defining \code{data} columns containing grouping 
#'                   variables. If \code{groups=NULL}, then do not group data.
#' @param    copy    name of the \code{data} column containing copy numbers for each 
#'                   sequence. If this value is specified, then total copy abundance
#'                   is determined by the sum of copy numbers within each clonal group.
#' @param    clone   name of the \code{data} column containing clone identifiers.
#' 
#' @return   A data.frame summarizing clone counts and frequencies with columns:
#'           \itemize{
#'             \item \code{clone_id}:    clone identifier. This is the default column
#'                                       name, specified with \code{clone='clone_id'}.
#'                                       If the function call uses Change-O 
#'                                       formatted data and \code{clone='CLONE'}, this
#'                                       column will have name \code{CLONE}.
#'             \item \code{seq_count}:   total number of sequences for the clone.
#'             \item \code{seq_freq}:    frequency of the clone as a fraction of the total
#'                                       number of sequences within each group.
#'             \item \code{copy_count}:  sum of the copy counts in the \code{copy} column.
#'                                       Only present if the \code{copy} argument is 
#'                                       specified.
#'             \item \code{copy_freq}:   frequency of the clone as a fraction of the total
#'                                       copy number within each group. Only present if 
#'                                       the \code{copy} argument is specified.
#'           }
#'           Also includes additional columns specified in the \code{groups} argument.
#' 
#' @examples
#' # Without copy numbers
#' clones <- countClones(ExampleDb, groups="sample_id")
#'
#' # With copy numbers and multiple groups
#' clones <- countClones(ExampleDb, groups=c("sample_id", "c_call"), copy="duplicate_count")
#' 
#' @export
countClones <- function(data, groups=NULL, copy=NULL, clone="clone_id") {
    # Check input
    check <- checkColumns(data, c(clone, copy, groups))
    if (check != TRUE) { stop(check) }
    
    # Tabulate clonal abundance
    if (is.null(copy)) {
        clone_tab <- data %>% 
            group_by(!!!rlang::syms(c(groups, clone))) %>%
            dplyr::summarize(seq_count=n()) %>%
            dplyr::mutate(seq_freq=!!rlang::sym("seq_count")/sum(!!rlang::sym("seq_count"), na.rm=TRUE)) %>%
            dplyr::arrange(desc(!!rlang::sym("seq_count")))
    } else {
        clone_tab <- data %>% 
            group_by(!!!rlang::syms(c(groups, clone))) %>%
            dplyr::summarize(seq_count=length(.data[[clone]]),
                              copy_count=sum(.data[[copy]], na.rm=TRUE)) %>%
            dplyr::mutate(seq_freq=!!rlang::sym("seq_count")/sum(!!rlang::sym("seq_count"), na.rm=TRUE),
                          copy_freq=!!rlang::sym("copy_count")/sum(!!rlang::sym("copy_count"), na.rm=TRUE)) %>%
            dplyr::arrange(desc(!!rlang::sym("copy_count"))) 
    }
    return(clone_tab)
}


# Perform boostrap abundance calculation
# 
# @param    x       named vector of observed abundance values.
# @param    n       number of samples to draw from the estimate complete abundance distribution.
# @param    nboot   number of bootstrap realizations.
# @param    method  complete abundance inferrence method. 
#                   One of "before", "after" or "none" for complete abundance distribution
#                   inferrence before sampling, after sampling, or uncorrected, respectively.
# 
# @return   A matrix of bootstrap results.
bootstrapAbundance <- function(x, n, nboot=200, method="before") {
    ## DEBUG
    # x=abund_obs; method="before"
    # Check argumets
    method <- match.arg(method)
  
    if (method == "before") {
        # Calculate estimated complete abundance distribution
        p <- inferCompleteAbundance(x)
        # Bootstrap abundance
        boot_mat <- rmultinom(nboot, n, p) / n
    } else if (method == "after") {
        # Calculate estimated complete abundance distribution
        p <- x / sum(x, na.rm=TRUE)
        boot_sam <- rmultinom(nboot, n, p)
        boot_list <- apply(boot_sam, 2, inferCompleteAbundance)
        
        # Convert to matrix
        boot_names <- unique(unlist(sapply(boot_list, names)))
        boot_mat <- matrix(0, nrow=length(boot_names), ncol=nboot)
        rownames(boot_mat) <- boot_names
        for (i in 1:nboot) {
            boot_mat[names(boot_list[[i]]), i] <- boot_list[[i]]
        } 
    } else if (method == "none") {
        # Raw sampling of input
        p <- x / sum(x, na.rm=TRUE)
        boot_sam <- rmultinom(nboot, n, p)
    } else {
        stop("Invalid method: ", method)
    }
      
    return(boot_mat)
}

#' Estimates the complete clonal relative abundance distribution
#' 
#' \code{estimateAbundance} estimates the complete clonal relative abundance distribution 
#' and confidence intervals on clone sizes using bootstrapping.
#' 
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    copy      name of the \code{data} column containing copy numbers for each 
#'                     sequence. If \code{copy=NULL} (the default), then clone abundance
#'                     is determined by the number of sequences. If a \code{copy} column
#'                     is specified, then clone abundances is determined by the sum of 
#'                     copy numbers within each clonal group.
#' @param    group     name of the \code{data} column containing group identifiers. 
#'                     If \code{NULL} then no grouping is performed and the \code{group} 
#'                     column of the output will contain the value \code{NA} for each row.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} then no 
#'                     maximum is set.
#' @param    uniform   if \code{TRUE} then uniformly resample each group to the same 
#'                     number of observations. If \code{FALSE} then allow each group to
#'                     be resampled to its original size or, if specified, \code{max_size}.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot     number of bootstrap realizations to generate.
#' @param    progress  if \code{TRUE} show a progress bar. 
#' 
#' @return   A \link{AbundanceCurve} object summarizing the abundances.
#'           
#' @references
#' \enumerate{
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
#'            A framework for sampling and estimation in species diversity studies. 
#'            Ecol Monogr. 2014 84:45-67.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by 
#'            generalizing the Good-Turing sample coverage theory. 
#'            Ecology. 2015 96, 11891201.
#' }
#' 
#' @seealso  
#' See \link{plotAbundanceCurve} for plotting of the abundance distribution.
#' See \link{alphaDiversity} for a similar application to clonal diversity.
#'           
#' @examples
#' abund <- estimateAbundance(ExampleDb, group="sample_id", nboot=100)
#'
#' @export
estimateAbundance <- function(data, clone="clone_id", copy=NULL, group=NULL, 
                              min_n=30, max_n=NULL, uniform=TRUE, ci=0.95, nboot=200,
                              progress=FALSE) {
    ## DEBUG
    # data=ExampleDb; group="sample_id"; clone="clone_id"; copy=NULL; min_n=1; max_n=NULL; ci=0.95; uniform=F; nboot=100
    # copy="duplicate_count"
    # group=NULL

    # Hack for visibility of dplyr variables
    . <- NULL
    
    # Check input
    if (!is.data.frame(data)) {
        stop("Input data is not a data.frame")
    }
    
    # Check columns that are reported are real columns (can be NULL)
    check <- checkColumns(data, c(clone, copy, group))
    if (check != TRUE) { stop(check) }
    
    # Set confidence interval
    ci_z <- ci + (1 - ci) / 2
    ci_x <- qnorm(ci_z)
    
    # Tabulate clonal abundance
    count_col <- if (!is.null(copy)) { "copy_count" } else { "seq_count" }
    clone_tab <- countClones(data, copy=copy, clone=clone, groups=group) %>%
        dplyr::mutate(clone_count=!!rlang::sym(count_col))

    # Tabulate group sizes
    if (!is.null(group)) {
        # Summarize groups
        group_tab <- clone_tab %>%
            group_by(!!rlang::sym(group)) %>%
            dplyr::summarize(count=sum(!!rlang::sym("clone_count"), na.rm=TRUE)) %>%
            rename(group=!!rlang::sym(group))
    } else {
        group_tab <- data.frame(v="All", count=sum(clone_tab$clone_count, na.rm=T))
        names(group_tab)[1] <- "group"
    }
    group_all <- as.character(group_tab$group)
    group_tab <- group_tab[group_tab$count >= min_n, ]
    group_keep <- as.character(group_tab$group)
    
    # Set number of sampled sequence
    if (uniform) {
        nsam <- min(group_tab$count, max_n)
        nsam <- setNames(rep(nsam, length(group_keep)), group_keep)
    } else {
        nsam <- if (is.null(max_n)) { group_tab$count } else { pmin(group_tab$count, max_n) }
        nsam <- setNames(nsam, group_keep)
    }
    
    # Warn if groups removed
    if (length(group_keep) < length(group_all)) {
        warning("Not all groups passed threshold min_n=", min_n, ".", 
                " Excluded: ", paste(setdiff(group_all, group_keep), collapse=", "))
    }
    
    # Generate abundance bootstrap
    if (progress) { 
        pb <- progressBar(length(group_keep))
    }
    boot_list <- list()
    abund_list <- list()
    for (g in group_keep) {
        n <- nsam[g]
        
        # Extract abundance vector
        if (!is.null(group)) {
            abund_obs <- clone_tab$clone_count[clone_tab[[group]] == g]
            names(abund_obs) <- clone_tab[[clone]][clone_tab[[group]] == g]
        } else {
            # Extract abundance vector
            abund_obs <- clone_tab$clone_count
            names(abund_obs) <- clone_tab[[clone]]
        } 
        
        # Infer complete abundance distribution
        boot_mat <- bootstrapAbundance(abund_obs, n, nboot=nboot, method="before")
        
        # Assign confidence intervals based on variance of bootstrap realizations
        p_mean <- apply(boot_mat, 1, mean)
        p_sd <- apply(boot_mat, 1, sd)
        p_err <- ci_x * p_sd
        p_lower <- pmax(p_mean - p_err, 0)
        p_upper <- p_mean + p_err
        
        # Assemble and sort abundance data.frame
	    abund_df <- tibble::tibble(!!clone := rownames(boot_mat), p=p_mean, p_sd=p_sd,
	                           lower=p_lower, upper=p_upper) %>%
	        dplyr::arrange(desc(!!rlang::sym("p"))) %>%
	        dplyr::mutate(rank=1:n())
			
        # Save summary
        abund_list[[g]] <- abund_df
        
        # Save bootstrap
        boot_list[[g]] <- as.data.frame(boot_mat) %>%
          tibble::rownames_to_column(clone)
        
        if (progress) { pb$tick() }
    }
    id_col <- if_else(is.null(group), "group", group)
    abundance_df <- as.data.frame(bind_rows(abund_list, .id=id_col))
    bootstrap_df <- as.data.frame(bind_rows(boot_list, .id=id_col))
    
    # Create a new diversity object with bootstrap
    abund_obj <- new("AbundanceCurve",
                     bootstrap=bootstrap_df, 
                     abundance=abundance_df,
                     clone_by=clone,
                     group_by=id_col,
                     #groups=if_else(is.null(group), as.character(NA), group_keep),
                     groups=group_keep,
                     n=nsam, 
                     nboot=nboot, 
                     ci=ci)
    
    return(abund_obj)
}

#### Diversity functions ####

#' Calculate the diversity index
#' 
#' \code{calcDiversity} calculates the clonal diversity index for a vector of diversity 
#' orders. 
#'
#' @param    p  numeric vector of clone (species) counts or proportions.
#' @param    q  numeric vector of diversity orders.
#' 
#' @return   A vector of diversity scores \eqn{D} for each \eqn{q}.
#' 
#' @details
#' This method, proposed by Hill (Hill, 1973), quantifies diversity as a smooth function 
#' (\eqn{D}) of a single parameter \eqn{q}. Special cases of the generalized diversity 
#' index correspond to the most popular diversity measures in ecology: species richness 
#' (\eqn{q = 0}), the exponential of the Shannon-Weiner index (\eqn{q} approaches \eqn{1}), the 
#' inverse of the Simpson index (\eqn{q = 2}), and the reciprocal abundance of the largest 
#' clone (\eqn{q} approaches \eqn{+\infty}). At \eqn{q = 0} different clones weight equally, 
#' regardless of their size. As the parameter \eqn{q} increase from \eqn{0} to \eqn{+\infty} 
#' the diversity index (\eqn{D}) depends less on rare clones and more on common (abundant) 
#' ones, thus encompassing a range of definitions that can be visualized as a single curve. 
#' 
#' Values of \eqn{q < 0} are valid, but are generally not meaningful. The value of \eqn{D} 
#' at \eqn{q=1} is estimated by \eqn{D} at \eqn{q=0.9999}. 
#'
#' @references
#' \enumerate{
#'   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#'            Ecology. 1973 54(2):427-32.
#' }
#' 
#' @seealso  Used by \link{alphaDiversity}.
#' 
#' @examples
#' # May define p as clonal member counts
#' p <- c(1, 1, 3, 10)
#' q <- c(0, 1, 2)
#' calcDiversity(p, q)
#'
#' # Or proportional abundance
#' p <- c(1/15, 1/15, 1/5, 2/3)
#' calcDiversity(p, q)
#' 
#' @export
calcDiversity <- function(p, q) {
    # Add jitter to q=1
    q[q == 1] <- 0.9999
    # Remove zeros
    p <- p[p > 0]
    # Convert p to proportional abundance
    p <- p / sum(p)
    # Calculate D for each q
    D <- sapply(q, function(x) sum(p^x)^(1 / (1 - x)))
    
    return(D)
}

# Calculate the inferred diversity index
# 
# \code{calcInferredDiversity} calculates the clonal diversity index for a vector of diversity 
# orders with a correction for the presence of unseen species. Does not take proportional abundance.
#
# @param    p  numeric vector of clone (species) counts.
# @param    q  numeric vector of diversity orders.
# 
# @return   A vector of diversity scores \eqn{D} for each \eqn{q}.
# 
# @details
# This method, proposed by Hill (Hill, 1973), quantifies diversity as a smooth function 
# (\eqn{D}) of a single parameter \eqn{q}. Special cases of the generalized diversity 
# index correspond to the most popular diversity measures in ecology: species richness 
# (\eqn{q = 0}), the exponential of the Shannon-Weiner index (\eqn{q} approaches \eqn{1}), the 
# inverse of the Simpson index (\eqn{q = 2}), and the reciprocal abundance of the largest 
# clone (\eqn{q} approaches \eqn{+\infty}). At \eqn{q = 0} different clones weight equally, 
# regardless of their size. As the parameter \eqn{q} increase from \eqn{0} to \eqn{+\infty} 
# the diversity index (\eqn{D}) depends less on rare clones and more on common (abundant) 
# ones, thus encompassing a range of definitions that can be visualized as a single curve. 
# 
# Values of \eqn{q < 0} are valid, but are generally not meaningful. The value of \eqn{D} 
# at \eqn{q=1} is estimated by \eqn{D} at \eqn{q=0.9999}. 
# 
# An adjusted detected species relative abundance distribution is applied before calculating diversity.
#
# @references
# \enumerate{
#   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#            Ecology. 1973 54(2):427-32.
# }
# 
# @seealso  Used by \link{alphaDiversity}
# 
# @examples
# # May define p as clonal member counts
# p <- c(1, 1, 3, 10)
# q <- c(0, 1, 2)
# calcInferredDiversity(p, q)
#
# 
# @export
# calcInferredDiversity <- function(p, q) {
#    # Correct abundance
#     .infer <- function(y) {
#         # Infer complete abundance distribution
#         p1 <- adjustObservedAbundance(y)
#         p2 <- inferUnseenAbundance(y)
#         names(p2) <- if (length(p2) > 0) { paste0("U", 1:length(p2)) } else { NULL }
#         return(c(p1, p2))
#     }
#    
#    # Correct abundance
#    p <- .infer(p)
#     # Add jitter to q=1
#     q[q == 1] <- 0.9999
#     # Remove zeros
#     p <- p[p > 0]
#     # Convert p to proportional abundance
#     p <- p / sum(p)
#     # Calculate D for each q
#     D <- sapply(q, function(x) sum(p^x)^(1 / (1 - x)))
#     
#     return(D)
# }


# Calculates diversity under rarefaction
# 
# Calculates Hill numbers under rarefaction
#
# @param    x  vector of observed abundance counts.
# @param    q  numeric vector of diversity orders.
# @param    m  the sequence count to rarefy to.
#
# @return   A vector of diversity scores \eqn{D} for each \eqn{q}.
inferRarefiedDiversity <- function(x, q, m) {
    x <- x[x >= 1]
    n <- sum(x)
    if (m > n) {
        stop("m must be <= the total count of observed sequences.")
    }
    q[q == 1] <- 0.9999
    
    # Tabulate frequency counts from 1:n
    fk_n <- tabulate(x, nbins=n)
    
    # Calculate estimated fk(m)
    fk_m <- sapply(1:m, function(k) sum(exp(lchoose(k:m, k) + 
                                            lchoose(n - k:m, m - k) - 
                                            lchoose(n, m))*fk_n[k:m]))
    
    # Calculate diversity
    D <- sapply(q, function(r) sum((1:m / m)^r * fk_m)^(1 / (1 - r)))
    
    return(D)
}


# Helper function for computing alpha diversity from bootrstrap outputs
#
# \code{helperAlpha} divides a set of bootstrapped clones by group annotation,
# and computes the diversity of each set. 
#
# @param    boot_output  data.frame from\link{AbundanceCurve} object containing bootstrapped clonal 
#                        abundance curves.
# @param    q            vector of Hill Diversity indices to test for diversity calculations.
# @param    clone        name of the \code{boot_output} column containing clone identifiers.
# @param    group        name of the \code{boot_output} column containing grouping information for 
#                        diversity calculation.
#
# @return   data.frame containing diversity calculations for each bootstrap iteration.
helperAlpha <- function(boot_output, q, clone="clone_id", group=NULL) {
    ## DEBUG
    # abundance <- estimateAbundance(ExampleDb, group="sample_id", nboot=100)
    # clone <- abundance@clone_by
    # group <- abundance@group_by
  
    # Compute diversity from a column of each bootstrap
    output <- boot_output %>% 
        dplyr::ungroup() %>%
        dplyr::select(-one_of(c(clone, group))) %>%
        as.matrix() %>% 
        apply(2, calcDiversity, q=q) %>%
        data.frame() %>% 
        mutate(q=q)

    return(output)
}


# Helper function for computing beta diversity from bootrstrap outputs
#
# \code{helperBeta} divides a set of bootstrapped clones by group annotation,
# and computes the alpha diversity. Group annotations are then ignored and 
# gamma diversity is computed. A multiplicative beta diversity is used corresponding
# to the gamma diversity divided by the average alpha diversity of each group.
#
# @param    boot_output   data.frame from\link{AbundanceCurve} object containing bootstrapped clonal abundance curves.
# @param    q             vector of Hill Diversity indices to test for diversity calculations.
# @param    ci_z          numeric value corresponding to confidence interval for calculating beta diversity.
# @param    clone         name of the \code{boot_output} column containing clone identifiers.
# @param    group         name of the \code{boot_output} column containing grouping information for diversity 
#                         calculation.
#
# @return   data.frame containing diversity calculations for each bootstrap iteration.
helperBeta <- function(boot_output, q, ci_x, clone="clone_id", group="group") { 
    # Hack for visibility of dplyr variables
    . <- NULL
        
    # Compute gamma diversity metrics
    gamma <- boot_output %>%
        dplyr::group_by(!!rlang::sym(clone)) %>%
        dplyr::select(-one_of(c(group))) %>%
        dplyr::summarize_all(sum) %>%
        dplyr::do(helperAlpha(., q=q, clone=clone)) %>%
        tidyr::gather(key="n", value="gamma", -!!rlang::sym("q")) %>%
        dplyr::mutate(gamma=as.numeric(!!rlang::sym("gamma")))

    # Compute alpha diversity metrics
    alpha <- boot_output %>%
        dplyr::group_by(!!rlang::sym(group)) %>%
        dplyr::do(helperAlpha(., q=q, clone=clone, group=group)) %>%
        dplyr::group_by(!!rlang::sym("q")) %>%
        dplyr::select(-one_of(c(group))) %>%
        dplyr::summarize_all(mean) %>%
        tidyr::gather(key="n", value="alpha", -!!rlang::sym("q")) %>%
        dplyr::mutate(alpha=as.numeric(!!rlang::sym("alpha")))

    # Perform comparisons of alpha and gamma to extract beta
    beta <- bind_cols(gamma, alpha) %>%
        dplyr::group_by(!!rlang::sym("q")) %>%
        dplyr::mutate(X=!!rlang::sym("gamma") / !!rlang::sym("alpha")) %>%
        dplyr::summarize(d=mean(!!rlang::sym("X"), na.rm=TRUE),
                         d_sd=sd(!!rlang::sym("X"), na.rm=TRUE)) %>%
        dplyr::mutate(d_lower=pmax(!!rlang::sym("d") - !!rlang::sym("d_sd") * ci_x, 0), 
                      d_upper=!!rlang::sym("d") + !!rlang::sym("d_sd") * ci_x)

    return(beta)
}

# Helper function for computing statistical significance
#
# \code{helperTest} computes the pairwise statistical significance of differences
# in bootstrapped diversity values between two sets defined by the group column. 
# A p-value is computed using the ECDF distribution as the frequency of bootstrap iterations
# for which no difference is observed. 
#
# @param    div_df  data.frame from\link{DiversityCurve} object containing bootstrapped 
#                   diversity curves.
# @param    group   name of the \code{boot_output} column containing grouping information 
#                   for diversity calculation.
# @param    q       vector of Hill Diversity indices to test for diversity calculations.
#
# @return   data.frame containing test results for each value of q.
helperTest <- function(div_df, q, group="group") {
    # Hack for visibility of dplyr variables
    . <- NULL
    
    # Pairwise test
    group_pairs <- combn(unique(div_df[[group]]), 2, simplify=F)
    pvalue_list <- list()
    for (group_pair in group_pairs) {
        pair_list <- list()
        for(q_i in q) {
            
            # Currently just testing for one diversity order
            mat1 <- div_df %>%
                dplyr::filter(!!rlang::sym(group) == group_pair[1], !!rlang::sym("q") == q_i) %>%
                dplyr::select(-one_of(c(group, "q"))) %>% unlist()
            mat2 <- div_df %>%
                dplyr::filter(!!rlang::sym(group) == group_pair[2], !!rlang::sym("q") == q_i) %>%
                dplyr::select(-one_of(c(group, "q"))) %>% unlist()

            if (mean(mat1) >= mean(mat2)) { 
                g_delta <- mat1 - mat2 
            } else { 
                g_delta <- mat2 - mat1 
            }  

            # Compute p-value from ecdf
            p <- ecdf(g_delta)(0)
            p <- ifelse(p <= 0.5, p * 2, (1 - p) * 2)
            
            pair_list[[as.character(q_i)]] <- list(delta_mean=mean(g_delta), 
                                                   delta_sd=sd(g_delta), 
                                                   pvalue=p)
        }
        pvalue_list[[paste(group_pair, collapse=" != ")]] <- bind_rows(pair_list, .id="q")

    }
    test_df <- bind_rows(pvalue_list, .id="test")

    return(test_df)
}


#' Calculate clonal alpha diversity
#'
#' \code{alphaDiversity} takes in a data.frame or \link{AbundanceCurve} and computes
#' diversity scores (\eqn{D}) over an interval of diversity orders (\eqn{q}).
#' 
#' @param    data      data.frame with Change-O style columns containing clonal assignments or
#'                     a \link{AbundanceCurve} generate by \link{estimateAbundance} object 
#'                     containing a previously calculated bootstrap distributions of clonal abundance.
#' @param    min_q     minimum value of \eqn{q}.
#' @param    max_q     maximum value of \eqn{q}.
#' @param    step_q    value by which to increment \eqn{q}.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    ...       additional arguments to pass to \link{estimateAbundance}. Additional arguments
#'                     are ignored if a \link{AbundanceCurve} is provided as input.
#' 
#' @return   A \link{DiversityCurve} object summarizing the diversity scores.
#' 
#' @references
#' \enumerate{
#'   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#'            Ecology. 1973 54(2):427-32.
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
#'            A framework for sampling and estimation in species diversity studies. 
#'            Ecol Monogr. 2014 84:45-67.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by 
#'            generalizing the Good-Turing sample coverage theory. 
#'            Ecology. 2015 96, 11891201.
#' }
#'  
#' @seealso  See \link{calcDiversity} for the basic calculation and 
#'           \link{DiversityCurve} for the return object. 
#'           See \link{plotDiversityCurve} for plotting the return object.
#' 
#' @details
#' Clonal diversity is calculated using the generalized diversity index (Hill numbers) 
#' proposed by Hill (Hill, 1973). See \link{calcDiversity} for further details.
#'
#' To generate a smooth curve, \eqn{D} is calculated for each value of \eqn{q} from
#' \code{min_q} to \code{max_q} incremented by \code{step_q}.  When \code{uniform=TRUE}
#' variability in total sequence counts across unique values in the \code{group} column 
#' is corrected by repeated resampling from the estimated complete clonal distribution to a 
#' common number of sequences. The complete clonal abundance distribution that is resampled 
#' from is inferred by using the Chao1 estimator to infer the number of unseen clones, 
#' followed by applying the relative abundance correction and unseen clone frequencies 
#' described in Chao et al, 2015.
#' 
#' The diversity index (\eqn{D}) for each group is the mean value of over all resampling 
#' realizations. Confidence intervals are derived using the standard deviation of the 
#' resampling realizations, as described in Chao et al, 2015.
#'
#' Significance of the difference in diversity index (\code{D}) between groups is tested by 
#' constructing a bootstrap delta distribution for each pair of unique values in the 
#' \code{group} column. The bootstrap delta distribution is built by subtracting the diversity 
#' index \code{Da} in group \code{a} from the corresponding value \eqn{Db} in group \code{b}, 
#' for all bootstrap realizations, yielding a distribution of \code{nboot} total deltas; where 
#' group \code{a} is the group with the greater mean \code{D}. The p-value for hypothesis 
#' \code{Da  !=  Db} is the value of \code{P(0)} from the empirical cumulative distribution 
#' function of the bootstrap delta distribution, multiplied by 2 for the two-tailed correction.
#' 
#' Note, this method may inflate statistical significance when clone sizes are uniformly small,
#' such as when most clones sizes are 1, sample size is small, and \code{max_n} is near
#' the total count of the smallest data group. Use caution when interpreting the results 
#' in such cases.
#'
#' @examples
#' # Group by sample identifier in two steps
#' abund <- estimateAbundance(ExampleDb, group="sample_id", nboot=100)
#' div <- alphaDiversity(abund, step_q=1, max_q=10)
#' plotDiversityCurve(div, legend_title="Sample")
#'                    
#' # Grouping by isotype rather than sample identifier in one step
#' div <- alphaDiversity(ExampleDb, group="c_call", min_n=40, step_q=1, max_q=10, 
#'                       nboot=100)
#' plotDiversityCurve(div, legend_title="Isotype")
#'
#' @export
alphaDiversity <- function(data, min_q=0, max_q=4, step_q=0.1, ci=0.95, ...) {
    # Hack for visibility of dplyr variables
    . <- NULL
    
    # Check input object and call estimateAbundance if required
    if (is(data, "AbundanceCurve")) {
        abundance <- data
    } else if (is(data, "data.frame")) {
        abundance <- estimateAbundance(data, ci=0.95, ...)
    } else {
        stop("Input must be either a data.frame or AbundanceCurve object.")
    }
    
    # Set diversity orders and confidence interval
    ci_z <- ci + (1 - ci) / 2
    ci_x <- qnorm(ci_z)
    q <- seq(min_q, max_q, step_q)
    if (!(0 %in% q)) { q <- c(0, q) }
    
    # Set grouping variables
    clone <- abundance@clone_by
    group <- abundance@group_by
        
    # Compute diversity metric for bootstrap instances
    boot_df <- abundance@bootstrap %>%
        dplyr::group_by(!!rlang::sym(group)) %>%
        dplyr::do(helperAlpha(., q=q, clone=clone, group=group)) %>%
        dplyr::ungroup()
    
    # Summarize diversity
    div_df <- boot_df %>%
        tidyr::gather(key="n", value="X", -one_of(c(group, "q"))) %>%
        dplyr::mutate(X=as.numeric(!!rlang::sym("X"))) %>%
        dplyr::group_by(!!!rlang::syms(c(group, "q"))) %>%
        dplyr::summarize(d=mean(!!rlang::sym("X"), na.rm=TRUE),
                         d_sd=sd(!!rlang::sym("X"), na.rm=TRUE)) %>%
        dplyr::mutate(d_lower=pmax(!!rlang::sym("d") - !!rlang::sym("d_sd") * ci_x, 0), 
                      d_upper=!!rlang::sym("d") + !!rlang::sym("d_sd") * ci_x)
    
    # Compute evenness
    div_qi <- div_df %>%
        filter(!!rlang::sym("q") == 0) %>%
        select(one_of(c(group, "d")))
    div_df <- div_df %>%
        dplyr::right_join(div_qi, by=group, suffix=c("", "_0")) %>%
        mutate(e=!!rlang::sym("d")/!!rlang::sym("d_0"), 
               e_lower=!!rlang::sym("d_lower")/!!rlang::sym("d_0"), 
               e_upper=!!rlang::sym("d_upper")/!!rlang::sym("d_0")) %>%
        select(-!!rlang::sym("d_0"))
    
    # Test
    if (length(abundance@groups) > 1) {
        test_df <- helperTest(boot_df, q=q, group=group)
    } else {
        test_df <- NULL
    }

    # Build return object
    group_set <- unique(div_df[[group]])
    div_obj <- new("DiversityCurve",
                   diversity=div_df, 
                   tests=test_df,
                   method="alpha",
                   group_by=group,
                   groups=group_set,
                   q=q,  
                   n=abundance@n,
                   ci=ci)
        
   return(div_obj) 
}

# Calculates the pairwise beta diversity
# 
# \code{betaDiversity} takes in a data.frame or \link{AbundanceCurve} and computes
# the multiplicative beta diversity across a range of Hill diversity indices.
# 
# @param    data         data.frame with Change-O style columns containing clonal assignments or
#                        an \link{AbundanceCurve} object generate by \link{estimateAbundance}.
#                        containing a previously calculated bootstrap distributions of clonal abundance.
# @param    comparisons  named list of comparisons between group members for computing beta diversity.
# @param    min_q        minimum value of \eqn{q}.
# @param    max_q        maximum value of \eqn{q}.
# @param    step_q       value by which to increment \eqn{q}.
# @param    ci           confidence interval to calculate; the value must be between 0 and 1.
# @param    ...          additional arguments to pass to \link{estimateAbundance}. Additional arguments
#                        are ignored if a \link{AbundanceCurve} is provided as input.
# 
# @return   A \link{DiversityCurve} object summarizing the diversity scores.
# 
# @details
# Beta diversity or the comparative difference between two samples as quantified using Hill
# diversity indices proposed by Jost (Jost, 2007).
# 
# Briefly, the alpha and gamma diversity components are calculated for each comparison.
# Alpha diversity is calculated as the average hill diversity across each independent sample
# while Gamma diversity is calculated as the total diversity without distinguishing between
# samples. Beta diversity is computed as Gamma/Alpha.
# 
# Diversity is calculated on the estimated clonal abundance distribution with a correction
# for unseen species much like the calculation for alpha diversity \link{alphaDiversity}.
# A smooth curve is generated in the same manner as in \link{alphaDiversity}.
# Confidence intervals are derived using the standard deviation of the resampling realizations.
# 
# \enumerate{
#   \item  Hill M. Diversity and evenness: a unifying notation and its consequences.
#            Ecology. 1973 54(2):427-32.
#   \item  Jost L. Partitioning Diversity Into Independent Alpha and Beta Components.
#            Ecology. 2007 88(10):2427–2439.
#   \item  Jost L, et al. Partitioning diversity for conservation analyses.
#            Diversity Distrib. 2010 16(1):65–76
# }
# 
# @examples
# div <- betaDiversity(ExampleDb, comparisons=list("TIME"=c("-1h", "+7d")), group="sample_id",
#                      min_n=40, step_q=1, max_q=10, nboot=100)
# 
# plotDiversityCurve(div, legend_title="Isotype")
# 
# @export
betaDiversity <- function(data, comparisons, min_q=0, max_q=4, step_q=0.1, ci=0.95, ...) {
    # Hack for visibility of dplyr variables
    . <- NULL
    
    if (!is.list(comparisons) || is.null(names(comparisons))) {
      stop("'comparisons' must be a named list")
    }
    
    # Check input object and call estimateAbundance if required
    if (is(data, "AbundanceCurve")) {
        abundance <- data
    } else if (is(data, "data.frame")) {
        abundance <- estimateAbundance(data, ci=0.95, ...)
    } else {
        stop("Input must be either a data.frame or AbundanceCurve object.")
    }
        
    # Set diversity orders and confidence interval
    ci_z <- ci + (1 - ci) / 2
    ci_x <- qnorm(ci_z)
    q <- seq(min_q, max_q, step_q)
    if (!(0 %in% q)) { q <- c(0, q) }
        
    # Compute pairwise beta diversity for bootstrap instances
    beta_diversity_list <- list()

    for (comparison in names(comparisons)) {
        beta_diversity_list[[comparison]] <- abundance@bootstrap %>%
            dplyr::ungroup() %>%
            dplyr::filter(.[[abundance@group_by]] %in% comparisons[[comparison]]) %>%
            dplyr::do(helperBeta(., q=q, clone=abundance@clone_by, group=abundance@group_by, ci_x=ci_x))
    }

    # Generate summary diversity output
    div_df <- bind_rows(beta_diversity_list, .id = "comparison")
    
    # Beta groups
    group_set <- unique(div_df[["comparison"]])
    
    # Compute evenness
    div_qi <- div_df %>%
        filter(!!rlang::sym("q") == 0) %>%
        select(one_of(c("comparison", "D")))

    div <- div_df %>%
        right_join(div_qi, by = "comparison", suffix = c("", "_0")) %>%
        mutate(d = !!rlang::sym("d")/!!rlang::sym("d_0"), 
            e_lower = !!rlang::sym("d_lower")/!!rlang::sym("d_0"), 
            e_upper = !!rlang::sym("d_upper")/!!rlang::sym("d_0")) %>%
        select(-!!rlang::sym("d_0"))
        
    # Test
    if (length(group_set) > 1) {
       test_df <- helperTest(div_df, q=q, group="comparison")
    } else {
       test_df <- NULL
    }
    
    # Build return object
    div_obj <- new("DiversityCurve",
                   diversity=div, 
                   tests=test_df,
                   method="beta",
                   group_by="comparison",
                   groups=group_set,
                   n=abundance@n,
                   q=q,  
                   ci=ci)
        
    return(div_obj)
}


#### Plotting functions ####

#' Plots a clonal abundance distribution
#' 
#' \code{plotAbundanceCurve} plots the results from estimating the complete clonal 
#' relative abundance distribution. The distribution is plotted as a log rank abundance 
#' distribution.
#' 
#' @param    data          \link{AbundanceCurve} object returned by \link{estimateAbundance}.
#' @param    colors        named character vector whose names are values in the 
#'                         \code{group} column of \code{data} and whose values are 
#'                         colors to assign to those group names.
#' @param    main_title    string specifying the plot title.
#' @param    legend_title  string specifying the legend title.
#' @param    xlim          numeric vector of two values specifying the 
#'                         \code{c(lower, upper)} x-axis limits.
#' @param    ylim          numeric vector of two values specifying the 
#'                         \code{c(lower, upper)} y-axis limits.
#' @param    annotate      string defining whether to added values to the group labels 
#'                         of the legend. When \code{"none"} (default) is specified no
#'                         annotations are added. Specifying (\code{"depth"}) adds 
#'                         sequence counts to the labels.
#' @param    silent        if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                         object; if \code{FALSE} draw the plot.
#' @param    ...           additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  
#' See \link{AbundanceCurve} for the input object and \link{estimateAbundance} for 
#' generating the input abundance distribution.
#' Plotting is performed with \link{ggplot}.
#'           
#' @examples
#' # Estimate abundance by sample and plot
#' abund <- estimateAbundance(ExampleDb, group="sample_id", nboot=100)
#' plotAbundanceCurve(abund, legend_title="Sample")
#' 
#' @export
plotAbundanceCurve <- function(data, colors=NULL, main_title="Rank Abundance", 
                               legend_title=NULL, xlim=NULL, ylim=NULL, 
                               annotate=c("none", "depth"),
                               silent=FALSE, ...) {
    
    # Check if abundance is in data
    if (is.null(data@abundance)) { stop("Missing abundance data.") }
    
    # Check arguments
    annotate <- match.arg(annotate)
    
    # Define group label annotations
    if (all(is.na(data@groups)) || length(data@groups) == 1) {
        group_labels <- NA  
    } else if (annotate == "none") {
        group_labels <- setNames(data@groups, data@groups)
    } else if (annotate == "depth") {
        group_labels <- setNames(paste0(data@groups, " (N=", data@n, ")"),
                                 data@groups)
    }

    # Stupid hack for check NOTE about `.x` in math_format
    .x <- NULL
    
    if (!all(is.na(group_labels))) {
        # Define grouped plot
        p1 <- ggplot(data@abundance, aes_string(x="rank", y="p", group=data@group_by)) + 
            ggtitle(main_title) + 
            baseTheme() + 
            xlab("Rank") +
            ylab("Abundance") +
            scale_x_log10(limits=xlim,
                          breaks=scales::trans_breaks("log10", function(x) 10^x),
                          labels=scales::trans_format("log10", scales::math_format(10^.x))) +
            scale_y_continuous(labels=scales::percent) +
            geom_ribbon(aes_string(ymin="lower", ymax="upper", fill=data@group_by), alpha=0.4) +
            geom_line(aes_string(color=data@group_by))
        
        # Set colors and legend
        if (!is.null(colors)) {
            p1 <- p1 + scale_color_manual(name=legend_title, labels=group_labels, values=colors) +
                scale_fill_manual(name=legend_title, labels=group_labels, values=colors)
        } else {
            p1 <- p1 + scale_color_discrete(name=legend_title, labels=group_labels) +
                scale_fill_discrete(name=legend_title, labels=group_labels)
        }
    } else {
        # Set color
        if (!is.null(colors) & length(colors) == 1) {
            line_color <- colors
        } else {
            line_color <- "black"
        }
        # Define plot
        p1 <- ggplot(data@abundance, aes_string(x="rank", y="p")) + 
            ggtitle(main_title) + 
            baseTheme() + 
            xlab("Rank") +
            ylab("Abundance") +
            scale_x_log10(limits=xlim,
                          breaks=scales::trans_breaks("log10", function(x) 10^x),
                          labels=scales::trans_format("log10", scales::math_format(10^.x))) +
            scale_y_continuous(labels=scales::percent) +
            geom_ribbon(aes_string(ymin="lower", ymax="upper"), fill=line_color, alpha=0.4) +
            geom_line(color=line_color)
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}


#' Plot the results of alphaDiversity
#' 
#' \code{plotDiversityCurve} plots a \code{DiversityCurve} object.
#'
#' @param    data            \link{DiversityCurve} object returned by 
#'                           \link{alphaDiversity}.
#' @param    colors          named character vector whose names are values in the 
#'                           \code{group} column of the \code{data} slot of \code{data},
#'                           and whose values are colors to assign to those group names.
#' @param    main_title      string specifying the plot title.
#' @param    legend_title    string specifying the legend title.
#' @param    log_x           if \code{TRUE} then plot \eqn{q} on a log scale;
#'                           if \code{FALSE} plot on a linear scale.
#' @param    log_y           if \code{TRUE} then plot the diversity/evenness scores 
#'                           on a log scale; if \code{FALSE} plot on a linear scale.
#' @param    xlim            numeric vector of two values specifying the 
#'                           \code{c(lower, upper)} x-axis limits.
#' @param    ylim            numeric vector of two values specifying the 
#'                           \code{c(lower, upper)} y-axis limits.
#' @param    annotate        string defining whether to added values to the group labels 
#'                           of the legend. When \code{"none"} (default) is specified no
#'                           annotations are added. Specifying (\code{"depth"}) adds 
#'                           sequence counts to the labels.
#' @param    score           one of \code{"diversity"} or \code{"evenness"} specifying which
#'                           score to plot on the y-asis.
#' @param    silent          if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                           object; if \code{FALSE} draw the plot.
#' @param    ...             additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  See \link{alphaDiversity} and \link{alphaDiversity} for generating 
#'           \link{DiversityCurve} objects for input. Plotting is performed with \link{ggplot}.
#' 
#' @examples
#' # Calculate diversity
#' div <- alphaDiversity(ExampleDb, group="sample_id", nboot=100)
#' 
#' # Plot diversity
#' plotDiversityCurve(div, legend_title="Sample")
#'
#' #' # Plot diversity
#' plotDiversityCurve(div, legend_title="Sample", score="evenness")
#' 
#' @export
plotDiversityCurve <- function(data, colors=NULL, main_title="Diversity", 
                               legend_title="Group", log_x=FALSE, log_y=FALSE,
                               xlim=NULL, ylim=NULL, annotate=c("none", "depth"), 
                               score=c("diversity", "evenness"),
                               silent=FALSE, ...) {
    # Check arguments
    annotate <- match.arg(annotate)
    score <- match.arg(score)
    
    # Define group label annotations
    if (all(is.na(data@groups)) || length(data@groups) == 1) {
        group_labels <- NA  
    } else if (annotate == "none") {
        group_labels <- setNames(data@groups, data@groups)
    } else if (annotate == "depth") {
        group_labels <- setNames(paste0(data@groups, " (n=", data@n, ")"),
                                 data@groups)
    }
    
    # Define y-axis scores
    if (score == "diversity") {
        y_value <- "d"
        y_min <- "d_lower"
        y_max <- "d_upper"
        y_label <- expression(''^q * D)
    } else if (score == "evenness") {
        y_value <- "e"
        y_min <- "e_lower"
        y_max <- "e_upper"
        y_label <- expression(''^q * e)
    }
    
    # Stupid hack for check NOTE about `.x` in math_format
    .x <- NULL
    
    if (!all(is.na(group_labels))) {
        # Define grouped plot
        p1 <- ggplot(data@diversity, aes_string(x="q", y=y_value, group=data@group_by)) + 
            ggtitle(main_title) + 
            baseTheme() + 
            xlab('q') +
            ylab(y_label) +
            geom_ribbon(aes_string(ymin=y_min, ymax=y_max, fill=data@group_by), alpha=0.4) +
            geom_line(aes_string(color=data@group_by))
    
        # Set colors and legend
        if (!is.null(colors)) {
            p1 <- p1 + scale_color_manual(name=legend_title, labels=group_labels, values=colors) +
                scale_fill_manual(name=legend_title, labels=group_labels, values=colors)
        } else {
            p1 <- p1 + scale_color_discrete(name=legend_title, labels=group_labels) +
                scale_fill_discrete(name=legend_title, labels=group_labels)
        }
    } else {
        # Set color
        if (!is.null(colors) & length(colors) == 1) {
          line_color <- colors
        } else {
          line_color <- "black"
        }
      
        # Define ungrouped plot
        p1 <- ggplot(data@diversity, aes_string(x="q", y=y_value)) + 
            ggtitle(main_title) + 
            baseTheme() + 
            xlab('q') +
            ylab(y_label) +
            geom_ribbon(aes_string(ymin=y_min, ymax=y_max), fill=line_color, alpha=0.4) +
            geom_line(color=line_color)
    }
    
    # Set x-axis style
    if (log_x) {
        p1 <- p1 + scale_x_continuous(trans=scales::log2_trans(), limits=xlim,
                                      breaks=scales::trans_breaks('log2', function(x) 2^x),
                                      labels=scales::trans_format('log2', scales::math_format(2^.x)))
    } else {
        p1 <- p1 + scale_x_continuous(limits=xlim)
    }
    
    # Set y-axis style
    if (log_y) {
        p1 <- p1 + scale_y_continuous(trans=scales::log2_trans(), limits=ylim,
                                      breaks=scales::trans_breaks('log2', function(x) 2^x),
                                      labels=scales::trans_format('log2', scales::math_format(2^.x)))
    } else {
        p1 <- p1 + scale_y_continuous(limits=ylim)
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))

    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}


#' Plot the results of diversity testing
#' 
#' \code{plotDiversityTest} plots summary data for a \code{DiversityCurve} object 
#' with mean and a line range indicating plus/minus one standard deviation.
#'
#' @param    data            \link{DiversityCurve} object returned by 
#'                           \link{alphaDiversity}.
#' @param    q               diversity order to plot the test for.
#' @param    colors          named character vector whose names are values in the 
#'                           \code{group} column of the \code{data} slot of \code{data},
#'                           and whose values are colors to assign to those group names.
#' @param    main_title      string specifying the plot title.
#' @param    legend_title    string specifying the legend title.
#' @param    log_d           if \code{TRUE} then plot the diversity scores \eqn{D} 
#'                           on a log scale; if \code{FALSE} plot on a linear scale.
#' @param    annotate        string defining whether to added values to the group labels 
#'                           of the legend. When \code{"none"} (default) is specified no
#'                           annotations are added. Specifying (\code{"depth"}) adds 
#'                           sequence counts to the labels.
#' @param    silent          if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                           object; if \code{FALSE} draw the plot.
#' @param    ...             additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  See \link{alphaDiversity} for generating input.
#'           Plotting is performed with \link{ggplot}.
#' 
#' @examples
#' # Calculate diversity
#' div <- alphaDiversity(ExampleDb, group="sample_id", min_q=0, max_q=2, step_q=1, nboot=100)
#' 
#' # Plot results at q=0 (equivalent to species richness)
#' plotDiversityTest(div, 0, legend_title="Sample")
#' 
#' # Plot results at q=2 (equivalent to Simpson's index)
#' plotDiversityTest(div, q=2, legend_title="Sample")
#' 
#' @export
plotDiversityTest <- function(data, q, colors=NULL, main_title="Diversity", legend_title="Group", 
                              log_d=FALSE, annotate=c("none", "depth"), silent=FALSE, ...) {
    # Stupid hack for check NOTE about `.x` in math_format
    .x <- NULL
    
    # Check arguments
    annotate <- match.arg(annotate)
    
    # Check if abundance is in data
    if (is.null(data@tests)) { 
        stop("Test data missing from input object.")
    }

    # Check if q is in data
    if (!(q %in% data@q)) { 
        stop("Test for order q=", q, " not found in input object.") 
    }

    # Define group label annotations
    if (annotate == "none") {
        group_labels <- setNames(data@groups, data@groups)
    } else if (annotate == "depth") {
        group_labels <- setNames(paste0(data@groups, " (N=", data@n, ")"),
                                 data@groups)
    }
    # Define plot values
    df <- data@diversity %>%
        dplyr::filter(!!rlang::sym("q") == !!rlang::enquo(q)) %>%
        dplyr::mutate(lower=!!rlang::sym("d") - !!rlang::sym("d_sd"), 
                      upper=!!rlang::sym("d") + !!rlang::sym("d_sd"))
    
    # Define base plot elements
    p1 <- ggplot(df, aes_string(x=data@group_by)) + 
        ggtitle(main_title) + 
        baseTheme() + 
        xlab("") +
        ylab(bquote("Mean " ^ .(q) * D %+-% "SD")) +
        geom_linerange(aes_string(ymin="lower", ymax="upper", color=data@group_by), alpha=0.8) +
        geom_point(aes_string(y="d", color=data@group_by))
    
    # Set colors and legend
    if (!is.null(colors)) {
        p1 <- p1 + scale_color_manual(name=legend_title, labels=group_labels, values=colors)
    } else {
        p1 <- p1 + scale_color_discrete(name=legend_title, labels=group_labels)
    }

    # Set x-axis style
    if (log_d) {
        p1 <- p1 + scale_y_continuous(trans=scales::log2_trans(),
                                      breaks=scales::trans_breaks('log2', function(x) 2^x),
                                      labels=scales::trans_format('log2', scales::math_format(2^.x)))
    } else {
        p1 <- p1 + scale_y_continuous()
    }

    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}
