#' Calculate expected leaves of a species tree
#'
#' @details Calculates the expected number of leaves for a birth-death
#'    simulation given a speciation and extinction rate and a time.
#'
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param t time to simulate until
#' @return The expected number of leaves
#' @references
#' Mooers, A., Gascuel, O., Stadler, T., Li, H., & Steel, M. (2012).
#'     Branch lengths on birth–death trees and the expected loss of phylogenetic
#'      diversity. Systematic biology, 61(2), 195-203.
#' @examples
#' spec_rate <- 1.0
#' ext_rate <- 0.5
#' time <- 2
#' calculate_expected_leaves_sptree(spec_rate, ext_rate, time)
#' @export
calculate_expected_leaves_sptree <- function(lambda,
                                      mu,
                                      t) {
    if (!is.numeric(lambda)) {
        stop("'lambda' needs to be a number")
    }
    if (lambda < 0) {
        stop("'lambda' needs to be positive")
    }
    if (!is.numeric(t)) {
        stop("'t' needs to be a number")
    }
    if (t <= 0) {
        stop("'t' needs to be positive")
    }
    if (!is.numeric(mu)) {
        stop("'mu' needs to be a number")
    }
    if (mu < 0) {
        stop("'mu' needs to be positive or 0.0")
    }
    2 * exp((lambda - mu) * t)
}
#' Calculate expected leaves of a locus tree
#'
#' @details Calculates the expected number of leaves for a birth-death
#'    simulation given a gene birth and death rate, a time, and the number of
#'    leaves on the species tree that the locus tree is to be simulated upon.
#'
#' @param t time to simulate until (the length of the species tree)
#' @param dup_rate gene birth rate
#' @param loss_rate gene death rate
#' @param num_species number of leaves on the species tree
#' @return The expected number of leaves
#' @references
#' Mallo, D., de Oliveira Martins, L., & Posada, D. (2016). SimPhy: phylogenomic
#'      simulation of gene, locus, and species trees. Systematic biology, 65(2),
#'       334-344.

#' @examples
#' gene_birth_rate <- 1.0
#' gene_death_rate <- 0.5
#' time <- 2
#' num_species <- 10
#' calculate_expected_leaves_locustree(time,
#'                                     gene_birth_rate,
#'                                     gene_death_rate,
#'                                     num_species)
#' @export
calculate_expected_leaves_locustree <- function(t,
                                                dup_rate,
                                                loss_rate,
                                                num_species) {
    if (!is.numeric(t)) {
        stop("'t' needs to be a number")
    }
    if (t <= 0) {
        stop("'t' needs to be positive")
    }
    if (!is.numeric(dup_rate)) {
        stop("'dup_rate' needs to be a number")
    }
    if (dup_rate < 0) {
        stop("'dup_rate' needs to be positive or 0.0.")
    }
    if (!is.numeric(loss_rate)) {
        stop("'loss_rate' needs to be a number")
    }
    if (loss_rate < 0) {
        stop("'loss_rate' needs to be positive or 0.0.")
    }
    f_numer <- (loss_rate * (1 - exp(- (dup_rate - loss_rate) * t)))
    f_denom <- dup_rate - loss_rate * (exp(- (dup_rate - loss_rate) * t))
    f <- f_numer / f_denom
    (num_species * exp(t * (dup_rate - loss_rate))) /
        (1 - (f ^ (num_species - 2)))
}
#' Calculate expected time to branching point of a species tree
#'
#' @description Calculates the expected time to branching point of a
#'    species tree for a birth-death simulation given a speciation
#'    and extinction rate and a number of leaves, and a branching point.
#'
#'
#' @details By default this branching point is 1 which corresponds
#' to the root, k = 2 corresponds to the first branching point after
#' the root, k = 3 the second, and so on. For more details see Gernhard 2008.
#'
#'
#'
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param n number of tips on tree
#' @param k branching point (k = 1 is the root and is the default)
#' @return The expected branching time
#' @references
#' Gernhard, T. (2008). The conditioned reconstructed process. Journal of
#'     theoretical biology, 253(4), 769-778.
#' @examples
#' spec_rate <- 1.0
#' ext_rate <- 0.5
#' nt <- 10
#' estimate_node_heights(lambda = spec_rate, mu = ext_rate, n = nt)
#'
#' estimate_node_heights(lambda = spec_rate, mu = ext_rate, n = nt, k = 2)
#' @export
estimate_node_heights <- function(lambda,
                                  mu,
                                  n,
                                  k = 1) {
    if (!is.numeric(lambda))
        stop("'lambda' needs to be a number")
    if (!is.numeric(mu))
        stop("'mu' needs to be a number")
    if (!is.numeric(n))
        stop("'n' needs to be a number")
    if (!is.numeric(k))
        stop("'k' needs to be a number")
    if (n < 1)
        stop("'n' must be greater than or equal to 1.")
    if (lambda <= 0.0)
        stop("'lambda' is less than or equal to 0, it must be greater than 0.")
    if (mu > lambda)
        stop("'mu' is greater than 'lambda', please correct this.")
    if (mu < 0.0)
        stop("'mu' cannot be less than 0.0.")
    if (mu == 0.0) {
        s <- 0
        for (i in k + 1:n) {
            s <- s + 1 / (lambda * i)
        }
        speciation_event_time <- s
    }
    else if (mu == lambda) {
        speciation_event_time <- (n - k) / (lambda * k)
    }
    else{
        rho <- mu / lambda
        first_part <- ((k + 1) / lambda) * exp(lchoose(n, k + 1)) * (-1)^k
        out_sum <- 0
        for (i in 0:(n - k - 1)) {
            s1 <- exp(lchoose(n - k - 1, i))
            s2 <- 1 / ((k + i + 1) * rho)
            s3 <- exp((k + i) * log(1 / rho - 1))
            in_sum <- 0
            for (j in 1:k + i) {
                ins1 <- exp(lchoose(k + i, j))
                ins2 <- (-1^j) / j
                ins3 <- 1 - (1 / (1 - rho))^j
                in_sum <- in_sum + ins1 * ins2 * ins3
            }
            s4 <- log(1 / (1 - rho)) - in_sum
            out_sum <- out_sum + s1 * s2 * s3 * s4
        }
        speciation_event_time <- first_part * out_sum
    }
    speciation_event_time
}