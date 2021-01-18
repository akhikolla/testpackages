outlier_model <- function(A,
                          adj,
                          nsim       = 10000,
                          ncores     = 1,
                          validate   = TRUE
                          ) {
  stopifnot(is.matrix(A))
  if (inherits(adj, "gengraph")) adj <- ess::adj_lst(adj)
  if (validate ) {
    if (any(is.na(A))) message("  Note: A has NA values. These have been treated as ordinay values.")
    if( !only_single_chars(A)) stop("All values in A must be represented as a single character. Use to_single_chars(A)")
  }
  RIP   <- ess::rip(adj) # the rip (or actually mcs) will check for decomposability here
  cms   <- a_marginals(A, RIP$C)
  sms   <- a_marginals(A, RIP$S)
  sims  <- .sim_internal(A, cms, sms, nsim = nsim, type = "deviance", ncores = ncores)
  mu    <- mean(sims)
  sigma <- stats::var(sims)
  cdf   <- stats::ecdf(sims)
  return(new_outlier_model(A, sims, mu, sigma, cdf, cms, sms))
}

#' Outlier detection
#'
#' Detecting outliers within a dataset or test if a new (novel) observation is an outlier.
#'
#' @param A Character matrix or data.frame. All values must be limited to a single character.
#' @param adj Adjacency list or \code{gengraph} object of a decomposable graph.
#' See package \code{ess} for \code{gengraph} objects.
#' @param z Named vector (same names as \code{colnames(A)}) or \code{NULL}. See details.
#' Values must be limited to a single character.
#' @param alpha Significance level
#' @param nsim Number of simulations
#' @param ncores Number of cores to use in parallelization
#' @param validate Logical. If true, it checks if \code{A} only has single character
#' values and converts it if not.
#' @return A \code{outlier_model} object with either \code{novelty} or \code{outlier}
#' as child classes. These are used for different purposes. See the details
#' @details If the goal is to detect outliers within \code{A} set \code{z} to \code{NULL};
#' this procedure is most often just referred to as outlier detection. Once \code{fit_outlier}
#' has been called in this situation, one can exploit the \code{outliers} function to get the
#' indicies for which observations in \code{A} that are outliers. See the examples.
#'
#' On the other hand, if the goal is test if the new unseen observation \code{z} is an outlier
#' in\code{A}, then supply a named vector to \code{z}.
#'
#' All values must be limited to a single character representation; if not, the function will
#' internally convert to one such representation. The reason for this, is a speedup in runtime
#' performance. One can also use the exported function \code{to_single_chars} on \code{A} in
#' advance and set \code{validate} to \code{FALSE}. 
#'
#' The \code{adj} object is most typically found using \code{fit_graph} from the \code{ess}
#' package. But the user can supply an adjacency list, just a named \code{list}, of their own
#' choice if needed.
#' @seealso \code{\link{fit_mixed_outlier}}, \code{\link{fit_multiple_models}},
#' \code{\link{outliers}}, \code{\link{pval}}, \code{\link{deviance}}
#' @examples
#'
#' library(dplyr)
#' library(ess)  # For the fit_graph function
#' set.seed(7)   # For reproducibility
#' 
#' # Psoriasis patients
#' d <- derma %>%
#'   filter(ES == "psoriasis") %>%
#'   select(1:20) %>% # only a subset of data is used to exemplify
#'   as_tibble()
#'
#' # Fitting the interaction graph
#' # see package ess for details
#' g <- fit_graph(d, trace = FALSE) 
#' plot(g)
#'
#' # -----------------------------------------------------------
#' #                        EXAMPLE 1
#' #    Testing which observations within d are outliers
#' # -----------------------------------------------------------
#'
#' # Only 500 simulations is used here to exeplify
#' # The default number of simulations is 10,000
#' m1 <- fit_outlier(d, g, nsim = 500)
#' print(m1)
#' outs  <- outliers(m1)
#' douts <- d[which(outs), ]
#' douts
#'
#' # Notice that m1 is of class 'outlier'. This means, that the procedure has tested which
#' # observations _within_ the data are outliers. This method is most often just referred to
#' # as outlier detection. The following plot is the distribution of the test statistic. Think
#' # of a simple t-test, where the distribution of the test statistic is a t-distribution.
#' # In order to conclude on the hypothesis, one finds the critical value and verify if the
#' # test statistic is greater or less than this.
#'
#' # Retrieving the test statistic for individual observations
#' x1 <- douts[1, ] %>% unlist()
#' x2 <- d[1, ] %>% unlist()
#' dev1 <- deviance(m1, x1) # falls within the critical region in the plot (the red area)
#' dev2 <- deviance(m1, x2) # falls within the acceptable region in the plot
#'
#' dev1
#' dev2
#'
#' # Retrieving the pvalues
#' pval(m1, dev1)
#' pval(m1, dev2)
#' 
#' # -----------------------------------------------------------
#' #                        EXAMPLE 2
#' #         Testing if a new observation is an outlier
#' # -----------------------------------------------------------
#' 
#' # An observation from class "chronic dermatitis"
#' z <- derma %>%
#'   filter(ES == "chronic dermatitis") %>%
#'   select(1:20) %>%
#'   slice(1) %>%
#'   unlist()
#' 
#' # Test if z is an outlier in class "psoriasis"
#' # Only 500 simulations is used here to exeplify
#' # The default number of simulations is 10,000
#' m2 <- fit_outlier(d, g, z, nsim = 500)
#' print(m2)
#' plot(m2) # Try using more simulations and the complete derma data
#'
#' # Notice that m2 is of class 'novelty'. The term novelty detection
#' # is sometimes used in the litterature when the goal is to verify
#' # if a new unseen observation is an outlier in a homogen dataset.
#'
#' # Retrieving the test statistic and pvalue for z
#' dz <- deviance(m2, z)
#' pval(m2, dz)
#'
#' @export
fit_outlier <- function(A,
                        adj,
                        z        = NULL,
                        alpha    = 0.05,
                        nsim     = 10000,
                        ncores   = 1,
                        validate = TRUE) {

  # TODO: z may now contain NAs since we are armed with the junction max-flow alg.
  # - in this case, we impute the NAs given the observed (evidence)

  if (!(is.data.frame(A) || is.matrix(A))) stop("A must be either a matrix or a data.frame", call. = FALSE)

  if (any(is.na(A))) stop("A has NA values.")
  
  novelty_detection <- !is.null(z)
  
  if (novelty_detection) {
    if (!identical(colnames(A), names(z))) {
      stop("Variables in A and the names of z is not in agreement!")
    }
  }

  if (inherits(adj, "gengraph")) adj <- ess::adj_lst(adj)
  
  if (is.data.frame(A)) A <- as.matrix(A)

  Az <- if (novelty_detection) rbind(A, z) else A

  if (validate) {
    if( !only_single_chars(Az) ) {
      message("A has values longer than a single character. to_single_chars() was used to correct this.")
      Az <- to_single_chars(Az)
      if (novelty_detection) z <- Az[nrow(Az), ]
    }    
  }
  
  m <- outlier_model(Az, adj, nsim = nsim, ncores = ncores, validate = FALSE)

  m <- if (novelty_detection) {
    dev_z <- deviance(m, z)
    new_novelty(m, dev_z, pval(m, dev_z), critval(m, alpha), alpha)
  } else {
    new_outlier(m, critval(m, alpha), alpha)
  }
  
  return(m)
}


#' Mixed Outlier Test
#'
#' A function for outlier detection with mixed, but independen, information
#'
#' @param m1 An object returned from \code{fit_outlier}
#' @param m2 An object returned from \code{fit_outlier}
#' @details It is assumed that the input data to \code{m1} and \code{m2}
#' holds information about the same observation in corresponding rows.
#' Thus, the two datasets must also be of same dimension.
#' @return An object of type \code{mixed_outlier} with \code{novelty} or \code{outlier}
#' as child classes. These are used for different purposes. See \code{fit_outlier}.
#' @seealso \code{\link{fit_outlier}}, \code{\link{fit_multiple_models}},
#' \code{\link{outliers}}, \code{\link{pval}}, \code{\link{deviance}}
#' @examples
#'
#' library(dplyr)
#' library(ess)  # for fit_components
#' set.seed(7)   # for reproducibility
#'
#' ## Data
#'
#' # The components - here microhaplotypes
#' haps <- tgp_haps[1:5] # only a subset of data is used to exemplify
#' dat <- tgp_dat %>%
#'        select(pop_meta, sample_name, all_of(unname(unlist(haps))))
#' 
#' # All the Europeans
#' eur <- dat %>%
#'   as_tibble() %>%
#'   filter(pop_meta == "EUR")
#' 
#' # Extracting the two databases for each copy of the chromosomes
#' eur_a <- eur %>%
#'   filter(grepl("a$", sample_name)) %>%
#'   select(-c(1:2))
#' 
#' eur_b <- eur %>%
#'   filter(grepl("b$", sample_name)) %>%
#'   select(-c(1:2))
#'
#'
#' # Fitting the interaction graphs on the EUR data 
#' ga <- fit_components(eur_a, comp = haps, trace = FALSE)
#' gb <- fit_components(eur_b, comp = haps, trace = FALSE)
#' print(ga)
#' plot(ga, vertex.size = 1)
#' 
#' ## ---------------------------------------------------------
#' ##                       EXAMPLE 1
#' ##   Testing which observations within data are outliers
#' ## ---------------------------------------------------------
#'
#' # Only 500 simulations is used here to exeplify
#' # The default number of simulations is 10,000
#' m1 <- fit_outlier(eur_a, ga, nsim = 500) # consider using more cores (ncores argument)
#' m2 <- fit_outlier(eur_b, gb, nsim = 500) # consider using more cores (ncores argument)
#' m  <- fit_mixed_outlier(m1, m2)
#' print(m)
#' plot(m)
#'
#' outs <- outliers(m)
#' eur_a_outs <- eur_a[which(outs), ]
#' eur_b_outs <- eur_b[which(outs), ]
#'
#' # Retrieving the test statistic for individual observations
#' x1 <- rbind(eur_a_outs[1, ], eur_b_outs[1, ])
#' x2 <- rbind(eur_a[1, ], eur_b[1, ])
#' dev1 <- deviance(m, x1) # falls within the critical region in the plot (the red area)
#' dev2 <- deviance(m, x2) # falls within the acceptable region in the plot
#'
#' dev1
#' dev2
#'
#' # Retrieving the pvalues
#' pval(m, dev1)
#' pval(m, dev2)
#' 
#'
#' \donttest{
#'  
#' ## ---------------------------------------------------------
#' ##                       EXAMPLE 2
#' ##      Testing if a new observation is an outlier
#' ## ---------------------------------------------------------
#' 
#' # Testing if an American is an outlier in Europe
#' amr <- dat %>%
#'   as_tibble() %>%
#'   filter(pop_meta == "AMR")
#' 
#' z1  <- amr %>%
#'   filter(grepl("a$", sample_name)) %>% 
#'   select(unname(unlist(haps))) %>%
#'   slice(1) %>%
#'   unlist()
#' 
#' z2  <- amr %>%
#'   filter(grepl("b$", sample_name)) %>% 
#'   select(unname(unlist(haps))) %>%
#'   slice(1) %>%
#'   unlist()
#'
#' # Only 500 simulations is used here to exemplify
#' # The default number of simulations is 10,000
#' m3 <- fit_outlier(eur_a, ga, z1, nsim = 500) # consider using more cores (ncores argument)
#' m4 <- fit_outlier(eur_b, gb, z2, nsim = 500) # consider using more cores (ncores argument)
#' m5 <- fit_mixed_outlier(m3, m4)
#' print(m5)
#' plot(m5)
#'
#' }
#' 
#' @export
fit_mixed_outlier <- function(m1, m2) {
  if (!identical(class(m1), class(m2))) stop("m1 and m2 are different models")
  if (!identical(dim(m1$A), dim(m2$A))) stop("m1 and m2 was generated from data with different dimensixons")
  m        <- convolute(m1, m2)
  type     <- ifelse(inherits(m1, "novelty"), "novelty", "outlier")
  dev      <- if (type == "novelty") m1$dev + m2$dev else NULL
  new_mixed_outlier(m, type, m1$alpha, dev)
}


#' Fit Multiple Models
#'
#' Conduct multiple novelty tests for a new observation
#'
#' @param A A character matrix or data.frame
#' @param z Named vector (same names as \code{colnames(A)} but without the class variable)
#' @param response A character with the name of the class variable of interest
#' @param alpha The significance level
#' @param type Character ("fwd", "bwd", "tree" or "tfwd") - the type of interaction graph to be used
#' @param q Penalty term in the stopping criterion when fitting the interaction graph (\code{0} = AIC and \code{1} = BIC)
#' @param comp A list with character vectors. Each element in the list is a component in the graph (using expert knowledge)
#' @param nsim Number of simulations
#' @param ncores Number of cores to use in parallelization
#' @param trace Logical indicating whether or not to trace the procedure
#' @param validate Logical. If true, it checks if \code{A} has only single character values and converts it if not.
#' @return An object of type \code{multiple_models}; a list of of \code{novely} objects from which one
#' can query pvalues etc. for outlierdetection.
#' @seealso \code{\link{fit_outlier}}, \code{\link{fit_mixed_outlier}}
#' @examples
#'
#' library(dplyr)
#' set.seed(1)
#' 
#' # A patient with psoriasis
#' z <- unlist(derma[2, 1:10])
#' 
#' d <- derma[, c(names(z), "ES")] %>%
#'      filter(ES %in% c("chronic dermatitis", "psoriasis"))
#' 
#' m <- fit_multiple_models(d, z, "ES", nsim = 1000, trace = FALSE, validate = FALSE)
#'
#' plot(m)
#' print(m)
#'
#' @export
fit_multiple_models <- function(A,
                                z,
                                response,
                                alpha      = 0.05,
                                type       = "fwd",
                                q          = 0.5,
                                comp       = NULL,
                                nsim       = 10000,
                                ncores     = 1,
                                trace      = TRUE,
                                validate   = TRUE) {


  if (!(is.data.frame(A) || is.matrix(A))) stop("A must be either a matrix or a data.frame", call. = FALSE)
  if (is.matrix(A)) A <- as.data.frame(A)

  res_vec  <- A[, response, drop = TRUE]
  res_lvls <- unique(res_vec)
  
  models <- lapply(seq_along(res_lvls), function(i) {

    if (trace) cat(i, "/", length(res_lvls), " ... \n")

    Ai <- A[res_vec == res_lvls[i], -which(colnames(A) == response)]

    if (!is.null(comp)) {
      gi <- ess::fit_components(Ai, comp = comp, type = type, q = q, trace = FALSE)
    } else {
      gi <- ess::fit_graph(Ai, type = type, q = q, trace = FALSE)
    }

    fit_outlier(A = Ai,
      adj         = gi,
      z           = z,
      alpha       = alpha,
      nsim        = nsim,
      ncores      = ncores,
      validate    = validate)
    
  })
  
  names(models) <- res_lvls
  structure(models, class = c("multiple_models", class(models)))
}
