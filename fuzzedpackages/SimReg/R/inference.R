#' Get full set of terms to use in inference procedure based on similarity function arguments
#'
#' @param args Named list of named arguments which gets passed to ontological similarity function by \code{sim_reg}.
#' @return Character vector of term IDs.
#' @export
get_terms <- function(args) {
	if (any(names(args) == "term_sim_mat")) colnames(args[["term_sim_mat"]])
	else if (any(names(args) == "information_content")) names(args[["information_content"]])
	else stop("Either 'information_content' or 'term_sim_mat' is required for similarity function arguments")
}

#' Similarity regression
#'
#' Performs Bayesian `similarity regression' on given \code{logical} response vector \code{y} against \code{list} of ontological term sets \code{x}. It returns an object of class \code{sim_reg_output}. Of particular interest are the probability of an association, which can be calculated with \code{\link{prob_association}}, and the characteristic ontological profile phi, which can be visualised using the functions \code{\link{plot_term_marginals}}, and \code{\link{term_marginals}}). The results can be summarised with \code{summary}.
#'
#' @template ontology
#' @param x \code{list} of \code{character} vectors of ontological terms.
#' @param y \code{logical} response vector.
#' @export
#' @param information_content Numeric vector of information contents of terms named by term ID. Defaults to information content based on frequencies of annotation in \code{x}.
#' @param sim_params List of arguments to pass to \code{get_asym_sim_grid}.
#' @param using_terms Character vector of term IDs giving the complete set of terms to include in the the \code{phi} parameter space.
#' @param term_weights Numeric vector of prior weights for individual terms.
#' @param prior Function for computing the unweighted prior probability of a \code{phi} value.
#' @param min_BF Bayes factor threshold below which to terminate computation, enabling faster execution time at the expense of accuracy and precision.
#' @param max_select Upper bound for number of \code{phi} values to sample.
#' @param max_phi_count Upper bound for number of \code{phi} values to include in final likelihood sum.
#' @param two_way Boolean value determining whether to calculate semantic similarity `in both directions' (i.e. compute \code{s_x} and \code{s_phi} or just \code{s_phi}).
#' @param selection_fn Function for selecting values of \code{phi} with high posterior mass.
#' @param lik_method Function for calculating marginal likelihood contional on values of \code{phi}.
#' @param lik_method_args List of additional arguments to pass to \code{lik_method}.
#' @param gamma0_ml Function for computing marginal likelihood of data under baseline model \code{gamma=0}.
#' @param min_ratio Lower bound on ratio below which to discard \code{phi} values.
#' @param ... Additional arguments to pass to \code{selection_fn}.
#' @importFrom ontologyIndex get_ancestors get_term_info_content
#' @importFrom ontologySimilarity get_asym_sim_grid
#' @examples
#' \dontrun{
#' set.seed(0)
#' data(hpo)
#' disease_terms <- c("HP:0005537", "HP:0000729", "HP:0001873")
#' all_terms <- get_ancestors(hpo, 
#'	c(disease_terms, sample(hpo$id, size=50)))
#' y <- c(rep(FALSE, 96), rep(TRUE, 3))
#' x <- lapply(y, function(.y) minimal_set(
#'	hpo, if (!.y) sample(all_terms, size=3) else 
#'		c(sample(all_terms, size=1), disease_terms[runif(n=3) < 0.8])))
#' sim_reg_out <- sim_reg(ontology=hpo, x=x, y=y)
#' }
sim_reg <- function(
	ontology,
	x,
	y,
	information_content=get_term_info_content(ontology, x),
	sim_params=list(ontology=ontology, information_content=information_content),
	using_terms=get_terms(sim_params),
	term_weights=rep(0, length(using_terms)),
	prior=discrete_gamma(using_terms),
	min_BF=-Inf,
	max_select=2000L,
	max_phi_count=200L,
	two_way=TRUE,
	selection_fn=fg_step,
	lik_method=NULL,
	lik_method_args=list(),
	gamma0_ml=bg_rate,
	min_ratio=1e-4,
	...
) {
	stopifnot(sum(y) > 0)
	stopifnot(all(unique(unlist(use.names=FALSE, x)) %in% using_terms))
	
	prior_weights <- if (!is.null(names(term_weights))) {
		stopifnot(all(using_terms %in% names(term_weights)))
		term_weights[using_terms]
	} else {
		stopifnot(length(using_terms) == length(term_weights))
		term_weights
	}

	bgml <- gamma0_ml(y)

	phi_roots <- Filter(f=function(x) length(x) > 0, x=get_phi_roots(ontology, information_content=information_content, min_intersect=1L, term_sets=x[y]))

	f_priors <- function(phis) if (length(phis) == 0) numeric(0) else { sapply(phis, prior) + sapply(phis, function(phi) mean(prior_weights[match(phi, using_terms)])-mean(prior_weights)) }

	ml <- function(phis) {
		if (length(phis) == 0) { numeric(0) } else {
			if (two_way) {
				s_phis <- split(do.call(what=get_asym_sim_grid, c(sim_params, list(B=x, A=phis))), seq(length(phis)))
				s_xs <- split(t(do.call(what=get_asym_sim_grid, c(sim_params, list(A=x, B=phis)))), seq(length(phis)))
				mapply(FUN=function(s_phi, s_x) selection_fn(s_phi=s_phi, s_x=s_x, y=y, ...), s_phis, s_xs)
			} else {
				mapply(FUN=function(x, pr) selection_fn(x=x, y=y, ...), split(do.call(what=get_asym_sim_grid, c(sim_params, list(B=x, A=phis))), seq(length(phis))))
			}
		}
	}

	pp_ml <- function(phis, phi_priors) {
		if (length(phis) == 0) { numeric(0) } else {
			best <- -Inf
			out <- numeric(length(phis))
			if (two_way) {
				s_phis <- split(do.call(what=get_asym_sim_grid, c(sim_params, list(B=x, A=phis))), seq(length(phis)))
				s_xs <- split(t(do.call(what=get_asym_sim_grid, c(sim_params, list(A=x, B=phis)))), seq(length(phis)))
				for (i in seq_along(phis)) { 
					out[i] <- do.call(what=lik_method, c(list(s_phi=s_phis[[i]], s_x=s_xs[[i]], y=y, min_log_ML=log(min_ratio) + best - phi_priors[i]), lik_method_args))
					best <- max(out[i] + phi_priors[i], best)
				}
			} else {
				s_phis <- split(do.call(what=get_asym_sim_grid, c(sim_params, list(B=x, A=phis))), seq(length(phis)))
				for (i in seq_along(phis)) {
					out[i] <- do.call(what=lik_method, c(list(x=s_phis[[i]], y=y, min_log_ML=log(min_ratio) + best - phi_priors[i]), lik_method_args))
					best <- max(out[i] + phi_priors[i], best)
				}
			}
			out
		}
	}

	priors <- f_priors(phi_roots)
	mls <- ml(phi_roots)

	best_root <- phi_roots[[which.max(mls)]]

	ont_roots <- ontology$id[!ontology$obsolete & sapply(ontology$parents, length) == 0]
	core_phis <- c(unname(as.list(ont_roots)), Filter(f=Negate(is.null), x=(function(terms_in, terms_out) {
		if (length(terms_out) == 0) { terms_in }
		else { 
			cand_phis <- lapply(terms_out, c, terms_in[[length(terms_in)]])
			best <- which.max(ml(cand_phis) + f_priors(cand_phis))
			sys.function(0)(c(terms_in, cand_phis[best]), terms_out[-best])
		}
	})(list(NULL), best_root)))

	core_phis <- c(core_phis, as.list(best_root))
	core_phis <- core_phis[!duplicated(term_set_names(core_phis))]

	phi_list <- core_phis
	result_mls <- ml(phi_list)
	result_priors <- f_priors(phi_list)
	result_totals <- result_mls + result_priors
	counted_names <- term_set_names(phi_list)
	best <- max(result_totals)

	result <- if (best-bgml+log(max_phi_count) < min_BF) {
		list(
			phis=phi_list,
			ml=result_totals-bgml,
			priors=result_priors,
			lik=result_mls
		)
	} else {

		restrict <- unique(unlist(use.names=FALSE, lapply(x[y], get_ancestors, ontology=ontology)))

		s <- core_phis
		while (length(s) > 0 & length(phi_list) < max_select) {
			d_prom <- lapply(do.call(what=c, lapply(s[result_totals[match(term_set_names(s), counted_names)] > log(min_ratio) + best], promotions, restrict=restrict, ontology=ontology)), minimal_set, ontology=ontology)
			d_dem <- lapply(do.call(what=c, lapply(s[result_totals[match(term_set_names(s), counted_names)] > log(min_ratio) + best], demotions, ontology=ontology)), minimal_set, ontology=ontology)
			d <- c(d_prom, d_dem)
			d_names <- term_set_names(d)
			s <- d[!duplicated(d_names) & !d_names %in% counted_names]

			phi_list <- c(phi_list, s)
			counted_names <- c(counted_names, term_set_names(s))
			result_priors <- c(result_priors, f_priors(s))
			result_mls <- c(result_mls, ml(s))
			result_totals <- result_mls + result_priors
			best <- max(result_totals)
		}


		if (!is.null(lik_method)) {
			best_total_lik <- -Inf
			selected <- order(result_totals, decreasing=TRUE)[seq(min(length(phi_list), max_phi_count))]
			pp_phis <- phi_list[selected]
			pp_priors <- result_priors[selected]
			liks <- pp_ml(pp_phis, pp_priors)
			list(
				phis=pp_phis,
				priors=pp_priors,
				ml=liks + pp_priors-bgml,
				lik=liks
			)
		} else {
			list(
				phis=phi_list,
				ml=result_totals-bgml,
				priors=result_priors,
				lik=result_mls
			)
		}
	}
	structure(class="sim_reg_output", c(list(term_names_used=ontology$name[unique(unlist(use.names=FALSE, result$phis))]), result))
}

#' Calculate log Bayes factor for similarity the model, \code{gamma=1} and baseline model, \code{gamma=0}.
#'
#' @param x \code{list} of term sets or \code{sim_reg_output} object.
#' @param ... If x is a \code{list} term sets, other arguments to pass to \code{\link{sim_reg}}, otherwise this is not used.
#' @return Numeric value.
#' @export
log_BF <- function(x, ...) {
	UseMethod("log_BF")
}

#' @rdname log_BF
#' @export
log_BF.default <- function(x, ...) {
	sum_log_probs(sim_reg(x=x, ...)$ml)
}

#' @rdname log_BF
#' @export
log_BF.sim_reg_output <- function(x, ...) {
	sum_log_probs(x$ml)
}

#' Calculate probability of association between \code{y} and \code{x}
#'
#' @param prior Numeric value determing prior probability that \code{gamma=1}.
#' @param ... Arguments to pass to \code{\link{log_BF}}.
#' @return Numeric value.
#' @export
prob_association <- function(..., prior=0.05) {
	bf <- log_BF(...)
	exp(bf) * prior / (1 - prior + exp(bf) * prior)
}

#' Calculate marginal probability of terms inclusion in \code{phi} from \code{sim_reg_out} object
#' 
#' @param sim_reg_out Object of class \code{sim_reg_output}.
#' @return Numeric vector of probabilities, named by term ID.
#' @export
get_term_marginals <- function(sim_reg_out) {
	stopifnot(class(sim_reg_out) == "sim_reg_output")
	all_terms <- unique(unlist(use.names=FALSE, sim_reg_out$phi))
	margs <- setNames(nm=all_terms, numeric(length(all_terms)))
	for (i in seq_along(sim_reg_out$phi)) {
		margs[match(sim_reg_out$phi[[i]], all_terms)] <- margs[match(sim_reg_out$phi[[i]], all_terms)] + exp(sim_reg_out$ml[i])
	}
	sort(margs/sum(exp(sim_reg_out$ml)), decreasing=TRUE)
}

#' Calculate marginal probability of terms inclusion in \code{phi} 
#'
#' @param ... Arguments to pass to \code{\link{sim_reg}}.
#' @return Numeric vector of probabilities, named by term ID.
#' @export
term_marginals <- function(...) {
	parts <- sim_reg(...)
	get_term_marginals(parts)
}
