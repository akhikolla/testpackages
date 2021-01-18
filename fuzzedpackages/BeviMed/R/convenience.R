subset_bevimed_m <- function(G, min_ac=1L, variant_weights=NULL, ...) {
	vars <- subset_variants(G=G, min_ac=min_ac, return_variants=TRUE)

	bevimed_m(
		G=G[,vars,drop=FALSE],
		min_ac=min_ac,
		variant_weights=if (is.null(variant_weights)) NULL else variant_weights[vars],
		...
	)
}

#' @title Bayesian Evaluation of Variant Involvement in Mendelian Disease 
#'
#' @description Infer probabilities of association between disease label and locus and posterior parameter values under BeviMed model.
#' 
#' @template y
#' @template G_matrix
#' @template ploidy
#' @template prior_prob_association
#' @template prior_prob_dominant
#' @param dominant_args Arguments to pass to \code{\link{bevimed_m}} conditioning on dominant inheritance.
#' @param recessive_args Arguments to pass to \code{\link{bevimed_m}} conditioning on recessive inheritance.
#' @param ... Arguments to be passed to \code{\link{bevimed_m}} for both modes of inheritance.
#' @return \code{BeviMed} object containing results of inference.
#' @export
#' @seealso \code{\link{prob_association}}, \code{\link{bevimed_m}}, \code{\link{summary.BeviMed}}, \code{\link{bevimed_polytomous}}
#' @template paper
bevimed <- function(
	y,
	G,
	ploidy=rep(2L, length(y)),
	prior_prob_association=0.01,
	prior_prob_dominant=0.5,
	dominant_args=NULL,	
	recessive_args=NULL,
	...
) {
	if (prior_prob_association > 1 | prior_prob_association < 0)
		stop("'prior_prob_association' must be between 0 and 1")
	if (prior_prob_dominant > 1 | prior_prob_dominant < 0)
		stop("'prior_prob_dominant' must be between 0 and 1")

	priors <- prior_prob_association * c(dominant=prior_prob_dominant, recessive=1-prior_prob_dominant)

	bevimed_polytomous(
		y=y,
		G=G,
		ploidy=ploidy,
		prior_prob_association=priors,
		variant_sets=list(dominant=seq(length.out=ncol(G)), recessive=seq(length.out=ncol(G))),
		moi=c("dominant","recessive"),
		model_specific_args=list(dominant=dominant_args, recessive=recessive_args),
		...
	)
}

#' @title Calculate marginal probability of observed case-control status y under model gamma = 0
#' @description Marginal probability calculated exactly by integration.
#' @template y
#' @param tau0_shape Beta shape hyper-priors for prior on rate of case labels
#' @return Log marginal likelihood.
#' @export
#' @seealso \code{\link{bevimed}}, \code{\link{gamma1_evidence}}
gamma0_evidence <- function(
	y,
	tau0_shape=c(1, 1)
) {
	lbeta(sum(y)+tau0_shape[1], length(y)-sum(y)+tau0_shape[2]) -
	lbeta(tau0_shape[1], tau0_shape[2])
}

#' @title Extract evidence for model gamma = 1
#' 
#' @description Extract evidence from \code{BeviMed_m} object.
#' @template x_BeviMed_m
#' @return Log marginal likelihood.
#' @export
#' @seealso \code{\link{gamma1_evidence}}, \code{\link{bevimed_m}}
extract_gamma1_evidence <- function(x) {
	stopifnot(class(x) == "BeviMed_m")
	sum_ML_over_PP(x[["traces"]][["y_log_lik_t_equals_1"]], x[["parameters"]][["temperatures"]])
}

#' @title Calculate evidence under model gamma = 1 
#'
#' @description Use \code{\link{bevimed_m}} to perform inference under model gamma = 1 and return only the log evidence/integrated likelihood.
#' 
#' @template dots_to_bevimed_m
#' @return Log marginal likelihood.
#' @export
#' @seealso \code{\link{bevimed_m}}, \code{\link{extract_gamma1_evidence}}
gamma1_evidence <- function(
	...
) {
	bv <- subset_bevimed_m(
		return_z_trace=FALSE,
		return_x_trace=FALSE,
		...
	)
	extract_gamma1_evidence(bv)
}

#' @title Extract expected number of explained cases
#'
#' @description Extract expected number of cases explained by pathogenic configurations of alleles from \code{BeviMed_m} object.
#' 
#' @template x_BeviMed_m
#' @return Numeric value.
#' @export
#' @seealso \code{\link{expected_explained}}, \code{\link{bevimed_m}}
extract_expected_explained <- function(x) {
	stopifnot(class(x) == "BeviMed_m")
	if (("x" %in% names(x[["traces"]]))*dim(x[["traces"]][["x"]])[1] == 0)
		stop("Must make sure to set 'return_x_trace=TRUE' in call to 'bevimed_m' to use this function")
	mean(apply(x$traces$x[,x$parameters$y,drop=FALSE], 1, sum))
}

#' @title Calculate expected number of explained cases
#'
#' @description Use \code{\link{bevimed_m}} to perform inference under model gamma = 1 and return only the expected number of cases explained by pathogenic allele configurations.
#'
#' @template dots_to_bevimed_m
#' @return Numeric value.
#' @export
#' @seealso \code{\link{bevimed_m}}, \code{\link{extract_expected_explained}}
expected_explained <- function(...) {
	extract_expected_explained(subset_bevimed_m(
		return_z_trace=FALSE,
		return_x_trace=TRUE,
		...
	))
}

#' @title Extract expected number of pathogenic variants in cases
#' @description Extract expected number of variants involved in cases explained by pathogenic conigurations of alleles from \code{BeviMed_m} object.
#' 
#' @template x_BeviMed_m
#' @return Numeric value.
#' @export
#' @seealso \code{\link{explaining_variants}}, \code{\link{bevimed_m}}
extract_explaining_variants <- function(x) {
	stopifnot(class(x) == "BeviMed_m")
	if (("z" %in% names(x[["traces"]]))*dim(x[["traces"]][["z"]])[1] == 0)
		stop("Must make sure to set 'return_z_trace=TRUE' and 'return_x_trace=TRUE' in call to 'bevimed_m' to use this function")
	G <- x$parameters$G
	logicalG <- matrix(G > 0, nrow=nrow(G), ncol=ncol(G))
	y <- x$parameters$y
	mean(mapply(
		SIMPLIFY=TRUE,
		FUN=function(z, x) sum(apply(logicalG[y & x, z, drop=FALSE], 2, any)),
		split(x$traces$z[,seq(to=ncol(x$traces$z), length.out=x$parameters$k),drop=FALSE], seq(length.out=nrow(x$traces$z))),
		split(x$traces$x, seq(length.out=nrow(x$traces$z)))
	))
}

#' @title Calculate expected number of pathogenic variants in cases
#'
#' @description Use \code{\link{bevimed_m}} to perform inference under model gamma = 1 and return only the expected number of pathogenic variants in cases.
#'
#' @template dots_to_bevimed_m
#' @return Numeric value.
#' @export
#' @seealso \code{\link{extract_explaining_variants}}, \code{\link{bevimed_m}}
explaining_variants <- function(...) {
	extract_explaining_variants(subset_bevimed_m(return_z_trace=TRUE, return_x_trace=TRUE, ...))
}

#' Calculate log Bayes factor between an association model with a given mode of inheritance and model gamma = 0
#' 
#' @description Compute log Bayes factor of an association model and model gamma = 0. 
#'
#' @template y
#' @template dots_to_bevimed_m
#' @template tau0_shape
#' @return Log Bayes factor.
#' @export
#' @seealso \code{\link{bevimed_m}}, \code{\link{prob_association_m}}
log_BF <- function(
	y,
	tau0_shape=c(1, 1),
	...
) {
	gamma1_evidence(
		y=y,
		...
	) - 
	gamma0_evidence(y, tau0_shape=tau0_shape)
}

#' @title Calculate probability of association for one mode of inheritance
#'
#' @description Equivalent to \code{\link{prob_association}} where the prior probability of one mode of inheritance is 1. This function is faster, as it only calls \code{\link{bevimed_m}} once.
#'
#' @template y
#' @template min_ac
#' @template prior_prob_association
#' @param ... Other arguments to pass to \code{\link{log_BF}}.
#' @return Probability value.
#' @export
#' @seealso \code{\link{log_BF}}, \code{\link{prob_association}}, \code{\link{bevimed_m}}
prob_association_m <- function(
	y,
	min_ac=1L,
	prior_prob_association=0.01,
	...
) {
	bf <- log_BF(y, min_ac=min_ac, ...)
	num <- prior_prob_association*exp(bf)
	num/(1-prior_prob_association+num)
}

#' @title Extract the posterior probability of association
#'
#' @description Get posterior probability of association as numeric value, or optionally as numeric vector of length two with probabilities broken down by mode of inheritance (by passing \code{by_model=TRUE}), from a \code{BeviMed} object.
#' @template x_BeviMed
#' @template by_model
#' @return Probability values.
#' @export
#' @seealso \code{\link{prob_association}}, \code{\link{bevimed}}
extract_prob_association <- function(x, by_model=FALSE) {
	stopifnot(class(x) == "BeviMed")
	priors_by_model <- x[["parameters"]][["prior_prob_association"]]
	evidence_by_model <- sapply(x[["models"]], extract_gamma1_evidence)
	g0_ev <- gamma0_evidence(y=x[["parameters"]][["y"]], tau0_shape=x[["parameters"]][["tau0_shape"]])

	modal_model_ev <- max(c(evidence_by_model, g0_ev))
	numerator <- priors_by_model * exp(evidence_by_model-modal_model_ev)
	numerator_g0 <- (1-sum(priors_by_model)) * exp(g0_ev-modal_model_ev)
	denom <- sum(c(numerator, numerator_g0))

	probs_by_model <- numerator/denom
	if (by_model) probs_by_model 
	else sum(probs_by_model)
}

#' @title Calculate probability of association
#'
#' @description Calculate probability of an association between case/control label and allele configuration, optionally broken down by mode of inheritance/model.
#' @template by_model
#' @param ... Arguments to pass to \code{\link{bevimed}}. 
#' @return Probability of association.
#' @export
#' @seealso \code{\link{bevimed}}, \code{\link{extract_prob_association}}
prob_association <- function(
	by_model=FALSE,
	...
) {
	bv <- bevimed(..., return_z_trace=FALSE, return_x_trace=FALSE)
	extract_prob_association(bv, by_model=by_model)
}

variant_marginals <- function(var_probs_by_model, variant_sets, k) {
	if (k == 0) numeric(0)
	else apply(matrix(nrow=k, ncol=length(variant_sets), data=mapply(
		SIMPLIFY=TRUE,
		FUN=function(probs, inds) {
			p <- rep(0, k)
			p[inds] <- probs
			p
		},
		var_probs_by_model,
		variant_sets
	)), 1, sum)
}

#' @title Extract variant marginal probabilities of pathogenicity
#'
#' @description Extract the marginal probability of pathogenicity for individual variants from \code{BeviMed} object, optionally broken down by mode of inheritance/model.
#'
#' @template x_BeviMed
#' @template by_model
#' @return A vector of probabilities of pathogenicity for individual variants, or if \code{by_model} is \code{TRUE}, then a matrix of probabilities, with rows corresponding to modes of inheritance and columns to variants.
#' @export
#' @seealso \code{\link{prob_pathogenic}}, \code{\link{bevimed}}
extract_prob_pathogenic <- function(x, by_model=TRUE) {
	probs <- extract_prob_association(x, by_model=TRUE)
	var_probs_by_model <- mapply(SIMPLIFY=FALSE, FUN="*", probs, lapply(x$models, extract_conditional_prob_pathogenic))
	if (by_model) {
		var_probs_by_model
	} else {
		k <- ncol(x[["parameters"]][["G"]])
		variant_marginals(var_probs_by_model, x[["parameters"]][["variant_sets"]], k)
	}
}

#' Calculate variant marginal probabilities of pathogencity 
#'
#' Calls \code{\link{bevimed}} and \code{\link{extract_prob_pathogenic}} to obtain marginal probabilities of pathogenicity.
#'
#' @template by_model
#' @param ... Arguments to pass to \code{\link{bevimed}}. 
#' @return If \code{by_model} is \code{FALSE}, a vector of probabilities of pathogenicity for each variant, otherwise a list of vectors of probabilities of pathogenicity conditional on each compared association model.
#' @export
#' @seealso \code{\link{extract_prob_pathogenic}}, \code{\link{bevimed}}
prob_pathogenic <- function(
	by_model=FALSE,
	...
) {
	bv <- bevimed(..., return_z_trace=TRUE, return_x_trace=FALSE)
	extract_prob_pathogenic(bv, by_model=by_model)
}

t1_z_trace <- function(x) {
	k <- x[["parameters"]][["k"]]
	x[["traces"]][["z"]][,k*(length(x[["parameters"]][["temperatures"]])-1L)+seq(length.out=k),drop=FALSE]
}

#' @title Extract probability of pathogenicity for variant conditional on a given association model
#'
#' @description Extract the probability of pathogenicity for individual variants from a \code{BeviMed_m} object.
#'
#' @template x_BeviMed_m
#' @return Vector of probabilities of pathogenicity for individual variants.
#' @export
#' @seealso \code{\link{conditional_prob_pathogenic}}, \code{\link{bevimed_m}}
extract_conditional_prob_pathogenic <- function(x) {
	if (x[["parameters"]][["k"]] == 0) numeric(0) else apply(t1_z_trace(x), 2, mean)
}

#' Calculate probability of pathogencity for variants conditional on mode of inheritance.
#'
#' Calls \code{\link{bevimed_m}} and \code{\link{extract_conditional_prob_pathogenic}} to obtain probabilities of pathogenicity.
#' @template dots_to_bevimed_m
#' @return Probabilities of pathogenicity.
#' @export
#' @seealso \code{\link{extract_conditional_prob_pathogenic}}, \code{\link{bevimed_m}}
conditional_prob_pathogenic <- function(
	...
) {
	bv <- bevimed_m(
		return_z_trace=TRUE,
		return_x_trace=FALSE,
		...
	)
	extract_conditional_prob_pathogenic(bv)
}

#' @title Model selection for multiple association models
#' @description Apply bevimed to the no association model (gamma = 0) and multiple association models for different sets of variants, for instance, corresponding to different functional consequences.
#' @template y
#' @template G_matrix
#' @template ploidy
#' @param variant_sets List of integer vectors corresponding to sets of indices of \code{G}, each of which is to be considered in a model explaining the phenotype, \code{y}.
#' @template prior_prob_association
#' @template tau0_shape
#' @param moi Character vector giving mode of inheritance for each model.
#' @param model_specific_args List of named lists of parameters to use in \code{\link{bevimed_m}} applications for specific models.
#' @param ... Other arguments to pass to \code{\link{bevimed_m}}.
#' @seealso \code{\link{bevimed_m}}, \code{\link{bevimed}}
#' @template paper
#' @export
bevimed_polytomous <- function(
	y,
	G,
	ploidy=rep(2L, length(y)),
	variant_sets, 
	prior_prob_association=rep(0.01/length(variant_sets), length(variant_sets)),
	tau0_shape=c(1, 1),
	moi=rep("dominant", length(variant_sets)),
	model_specific_args=vector(mode="list", length=length(variant_sets)),
	...
) {
	if (sum(prior_prob_association) > 1 | sum(prior_prob_association) < 0)
		stop("The sum of 'prior_prob_association' must be between 0 and 1")
	if (any(prior_prob_association > 1) | any(prior_prob_association < 0))
		stop("Each element of 'prior_prob_association' must be between 0 and 1")
	n_models <- length(variant_sets)
	if (length(model_specific_args) != n_models)
		stop("The length of 'variant_sets' must be the same as 'model_specific_args'")
	if (length(prior_prob_association) != n_models)
		stop("The length of 'variant_sets' must be the same as 'prior_prob_association'")
	if (length(moi) != n_models)
		stop("The length of 'variant_sets' must be the same as 'moi'")
	if (!all(unlist(use.names=FALSE, variant_sets) %in% seq(length.out=ncol(G))))
		stop("All indices given in 'variant_sets' must be between one and the number of variants in 'G'")

	shared <- list(...)

	structure(
		class="BeviMed",
		list(
			parameters=list(
				tau0_shape=tau0_shape,
				prior_prob_association=setNames(nm=names(variant_sets), prior_prob_association),
				ploidy=ploidy,
				y=y,
				G=G,
				variant_table=to_var_tab(get_G_args(G)),
				variant_sets=variant_sets,
				moi=moi
			),
			models=Map(
				f=function(var_inds, min_ac, args) { do.call(what=bevimed_m, c(list(y=y, G=G[,var_inds,drop=FALSE], min_ac=min_ac), args, shared)) },
				variant_sets,
				lapply(moi, function(m) if (m=="dominant") { 1L } else { if (length(ploidy) == 0) 2L else ploidy }),
				model_specific_args
			)
		)
	)
}
