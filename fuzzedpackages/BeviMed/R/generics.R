#' Summarise a \code{BeviMed_m} object
#'
#' Create a summary of inference conditional on mode of inheritance. 
#'
#' @param object Object of class \code{BeviMed_m}. See function \code{\link{bevimed_m}}.
#' @template confidence
#' @template simulations
#' @param ... Unused arguments.
#' @details Returns a \code{BeviMed_m_summary} object, which is a list containing elements:
#' \itemize{
#' \item `gamma1_evidence': the log evidence under model gamma = 1,
#' \item `gamma1_evidence_confidence_interval': a confidence interval for the log evidence under model gamma = 1,
#' \item `conditional_prob_pathogenic': vector of marginal probabilities of pathogenicity for individual variants,
#' \item `expected_explained': the expected number of cases with a pathogenic configuration of alleles,
#' \item `explaining_variants': the expected number of variants present for which cases harbour a rare allele,
#' \item `number_of_posterior_samples': the number of samples from the posterior distribution of the model parameters which upon which the summary is based,
#' \item `omega_estimated': logical value indicating whether the parameter omega was estimated,
#' \item `omega': the posterior mean of omega,
#' \item `omega_acceptance_rate': if omega was estimated, the rate of acceptance of proposed omega values in the Metropolis-Hastings sampling routine,  
#' \item `phi_estimated': logical value indicating whether the parameter phi was estimated,
#' \item `phi': the posterior mean of phi,
#' \item `phi_acceptance_rate': if phi was estimated, the rate of acceptance of proposed phi values in the Metropolis-Hastings sampling routine, 
#' \item `N`: number of samples in the analysis,
#' \item `k`: number of variants in the analysis,
#' \item `variant_counts': list of counts of each variant for cases and controls,
#' \item `temperatures': numeric vector of temperatures used as temperatures for tempered MCMC chains 
#' }
#' @return Object of class \code{BeviMed_m_summary}.
#' @method summary BeviMed_m
#' @export
#' @seealso \code{\link{summary.BeviMed}}
summary.BeviMed_m <- function(object, confidence=0.95, simulations=1000, ...) {
	vt <- object[["parameters"]][["variant_table"]]
	y <- object[["parameters"]][["y"]]
	k <- object[["parameters"]][["k"]]
	variant_counts <- lapply(setNames(nm=c(F,T)), function(y_is) as.integer(table(factor(vt$variant[y[vt$case]==y_is], levels=seq(length.out=k)))))

	gamma1_evidence <- extract_gamma1_evidence(object)
	temps <- object[["parameters"]][["temperatures"]]
	num_temps <- length(temps)
	
	phi_estimated <- object[["parameters"]][["estimate_phi"]]
	omega_estimated <- object[["parameters"]][["estimate_omega"]]

	has_z <- dim(object[["traces"]][["z"]])[1] > 0
	has_x <- dim(object[["traces"]][["x"]])[1] > 0

	structure(list(
		gamma1_evidence=gamma1_evidence,
		gamma1_evidence_confidence_interval=CI_gamma1_evidence(
			temperatures=temps,
			y_log_lik_t_equals_1_traces=object[["traces"]][["y_log_lik_t_equals_1"]],
			confidence=confidence,
			simulations=simulations
		),
		conditional_prob_pathogenic=if (has_z) setNames(nm=colnames(object[["parameters"]][["G"]]), extract_conditional_prob_pathogenic(object)),
		expected_explained=if (has_x) extract_expected_explained(object) else NULL,
		explaining_variants=if (has_z & has_x) extract_explaining_variants(object) else NULL,
		number_of_posterior_samples=nrow(object[["traces"]][["y_log_lik_t_equals_1"]]),
		omega_estimated=omega_estimated,
		phi_estimated=phi_estimated,
		phi=if (phi_estimated) mean(exp(object[["traces"]][["log_phi"]][,num_temps])) else NA,
		omega=if (omega_estimated) { 
			mean(1-1/(1+exp(object[["traces"]][["logit_omega"]][,num_temps]))) 
		} else { 
			if (has_z & all(object[["parameters"]][["c_weights"]] == 0)) {
				sumZ <- apply(t1_z_trace(object),1,sum)
				(object[["parameters"]][["omega_shape"]][1] + mean(sumZ))/sum(c(k,object[["parameters"]][["omega_shape"]]))
			} else {
				NA 
			}
		},
		phi_acceptance_rate=if (phi_estimated) apply(object[["traces"]][["log_phi"]], 2, function(log_phis) mean(log_phis[-length(log_phis)] != log_phis[-1])) else NA, 
		omega_acceptance_rate=if (omega_estimated) apply(object[["traces"]][["logit_omega"]], 2, function(logit_omegas) mean(logit_omegas[-length(logit_omegas)] != logit_omegas[-1])) else NA, 
		N=length(object[["parameters"]][["y"]]),
		k=object[["parameters"]][["k"]],
		variant_counts=variant_counts,
		temperatures=temps
	), class="BeviMed_m_summary")
}

#' Summarise a \code{BeviMed} object
#'
#' Create a summary of inference over model gamma = 0 and association models. 
#'
#' @param object Object of class \code{BeviMed}.
#' @param ... Arguments passed to \code{summary.BeviMed_m}.
#' @details Returns a \code{BeviMed_summary} object, which is a list containing elements:
#' \itemize{
#' \item `prob_association`: the probability of association under each association model,
#' \item `prior_prob_association`: the prior probability of association for each association model,
#' \item `gamma0_evidence': the log evidence under model gamma = 0,
#' \item `models': a list of summaries of model conditional inferences, i.e. objects of class \code{BeviMed_m_summary}. See \code{\link{summary.BeviMed_m}} for more details.
#' }
#' @return Object of class \code{BeviMed_summary}.
#' @method summary BeviMed
#' @seealso \code{\link{summary.BeviMed_m}}
#' @export
summary.BeviMed <- function(object, ...) {
	structure(
		class="BeviMed_summary",
		c(
			object[["parameters"]][c(
				"prior_prob_association",
				"variant_sets",
				"moi"
			)],
		  	list(
				gamma0_evidence=gamma0_evidence(y=object[["parameters"]][["y"]], tau0_shape=object[["parameters"]][["tau0_shape"]]),
				prob_association=extract_prob_association(object, by_model=TRUE),
				models=lapply(object[["models"]], summary.BeviMed_m, ...),
				N=length(object[["parameters"]][["y"]]),
				k=ncol(object[["parameters"]][["G"]]),
				variant_names=colnames(object[["parameters"]][["G"]])
			)
		)
	)
}

#' Print readable summary of \code{BeviMed_summary} object.
#'
#' @template print_description
#' @param x \code{BeviMed_summary} object.
#' @param print_prob_pathogenic Logical value indicating whether to print list of marginal probabilities of \code{z_j = 1} for all variants \code{j} under each mode of inheritance.
#' @param ... Unused arguments 
#' @return Prints a summary
#' @method print BeviMed_summary
#' @export
print.BeviMed_summary <- function(x, print_prob_pathogenic=TRUE, ...) {
	stopifnot(class(x) == "BeviMed_summary")
	dashed <- paste0(rep("-", getOption("width")), collapse="")
	cat(dashed, "\n")
	cat("Posterior probability of association: \n\t", round(sum(x[["prob_association"]]), digits=3), " [prior: ", round(sum(x[["prior_prob_association"]]), digits=3), "]\n", sep="")

	cat(dashed, "\n")

	model_names <- if (!is.null(names(x[["prob_association"]]))) names(x[["prob_association"]]) else seq(length.out=length(x[["prob_association"]]))

	summary_mat <- data.frame(
		check.names=FALSE,
		stringsAsFactors=FALSE,
		Model=model_names,
		`MOI`=substr(x[["moi"]], 1, 3),
		`Prior`=x[["prior_prob_association"]]/sum(x[["prior_prob_association"]]),
		`Post`=x[["prob_association"]]/sum(x[["prob_association"]]),
		`Cases`=sapply(x[["models"]], function(m) { ee <- "expected_explained"; if (ee %in% names(m)) { if (is.null(m[[ee]])) NA else m[[ee]] } else { NA } }),
		`Variants`=sapply(x[["models"]], function(m) { ee <- "explaining_variants"; if (ee %in% names(m)) { if (is.null(m[[ee]])) NA else m[[ee]] } else { NA } })
	)[order(x[["prob_association"]], decreasing=TRUE),]

	print(row.names=FALSE, digits=3, summary_mat)

	cat("\n")
	cat("MOI: mode of inheritance, dominant (dom) or recessive (rec)\n")
	cat("Prior: prior probability of model given association\n")
	cat("Post: posterior probability of model given association\n")
	cat("Cases: posterior expected number of cases explained\n")
	cat("Variants: posterior expected number of variants involved in explained cases\n")
	cat(dashed, "\n")
	if (print_prob_pathogenic) {
		if (x[["k"]] > 0) {
			cat("Probabilities of pathogenicity for individual variants given association\n\n")

			patho <- variant_marginals(Map(f="*", lapply(x[["models"]], "[[", "conditional_prob_pathogenic"), x[["prob_association"]]), x[["variant_sets"]], x[["k"]])/sum(x[["prob_association"]])
			patho_names <- if (!is.null(x[["variant_names"]])) x[["variant_names"]] else seq_along(patho)

			patho_rounded <- lapply(patho, function(vals) sprintf("%.2f", vals))

			bar_width <- 17
			print(row.names=FALSE, data.frame(
				check.names=FALSE,
				stringsAsFactors=FALSE,
				Var=substr(patho_names, 1, 22),
				`Probability pathogenic`=sapply(seq(length.out=length(patho)), function(j) paste0("[", patho_rounded[j], " ", paste0(collapse="", rep("=", as.integer(round(patho[j]*bar_width, digits=0L)))), paste0(collapse="", rep(" ", bar_width-as.integer(round(patho[j]*bar_width, digits=0L)))), "]"))
			))

		} else {
			cat("Specified models contain no variants\n")
		}
		cat(dashed, "\n")
	}
}

#' @title Print readable summary of \code{BeviMed} object
#'
#' @template print_description
#' @param x \code{BeviMed} object.
#' @param ... Arguments passed to \code{\link{summary.BeviMed}} 
#' @return Prints a summary.
#' @method print BeviMed
#' @export
#' @seealso \code{\link{summary.BeviMed}}
print.BeviMed <- function(x, ...) {
	stopifnot(class(x) == "BeviMed")
	print(summary(x, ...))
}

#' @title Print \code{BeviMed_m} object
#'
#' @description Print summary statistics for \code{BeviMed_m} object.
#' @template x_BeviMed_m
#' @param ... Unused arguments.
#' @return Prints a summary.
#' @seealso \code{\link{summary.BeviMed_m}}
#' @export
#' @method print BeviMed_m
print.BeviMed_m <- function(x, ...) {
	stopifnot(class(x) == "BeviMed_m")
	print(summary(x, ...))
}
