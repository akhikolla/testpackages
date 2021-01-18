#' @title R interface to BeviMed c++ MCMC procedure 
#' @description Allows other functions in the package to call the c++ function passing arguments more succinctly and by name.
#' @template samples_per_chain
#' @param y Logical vector of subject affectedness status.
#' @param block_starts Integer vector of k 0-indexed start positions (with respect to \code{cases} and \code{counts}) for contiguous blocks relating to the k variants.
#' @param block_ends Integer vector of (exclusive) k 0-indexed end positions.
#' @param cases 0 based vector of case indices with respect to y.
#' @param counts Vector of variant counts.
#' @template min_ac
#' @param tau_shape Beta distribution parameterisation of benign variant configuration rate of affection, q.
#' @param pi_shape Beta distribution parameterisation of pathogenic variant configuration rate of affection, p.
#' @param omega_shape Beta distribution of global rate of pathogenicty of variants in gene given pathogenicity of gene, omega.
#' @template temperatures
#' @param z0_matrix Matrix of logicals, where the rows are used as an initial zs for the chains.
#' @param estimate_omega Logical value determining whether to estimate the parameter omega.
#' @param logit_omegas Numeric vector of logit omega values, one value per chain.
#' @param logit_omega_proposal_sds Numeric vector of proposal standard deviations for Metropolis-Hastings sampling of logit omega parameter, one value per chain.
#' @template variant_weights
#' @template estimate_phi
#' @param log_phis Numeric vector of log phi values, one value per chain.
#' @template log_phi_mean
#' @template log_phi_sd
#' @param log_phi_proposal_sds Numeric vector of proposal standard deviations for Metropolis-Hastings sampling of log phi parameter, one value per chain.
#' @param chain_swaps_per_cycle Number of chain swaps to propose per update cycle.
#' @param annealing Logical value determining whether to anneal the chains, e.g. for optimisation.
#' @template tandem_variant_updates
#' @param comphet_variant_block_starts 0-indexed start positions for contiguous blocks of variants in \code{comphet_variants}.
#' @param comphet_variant_block_ends As \code{comphet_variant_block_starts} for (exclusive) stop positions.
#' @param comphet_variants Integer vector giving variant numbers (0-based, i.e. between 0 and k-1). Used to pick pairs of variants for tandem updates from.
#' @template return_z_trace
#' @template return_x_trace
#' @template burn
#' @param check Logical value indicating whether to perform validation on the arguments before calling the c++ function. 
#' @return Object of class \code{BeviMed_raw}, containing the output of the MCMC sampling.
#' @importFrom Rcpp evalCpp
#' @useDynLib BeviMed
call_cpp <- function(
	samples_per_chain,
	y,
	block_starts,
	block_ends,
	cases,
	counts,
	min_ac,
	tau_shape,
	pi_shape,
	omega_shape,
	temperatures,
	z0_matrix,
	estimate_omega,
	logit_omegas,
	logit_omega_proposal_sds,
	variant_weights,
	estimate_phi,
	log_phis,
	log_phi_mean,
	log_phi_sd,
	log_phi_proposal_sds,
	chain_swaps_per_cycle,
	annealing,
	tandem_variant_updates,
	comphet_variant_block_starts,
	comphet_variant_block_ends,
	comphet_variants,
	return_z_trace,
	return_x_trace,
	burn=0,
	check=TRUE
) {
	if (check) {
		stopifnot(length(omega_shape) == 2 & min(omega_shape) > 0)
		stopifnot(length(tau_shape) == 2 & min(tau_shape) > 0)
		stopifnot(length(pi_shape) == 2 & min(pi_shape) > 0)
		stopifnot(is.matrix(z0_matrix))
		stopifnot(ncol(z0_matrix) == length(variant_weights))
		stopifnot(length(temperatures) == nrow(z0_matrix))
		stopifnot(length(temperatures) == nrow(logit_omegas))
		stopifnot(length(temperatures) == nrow(log_phis))
		stopifnot(length(temperatures) == length(logit_omega_proposal_sds))
		stopifnot(length(temperatures) == length(log_phi_proposal_sds))
		stopifnot(estimate_phi | all(log_phis == 0))
		stopifnot(length(block_ends) == length(block_starts))
		stopifnot(block_ends[length(block_ends)] == length(cases))
		stopifnot(length(cases) == length(counts))
		if (tandem_variant_updates > 0 & length(unique(comphet_variants)) < 2) stop("Must have more than 1 variant to select from if making tandem updates")
	}

	raw <- bevimed_mc(
		samples_per_chain,
		y,
		block_starts,
		block_ends,
		cases,
		counts,
		if (length(min_ac) == length(y)) min_ac else rep(min_ac, length(y)),
		tau_shape[1],
		tau_shape[2],
		pi_shape[1],
		pi_shape[2],
		omega_shape[1],
		omega_shape[2],
		z0_matrix,
		estimate_omega,
		logit_omegas,
		logit_omega_proposal_sds,
		variant_weights,
		estimate_phi,
		log_phis,
		log_phi_mean,
		log_phi_sd,
		log_phi_proposal_sds,
		temperatures,
		chain_swaps_per_cycle,
		annealing,
		tandem_variant_updates,
		comphet_variant_block_starts,
		comphet_variant_block_ends,
		comphet_variants,
		return_z_trace,
		return_x_trace
	)

	if (burn > 0) {
		raw$traces <- lapply(
			raw$traces,
			function(x) x[-seq(length.out=burn),,drop=FALSE]
		)
		raw$swaps <- lapply(
			raw$swaps,
			function(x) x[-seq(length.out=chain_swaps_per_cycle * burn)]
		)
	}

	structure(
		class="BeviMed_raw",
		raw
	)
}

#' @title Calculate marginal likelihood from power posteriors output
#' @description Calculate the Marginal Likelihood by summation over power posterior likelihood exptectances
#' @template y_log_lik_t_equals_1_traces
#' @param temperatures Numeric vector of temperatures used to produce \code{y_log_lik_t_equals_1_traces}.
#' @return Numeric value of estimated log marginal likelihood.
sum_ML_over_PP <- function(y_log_lik_t_equals_1_traces, temperatures) {
	sum(mapply(FUN=function(y_lik, t_diff) { log(mean(exp(t_diff*y_lik-max(t_diff*y_lik))))+max(t_diff*y_lik) }, split(t(y_log_lik_t_equals_1_traces)[-length(temperatures),], seq(length.out=length(temperatures)-1)), diff(temperatures)))
}

#' @title Concatenate objects of class \code{BeviMed_raw}
#' @description This function could be used to stitch together consecutive chains to create one larger sampled set of states from the MCMC procedure.
#' @param objects \code{list} of \code{BeviMed_raw} objects.
#' @return \code{BeviMed} object.
#' @importFrom stats setNames
stack_BeviMeds <- function(objects) {
	stopifnot(all(sapply(objects, class) == "BeviMed_raw"))

	structure(
		list(
			traces=do.call(what=Map, c(list(f=rbind), lapply(objects, "[[", "traces"))),
			swaps=do.call(what=c, c(list(f=rbind), lapply(objects, "[[", "swaps"))),
			final=objects[[length(objects)]]$final
		),
		class="BeviMed_raw"
	)
}

#' @title Estimate confidence interval for estimated marginal likelihood
#' @description Central limit theorem not applicable so use simulation to estimate confidence interval for evidence.
#' @template temperatures
#' @template y_log_lik_t_equals_1_traces
#' @template confidence
#' @template simulations
#' @return Confidence interval as numeric vector of length 2.
#' @importFrom stats rt quantile
CI_gamma1_evidence <- function(
	temperatures,
	y_log_lik_t_equals_1_traces,
	confidence=0.95,
	simulations=1000
) {
	#we want all blocks to have the same length
	a <- b <- as.integer(sqrt(nrow(y_log_lik_t_equals_1_traces)))
	if (a < 2) stop("Longer sample-block lengths are required to make confidence interval estimation")

	#with columns for batches and rows for temperatures...
	batch_means <- simplify2array(tapply(X=(nrow(y_log_lik_t_equals_1_traces)-a*b+1):nrow(y_log_lik_t_equals_1_traces), INDEX=gl(n=a, k=b), FUN=function(rows) apply(MARGIN=1, FUN=mean, X=exp(t(y_log_lik_t_equals_1_traces[rows,-ncol(y_log_lik_t_equals_1_traces),drop=FALSE]) * diff(temperatures)))))

	overall_mean <- apply(MARGIN=1, FUN=mean, X=exp(t(y_log_lik_t_equals_1_traces[,-ncol(y_log_lik_t_equals_1_traces)]) * diff(temperatures)))

	estimate_var_by_temp <- mapply(FUN=function(overall, batch) b / (a-1) * sum((batch-overall)^2), overall_mean, split(batch_means, seq(length.out=nrow(batch_means))))

	samples <- mapply(SIMPLIFY=TRUE, FUN=function(est_mean, est_var) est_mean + (rt(df=a-1, n=simulations) * sqrt(est_var) / sqrt(a)), overall_mean, estimate_var_by_temp)

	simulated_MLs <- apply(samples, 1, function(expectance_at_temp_i) sum(log(ifelse(expectance_at_temp_i < 0, 0, ifelse(expectance_at_temp_i > 1, 1, expectance_at_temp_i)))))

	quantile(probs=c((1-confidence)/2,1-(1-confidence)/2), simulated_MLs)
}

#' @title Remove variants with no data for pathogenicity
#' @description Subset an allele count matrix given a minimum allele count threshold for pathogenicity per individual so that only variants for which data relevant to pathogencity are retained. This is useful to apply before running \code{\link{bevimed}} as it reduces the size of the parameter space used in the inference. 
#' @template G_matrix
#' @template min_ac 
#' @param return_variants Logical value determining whether to return an integer vector of indices of retained variants or the subsetted allele count matrix
#' @export
subset_variants <- function(G, min_ac=1L, return_variants=FALSE) {
	vars <- which(apply(G[apply(G, 1, sum) >= min_ac,,drop=FALSE], 2, sum) > 0L)
	if (return_variants)
		vars
	else
		G[,vars,drop=FALSE]
}

#' @title Apply the MCMC algorithm in blocks until conditions are met
#' @description Sample blocks of a given size until either the estimated log marginal likelihood falls within a given confidence interval, there is sufficient confidence that the evidence model gamma = 1 is at most a certain quantity, or a certain number of blocks have been sampled.
#' @template y
#' @param blocks_remaining Maximum number of blocks left before termination.
#' @param start_zs Initial (logical) z-matrix.
#' @param start_logit_omegas Initial values of logit_omega (numeric vector - one value per chain).
#' @param start_log_phis Initial values of log_phi (numeric vector - one value per chain).
#' @template temperatures
#' @param tolerance Maximum width for confidence_interval of log marginal likelihood to allow before stopping the chain.
#' @template confidence
#' @template simulations
#' @param log_evidence_threshold Numeric value used to determine whether to stop the sampling procedure after successive blocks. If we are confident (to the level of \code{confidence}) that the evidence for model gamma = 1 is under this value, sampling is halted.
#' @template y_log_lik_t_equals_1_traces
#' @param full_block_traces List of outputs of calls to MCMC routine. 
#' @template verbose
#' @param ... Other arguments passed to \code{\link{call_cpp}}
#' @return An object of class \code{BeviMed}.
stop_chain <- function(
	y,
	blocks_remaining,
	start_zs,
	start_logit_omegas,
	start_log_phis,
	temperatures,
	tolerance=1,
	confidence=0.95,
	simulations=1000,
	log_evidence_threshold=-Inf,
	y_log_lik_t_equals_1_traces=matrix(ncol=length(temperatures),nrow=0),
	full_block_traces=list(),
	verbose=FALSE,
	...
) {
	if (verbose) cat("Sampling up to ", blocks_remaining, " more blocks to get marginal likelihood within tolerance of ", tolerance, "\n", sep="")

	mc <- call_cpp(
		z0_matrix=start_zs,
		logit_omegas=start_logit_omegas,
		log_phis=start_log_phis,
		temperatures=temperatures,
		y=y,
		...
	)

	y_ll <- rbind(
		y_log_lik_t_equals_1_traces,
		mc[["traces"]][["y_log_lik_t_equals_1"]]
	)
	
	confidence_interval <- CI_gamma1_evidence(temperatures, y_ll, confidence=confidence, simulations=simulations)

	if (verbose) {
		cat(paste0(round(confidence * 100), "% confidence interval, width = ", round(digits=2, diff(confidence_interval)), ":\n"))
		print(round(digits=2, confidence_interval))
	}

	if (diff(confidence_interval) < tolerance | confidence_interval[2] < log_evidence_threshold | blocks_remaining <= 1) {
		stack_BeviMeds(c(full_block_traces, list(mc)))
	} else {
		stop_chain(
			y=y,
			blocks_remaining=blocks_remaining - 1,
			tolerance=tolerance,
			confidence=confidence,
			start_zs=mc[["final"]][["z"]],
			start_logit_omegas=mc[["final"]][["logit_omega"]],
			start_log_phis=mc[["final"]][["log_phi"]],
			temperatures=temperatures,
			log_evidence_threshold=log_evidence_threshold,
			y_log_lik_t_equals_1_traces=y_ll,
			full_block_traces=c(full_block_traces, list(mc)),
			verbose=verbose,
			...
		)
	}
}

#' @title Tune proposal standard deviation for MH sampled parameters
#' @description Tune the proposal standard deviations for the Metropolis-Hastings updates of either phi or omega
#' @param tune_for Character vector of length one, naming which variable to tune the proposal SDs for: either \code{"logit_omega"} or \code{"log_phi"}.
#' @param initial_proposal_sds Numeric vector with the initial values of the proposal SDs.
#' @param target_acceptance_range Numeric vector of length 2 where the first element is the lower bound for the acceptance interval and the second is the upper bound.
#' @param other_param_proposal_sd The proposal SD to use for \code{log_phi} when tuning \code{logit_omega} or vice versa.
#' @param max_tuning_cycles Maximum number of tuning cycles to perform before returning the proposal SDs as they are.
#' @param initial_rate Initial rate at which to mutate the proposal SDs.
#' @param rate_decay Geometric rate of decay for size of proposal SD mutation with each successive tuning cycle.
#' @template verbose
#' @param ... Other arguments to be passed to \code{\link{call_cpp}}.
#' @return Numeric vector of proposal SDs for the different temperature chains.
tune_proposal_sds <- function(tune_for=c("logit_omega"), initial_proposal_sds, target_acceptance_range=c(0.3,0.7), other_param_proposal_sd=0.7, max_tuning_cycles=10, initial_rate=1, rate_decay=1.2, verbose=FALSE, ...) {
	stopifnot(all(tune_for %in% c("logit_omega", "log_phi")))

	if (verbose) {
		cat("Tuning proposal standard deviations for ", tune_for[1], " targeting acceptance rate range (", target_acceptance_range[1], ",", target_acceptance_range[2], ")\n", sep="")
	}

	current_proposal_sds <- initial_proposal_sds
	acceptances <- rep(-Inf, length(current_proposal_sds))
	cycle <- 0
	while (
		cycle <= max_tuning_cycles 
		& ( any(acceptances < target_acceptance_range[1])
		   |any(acceptances > target_acceptance_range[2]))
	) {
		if (verbose) cat("Tuning cycle ", cycle, "\n\tTest proposal SDs:\n\t\t", paste0(collapse=" : ", round(digits=2, current_proposal_sds)), "\n", sep="")
		out <- if (tune_for[1] == "logit_omega") call_cpp(logit_omega_proposal_sds=current_proposal_sds, log_phi_proposal_sds=rep(other_param_proposal_sd, length(current_proposal_sds)), ...) else call_cpp(log_phi_proposal_sds=current_proposal_sds, logit_omega_proposal_sds=rep(other_param_proposal_sd, length(current_proposal_sds)), ...)

		acceptances <- apply(out[["traces"]][[tune_for[1]]], 2, function(var_vals) mean(var_vals[-length(var_vals)] != var_vals[-1]))
		current_proposal_sds <- mapply(SIMPLIFY=TRUE, FUN=function(prop_sd, acc_rate) if (acc_rate < target_acceptance_range[1] | acc_rate > target_acceptance_range[2]) { match.fun(if (acc_rate < target_acceptance_range[1]) "/" else "*")(prop_sd, (1 + initial_rate * rate_decay ^ (-cycle+1))) } else { prop_sd }, current_proposal_sds, acceptances)

		if (verbose) cat("\tAcceptance rates:\n\t\t", paste0(collapse=" : ", round(digits=2, acceptances)), "\n", sep="")
		cycle <- cycle + 1
	}
	if (verbose) cat("Terminating\n")
	current_proposal_sds
}

#' @title Tune temperatures
#' @description Tune temperatures using interval bisection to minimimise Kullback-Liebler divergence between adjacent power posteriors
#' @param number_of_temperatures Integer value giving number of tuned temperatures (including 0 and 1) to obtain.
#' @param return_temperatures Logical value determining whether to return just the numeric vector of tuned temperatures or to return the \code{BeviMed_m}-classed object containing the output of the MCMC sampling.
#' @param ... Other arguments to pass to \code{call_cpp}.
#' @return If \code{return_temperatures == TRUE}, a numeric vector of tuned temperatures, otherwise an object of class \code{BeviMed_m}.
#' @importFrom stats var
tune_temperatures <- function(
	number_of_temperatures,
	return_temperatures=FALSE,
	...
) {
	temperatures <- 0:1
	chains <- lapply(temperatures, function(t) call_cpp(
		temperatures=t,
		...
	))

	E <- sapply(chains, function(ch) mean(ch[["traces"]][["y_log_lik_t_equals_1"]]))
	V <- sapply(chains, function(ch) var(ch[["traces"]][["y_log_lik_t_equals_1"]]))

	while (length(temperatures) <= number_of_temperatures) {
		areas <- diff(temperatures) * diff(E)
		largest <- which.max(areas)
		t_pair <- c(largest, largest+1)
		t_intersect <- mean(temperatures[t_pair])
		temperatures <- c(temperatures[seq(length.out=largest)], t_intersect, temperatures[(largest+1):length(temperatures)])
		chains <- c(chains[seq(length.out=largest)], list(call_cpp(temperatures=t_intersect, ...)), chains[(largest+1):length(chains)])

		E <- sapply(chains, function(ch) mean(ch[["traces"]][["y_log_lik_t_equals_1"]]))
		V <- sapply(chains, function(ch) var(ch[["traces"]][["y_log_lik_t_equals_1"]]))
	}
	
	if (return_temperatures) temperatures else list(chains=chains, E=E, V=V, temperatures=temperatures)
}

#' @importFrom methods is slot
get_G_args <- function(G) {
	if (is.matrix(G)) {	
		counts <- as.integer(G)
		variants <- rep(seq(length.out=ncol(G)), each=nrow(G))
		cases <- rep(seq(length.out=nrow(G)), times=ncol(G))

		block_ends <- unname(cumsum(lapply(split(counts, factor(variants, levels=seq(length.out=ncol(G)))), function(cnts) sum(cnts > 0))))
		block_starts <- if (length(block_ends) > 0) c(0, block_ends[-length(block_ends)]) else integer(0)
		list(
			cases=cases[counts > 0],
			counts=counts[counts > 0],
			block_ends=as.integer(block_ends),
			block_starts=as.integer(block_starts)
		)
	} else if (is(G, "sparseMatrix")) {
		p <- slot(G, "p")
		G_i <- slot(G, "i")
		return(list(
			cases=rep(1L, length(G_i))+G_i,
			counts=as.integer(slot(G, "x")),
			block_ends=p[-1L],
			block_starts=p[-length(p)]
		))
	} else {
		stop("'G' must be a matrix!")
	}
}

to_var_tab <- function(G_args) {
	data.frame(
		variant=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, seq(length.out=length(G_args$block_starts)), G_args$block_ends-G_args$block_starts)),
		case=G_args$cases,
		count=G_args$counts
	)
}

#' @title Perform inference under model gamma = 1 conditional on mode of inheritance
#'
#' @description Sample from posterior distribution of parameters under model gamma = 1 and conditional on mode of inheritance, set via the \code{min_ac} argument.
#' @template y
#' @template G_matrix
#' @template min_ac
#' @template tau_shape
#' @template pi_shape
#' @template omega_shape 
#' @template samples_per_chain
#' @param stop_early Logical value determining whether to attempt to stop the sampling as soon as certain conditions are met (i.e. either the estimated marginal log likelihood lies within a certain confidence interval, or we are sufficiently confidence that the log Bayes factor against of model gamma = 1 over model gamma = 0 is sufficiently low).
#' @param blocks Maximum number of blocks of \code{samples_per_chain} samples to draw before either the confidence interval for the marginal likelihood under the model gamma = 1 is sufficiently small or terminating the sampling. This parameter is ignored if unless \code{stop_early==TRUE}.
#' @template burn
#' @template temperatures
#' @param tune_temps Integer value - if greater than 0, the \code{temperatures} argument is ignored, and instead \code{tune_temps} tuned temperatures are used instead.
#' @template return_z_trace
#' @template return_x_trace
#' @param raw_only Logical value determining whether to return raw output of MCMC routine only.
#' @param swaps Number of swaps between adjacent tempered chains to perform per update cycle.
#' @param optimise_z0 Logical value determining whether to use a simulated annealing optimisation run to tune the initial values of \code{z}.
#' @param tune_omega_and_phi_proposal_sd Logical value determining whether the proposal SDs of the Metropolis-Hastings estimated parameters should be tuned for a target acceptance range.
#' @param tune_block_size Integer value giving number of samples to draw when estimatating the acceptance rate of the omega/phi proposals.
#' @template variant_weights
#' @param standardise_weights Boolean value determining whether weights should be standardised by subtracting their mean and dividing by their sample standard deviation. If \code{FALSE}, weights are untransformed.
#' @template log_phi_mean
#' @param log_phi_sd SD for normal prior on scaling factor phi. Setting to 0 causes the weights to be fixed and not estimated.
#' @template tandem_variant_updates
#' @param ... Other arguments to be passed to \code{\link{stop_chain}} and/or \code{\link{tune_proposal_sds}}.
#' @return An object of class \code{BeviMed_m}.
#' @details A \code{BeviMed_m} object is a list containing elements:
#' \itemize{
#' \item `parameters': a list containing arguments used in the function call, including the adjusted weights used in the inference in the `c_weights' slot,
#' \item `traces': a list of traces of model parameters from all MCMC chains for each parameter. Parameters sampled are z, omega, phi and x (the indicator of having a pathogenic configuration of alleles). The list of traces is named by parameter name, and each is a matrix where the rows correspond to samples. $z has k columns for each temperature, with the samples from the true posterior (i.e. with temperature equal to 1) of z corresponding to the final k columns. Likewise, the true posterior is given by the final column for the traces of phi and omega. The trace of x is only given for temperature equal to 1 to reduce memory usage.
#' \item `final': a list named by model parameter giving the final sample of each,
#' \item `swaps': a list with an element named `accept' which is a logical vector whose ith element indicates whether the ith swap between adjacent tempered chains was accepted or not, and an element named `at_temperature`, an integer vector whose ith element indicates which pair of consecutive temperatures was the ith to be proposed for swapping (giving the lowest one). 
#' }
#' @export
#' @importFrom stats rnorm runif rbeta sd
#' @importFrom methods is
#' @seealso \code{\link{bevimed_m}}, \code{\link{prob_association_m}}
#' @template paper
bevimed_m <- function(
	y,
	G,
	min_ac=1L,
	tau_shape=c(1, 1),
	pi_shape=c(6, 1),
	omega_shape=if (max(min_ac) == 1L) c(2, 8) else c(2, 2),
	samples_per_chain=1000,
	stop_early=FALSE,
	blocks=5,
	burn=as.integer(samples_per_chain/10),
	temperatures=(0:6/6)^2,
	tune_temps=0,
	return_z_trace=TRUE,
	return_x_trace=TRUE,
	raw_only=FALSE,
	swaps=as.integer(length(temperatures)/2),
	optimise_z0=FALSE,
	tune_omega_and_phi_proposal_sd=FALSE,
	tune_block_size=100,
	variant_weights=NULL,
	standardise_weights=TRUE,
	log_phi_mean=-0.15,
	log_phi_sd=sqrt(0.3),
	tandem_variant_updates=if (max(min_ac) == 1) 0 else min(sum(y), ncol(G)),
	...
) {
	stopifnot(is.matrix(G)||is(G, "sparseMatrix"))
	stopifnot(nrow(G)==length(y))
	stopifnot(!is.matrix(G)||is.numeric(G))
	stopifnot(is.logical(y))
	stopifnot(min_ac > 0)
	stopifnot(identical(as.numeric(range(temperatures)), as.numeric(c(0, 1))))

	estimate_phi <- !((log_phi_sd == 0) | is.null(variant_weights))

	c_weights <- 
		if (is.null(variant_weights)) { rep(0, ncol(G)) }
		else {
			stopifnot(is.numeric(variant_weights))
			stopifnot(length(variant_weights) == ncol(G))
			if (standardise_weights) {
				if (length(variant_weights) == 1L | sd(variant_weights) == 0)
					rep(0, ncol(G))
				else
					(variant_weights-mean(variant_weights))/sd(variant_weights)
			} else {
				variant_weights
			}
		}
		
	G_args <- get_G_args(G)

	G_logical <- matrix(data=G > 0, nrow=nrow(G), ncol=ncol(G))
	comphet_cases <- apply(G_logical, 1, sum) > 1
	comphet_variants <- lapply(split(G_logical[comphet_cases,,drop=FALSE], seq(length.out=sum(comphet_cases))), which)
	comphet_block_ends <- unname(cumsum(lapply(comphet_variants, length)))
	comphet_block_starts <- if (length(comphet_block_ends) > 0) c(0, comphet_block_ends[-length(comphet_block_ends)]) else integer(0)
	adjusted_tvu <- if (sum(comphet_cases) > 0) tandem_variant_updates else 0

	initial_log_phi_proposal_sd <- 0.5
	initial_logit_omega_proposal_sd <- 1
	estimate_omega <- !is.null(variant_weights)

	reused_arguments <- list(
		y=y,
		block_starts=G_args$block_starts,
		block_ends=G_args$block_ends,
		cases=G_args$cases-1L,
		counts=G_args$counts,
		min_ac=min_ac,
		tau_shape=tau_shape,
		pi_shape=pi_shape,
		omega_shape=omega_shape,
		estimate_omega=estimate_omega,
		variant_weights=c_weights,
		estimate_phi=estimate_phi,
		log_phi_mean=log_phi_mean,
		log_phi_sd=log_phi_sd,
		chain_swaps_per_cycle=swaps,
		tandem_variant_updates=adjusted_tvu,
		comphet_variant_block_starts=comphet_block_starts,
		comphet_variant_block_ends=comphet_block_ends,
		comphet_variants=unlist(use.names=FALSE, comphet_variants)-1
	)

	if (tune_temps > 0) {
		temperatures <- do.call(
			what=tune_temperatures,
			c(
				reused_arguments,
				list(
					samples_per_chain=samples_per_chain,
					number_of_temperatures=tune_temps,
					return_temperatures=TRUE,
					z0_matrix=matrix(runif(ncol(G)) < omega_shape[1]/sum(omega_shape), nrow=1, ncol=ncol(G)),
					logit_omegas=local({ w <- rbeta(n=1, shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
					logit_omega_proposal_sds=rep(initial_logit_omega_proposal_sd, 1),
					log_phis=if (estimate_phi) rnorm(n=1, mean=log_phi_mean, sd=log_phi_sd) else rep(0, 1),
					log_phi_proposal_sds=rep(initial_log_phi_proposal_sd, 1), 
					annealing=FALSE,
					return_z_trace=FALSE,
					return_x_trace=FALSE
				)
			)
		)
	}

	initial_z <- if (optimise_z0) {
		do.call(what=call_cpp, c(
			reused_arguments,
			list(
				samples_per_chain=samples_per_chain,
				z0_matrix=matrix(runif(ncol(G)*length(temperatures)) < omega_shape[1]/sum(omega_shape), nrow=length(temperatures), ncol=ncol(G)),
				logit_omegas=local({ w <- rbeta(n=temperatures, shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
				logit_omega_proposal_sds=rep(initial_logit_omega_proposal_sd, length(temperatures)),
				log_phis=if (estimate_phi) rnorm(n=length(temperatures), mean=log_phi_mean, sd=log_phi_sd) else rep(0, length(temperatures)),
				log_phi_proposal_sds=rep(initial_log_phi_proposal_sd, length(temperatures)), 
				temperatures=rep(1, length(temperatures)),
				chain_swaps_per_cycle=swaps,
				annealing=TRUE,
				return_z_trace=FALSE,
				return_x_trace=FALSE
			)
		))[["final"]][["z"]]
	} else {
		matrix(runif(ncol(G) * length(temperatures)) < omega_shape[1]/sum(omega_shape), nrow=length(temperatures), ncol=ncol(G))
	}

	proposal_sds <- lapply(
		setNames(nm=c("logit_omega", "log_phi")), 
		FUN=if (tune_omega_and_phi_proposal_sd & (estimate_phi | estimate_omega)) { 
			function(tune_for) do.call(what=tune_proposal_sds, c(
				reused_arguments, 
				list(
					tune_for=tune_for,
					initial_proposal_sds=rep(if (tune_for == "log_phi") initial_log_phi_proposal_sd else initial_logit_omega_proposal_sd, length(temperatures)),
					samples_per_chain=tune_block_size,
					burn=0,
					z0_matrix=initial_z,
					logit_omegas=local({ w <- rbeta(n=length(temperatures), shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
					log_phis=if (estimate_phi) rnorm(n=length(temperatures), mean=log_phi_mean, sd=log_phi_sd) else rep(0, length(temperatures)),
					temperatures=temperatures,
					annealing=FALSE,
					return_z_trace=FALSE,
					return_x_trace=FALSE
				),
				list(...)[intersect(names(list(...)), names(formals(tune_proposal_sds)))]
			))
		} else { 
			function(tune_for) rep(if (tune_for == "log_phi") initial_log_phi_proposal_sd else initial_logit_omega_proposal_sd, length(temperatures))
		})

	result <- if (stop_early) {
		burn <- do.call(what=call_cpp, c(
			reused_arguments,
			list(
				samples_per_chain=samples_per_chain,
				burn=0,
				z0_matrix=initial_z,
				logit_omegas=local({ w <- rbeta(n=length(temperatures), shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
				logit_omega_proposal_sds=proposal_sds[["logit_omega"]],
				log_phis=if (estimate_phi) rnorm(n=length(temperatures), mean=log_phi_mean, sd=log_phi_sd) else rep(0, length(temperatures)),
				log_phi_proposal_sds=proposal_sds[["log_phi"]], 
				temperatures=temperatures,
				annealing=FALSE,
				return_z_trace=FALSE,
				return_x_trace=FALSE
			)
		))

		do.call(what=stop_chain, c(
			reused_arguments,
			list(
				samples_per_chain=samples_per_chain,
				y_log_lik_t_equals_1=matrix(ncol=length(temperatures),nrow=0),
				burn=0,
				start_zs=burn[["final"]][["z"]],
				start_logit_omegas=burn[["final"]][["logit_omega"]],
				start_log_phis=burn[["final"]][["log_phi"]],
				temperatures=temperatures,
				blocks_remaining=max(1, blocks),
				logit_omega_proposal_sds=proposal_sds[["logit_omega"]],
				log_phi_proposal_sds=proposal_sds[["log_phi"]], 
				annealing=FALSE,
				return_z_trace=return_z_trace,
				return_x_trace=return_x_trace
			),
			list(...)[intersect(names(list(...)), names(formals(stop_chain)))]
		))
	} else {
		do.call(
			what=call_cpp,
			c(
				reused_arguments,
				list(
					samples_per_chain=samples_per_chain,
					burn=burn,
					z0_matrix=initial_z,
					logit_omegas=local({ w <- rbeta(n=length(temperatures), shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
					logit_omega_proposal_sds=proposal_sds[["logit_omega"]],
					log_phis=if (estimate_phi) rnorm(n=length(temperatures), mean=log_phi_mean, sd=log_phi_sd) else rep(0, length(temperatures)),
					log_phi_proposal_sds=proposal_sds[["log_phi"]], 
					temperatures=temperatures,
					annealing=FALSE,
					return_z_trace=return_z_trace,
					return_x_trace=return_x_trace
				)
			)
		)
	}

	if (raw_only)
		result
	else
		structure(
			class="BeviMed_m",
			c(
				result,
				list(parameters=list(
					omega_shape=omega_shape,
					pi_shape=pi_shape,
					tau_shape=tau_shape,
					variant_weights=variant_weights,
					temperatures=temperatures,
					estimate_phi=estimate_phi,
					estimate_omega=estimate_omega,
					y=y,
					min_ac=min_ac,
					variant_table=to_var_tab(G_args),
					c_weights=c_weights,
					G=G,
					N=length(y),
					k=ncol(G)
				))
			)
		)
}
