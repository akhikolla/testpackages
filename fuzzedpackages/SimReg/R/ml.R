#' Calculate sum of log probabilities on log scale without over/under-flow
#'
#' @param log_probs Numeric vector of probabilities on log scale.
#' @return Numeric value on log scale.
#' @export
sum_log_probs <- function(log_probs) log(sum(exp(log_probs-max(log_probs))))+max(log_probs)

cumsumgrid <- function(x) apply(apply(x, 1, cumsum), 1, cumsum)

#' @importFrom stats dbeta
fg_step <- function( s_phi, s_x, y, x_values=seq(from=0, to=1, length.out=30), tf_s1=3, tf_s2=1, tg_s1=1, tg_s2=1, q_s1=mean(y), q_s2=1, p_s1=3, p_s2=1) {
	tf_probs <- local({ h0 <- dbeta(x_values[-1], shape1=tf_s1, shape2=tf_s2); h <- c(h0[1], h0); h/sum(h) })
	tg_probs <- local({ h0 <- dbeta(x_values[-1], shape1=tg_s1, shape2=tg_s2); h <- c(h0[1], h0); h/sum(h) })

	y1_tab <- table(
		data.frame(
			phi=cut(s_phi[y], breaks=c(-Inf, x_values)),
			x=cut(s_x[y], breaks=c(-Inf, x_values))))
	y0_tab <- table(
		data.frame(
			phi=cut(s_phi[!y], breaks=c(-Inf, x_values)),
			x=cut(s_x[!y], breaks=c(-Inf, x_values))))

	n <- length(x_values)
	sq <- seq(from=n, to=1, by=-1)
	y1 <- cumsumgrid(y1_tab[sq,sq])[sq,sq]
	y0 <- cumsumgrid(y0_tab[sq,sq])[sq,sq]
	liks <- outer(FUN="+", log(tf_probs), log(tg_probs)) + lbeta(y1 + p_s1, y0 + p_s2) - lbeta(p_s1, p_s2) + lbeta(y1[1]-y1 + q_s1, y0[1]-y0 + q_s2) - lbeta(q_s1, q_s2)
	sum_log_probs(liks)
}

#' @importFrom stats dbeta
f_step <- function(x, y, x_values=seq(from=0, to=1, length.out=30), t_s1=3, t_s2=1, q_s1=mean(y), q_s2=1, p_s1=3, p_s2=1) {
	t_probs <- local({ h0 <- dbeta(x_values[-1], shape1=t_s1, shape2=t_s2); h <- c(h0[1], h0); h/sum(h) })
	y1 <- rev(cumsum(rev(as.integer(table(cut(x[y], breaks=c(-Inf, x_values)))))))
	y0 <- rev(cumsum(rev(as.integer(table(cut(x[!y], breaks=c(-Inf, x_values)))))))
	liks <- log(t_probs) + lbeta(y1 + p_s1, y0 + p_s2) - lbeta(p_s1, p_s2) + lbeta(y1[1]-y1 + q_s1, y0[1]-y0 + q_s2) - lbeta(q_s1, q_s2)
	sum_log_probs(liks)
}

#' @useDynLib SimReg 
fg <- function(
	s_phi, 
	s_x, 
	y, 
	x_values=seq(from=0, to=1, length.out=20), 
	t=(0:100/100)^2, 
	max_samples=30, 
	min_samples=7, 
	min_log_ML=-Inf,
	log_scale_tolerance=-2.3,
	alpha_mean=log(mean(y)), 
	alpha_sd=1, 
	log_beta_mean=1, 
	log_beta_sd=1, 
	logit_f_mean=1.5, 
	logit_f_sd=0.5, 
	log_f_a_plus_b_mean=1, 
	log_f_a_plus_b_sd=0.5, 
	logit_g_mean=0, 
	logit_g_sd=1, 
	log_g_a_plus_b_mean=1, 
	log_g_a_plus_b_sd=1, 
	alpha_prop_sd=0.8, 
	log_beta_prop_sd=0.8, 
	logit_f_mean_prop_sd=0.8, 
	log_f_a_plus_b_prop_sd=0.8, 
	logit_g_mean_prop_sd=0.8, 
	log_g_a_plus_b_prop_sd=0.8
) { 

	y1_tab <- table(
		data.frame(
			phi=cut(s_phi[y], breaks=c(-Inf, x_values)),
			x=cut(s_x[y], breaks=c(-Inf, x_values))))
	y0_tab <- table(
		data.frame(
			phi=cut(s_phi[!y], breaks=c(-Inf, x_values)),
			x=cut(s_x[!y], breaks=c(-Inf, x_values))))

	tab <- y1_tab + y0_tab

	who <- which(arr.ind=TRUE, tab > 0)

	ML(
		s_phi_values=x_values[who[,1]],
		s_x_values=x_values[who[,2]],
		num_y0_phi=y0_tab[who],
		num_y1_phi=y1_tab[who],
		t=t,
		log_scale_tolerance=log_scale_tolerance,
		max_samples=max_samples,
		min_samples=min_samples,
		min_log_ML=min_log_ML,
		alpha_mean=alpha_mean,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_f_mean=logit_f_mean,
		logit_f_sd=logit_f_sd,
		log_f_a_plus_b_mean=log_f_a_plus_b_mean,
		log_f_a_plus_b_sd=log_f_a_plus_b_sd,
		logit_g_mean=logit_g_mean,
		logit_g_sd=logit_g_sd,
		log_g_a_plus_b_mean=log_g_a_plus_b_mean,
		log_g_a_plus_b_sd=log_g_a_plus_b_sd,
		alpha_prop_sd=alpha_prop_sd,
		log_beta_prop_sd=log_beta_prop_sd,
		logit_f_mean_prop_sd=logit_f_mean_prop_sd,
		log_f_a_plus_b_prop_sd=log_f_a_plus_b_prop_sd,
		logit_g_mean_prop_sd=logit_g_mean_prop_sd,
		log_g_a_plus_b_prop_sd=log_g_a_plus_b_prop_sd
	)
}

#' @useDynLib SimReg 
f <- function(
	x, 
	y, 
	x_values=seq(from=0, to=1, length.out=20), 
	t=(0:100/100)^2, 
	max_samples=50, 
	min_samples=10, 
	min_log_ML=-Inf,
	log_scale_tolerance=-2.3,
	alpha_mean=log(mean(y)), 
	alpha_sd=1, 
	log_beta_mean=1, 
	log_beta_sd=1, 
	logit_f_mean=1.5, 
	logit_f_sd=0.5, 
	log_f_a_plus_b_mean=1, 
	log_f_a_plus_b_sd=0.5, 
	alpha_prop_sd=0.8, 
	log_beta_prop_sd=0.8, 
	logit_f_mean_prop_sd=0.8, 
	log_f_a_plus_b_prop_sd=0.8
) { 
	y1 <- as.integer(table(cut(x[y], breaks=c(-Inf, x_values))))
	y0 <- as.integer(table(cut(x[!y], breaks=c(-Inf, x_values))))
	f_ML(
		x_values=x_values,
		num_y0_phi=y0,
		num_y1_phi=y1,
		t=t,
		log_scale_tolerance=log_scale_tolerance,
		max_samples=max_samples,
		min_samples=min_samples,
		min_log_ML=min_log_ML,
		alpha_mean=alpha_mean,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_f_mean=logit_f_mean,
		logit_f_sd=logit_f_sd,
		log_f_a_plus_b_mean=log_f_a_plus_b_mean,
		log_f_a_plus_b_sd=log_f_a_plus_b_sd,
		alpha_prop_sd=alpha_prop_sd,
		log_beta_prop_sd=log_beta_prop_sd,
		logit_f_mean_prop_sd=logit_f_mean_prop_sd,
		log_f_a_plus_b_prop_sd=log_f_a_plus_b_prop_sd
	)
}

bg_rate <- function(y, q_s1=mean(y), q_s2=1) {
	liks <- lbeta(sum(y) + q_s1, sum(!y) + q_s2) - lbeta(q_s1, q_s2)
	log(sum(exp(liks-min(liks))))+min(liks)
}

#' @useDynLib SimReg 
bg <- function(y, alpha_mean=log(mean(y)), alpha_sd=1, t=(0:100/100)^2, n_samples=300, alpha_prop_sd=1) bg_ML(y0=sum(!y), y1=sum(y), t=t, n_samples=n_samples, alpha_mean=alpha_mean, alpha_sd=alpha_sd, alpha_prop_sd=alpha_prop_sd)
