#include "sampling.h"

List bevimed_mc(
	int its,
	LogicalVector y,
	IntegerVector var_block_start_index,
	IntegerVector var_block_stop_index,
	IntegerVector cases,
	IntegerVector counts,
	IntegerVector min_ac,
	double tau_shape1,
	double tau_shape2,
	double pi_shape1,
	double pi_shape2,
	double z_shape1,
	double z_shape2,
	LogicalMatrix z0,
	bool estimate_logit_z_rate,
	NumericVector logit_z_rates,
	NumericVector logit_z_rate_proposal_sds,
	NumericVector z_weights,
	bool estimate_phi,
	NumericVector log_phis,
	double log_phi_mean,
	double log_phi_sd,
	NumericVector log_phi_proposal_sds,
	NumericVector t,
	int swaps,
	bool annealing,
	int tandem_variant_updates,
	IntegerVector tandem_case_block_start_index,
	IntegerVector tandem_case_block_stop_index,
	IntegerVector tandem_variants,
	bool return_z_trace,
	bool return_x_trace
) {
	double logit_z_rate_mean = z_shape1;
	double logit_z_rate_sd = z_shape2;

	int n = y.length();

	int k = var_block_start_index.length();
	int num_temps = t.length();


	IntegerVector chain_temperature_reference(num_temps);
	IntegerVector temperature_chain_reference(num_temps);
	for (int temp = 0; temp < num_temps; temp++) {
		chain_temperature_reference[temp] = temp;
		temperature_chain_reference[temp] = temp;
	}

	LogicalMatrix z_trace(return_z_trace ? its : 0, return_z_trace ? (k * num_temps) : 0);
	LogicalMatrix x_trace(return_x_trace ? its : 0, return_x_trace ? n : 0);
	NumericMatrix y_log_lik_trace(its, num_temps);
	NumericMatrix y_log_lik_t_equals_1_trace(its, num_temps);
	LogicalMatrix z(num_temps, k);

	NumericMatrix logit_z_rates_trace(its, num_temps);
	NumericMatrix log_phis_trace(its, num_temps);

	LogicalVector swap_accept_trace(swaps * its);
	IntegerVector swap_temp1_trace(swaps * its);

	IntegerVector count_z1(num_temps, 0);
	for (int temp = 0; temp < num_temps; temp++) {
		for (int v = 0; v < k; v++) {
			z(temp, v) = z0(temp, v); 
			count_z1[temp] += (int)z(temp, v);	
		}
	}

	IntegerMatrix pathogenic_var_count(num_temps, n);
	IntegerVector temporary_pathogenic_var_count(n, 0);
	LogicalVector temporary_counted_indicator(n, false);

	LogicalMatrix x(num_temps, n);
	for (int temp = 0; temp < num_temps; temp++) {
		for (int i = 0; i < n; i++) {
			pathogenic_var_count(temp, i) = 0;
			x(temp, i) = false;
		}
	}

	IntegerVector count_y1(num_temps, 0);
	IntegerVector count_x1(num_temps, 0);
	IntegerVector count_y1x1(num_temps, 0);

	for (int temp = 0; temp < num_temps; temp++) {
		for (int v = 0; v < k; v++) {
			for (int case_with_v = var_block_start_index[v]; case_with_v < var_block_stop_index[v]; case_with_v++) {
				pathogenic_var_count(temp, cases[case_with_v]) += (int)z(temp, v) * counts[case_with_v];
			}
		}

		for (int i = 0; i < n; i++) {
			x(temp, i) = pathogenic_var_count(temp, i) >= min_ac[i];
			count_x1[temp] += (int)x(temp, i);
			if (y[i]) { 
				count_y1x1[temp] += (int)x(temp, i);
				count_y1[temp] += 1;
			}
		}
	}

	NumericVector q1_tab(n+1);
	NumericVector q2_tab(n+1);
	NumericVector p1_tab(n+1);
	NumericVector p2_tab(n+1);
	NumericVector q_tot_tab(n+1);
	NumericVector p_tot_tab(n+1);
	q1_tab[0] = 0.0;
	q2_tab[0] = 0.0;
	p1_tab[0] = 0.0;
	p2_tab[0] = 0.0;
	q_tot_tab[0] = 0.0;
	p_tot_tab[0] = 0.0;
	for (int i = 1; i <= n; i++) {
		q1_tab[i] = q1_tab[i-1] + log(i-1 + tau_shape1);
		q2_tab[i] = q2_tab[i-1] + log(i-1 + tau_shape2);
		p1_tab[i] = p1_tab[i-1] + log(i-1 + pi_shape1);
		p2_tab[i] = p2_tab[i-1] + log(i-1 + pi_shape2);
		q_tot_tab[i] = q_tot_tab[i-1] + log(i-1 + tau_shape1 + tau_shape2);
		p_tot_tab[i] = p_tot_tab[i-1] + log(i-1 + pi_shape1 + pi_shape2);
	}

	for (int it = 0; it < its; it++) {
		double annealing_factor = annealing ? (double)(its-it) : 1.0;
		for (int chain_number = 0; chain_number < num_temps; chain_number++) {
			if (estimate_logit_z_rate) {
				double proposal = logit_z_rates[chain_number] + norm_rand() * logit_z_rate_proposal_sds[chain_number];

				double ll_cur = logit_beta(logit_z_rate_mean, logit_z_rate_sd, logit_z_rates[chain_number]);
				double ll_upd = logit_beta(logit_z_rate_mean, logit_z_rate_sd, proposal);
				for (int v = 0; v < k; v++) {
					double pv_cur = expit(exp(log_phis[chain_number]) * z_weights[v] + logit_z_rates[chain_number]);
					double pv_upd = expit(exp(log_phis[chain_number]) * z_weights[v] + proposal);
					ll_cur += z(chain_number, v) ? log(pv_cur) : log(1.0-pv_cur);
					ll_upd += z(chain_number, v) ? log(pv_upd) : log(1.0-pv_upd);
				}
				if (log(unif_rand()) < (ll_upd-ll_cur)) {
					logit_z_rates[chain_number] = proposal;
				}

				if (estimate_phi) {
					double phi_proposal = log_phis[chain_number] + norm_rand() * log_phi_proposal_sds[chain_number];

					double ll_phi_cur = log_likelihood_normal(log_phi_mean, log_phi_sd, log_phis[chain_number]);
					double ll_phi_upd = log_likelihood_normal(log_phi_mean, log_phi_sd, phi_proposal);
					for (int v = 0; v < k; v++) {
						double pv_cur = expit(exp(log_phis[chain_number]) * z_weights[v] + logit_z_rates[chain_number]);
						double pv_upd = expit(exp(phi_proposal) * z_weights[v] + logit_z_rates[chain_number]);
						ll_phi_cur += z(chain_number, v) ? log(pv_cur) : log(1.0-pv_cur);
						ll_phi_upd += z(chain_number, v) ? log(pv_upd) : log(1.0-pv_upd);
					}
					if (log(unif_rand()) < (ll_phi_upd-ll_phi_cur)) {
						log_phis[chain_number] = phi_proposal;
					}

				}
			}
			for (int v = 0; v < k; v++) {
				int dx1 = 0;
				int dy1x1 = 0;

				int c_s01 = count_y1[chain_number] - count_y1x1[chain_number];
				int c_s00 = (n - count_x1[chain_number]) - (count_y1[chain_number] - count_y1x1[chain_number]);
				int c_s11 = count_y1x1[chain_number];
				int c_s10 = count_x1[chain_number] - count_y1x1[chain_number];

				for (int i = var_block_start_index[v]; i < var_block_stop_index[v]; i++) {
					int alt_pathogenic_var_count = pathogenic_var_count(chain_number, cases[i]) + (z(chain_number, v) ? (-counts[i]) : counts[i]);
					bool alt_x = alt_pathogenic_var_count >= min_ac[cases[i]];
					if (alt_x != x(chain_number, cases[i])) {
						dx1 += (int)alt_x - (int)x(chain_number, cases[i]);
						if (y[cases[i]]) dy1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
					}
				}

				int alt_count_x1 = count_x1[chain_number] + dx1;
				int alt_count_y1x1 = count_y1x1[chain_number] + dy1x1;

				int u_s01 = count_y1[chain_number] - alt_count_y1x1;
				int u_s00 = (n - alt_count_x1) - (count_y1[chain_number] - alt_count_y1x1);
				int u_s11 = alt_count_y1x1;
				int u_s10 = alt_count_x1 - alt_count_y1x1;

				double z_contribution;
				if (estimate_logit_z_rate) {
					int dz1 = z(chain_number, v) ? -1 : 1;
					double vp = expit(logit_z_rates[chain_number] + exp(log_phis[chain_number]) * z_weights[v]);
					z_contribution = + dz1 * log(1.0-vp) - dz1 * log(vp); 
				} 
				else {
					z_contribution = + ((z(chain_number, v)) ? (log(count_z1[chain_number] + z_shape1 - 1.0) - log(k - count_z1[chain_number] + z_shape2)) : (-log(count_z1[chain_number] + z_shape1) + log(k - count_z1[chain_number] + z_shape2 - 1.0)));
				}

				double current_zv_log_odds = 
					+ z_contribution
					+ (
						+ (q1_tab[c_s01]-q1_tab[u_s01])
						+ (q2_tab[c_s00]-q2_tab[u_s00])
						+ (p1_tab[c_s11]-p1_tab[u_s11])
						+ (p2_tab[c_s10]-p2_tab[u_s10])
						+ (q_tot_tab[u_s01 + u_s00]-q_tot_tab[c_s01 + c_s00])
						+ (p_tot_tab[u_s11 + u_s10]-p_tot_tab[c_s11 + c_s10])
					) * t[chain_temperature_reference[chain_number]] * annealing_factor
				;

				if (unif_rand() < (1.0-1.0/(1.0+exp(-current_zv_log_odds)))) {
					count_x1[chain_number] += dx1;
					count_y1x1[chain_number] += dy1x1;
					for (int i = var_block_start_index[v]; i < var_block_stop_index[v]; i++) {
						pathogenic_var_count(chain_number, cases[i]) += (z(chain_number, v) ? (-counts[i]) : counts[i]);
						x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
					}
					count_z1[chain_number] += z(chain_number, v) ? (-1) : 1;
					z(chain_number, v) = !z(chain_number, v);
				}
			}
		}

		if (tandem_variant_updates > 0) {
			for (int chain_number = 0; chain_number < num_temps; chain_number++) {
				for (int double_variant_update_number = 0; double_variant_update_number < tandem_variant_updates; double_variant_update_number ++) {

					int select_case = random_integer(tandem_case_block_start_index.length());
					int num_vars_in_case = tandem_case_block_stop_index[select_case]-tandem_case_block_start_index[select_case];
					int v1p = random_integer(num_vars_in_case);
					int v2p = random_integer(num_vars_in_case-1) + 1;
					if (v1p >= v2p) {
						v1p = num_vars_in_case - v1p - 1;
						v2p = num_vars_in_case - v2p;
					}
					int v1 = tandem_variants[tandem_case_block_start_index[select_case] + v1p];
					int v2 = tandem_variants[tandem_case_block_start_index[select_case] + v2p];

					int no_change_s01 = count_y1[chain_number] - count_y1x1[chain_number];
					int no_change_s00 = (n - count_x1[chain_number]) - (count_y1[chain_number] - count_y1x1[chain_number]);
					int no_change_s11 = count_y1x1[chain_number];
					int no_change_s10 = count_x1[chain_number] - count_y1x1[chain_number];


					int v1_change_count_x1 = count_x1[chain_number];
					int v1_change_count_y1x1 = count_y1x1[chain_number];

					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
						int alt_pathogenic_var_count = pathogenic_var_count(chain_number, cases[i]) + (z(chain_number, v1) ? (-counts[i]) : counts[i]);
						bool alt_x = alt_pathogenic_var_count >= min_ac[cases[i]];
						if (alt_x != x(chain_number, cases[i])) {
							v1_change_count_x1 += (int)alt_x - (int)x(chain_number, cases[i]);
							if (y[cases[i]]) v1_change_count_y1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
						}
					}

					int v1_change_s01 = count_y1[chain_number] - v1_change_count_y1x1;
					int v1_change_s00 = (n - v1_change_count_x1) - (count_y1[chain_number] - v1_change_count_y1x1);
					int v1_change_s11 = v1_change_count_y1x1;
					int v1_change_s10 = v1_change_count_x1 - v1_change_count_y1x1;

					int v2_change_count_x1 = count_x1[chain_number];
					int v2_change_count_y1x1 = count_y1x1[chain_number];

					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
						int alt_pathogenic_var_count = pathogenic_var_count(chain_number, cases[i]) + (z(chain_number, v2) ? (-counts[i]) : counts[i]);
						bool alt_x = alt_pathogenic_var_count >= min_ac[cases[i]];
						if (alt_x != x(chain_number, cases[i])) {
							v2_change_count_x1 += (int)alt_x - (int)x(chain_number, cases[i]);
							if (y[cases[i]]) v2_change_count_y1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
						}
					}

					int v2_change_s01 = count_y1[chain_number] - v2_change_count_y1x1;
					int v2_change_s00 = (n - v2_change_count_x1) - (count_y1[chain_number] - v2_change_count_y1x1);
					int v2_change_s11 = v2_change_count_y1x1;
					int v2_change_s10 = v2_change_count_x1 - v2_change_count_y1x1;

					int v1_and_v2_change_count_x1 = count_x1[chain_number];
					int v1_and_v2_change_count_y1x1 = count_y1x1[chain_number];
					
					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++)
						temporary_pathogenic_var_count[cases[i]] = pathogenic_var_count(chain_number, cases[i]);
					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++)
						temporary_pathogenic_var_count[cases[i]] = pathogenic_var_count(chain_number, cases[i]);
					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++)
						temporary_pathogenic_var_count[cases[i]] += (z(chain_number, v1) ? (-counts[i]) : counts[i]);
					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++)
						temporary_pathogenic_var_count[cases[i]] += (z(chain_number, v2) ? (-counts[i]) : counts[i]);

					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
						bool alt_x = temporary_pathogenic_var_count[cases[i]] >= min_ac[cases[i]];
						if ((!temporary_counted_indicator[cases[i]]) && (alt_x != x(chain_number, cases[i]))) {
							v1_and_v2_change_count_x1 += (int)alt_x - (int)x(chain_number, cases[i]);
							if (y[cases[i]]) v1_and_v2_change_count_y1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
						}
						temporary_counted_indicator[cases[i]] = true;
					}

					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
						bool alt_x = temporary_pathogenic_var_count[cases[i]] >= min_ac[cases[i]];
						if ((!temporary_counted_indicator[cases[i]]) && (alt_x != x(chain_number, cases[i]))) {
							v1_and_v2_change_count_x1 += (int)alt_x - (int)x(chain_number, cases[i]);
							if (y[cases[i]]) v1_and_v2_change_count_y1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
						}
						temporary_counted_indicator[cases[i]] = true;
					}

					int v1_and_v2_change_s01 = count_y1[chain_number] - v1_and_v2_change_count_y1x1;
					int v1_and_v2_change_s00 = (n - v1_and_v2_change_count_x1) - (count_y1[chain_number] - v1_and_v2_change_count_y1x1);
					int v1_and_v2_change_s11 = v1_and_v2_change_count_y1x1;
					int v1_and_v2_change_s10 = v1_and_v2_change_count_x1 - v1_and_v2_change_count_y1x1;

					double v1_z_log_odds_favour_current;
					if (estimate_logit_z_rate) {
						int dz1 = z(chain_number, v1) ? -1 : 1;
						double vp = expit(logit_z_rates[chain_number] + exp(log_phis[chain_number]) * z_weights[v1]);
						v1_z_log_odds_favour_current = + dz1 * log(1.0-vp) - dz1 * log(vp); 
					}
					else {
						v1_z_log_odds_favour_current = ((z(chain_number, v1)) ? (log(count_z1[chain_number] + z_shape1 - 1.0) - log(k - count_z1[chain_number] + z_shape2)) : (-log(count_z1[chain_number] + z_shape1) + log(k - count_z1[chain_number] + z_shape2 - 1.0)));
					}
					double v1_change_odds_ratio = exp((
						+ v1_z_log_odds_favour_current 
						+ (
							+ (q1_tab[no_change_s01]-q1_tab[v1_change_s01])
							+ (q2_tab[no_change_s00]-q2_tab[v1_change_s00])
							+ (p1_tab[no_change_s11]-p1_tab[v1_change_s11])
							+ (p2_tab[no_change_s10]-p2_tab[v1_change_s10])
							+ (q_tot_tab[v1_change_s01 + v1_change_s00]-q_tot_tab[no_change_s01 + no_change_s00])
							+ (p_tot_tab[v1_change_s11 + v1_change_s10]-p_tot_tab[no_change_s11 + no_change_s10])
						) * t[chain_temperature_reference[chain_number]] * annealing_factor
					) * (-1.0));

					double v2_z_log_odds_favour_current;
					if (estimate_logit_z_rate) {
						int dz1 = z(chain_number, v2) ? -1 : 1;
						double vp = expit(logit_z_rates[chain_number] + exp(log_phis[chain_number]) * z_weights[v2]);
						v2_z_log_odds_favour_current = + dz1 * log(1.0-vp) - dz1 * log(vp); 
					}
					else {
						v2_z_log_odds_favour_current = ((z(chain_number, v2)) ? (log(count_z1[chain_number] + z_shape1 - 1.0) - log(k - count_z1[chain_number] + z_shape2)) : (-log(count_z1[chain_number] + z_shape1) + log(k - count_z1[chain_number] + z_shape2 - 1.0)));
					}
					double v2_change_odds_ratio = exp((
						+ v2_z_log_odds_favour_current
						+ (
							+ (q1_tab[no_change_s01]-q1_tab[v2_change_s01])
							+ (q2_tab[no_change_s00]-q2_tab[v2_change_s00])
							+ (p1_tab[no_change_s11]-p1_tab[v2_change_s11])
							+ (p2_tab[no_change_s10]-p2_tab[v2_change_s10])
							+ (q_tot_tab[v2_change_s01 + v2_change_s00]-q_tot_tab[no_change_s01 + no_change_s00])
							+ (p_tot_tab[v2_change_s11 + v2_change_s10]-p_tot_tab[no_change_s11 + no_change_s10])
						) * t[chain_temperature_reference[chain_number]] * annealing_factor
					) * (-1.0));



					double v1_and_v2_z_log_odds_favour_current;
					if (estimate_logit_z_rate) {
						int dz1_v1 = z(chain_number, v1) ? -1 : 1;
						int dz1_v2 = z(chain_number, v2) ? -1 : 1;
						double v1p = expit(logit_z_rates[chain_number] + exp(log_phis[chain_number]) * z_weights[v1]);
						double v2p = expit(logit_z_rates[chain_number] + exp(log_phis[chain_number]) * z_weights[v2]);
						v1_and_v2_z_log_odds_favour_current = 
							+ dz1_v1 * log(1.0-v1p) - dz1_v1 * log(v1p)
							+ dz1_v2 * log(1.0-v2p) - dz1_v2 * log(v2p)
						;
					}
					else {
						if (z(chain_number, v1) != z(chain_number, v2)) {
							v1_and_v2_z_log_odds_favour_current = 0;
						}
						else {
							if (z(chain_number, v1)) {
								v1_and_v2_z_log_odds_favour_current = (
									+ log(z_shape1 + count_z1[chain_number] - 2.0)
									+ log(z_shape1 + count_z1[chain_number] - 1.0)
									- log(z_shape2 + k - count_z1[chain_number] + 1.0)
									- log(z_shape2 + k - count_z1[chain_number])
								);
							}
							else {
								v1_and_v2_z_log_odds_favour_current = (
									- log(z_shape1 + count_z1[chain_number] + 1.0) 
									- log(z_shape1 + count_z1[chain_number])
									+ log(z_shape2 + k - count_z1[chain_number] - 1.0)
									+ log(z_shape2 + k - count_z1[chain_number] - 2.0)
								);
							}
						}
					}
					double v1_and_v2_change_odds_ratio = exp((
						+ v1_and_v2_z_log_odds_favour_current
						+ (
							+ (q1_tab[no_change_s01]-q1_tab[v1_and_v2_change_s01])
							+ (q2_tab[no_change_s00]-q2_tab[v1_and_v2_change_s00])
							+ (p1_tab[no_change_s11]-p1_tab[v1_and_v2_change_s11])
							+ (p2_tab[no_change_s10]-p2_tab[v1_and_v2_change_s10])
							+ (q_tot_tab[v1_and_v2_change_s01 + v1_and_v2_change_s00]-q_tot_tab[no_change_s01 + no_change_s00])
							+ (p_tot_tab[v1_and_v2_change_s11 + v1_and_v2_change_s10]-p_tot_tab[no_change_s11 + no_change_s10])
						) * t[chain_temperature_reference[chain_number]] * annealing_factor
					) * (-1.0));

					double sum_ratios = 
						+ v1_change_odds_ratio
						+ v2_change_odds_ratio
						+ v1_and_v2_change_odds_ratio
					;

					double p_no_change = 1.0/(1.0+sum_ratios);
					double p_v1_change = p_no_change * v1_change_odds_ratio;
					double p_v2_change = p_no_change * v2_change_odds_ratio;


					double p_v2_threshold = p_no_change + p_v1_change + p_v2_change;
					double random_draw = unif_rand();



					if (random_draw >= p_no_change) {
						if (random_draw < (p_no_change + p_v1_change)) {
							count_x1[chain_number] = v1_change_count_x1;
							count_y1x1[chain_number] = v1_change_count_y1x1;
							for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
								pathogenic_var_count(chain_number, cases[i]) += (z(chain_number, v1) ? (-counts[i]) : counts[i]);
								x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
							}
							count_z1[chain_number] += z(chain_number, v1) ? (-1) : 1;
							z(chain_number, v1) = !z(chain_number, v1);
						}
						else if (random_draw < p_v2_threshold) {
							count_x1[chain_number] = v2_change_count_x1;
							count_y1x1[chain_number] = v2_change_count_y1x1;
							for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
								pathogenic_var_count(chain_number, cases[i]) += (z(chain_number, v2) ? (-counts[i]) : counts[i]);
								x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
							}
							count_z1[chain_number] += z(chain_number, v2) ? (-1) : 1;
							z(chain_number, v2) = !z(chain_number, v2);
						}
						else {
							count_x1[chain_number] = v1_and_v2_change_count_x1;
							count_y1x1[chain_number] = v1_and_v2_change_count_y1x1;
							for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
								pathogenic_var_count(chain_number, cases[i]) = temporary_pathogenic_var_count[cases[i]];
								x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
							}
							count_z1[chain_number] += z(chain_number, v1) ? (-1) : 1;
							z(chain_number, v1) = !z(chain_number, v1);
							for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
								pathogenic_var_count(chain_number, cases[i]) = temporary_pathogenic_var_count[cases[i]];
								x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
							}
							count_z1[chain_number] += z(chain_number, v2) ? (-1) : 1;
							z(chain_number, v2) = !z(chain_number, v2);
						}
					}

					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
						temporary_counted_indicator[cases[i]] = false;
						temporary_pathogenic_var_count[cases[i]] = 0;
					}
					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
						temporary_counted_indicator[cases[i]] = false;
						temporary_pathogenic_var_count[cases[i]] = 0;
					}
				}
			}
		}

		if ((swaps > 0) && (num_temps > 1)) {
			for (int swap_no = 0; swap_no < swaps; swap_no++) {
				int temperature_number1 = random_integer(num_temps-1);
				int chain1 = temperature_chain_reference[temperature_number1];
				int chain2 = temperature_chain_reference[temperature_number1+1];

				int s01_1 = count_y1[chain1] - count_y1x1[chain1];
				int s00_1 = (n - count_x1[chain1]) - (count_y1[chain1] - count_y1x1[chain1]);
				int s11_1 = count_y1x1[chain1];
				int s10_1 = count_x1[chain1] - count_y1x1[chain1];

				int s01_2 = count_y1[chain2] - count_y1x1[chain2];
				int s00_2 = (n - count_x1[chain2]) - (count_y1[chain2] - count_y1x1[chain2]);
				int s11_2 = count_y1x1[chain2];
				int s10_2 = count_x1[chain2] - count_y1x1[chain2];


				double chain1_y_log_lik_t_equals_1 =
					+ (q1_tab[s01_1]-q1_tab[0])
					+ (q2_tab[s00_1]-q2_tab[0])
					+ (p1_tab[s11_1]-p1_tab[0])
					+ (p2_tab[s10_1]-p2_tab[0])
					+ (q_tot_tab[0 + 0]-q_tot_tab[s01_1 + s00_1])
					+ (p_tot_tab[0 + 0]-p_tot_tab[s11_1 + s10_1])
				;
				
				double chain2_y_log_lik_t_equals_1 =
					+ (q1_tab[s01_2]-q1_tab[0])
					+ (q2_tab[s00_2]-q2_tab[0])
					+ (p1_tab[s11_2]-p1_tab[0])
					+ (p2_tab[s10_2]-p2_tab[0])
					+ (q_tot_tab[0 + 0]-q_tot_tab[s01_2 + s00_2])
					+ (p_tot_tab[0 + 0]-p_tot_tab[s11_2 + s10_2])
				;

				
				double logA = (
					- (t[chain_temperature_reference[chain1]] * chain1_y_log_lik_t_equals_1 + t[chain_temperature_reference[chain2]] * chain2_y_log_lik_t_equals_1)
					+ (t[chain_temperature_reference[chain2]] * chain1_y_log_lik_t_equals_1 + t[chain_temperature_reference[chain1]] * chain2_y_log_lik_t_equals_1)
				) * annealing_factor;

				swap_temp1_trace(it * swaps + swap_no) = chain_temperature_reference[chain1]+1;

				if (log(unif_rand()) < logA) {
					int temp_ref1 = chain_temperature_reference[chain1];
					chain_temperature_reference[chain1] = chain_temperature_reference[chain2];
					chain_temperature_reference[chain2] = temp_ref1;

					temperature_chain_reference[chain_temperature_reference[chain2]] = chain2;
					temperature_chain_reference[chain_temperature_reference[chain1]] = chain1;

					swap_accept_trace(it * swaps + swap_no) = true;

				} else {
					swap_accept_trace(it * swaps + swap_no) = false;
				}

			}
		}

		for (int chain_number = 0; chain_number < num_temps; chain_number++) {
			double y_log_lik_t_equals_1;
			
			int s01 = count_y1[chain_number] - count_y1x1[chain_number];
			int s00 = (n - count_x1[chain_number]) - (count_y1[chain_number] - count_y1x1[chain_number]);
			int s11 = count_y1x1[chain_number];
			int s10 = count_x1[chain_number] - count_y1x1[chain_number];

			y_log_lik_t_equals_1 =
				+ (q1_tab[s01]-q1_tab[0])
				+ (q2_tab[s00]-q2_tab[0])
				+ (p1_tab[s11]-p1_tab[0])
				+ (p2_tab[s10]-p2_tab[0])
				+ (q_tot_tab[0 + 0]-q_tot_tab[s01 + s00])
				+ (p_tot_tab[0 + 0]-p_tot_tab[s11 + s10])
			;

			double y_log_lik = t[chain_temperature_reference[chain_number]] * y_log_lik_t_equals_1 * annealing_factor;

			y_log_lik_trace(it, chain_temperature_reference[chain_number]) = y_log_lik;
			y_log_lik_t_equals_1_trace(it, chain_temperature_reference[chain_number]) = y_log_lik_t_equals_1;
			logit_z_rates_trace(it, chain_temperature_reference[chain_number]) = logit_z_rates[chain_number];
			log_phis_trace(it, chain_temperature_reference[chain_number]) = log_phis[chain_number];

			if (return_z_trace) {
				for (int v = 0; v < k; v++) {
					z_trace(it, chain_temperature_reference[chain_number] * k + v) = z(chain_number, v);
				}
			}
		}

		if (return_x_trace) {
			for (int samp = 0; samp < n; samp++) {
				x_trace(it, samp) = x(temperature_chain_reference[num_temps-1], samp);
			}
		}
	}

	LogicalMatrix terminal_z(num_temps, k);
	for (int chain_number = 0; chain_number < num_temps; chain_number++)
		for (int j = 0; j < k; j++)
			terminal_z(chain_temperature_reference[chain_number], j) = z(chain_number, j);

	NumericVector terminal_log_phi(num_temps);
	NumericVector terminal_logit_omega(num_temps);
	for (int chain_number = 0; chain_number < num_temps; chain_number++) {
		terminal_log_phi[chain_number] = log_phis_trace(its-1, chain_number);
		terminal_logit_omega[chain_number] = logit_z_rates_trace(its-1, chain_number);
	}

	return List::create(
		Named("traces")=List::create(
			Named("y_log_lik")=y_log_lik_trace,
			Named("y_log_lik_t_equals_1")=y_log_lik_t_equals_1_trace,
			Named("z")=z_trace,
			Named("x")=x_trace,
			Named("logit_omega")=logit_z_rates_trace,
			Named("log_phi")=log_phis_trace
		),
		Named("swaps")=List::create(
			Named("accept")=swap_accept_trace,
			Named("at_temperature")=swap_temp1_trace
		),
		Named("final")=List::create(
			Named("z")=terminal_z,
			Named("log_phi")=terminal_log_phi,
			Named("logit_omega")=terminal_logit_omega
		)
	);
}
