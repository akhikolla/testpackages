## obsolete

## converts inputs into list format for a given draw, needed for the
## simulation functions
## reshape_model_inputs <- function(X_comp, X_inter, beta, eta, beta_map, eta_map) {
##     num_groups <- dim(beta)[1]
##     num_comp <- dim(X_comp)[1]
##     num_strata <- dim(beta_map)[1]
##     eta_list <- lapply(1:num_groups, function(d) as.vector(eta[d,]))
##     beta_list <- lapply(1:num_groups, function(g) { lapply(1:num_comp, function(d) as.vector(beta[g,d,])) })
##     eta_map_list <- lapply(1:num_strata, function(s) as.vector(eta_map[s,]))
##     beta_map_list <- lapply(1:num_strata, function(s) { lapply(1:num_comp, function(d) as.vector(beta_map[s,d,])) })
##     X_comp_list <- lapply(1:num_comp, function(d) as.matrix(X_comp[d,,]))
##     list(X_comp=X_comp_list, X_inter=X_inter, beta=beta_list, eta=eta_list, beta_map=beta_map_list, eta_map=eta_map_list)
## }

##blrm_grouped_rng_auto <- function(group, n, X_comp, X_inter, beta, eta) {
##    do.call(blrm_grouped_rng, c(list(group=group, n=n), reshape_model_inputs(X_comp, X_inter, beta, eta)))
##}

##blrm_logit_grouped_auto <- function(group, stratum, X_comp, X_inter, beta, eta, beta_map, eta_map) {
##    do.call(blrm_logit_grouped, c(list(group=group, stratum=stratum), reshape_model_inputs(X_comp, X_inter, beta, eta, beta_map, eta_map)))
##}
