CombIncrease_sim = function(ndose_a1, ndose_a2, p_tox, target, target_min, target_max, prior_tox_a1, prior_tox_a2, n_cohort,
                        cohort, tite=FALSE, time_full=0, poisson_rate=0, nsim, c_e=0.85, c_d=0.45, c_stop=0.95, c_t=0.5,
                        c_over=0.25, cmin_overunder=2, cmin_mtd=3, cmin_recom=1, startup=1, alloc_rule=1, early_stop=1,
                        nburn=2000, niter=5000, seed=14061991){
  c_d = 1-c_d
  dim_ptox = dim(p_tox)

  if(dim_ptox[1] != ndose_a1 || dim_ptox[2] != ndose_a2){
    stop("Wrong dimension of the matrix for true toxicity probabilities.")
  }
  n_prior_tox_a1 = length(prior_tox_a1)
  if(n_prior_tox_a1 != ndose_a1){
    stop("The entered vector of initial guessed toxicity probabities for agent 1 is of wrong length.")
  }
  n_prior_tox_a2 = length(prior_tox_a2)
  if(n_prior_tox_a2 != ndose_a2){
    stop("The entered vector of initial guessed toxicity probabities for agent 2 is of wrong length.")
  }


  ndose_a1 = as.integer(ndose_a1)[1]
  ndose_a2 = as.integer(ndose_a2)[1]
  target = as.double(target)[1]
  target_min = as.double(target_min)[1]
  target_max = as.double(target_max)[1]
  prior_tox_a1 = as.double(prior_tox_a1)
  prior_tox_a2 = as.double(prior_tox_a2)
  n_cohort = as.integer(n_cohort)[1]
  cohort = as.integer(cohort)[1]
  tite = as.logical(tite)[1]
  time_full = as.double(time_full)[1]
  poisson_rate = as.double(poisson_rate)[1]
  nsim = as.integer(nsim)[1]
  c_e = as.double(c_e)[1]
  c_d = as.double(c_d)[1]
  c_stop = as.double(c_stop)[1]
  c_t = as.double(c_t)[1]
  c_over = as.double(c_over)[1]
  cmin_overunder = as.integer(cmin_overunder)[1]
  cmin_mtd = as.integer(cmin_mtd)[1]
  cmin_recom = as.integer(cmin_recom)[1]
  startup = as.integer(startup)[1]
  alloc_rule = as.integer(alloc_rule)[1]
  early_stop = as.integer(early_stop)[1]
  seed = as.integer(seed)[1]
  nburn = as.integer(nburn)[1]
  niter = as.integer(niter)[1]

  if(startup < 0 || startup > 3){
    stop("Unknown start-up id.")
  }
  if(alloc_rule != 1 && alloc_rule != 2 && alloc_rule != 3){
    stop("Unknown allocation rule id.")
  }
  if(early_stop != 1 && early_stop != 2 && early_stop != 3){
    stop("Unknown early stopping rule id.")
  }
  if(target < 0 || target > 1){
    stop("Targeted toxicity probability is not comprised between 0 and 1.")
  }
  if(target_max < 0 || target_max > 1){
    stop("Maximum targeted toxicity probability is not comprised between 0 and 1.")
  }
  if(target_min < 0 || target_min > 1){
    stop("Minimum targeted toxicity probability is not comprised between 0 and 1.")
  }
  if(n_cohort <= 0){
    stop("Number of cohorts must be positive.")
  }
  if(cohort <= 0){
    stop("Cohort size must be positive.")
  }
  if(time_full < 0){
    stop("Full follow-up time must be positive.")
  }
  if(poisson_rate < 0){
    stop("Parameter for Poisson process accrual must be positive.")
  }
  if(nsim <= 0){
    stop("Number of simulations must be positive.")
  }
  if(c_e < 0 || c_e > 1 || c_d < 0 || c_d > 1 || c_stop < 0 || c_stop > 1 || c_t < 0 || c_t > 1 || c_over < 0 || c_over > 1){
    stop("Probability thresholds are not comprised between 0 and 1.")
  }
  if(cmin_overunder < 0 || cmin_mtd < 0 || cmin_recom < 0){
    stop("Minimum number of cohorts for stopping or recommendation rule must be positive.")
  }
  if(nburn <= 0 || niter <= 0){
    stop("Number of iterations and burn-in for MCMC must be positive.")
  }
  for(a1 in 1:ndose_a1){
    if(prior_tox_a1[a1] < 0 || prior_tox_a1[a1] > 1){
      stop("At least one of the initial guessed toxicity probability for agent 1 is not comprised between 0 and 1.")
    }
  }
  for(a2 in 1:ndose_a2){
    if(prior_tox_a2[a2] < 0 || prior_tox_a2[a2] > 1){
      stop("At least one of the initial guessed toxicity probability for agent 2 is not comprised between 0 and 1.")
    }
  }
  for(a1 in 1:ndose_a1){
    for(a2 in 1:ndose_a2){
      if(p_tox[a1,a2] < 0 || p_tox[a1,a2] > 1){
        stop("At least one of the true toxicity probability is not comprised between 0 and 1.")
      }
    }
  }
  p_tox_na = matrix(NA, nrow=ndose_a1+1, ncol=ndose_a2+1)
  p_tox_na[1:ndose_a1, 1:ndose_a2] = p_tox
  for(a1 in 1:ndose_a1){
    for(a2 in 1:ndose_a2){
      if(p_tox[a1,a2] >
         min(1,p_tox_na[a1+1,a2],p_tox_na[a1,a2+1],p_tox_na[a1+1,a2+1],na.rm=TRUE)){
        stop("The partial ordering between true toxicity probabilities is not satisfied.")
      }
    }
  }

  p_tox = as.double(p_tox)
  inconc = as.double(numeric(1))
  n_pat_dose = as.double(numeric(ndose_a1*ndose_a2))
  rec_dose = as.double(numeric(ndose_a1*ndose_a2))
  n_tox_dose = as.double(numeric(ndose_a1*ndose_a2))
  early_conc = as.double(numeric(1))
  conc_max = as.double(numeric(1))
  tab_pat = as.double(numeric(nsim))

  # Appeler fonction C
  logistic = .C(C_logistic_sim, tite, ndose_a1, ndose_a2, time_full, poisson_rate, p_tox, target,
    target_max, target_min, prior_tox_a1, prior_tox_a2, n_cohort, cohort, nsim, c_e, c_d, c_stop,
    c_t, c_over, cmin_overunder, cmin_mtd, cmin_recom, seed, startup, alloc_rule, early_stop,
    nburn, niter,

    rec_dose=rec_dose, n_pat_dose=n_pat_dose, n_tox_dose=n_tox_dose, inconc=inconc, early_conc=early_conc, tab_pat=tab_pat)

  nsim = logistic$nsim
  inconc=logistic$inconc*100
  early_conc=logistic$early_conc*100
  conc_max=100-early_conc-inconc
  tab_pat=logistic$tab_pat
  rec_dose=logistic$rec_dose*100
  n_pat_dose=logistic$n_pat_dose
  n_tox_dose=logistic$n_tox_dose

  # Reformat outputs
  p_tox= matrix(p_tox,nrow=ndose_a1)
  rec_dose=matrix(rec_dose,nrow=ndose_a1)
  n_pat_dose=matrix(n_pat_dose,nrow=ndose_a1)
  n_tox_dose=matrix(n_tox_dose,nrow=ndose_a1)
  p_tox_p = t(p_tox)[ndose_a2:1,]
  rec_dose_p = t(rec_dose)[ndose_a2:1,]
  n_pat_dose_p = t(n_pat_dose)[ndose_a2:1,]
  n_tox_dose_p = t(n_tox_dose)[ndose_a2:1,]
  dimnames(p_tox_p) = list("Agent 2" = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(rec_dose_p) = list("Agent 2 " = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(n_pat_dose_p) = list("Agent 2"=ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(n_tox_dose_p) = list("Agent 2" = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  pat_tot = round(sum(n_pat_dose),1)

  res = list(call = match.call(),
             tite=tite,
             ndose_a1=ndose_a1,
             ndose_a2=ndose_a2,
             time_full=time_full,
             poisson_rate=poisson_rate,
             startup=startup,
             alloc_rule=alloc_rule,
             early_stop=early_stop,
             p_tox=p_tox,
             p_tox_p=p_tox_p,
             target=target,
             target_min=target_min,
             target_max=target_max,
             prior_tox_a1=prior_tox_a1,
             prior_tox_a2=prior_tox_a2,
             n_cohort=n_cohort,
             cohort=cohort,
             pat_tot=pat_tot,
             nsim=nsim,
             c_e=c_e,
             c_d=c_d,
             c_stop=c_stop,
             c_t=c_t,
             c_over=c_over,
             cmin_overunder=cmin_overunder,
             cmin_mtd=cmin_mtd,
             cmin_recom=cmin_recom,
             nburn=nburn,
             niter=niter,
             seed=seed,
             rec_dose=rec_dose,
             n_pat_dose=n_pat_dose,
             n_tox_dose=n_tox_dose,
             rec_dose_p=rec_dose_p,
             n_pat_dose_p=n_pat_dose_p,
             n_tox_dose_p=n_tox_dose_p,
             inconc=inconc,
             early_conc=early_conc,
             conc_max=conc_max,
             tab_pat=tab_pat)

  class(res) = "CombIncrease_sim"

  return(res)
}



print.CombIncrease_sim = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print_rnd= function (hd, x) {cat(hd, "\n"); print(round(x, digits = dgt)); cat("\n")}
  print_rnd("True toxicities:", x$p_tox_p)
  print_rnd("Percentage of Selection:", x$rec_dose_p)
  print_rnd("Mean number of patients:" , x$n_pat_dose_p)
  print_rnd("Mean number of toxicities:", x$n_tox_dose_p)
  cat(paste("Percentage of inconclusive trials:\t",x$inconc,"\n",sep=""), sep="")
  cat(paste("Percentage of trials stopping with criterion for finding MTD:\t",x$early_conc,"\n",sep=""), sep="")
  cat(paste("Percentage of trials stopping with recommendation based on maximum sample size:\t",x$conc_max,"\n",sep=""), sep="")
  cat("\n")
  cat("Number of simulations:\t", x$nsim, "\n")
  cat("Total mean number of patients accrued:\t", x$pat_tot, "\n")
  cat("Quantiles for number of patients accrued:\t", "\n", quantile(x$tab_pat), "\n")
}




CombIncrease_next = function(ndose_a1, ndose_a2, target, target_min, target_max, prior_tox_a1, prior_tox_a2, cohort, final,
                         pat_incl, dose_adm1, dose_adm2, tite=FALSE, toxicity, time_full=0, time_tox=0, time_follow=0,
                         c_e=0.85, c_d=0.45, c_stop=0.95, c_t=0.5, c_over=0.25, cmin_overunder=2, cmin_mtd=3, cmin_recom=1,
                         early_stop=1, alloc_rule=1, nburn=2000, niter=5000){
  if(tite == TRUE) {
    toxicity = as.numeric(time_tox < time_follow)
  }
  if(pat_incl > 0) {
    cdose1 = dose_adm1[pat_incl]
    cdose2 = dose_adm2[pat_incl]
  }
  else {
    cdose1 = 0
    cdose2 = 0
  }

  n_prior_tox_a1 = length(prior_tox_a1)
  if(n_prior_tox_a1 != ndose_a1){
    stop("The entered vector of initial guessed toxicity probabities for agent 1 is of wrong length.")
  }
  n_prior_tox_a2 = length(prior_tox_a2)
  if(n_prior_tox_a2 != ndose_a2){
    stop("The entered vector of initial guessed toxicity probabities for agent 2 is of wrong length.")
  }

  n_toxicity = length(toxicity)
  n_time_follow = length(time_follow)
  n_time_tox = length(time_tox)
  n_dose_adm1 = length(dose_adm1)
  n_dose_adm2 = length(dose_adm2)
  if(tite==FALSE && n_toxicity != pat_incl){
    stop("The entered vector of observed toxicities is of wrong length.")
  }
  if(tite==TRUE && n_time_follow != pat_incl){
    stop("The entered vector for patients' follow-up time is of wrong length.")
  }
  if(tite==TRUE && n_time_tox != pat_incl){
    stop("The entered vector for patients' time-to-toxicity is of wrong length.")
  }
  if(n_dose_adm1 != pat_incl){
    stop("The entered vector for patients' dose of agent 1 is of wrong length.")
  }
  if(n_dose_adm2 != pat_incl){
    stop("The entered vector for patients' dose of agent 2 is of wrong length.")
  }

  tite = as.logical(tite)
  ndose_a1 = as.integer(ndose_a1)[1]
  ndose_a2 = as.integer(ndose_a2)[1]
  time_full = as.double(time_full)[1]
  target = as.double(target)[1]
  target_max = as.double(target_max)[1]
  target_min = as.double(target_min)[1]
  prior_tox_a1 = as.double(prior_tox_a1)
  prior_tox_a2 = as.double(prior_tox_a2)
  cohort = as.integer(cohort)[1]
  final = as.logical(final)
  c_e = as.double(c_e)[1]
  c_d = as.double(c_d)[1]
  c_stop = as.double(c_stop)[1]
  c_t = as.double(c_t)[1]
  c_over = as.double(c_over)[1]
  cmin_overunder = as.integer(cmin_overunder)[1]
  cmin_mtd = as.integer(cmin_mtd)[1]
  cmin_recom = as.integer(cmin_recom)[1]
  pat_incl = as.integer(pat_incl)[1]
  cdose1 = as.integer(cdose1-1)
  cdose2 = as.integer(cdose2-1)
  dose_adm1 = as.integer(dose_adm1-1)
  dose_adm2 = as.integer(dose_adm2-1)
  time_tox = as.double(time_tox)
  time_follow = as.double(time_follow)
  toxicity = as.logical(toxicity)
  alloc_rule = as.integer(alloc_rule)[1]
  early_stop = as.integer(early_stop)[1]
  nburn = as.integer(nburn)[1]
  niter = as.integer(niter)[1]

  if(alloc_rule != 1 && alloc_rule != 2 && alloc_rule != 3){
    stop("Unknown allocation rule id.")
  }
  if(early_stop != 1 && early_stop != 2 && early_stop != 3){
    stop("Unknown early stopping rule id.")
  }
  if(cohort <= 0){
    stop("Cohort size must be positive.")
  }
  if(time_full < 0){
    stop("Full follow-up time must be positive.")
  }
  if(c_e < 0 || c_e > 1 || c_d < 0 || c_d > 1 || c_stop < 0 || c_stop > 1 || c_t < 0 || c_t > 1 || c_over < 0 || c_over > 1){
    stop("Probability thresholds are not comprised between 0 and 1.")
  }
  if(cmin_overunder < 0 || cmin_mtd < 0 || cmin_recom < 0){
    stop("Minimum number of cohorts for stopping or recommendation rule must be positive.")
  }
  if(nburn <= 0 || niter <= 0){
    stop("Number of iterations and burn-in for MCMC must be positive.")
  }
  for(i in 1:ndose_a1){
    if(prior_tox_a1[i] < 0 || prior_tox_a1[i] > 1){
      stop("At least one of the initial guessed toxicity for agent 1 is not comprised between 0 and 1.")
    }
  }
  for(i in 1:ndose_a2){
    if(prior_tox_a2[i] < 0 || prior_tox_a2[i] > 1){
      stop("At least one of the initial guessed toxicity for agent 2 is not comprised between 0 and 1.")
    }
  }
  if(target < 0 || target > 1){
    stop("Targeted toxicity probability is not comprised between 0 and 1.")
  }
  if(target_max < 0 || target_max > 1){
    stop("Maximum targeted toxicity probability is not comprised between 0 and 1.")
  }
  if(target_min < 0 || target_min > 1){
    stop("Minimum targeted toxicity probability is not comprised between 0 and 1.")
  }


  inconc = as.logical(numeric(1))
  early_conc = as.logical(numeric(1))
  pi = as.double(numeric(ndose_a1*ndose_a2))
  ptox_inf = as.double(numeric(ndose_a1*ndose_a2))
  ptox_inf_targ = as.double(numeric(ndose_a1*ndose_a2))
  ptox_targ = as.double(numeric(ndose_a1*ndose_a2))
  ptox_sup_targ = as.double(numeric(ndose_a1*ndose_a2))


  logistic = .C(C_logistic_next, tite, ndose_a1, ndose_a2, time_full, target, target_max, target_min, prior_tox_a1, prior_tox_a2,
                cohort, final, c_e, c_d, c_stop, c_t, c_over, cmin_overunder, cmin_mtd, cmin_recom, early_stop, alloc_rule,
                nburn, niter, pat_incl, cdose1=cdose1, cdose2=cdose2, dose_adm1, dose_adm2, time_tox, time_follow, toxicity,
                inconc=inconc, early_conc=early_conc, pi=pi, ptox_inf=ptox_inf, ptox_inf_targ=ptox_inf_targ, ptox_targ=ptox_targ, ptox_sup_targ=ptox_sup_targ, NAOK=TRUE)

  # Reformat outputs
  cdose1=logistic$cdose1+1
  cdose2=logistic$cdose2+1
  dose_adm1=dose_adm1+1
  dose_adm2=dose_adm2+1
  pi=matrix(logistic$pi, nrow=ndose_a1)
  ptox_inf=matrix(logistic$ptox_inf, nrow=ndose_a1)
  ptox_inf_targ=matrix(logistic$ptox_inf_targ, nrow=ndose_a1)
  ptox_targ=matrix(logistic$ptox_targ, nrow=ndose_a1)
  ptox_sup_targ=matrix(logistic$ptox_sup_targ, nrow=ndose_a1)

  n_pat_comb = matrix(0, nrow=ndose_a1, ncol=ndose_a2)
  n_tox_comb = matrix(0, nrow=ndose_a1, ncol=ndose_a2)
  for(i in 1:pat_incl){
    n_pat_comb[dose_adm1[i],dose_adm2[i]] = n_pat_comb[dose_adm1[i],dose_adm2[i]]+1
    n_tox_comb[dose_adm1[i],dose_adm2[i]] = n_tox_comb[dose_adm1[i],dose_adm2[i]]+toxicity[i]
  }
  n_pat_comb_p = t(n_pat_comb)[ndose_a2:1,]
  n_tox_comb_p = t(n_tox_comb)[ndose_a2:1,]
  pi_p = t(pi)[ndose_a2:1,]
  ptox_inf_p = t(ptox_inf)[ndose_a2:1,]
  ptox_inf_targ_p = t(ptox_inf_targ)[ndose_a2:1,]
  ptox_targ_p = t(ptox_targ)[ndose_a2:1,]
  ptox_sup_targ_p = t(ptox_sup_targ)[ndose_a2:1,]
  dimnames(n_pat_comb_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(n_tox_comb_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(pi_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(ptox_inf_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(ptox_inf_targ_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(ptox_targ_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(ptox_sup_targ_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)

  res = list(call = match.call(),
             tite=tite,
             ndose_a1=ndose_a1,
             ndose_a2=ndose_a2,
             time_full=time_full,
             target=target,
             target_max=target_max,
             target_min=target_min,
             prior_tox_a1=prior_tox_a1,
             prior_tox_a2=prior_tox_a2,
             cohort=cohort,
             final=final,
             c_e=c_e,
             c_d=c_d,
             c_stop=c_stop,
             c_t=c_t,
             c_over=c_over,
             cmin_overunder=cmin_overunder,
             cmin_mtd=cmin_mtd,
             cmin_recom=cmin_recom,
             early_stop=early_stop,
             alloc_rule=alloc_rule,
             nburn=nburn,
             niter=niter,
             pat_incl=pat_incl,
             cdose1=cdose1,
             cdose2=cdose2,
             dose_adm1=dose_adm1,
             dose_adm2=dose_adm2,
             time_tox=time_tox,
             time_follow=time_follow,
             toxicity=toxicity,
             inconc=logistic$inconc,
             early_conc=logistic$early_conc,
             n_pat_comb=n_pat_comb,
             n_tox_comb=n_tox_comb,
             pi=pi,
             ptox_inf=ptox_inf,
             ptox_inf_targ=ptox_inf_targ,
             ptox_targ=ptox_targ,
             ptox_sup_targ=ptox_sup_targ,
             n_pat_comb_p=n_pat_comb_p,
             n_tox_comb_p=n_tox_comb_p,
             pi_p=pi_p,
             ptox_inf_p=ptox_inf_p,
             ptox_inf_targ_p=ptox_inf_targ_p,
             ptox_targ_p=ptox_targ_p,
             ptox_sup_targ_p=ptox_sup_targ_p)

  class(res) = "CombIncrease_next"

  return(res)
}




print.CombIncrease_next = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print_rnd= function (hd, x) {cat(hd, "\n"); print(round(x, digits = dgt)); cat("\n")}
  print_rnd("Number of patients:" , x$n_pat_comb_p)
  print_rnd("Number of toxicities:", x$n_tox_comb_p)
  print_rnd("Estimated toxicity probabilities:", x$pi_p)
  print_rnd("P(toxicity probability < target):", x$ptox_inf_p)
  print_rnd("Probabilities of underdosing:", x$ptox_inf_targ_p)
  print_rnd("Probabilities being in targeted interval:", x$ptox_targ_p)
  print_rnd("Probabilities of overdosing:", x$ptox_sup_targ_p)
  cat("Warning: recommendation for model-based phase (start-up phase must be ended).\n")
  if(!x$inconc){
    if(!x$early_conc){
      if(x$final){
        cat(paste("The RECOMMENDED COMBINATION at the end of the trial is:\t (",x$cdose1, ",", x$cdose2, ")\n",sep=""), sep="")
      }
      else{
        cat(paste("The next RECOMMENDED COMBINATION is:\t (",x$cdose1, ",", x$cdose2, ")\n",sep=""), sep="")
      }
    }
    else{
      cat(paste("The dose-finding process should be STOPPED (criterion for identifying MTD met) and the RECOMMENDED COMBINATION is:\t (",x$cdose1, ",", x$cdose2, ")\n",sep=""), sep="")
    }
  }
  else{
    cat(paste("The dose-finding process should be STOPPED WITHOUT COMBINATION RECOMMENDATION (criterion for over or under toxicity met)\n",sep=""), sep="")
  }
}
