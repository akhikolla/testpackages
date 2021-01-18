#' Internal function to simulate from the posterior new parameter draws
#'
#' @keywords internal
pp_data <- function(object, newdata, draws, re.form) {
    data <- object$data
    if(!missing(newdata))
        data <- newdata
    if(!missing(re.form))
        stop("ERROR: re.form not yet supported")

    strata_group_fct <- .get_strata_group_fct(object, data)

    f <- object$formula
    orig_mf <- object$model
    idx_group_term <- object$idx_group_term
    idx_inter_term <- object$idx_inter_term
    has_inter <- object$has_inter
    tt <- terms(f, data=orig_mf, lhs=1, rhs=1:idx_group_term)
    Terms <- delete.response(tt)
    mf <- model.frame(Terms, data)

    num_comp <- object$standata$num_comp

    ## setup design matrices which must have intercept and slope
    X_comp <- list()
    for (i in 1:num_comp) {
        X_comp[[i]] <- model.matrix(f, mf, rhs=i)
        assert_matrix(X_comp[[i]], ncols=2, any.missing=FALSE)
    }
    X_comp <- do.call(abind, c(X_comp, list(along=0)))

    if(has_inter)
        X_inter <- model.matrix(f, mf, rhs=idx_inter_term)
    else
        X_inter <- model.matrix(~0, mf)

    num_obs <- dim(X_comp)[2]


    group_idx  <- as.integer(unclass(strata_group_fct$group_fct))
    stratum_idx  <- as.integer(unclass(strata_group_fct$strata_fct))

    pred_data <- list(X_comp=X_comp,
                      X_inter=X_inter,
                      group=group_idx,
                      stratum=stratum_idx)

    post <- rstan::extract(object$stanfit, c("beta_group", "eta_group"))
    names(post) <- sub("_group", "", names(post))

    num_post_draws <- dim(post$beta)[1]
    if(!has_inter) {
        num_groups_fitted <- dim(post$beta)[2]
        post$eta <- array(NA, dim=c(num_post_draws, num_groups_fitted, 0))
        post$eta_map <- array(NA, dim=c(num_post_draws, num_groups_fitted, 0))
    }

    if(!missing(draws))
        post <- extract_draw(post, 1:draws)

    ##pp <- posterior_simulate(post, blrm_logit_grouped, envir=list2env(pred_data))
    pp <- posterior_simulate(post, blrm_logit_grouped_vec, vectorized=TRUE, envir=list2env(pred_data))
    colnames(pp) <- rownames(data)
    pp
}

.validate_factor <- function(test, expected, name) {
    expected_levels  <- levels(expected)
    if(is.factor(test)) {
        test_levels  <- levels(test)
        unseen_levels <- setdiff(test_levels, expected_levels)
        assert_that(length(unseen_levels) == 0,
                    msg=paste0("Found unkown factor levels in ", name, ": ", paste(unseen_levels, collapse=", ")))
        assert_that(all(expected_levels == test_levels), msg=paste0("Mismatch in factor level defintion of ", name))
        return(test)
    }

    unseen_levels <- setdiff(unique(test), expected_levels)
    assert_that(length(unseen_levels) == 0,
                msg=paste0("Found unkown factor levels in ", name, ": ", paste(unseen_levels, collapse=", ")))
    factor(test, levels=expected_levels)
}



#' Numerically stable summation of log inv logit
#' @keywords internal
log_inv_logit <- function(mat) {
    - ifelse(is.finite(mat) & (mat < 0), log1p(exp(mat)) - mat, log1p(exp(-mat)))
}



## blrm_logit_grouped <- function(group, stratum, X_comp, X_inter, beta, eta) {
##     num_comp <- dim(X_comp)[1]
##     num_inter <- dim(X_inter)[2]
##     num_obs <- length(group)
##     mu <- rep(NA, num_obs)
##     for(i in 1:num_obs) {
##         g <- group[i]
##         s <- stratum[i]
##         log_p0_nr <- 0
##         for(j in 1:num_comp) {
##             ##log_p0_nr <- log_p0_nr + log(1 - inv_logit(X_comp[j,i,] %*% beta[g,j,]))
##             log_p0_nr <- log_p0_nr + log_inv_logit(-1 * adrop(X_comp[j,i,,drop=FALSE], drop=1:2) %*% beta[g,j,,drop=TRUE])
##         }
        ##mu[i] <- log(1 - exp(log_p0_nr)) - log_p0_nr
##         mu[i] <- log1p(-exp(log_p0_nr)) - log_p0_nr
##         if(num_inter > 0) {
##             mu[i] <- mu[i] + adrop(X_inter[i,,drop=FALSE], drop=1) %*% eta[g,,drop=TRUE]
##         }
##     }
##     mu
## }

## numerically stable version of log1m_exp(x) = log(1-exp(x)) for x < 0
log1m_exp_max0  <- function(x) {
    ## qlogis(x) = logit(exp(x)) = x - log(1-exp(x))
    x - qlogis(x, log.p=TRUE)
}

## vectorized version
blrm_logit_grouped_vec <- function(group, stratum, X_comp, X_inter, beta, eta) {
    num_comp <- dim(X_comp)[1]
    num_inter <- dim(X_inter)[2]
    num_obs <- length(group)
    S <- dim(beta)[1]
    mu <- matrix(NA, S, num_obs)
    for(i in 1:num_obs) {
        g <- group[i]
        s <- stratum[i]
        log_p0_nr <- 0
        for(j in 1:num_comp) {
            ##log_p0_nr <- log_p0_nr + log(1 - inv_logit(X_comp[j,i,] %*% beta[g,j,]))
            log_p0_nr <- log_p0_nr + log_inv_logit(-1 * adrop(X_comp[j,i,,drop=FALSE], drop=1) %*% t(matrix(beta[,g,j,,drop=FALSE], S, 2)))
        }
        ##mu[i] <- log(1 - exp(log_p0_nr)) - log_p0_nr
        ##mu[,i] <- log1p(-exp(log_p0_nr)) - log_p0_nr
        mu[,i] <- log1m_exp_max0(log_p0_nr) - log_p0_nr
        if(num_inter > 0) {
            mu[,i] <- mu[,i] + X_inter[i,,drop=FALSE] %*% t(matrix(eta[,g,,drop=FALSE], S, num_inter))
        }
    }
    mu
}



pp_binomial_trials <- function(object, newdata) {
    data <- object$data
    if(!missing(newdata))
        data <- newdata

    f <- object$formula
    orig_mf <- object$model
    idx_group_term <- object$idx_group_term
    tt <- terms(f, data=orig_mf, lhs=1, rhs=1:idx_group_term)
    mf <- model.frame(tt, data)
    y <- model.response(mf)
    return(rowSums(y))
}


#' extracts from a blrmfit object and a given data-set the group and
#' stratum factor
#'
#' @keywords internal
.get_strata_group_fct  <- function(object, data) {
    f <- object$formula
    orig_mf <- object$model
    idx_group_term <- object$idx_group_term
    tt <- terms(f, data=orig_mf, lhs=1, rhs=1:idx_group_term)
    Terms <- delete.response(tt)
    mf <- model.frame(Terms, data)

    group_index_term <- model.part(f, data = mf, rhs = idx_group_term)
    if(ncol(group_index_term) == 2) {
        idx_group_index <- 2
        idx_strata_index <- 1
    } else {
        idx_group_index <- 1
        idx_strata_index <- NA
    }

    model_group_fct <- object$group_fct
    group_fct <- model.part(f, data = mf, rhs = idx_group_term)
    group_fct <- group_fct[,idx_group_index]
    group_fct <- .validate_factor(group_fct, model_group_fct, "grouping")

    ## enforce that all strata labels are known and match what has
    ## been defined at model fit
    if(!is.na(idx_strata_index)) {
        model_strata_fct <- object$strata_fct
        strata_fct <- model.part(f, data = mf, rhs = idx_group_term)[,idx_strata_index]
        strata_fct <- .validate_factor(strata_fct, model_strata_fct, "stratum")
        strata_group <- data.frame(strata_fct=strata_fct, group_fct=group_fct)
    } else {
        strata_group <- data.frame(strata_fct=1, group_fct=group_fct)
    }

    strata_group
}
