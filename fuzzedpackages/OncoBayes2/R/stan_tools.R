## that the first dimension of each list entry is the iteration.
extract_draw <- function(sims, draw) lapply(sims, asub, idx=draw, dim=1, drop=FALSE)

## applies a function over each entry of the posterior if
## vectorized=FALSE; for vectorized=TRUE the function is assumed to
## perform the simulation in a single sweep. Note that all arguments
## to the function are automatically deduced from it's formals and
## that all arguments which are not in the sims list are searched in
## the global environment.
posterior_simulate <- function(sims, fun, vectorized=FALSE, res_type, envir) {
    args <- setdiff(names(formals(fun)), "seed")

    from_draw <- intersect(args, names(sims))
    from_env  <- setdiff(args, names(sims))

    if (missing(envir))
        envir <- parent.frame()

    sims <- sims[from_draw]
    aux <- mget(from_env, envir=envir)

    if(!vectorized) {
        S <- NROW(sims[[1]])
        calc_draw <- function(i) do.call(fun, c(aux, extract_draw(sims, i)))
        if(missing(res_type))
            res_type <- calc_draw(1)
        res <- vapply(1:S, calc_draw, res_type)
        nd <- length(dim(res))
        if(nd == 0) {
            return(array(res, dim=c(length(res),1)))
        }
        if(2 > (nd - 1))
            aux_ind <- c()
        else
            aux_ind <- 2:(nd-1)
        return(aperm(res, c( nd, aux_ind, 1)))
    } else {
        return(do.call(fun, c(sims, aux)))
    }
}

stan_generate_model <- function(file, ..., postfix="_generated") {
    generated_file <- gsub(".stan$", paste0(postfix, ".stan"), file)
    stan_model <- stanc_builder(file, ...)
    cat(stan_model$model_code, file=generated_file)
    cat("\n", file=generated_file, append=TRUE)
    generated_file
}

