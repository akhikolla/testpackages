#' @title Plots \code{ABCSMC} objects
#'
#' @description Plot method for \code{ABCSMC} objects.
#'
#' @param x     An \code{ABCSMC} object.
#' @param type  Takes the value \code{"post"} if you want to plot posterior distributions.
#'              Takes the value \code{"output"} if you want to plot the simulated outputs.
#' @param gen   A vector of generations to plot. If left missing then defaults to all generations.
#' @param joint A logical describing whether joint or marginal distributions are wanted.
#' @param transfunc Is a \code{function} object where the arguments to the function must
#'                  match all or a subset of the parameters in the model. This function needs 
#'                  to return a \code{data.frame} object with columns containing the transformed
#'                  parameters.
#' @param ... Not used here.
#'
#' @return A plot of the ABC posterior distributions for different generations, or the distributions
#'         of the simulated summary measures for different generations.
#' 
#' @method plot ABCSMC
#'         
#' @export
#' 
#' @seealso \code{\link{ABCSMC}}, \code{\link{print.ABCSMC}}, \code{\link{summary.ABCSMC}}
#'     
#' @examples 
#' \donttest{
#' 
#' ## set up SIR simulation model
#' transitions <- c(
#'     "S -> beta * S * I -> I", 
#'     "I -> gamma * I -> R"
#' )
#' compartments <- c("S", "I", "R")
#' pars <- c("beta", "gamma")
#' model <- mparseRcpp(
#'     transitions = transitions, 
#'     compartments = compartments,
#'     pars = pars
#' )
#' model <- compileRcpp(model)
#' 
#' ## generate function to run simulators
#' ## and return summary statistics
#' simSIR <- function(pars, data, tols, u, model) {
#' 
#'     ## run model
#'     sims <- model(pars, 0, data[2] + tols[2], u)
#'     
#'     ## this returns a vector of the form:
#'     ## completed (1/0), t, S, I, R (here)
#'     if(sims[1] == 0) {
#'         ## if simulation rejected
#'         return(NA)
#'     } else {
#'         ## extract finaltime and finalsize
#'         finaltime <- sims[2]
#'         finalsize <- sims[5]
#'     }
#'     
#'     ## return vector if match, else return NA
#'     if(all(abs(c(finalsize, finaltime) - data) <= tols)){
#'         return(c(finalsize, finaltime))
#'     } else {
#'         return(NA)
#'     }
#' }
#' 
#' ## set priors
#' priors <- data.frame(
#'     parnames = c("beta", "gamma"), 
#'     dist = rep("gamma", 2), 
#'     stringsAsFactors = FALSE
#' )
#' priors$p1 <- c(10, 10)
#' priors$p2 <- c(10^4, 10^2)
#' 
#' ## define the targeted summary statistics
#' data <- c(
#'     finalsize = 30, 
#'     finaltime = 76
#' )
#' 
#' ## set initial states (1 initial infection 
#' ## in population of 120)
#' iniStates <- c(S = 119, I = 1, R = 0)
#' 
#' ## set initial tolerances
#' tols <- c(
#'     finalsize = 50,
#'     finaltime = 50
#' )
#' 
#' ## run 2 generations of ABC-SMC
#' ## setting tolerance to be 50th
#' ## percentile of the accepted 
#' ## tolerances at each generation
#' post <- ABCSMC(
#'     x = data, 
#'     priors = priors, 
#'     func = simSIR, 
#'     u = iniStates, 
#'     tols = tols, 
#'     ptol = 0.2, 
#'     ngen = 2, 
#'     npart = 50,
#'     model = model
#' )
#' post
#' 
#' ## run one further generation
#' post <- ABCSMC(post, ptols = 0.5, ngen = 1)
#' post
#' summary(post)
#' 
#' ## plot posteriors
#' plot(post)
#' 
#' ## plot outputs
#' plot(post, "output")
#' }
#' 

plot.ABCSMC <- function(x, type = c("post", "output"), gen = NA, joint = FALSE, transfunc = NA, ...) {
    
    ## check x
    if(class(x) != "ABCSMC"){
        stop("'x' is not a ABCSMC object")
    }
    
    ## check type
    type <- type[1]
    checkInput(type, "character", inSet = c("post", "output"))
    
    ## check gens
    if(is.na(gen[1])) {
        gen <- 1:length(x$pars)
    }
    checkInput(gen, c("vector", "numeric"), int = TRUE, inSet = 1:length(x$pars))
    gen <- as.character(sort(gen))
    
    ## check joint
    checkInput(joint, c("vector", "logical"), 1)
    if(joint) {
        if(!requireNamespace("GGally", quietly = TRUE)) {
            stop("'GGally' package required for joint distribution plots")
        }
    }
    
    if(type == "post") {
        ## generate colorRamp
        getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))
        fillCols <- getPalette(length(x$pars))
        fillCols <- fillCols[as.numeric(gen)]
        
        ## get parameter names
        pnames <- x$priors$parnames
        
        p <- x$pars %>%
                map(~{
                    as_tibble(.) %>%
                    set_names(pnames)
                }, pnames = pnames) %>%
                bind_rows(.id = "Generation") %>%
                dplyr::filter(Generation %in% gen) %>%
                mutate(Generation = factor(Generation, levels = gen))
                
        ## check for transformations if required
        if(length(transfunc) != 1){
            stop("'transfunc' not of length 1")
        }
        if(is.function(transfunc)) {
        
            ## check function arguments
            fargs <- formals(transfunc)
            checkInput(names(fargs), inSet = colnames(p))
            
            ## perform transformations if required
            temppars <- p[, match(names(fargs), colnames(p))]
            temppars <- as.list(temppars)
            names(temppars) <- names(fargs)
            temp <- do.call("transfunc", temppars)
            checkInput(temp, "data.frame", nrow = nrow(p))
            
            ## bind to current posterior samples
            p <- cbind(p, temp)
        } else {
            if(!is.na(transfunc)) {
                stop("'transfunc' incorrectly specified")
            }
        }
        
        ## plot posteriors
        if(!joint) {
            p <- p %>%
                gather(Parameter, value, -Generation) %>%
                ggplot(aes(x = value, fill = Generation)) +
                    geom_density(alpha = 0.8) +
                    xlab("Parameter value") + ylab("Density") +
                    scale_fill_manual(values = fillCols) +
                    facet_wrap(~ Parameter, scales = "free")
         } else {
            p <- GGally::ggpairs(p, mapping = aes(colour = Generation), 
                columns = 2:ncol(p),
                diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.8)),
                lower = list(continuous = "density"),
                upper = list(continuous = "blank"),
                legend = c(1, 1))
            ## amend colours
            for(i in 1:p$nrow) {
                for(j in 1:p$ncol){
                    p[i, j] <- p[i, j] + scale_fill_manual(values = fillCols)
                    p[i, j] <- p[i, j] + scale_colour_manual(values = fillCols)
                }
            }
         }
     } else {
        ## generate colorRamp
        getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))
        fillCols <- getPalette(length(x$pars))
        fillCols <- fillCols[as.numeric(gen)]
        
        ## get output names
        onames <- colnames(x$output[[1]])
        
        ## plot outputs
        p <- x$output %>%
            map(~{
                as_tibble(.) %>%
                set_names(onames)
            }, onames = onames) %>%
            bind_rows(.id = "Generation") %>%
            dplyr::filter(Generation %in% gen) %>%
            mutate(Generation = factor(Generation, levels = gen))
            
        if(length(transfunc) != 1){
            stop("'transfunc' not of length 1")
        }
        if(is.function(transfunc)) {
            stop("'transfunc' can't be used for output plots")
        } else {
            if(!is.na(transfunc)) {
                stop("'transfunc' incorrectly specified")
            }
        }
        
        ## extract data
        dat <- x$data %>% as.list() %>% data.frame() %>%
            gather(Output, value)
            
        if(!joint) {
            p <- p %>%
                gather(Output, value, -Generation)
            ## add noise if all values identical, else plot
            ## doesn't render properly 
            ident <- p %>% 
                group_by(Generation, Output) %>%
                summarise(nident = length(unique(value))) %>%
                dplyr::filter(nident == 1)
            if(nrow(ident) > 0){
                p <- full_join(p, ident, by = c("Generation", "Output"))
            } else {
                p <- mutate(p, nident = NA)
            }
            p <- dplyr::filter(p, is.na(nident)) %>% 
                ggplot(aes(x = value, fill = Generation)) +
                    geom_density(alpha = 0.8) +
                    geom_density(data = dplyr::filter(p, !is.na(nident)), alpha = 0.8, adjust = 0.0001) +
                    xlab("Output") + ylab("Density") +
                    scale_fill_manual(values = fillCols) +
                    facet_wrap(~ Output, scales = "free") +
                    geom_vline(aes(xintercept = value), data = dat, linetype = 2, colour = "red", size = 1.5)
        } else {
            p <- GGally::ggpairs(p, mapping = aes(colour = Generation), 
                columns = 2:ncol(p),
                diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.8)),
                lower = list(continuous = "density"),
                upper = list(continuous = "blank"),
                legend = c(1, 1))
            ## amend colours
            for(i in 1:p$nrow) {
                for(j in 1:p$ncol){
                    p[i, j] <- p[i, j] + scale_fill_manual(values = fillCols)
                    p[i, j] <- p[i, j] + scale_colour_manual(values = fillCols)
                }
            }
        }
            
     }
     p
}
    
