#' @title Plots \code{SimBIID_runs} objects
#'
#' @description Plot method for \code{SimBIID_runs} objects.
#'
#' @param x     An \code{SimBIID_runs} object.
#' @param which A character vector of states to plot. Can be \code{"all"} to plot all
#'              states (and final event times), or \code{"t"} to plot final event times.
#' @param type Character stating whether to plot full simulations over time (\code{"runs"}) or
#'             summaries (\code{"sums"}).
#' @param rep An integer vector of simulation runs to plot.
#' @param quant A vector of quantiles (> 0.5) to plot if \code{type == "runs"}.
#' @param data A \code{data.frame} containing time series count data, 
#'          with the first column called \code{t}, followed by columns of time-series counts.
#' @param matchData A character vector containing matches between the columns of \code{data} and
#'                  the columns of the model runs. Each entry must be of the form e.g. \code{"SD = SR"},
#'                  where \code{SD} is the name of the column in \code{data}, and \code{SR} is the name
#'                  of the column in \code{x}.
#' @param ... Not used here.
#'
#' @return A plot of individual simulations and/or summaries of repeated simulations 
#'         extracted from \code{SimBIID_runs} object.
#'  
#' @method plot SimBIID_runs 
#' 
#' @seealso \code{\link{mparseRcpp}}, \code{\link{print.SimBIID_runs}}, \code{\link{run}} 
#'      
#' @export
#' 
#' @examples 
#' \donttest{
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
#'     pars = pars,
#'     tspan = TRUE
#' )
#' 
#' ## run 100 replicate simulations and
#' ## plot outputs
#' sims <- run(
#'     model = model,
#'     pars = c(beta = 0.001, gamma = 0.1),
#'     tstart = 0,
#'     tstop = 100,
#'     u = c(S = 119, I = 1, R = 0),
#'     tspan = seq(1, 100, length.out = 10),
#'     nrep = 100
#' )
#' plot(sims, quant = c(0.55, 0.75, 0.9))
#' 
#' ## add replicate 1 to plot
#' plot(sims, quant = c(0.55, 0.75, 0.9), rep = 1)
#' }
#' 

plot.SimBIID_runs <- function(x, which = c("all", "t"), type = c("runs", "sums"), 
                              rep = NA, quant = 0.9, data = NULL, matchData = NULL, ...) {
    ## check x
    if(class(x) != "SimBIID_runs"){
        stop("'x' is not a SimBIID_runs object")
    }
    ## check type
    checkInput(type, c("vector", "character"))
    type <- type[1]
    checkInput(type, inSet = c("runs", "sums"))
    ## check which
    checkInput(which, c("vector", "character"))
    if(which[1] == "t" & type == "runs") {
        stop("Can't plot 't' variable with 'type == \"runs\"'")
    }
    whichset <- c("t", colnames(x$sums)[-c(1:3)])
    if(which[1] != "all") {
        checkInput(which, c("vector", "character"), inSet = whichset)
    } else {
        which <- whichset
    }
    ## check rep
    if(!is.na(rep[1])){
        checkInput(rep, c("vector", "numeric"), inSet = x$sums$rep, int = TRUE)
    }
    ## check quant
    checkInput(quant, c("vector", "numeric"))
    checkInput(format(quant, drop0trailing = TRUE), inSet = format(seq(0.55, 0.95, by = 0.05), drop0trailing = TRUE))
    quant <- sort(quant)
    quant <- cbind(rev(1 - quant), rev(quant))
    quant1 <- sort(as.vector(quant))
    
    ## check data
    if(!is.null(data[1])){
        if(is.null(matchData[1])){
            stop("'matchData' must be provided if 'data' is provided")
        }
        checkInput(data, "data.frame")
        if(colnames(data)[1] != "t"){
            stop("First column of 'data' must be 't'")
        }
        if(ncol(data) < 2) {
            stop("Must have at least one count column in 'data'")
        }
        checkInput(data$t, "numeric")
        for(j in 2:ncol(data)) {
            checkInput(data[, j, drop = TRUE], "numeric", int = TRUE)
        }
        ## check matchData
        matchData <- strsplit(matchData, "=")
        if(!all(sapply(matchData, length) == 2)){
            stop("'matchData' not correctly specified")
        }
        ## trim whitespace
        matchData <- lapply(matchData, trimws)
        ## extract colnames of data
        datNames <- sapply(matchData, function(x){
            x[[1]]
        })
        if(!all(datNames %in% colnames(data))){
            stop("Can't match data names in 'matchData' to 'data'")
        }
        ## extract colnames of sims
        simNames <- sapply(matchData, function(x){
            x[[2]]
        })
        if(!all(simNames %in% colnames(x$sums))){
            stop("Can't match simulated data names in 'matchData' to 'x$runs'")
        }
        ## match to simulations
        data <- data %>%
            arrange(t) 
        data1 <- dplyr::select(data, t)
        for(j in 1:length(datNames)) {
            data1 <- cbind(data1, data[, match(datNames[j], colnames(data))]) %>%
                set_names(c("t", simNames[1:j])) 
        }
        data <- data1 %>%
            gather(output, value, -t) %>%
            mutate(output = factor(output, levels = which))
    } else {
        if(!is.null(matchData[1])){
            stop("'data' must be provided if 'matchData' is provided")
        }
    }
    
    ## plot final times and final epidemic sizes
    if(!is.data.frame(x$runs) | type == "sums"){
        if(nrow(x$sums) < 20){
            stop("Can't plot summary of final sizes and times for n < 20 replicates")
        }
        p <- x$sums %>%
            dplyr::select(!!which) %>%
            gather(output, value) %>%
            mutate(output = factor(output, levels = which)) %>%
            ggplot(aes(x = value)) +
                geom_histogram() +
                facet_wrap(~ output, scales = "free") +
                xlab("")
        if(!is.na(rep[1])){
            rep1 <- rep
            repSums <- x$sums %>%
                slice(rep1) %>%
                dplyr::select(!!which) %>%
                gather(output, value) %>%
                mutate(output = factor(output, levels = which))
            p <- p + geom_point(aes(x = value), data = repSums, y = 0, colour = "red", shape = 16)
        }
        ## add observed data if specified
        if(!is.null(data[1])){
            if(nrow(data) > 1){
                stop("Must have only 1 data point to plot if 'type = \"sums\"")
            }
            if(any(simNames %in% which)) {
                p <- p + geom_point(aes(x = value), data = data, y = 0, shape = 16)
            }
        }   
    } else {
        ## produce plot
        which <- unique(c("rep", "t", which))
        if(nrow(x$sums) < 20) {
            ## produce plot
            p <- x$runs %>%
                gather(output, value, -rep, -t) %>%
                mutate(output = factor(output, levels = which[-match(c("rep", "t"), which)])) %>%
                ggplot(aes(x = t, y = value)) +
                xlab("Time") + ylab("Counts") +
                facet_wrap(~ output)
            ## add individual simulations if required
            rep1 <- unique(x$runs$rep)
            repSums <- x$runs %>%
                dplyr::select(!!which) %>%
                gather(output, value, -rep, -t) %>%
                mutate(output = factor(output, levels = which[-match(c("rep", "t"), which)]))
            for(i in rep1) {
                temp <- dplyr::filter(repSums, rep == i)
                p <- p + geom_line(aes(x = t, y = value), data = temp)
            }
            ## add observed data if specified
            if(!is.null(data[1])){
                if(any(simNames %in% which)) {
                    p <- p + geom_line(aes(x = t, y = value), data = data, linetype = "dashed")
                }
            }  
        } else {
            p <- x$runs %>%
                dplyr::select(!!which) %>%
                gather(output, value, -rep, -t) %>%
                group_by(output, t) %>%
                summarise(med = stats::median(value), value = list(enframe(stats::quantile(value, probs = quant1)))) %>%
                unnest() 
            p1 <- quant %>%
                as.data.frame() %>%
                rename(lci = V1, uci = V2) %>%
                mutate(pair = 1:n()) %>%
                gather(output, name, -pair) %>%
                mutate(name = 100 * name) %>%
                mutate(name = paste0(name, "%")) %>%
                inner_join(p, by = "name") %>%
                dplyr::select(-name) %>%
                spread(output.x, value) %>%
                mutate(output.y = factor(output.y, levels = which[-match(c("rep", "t"), which)])) %>%
                mutate(pair = as.character(quant[pair, 2]))
           
           ## produce plot
           p <- p1 %>%
               ggplot(aes(x = t)) +
                   xlab("Time") + ylab("Counts")
           for(i in unique(p1$pair)){
                temp <- dplyr::filter(p1, pair == i)
                p <- p + geom_ribbon(aes(ymin = lci, ymax = uci), 
                     data = temp, alpha = 0.2)
           }
           
           ## add individual simulations if required
           if(!is.na(rep[1])){
               rep1 <- rep
               repSums <- x$runs %>%
                   dplyr::filter(rep %in% rep1) %>%
                   dplyr::select(!!which) %>%
                   gather(output, value, -rep, -t) %>%
                   mutate(output.y = factor(output, levels = levels(p1$output.y)))
               for(i in rep1){
                   temp <- dplyr::filter(repSums, rep == i)
                   p <- p + geom_line(aes(x = t, y = value), data = temp)
               }
           }
           p <- p + geom_line(aes(y = med), colour = "red") +
               facet_wrap(~ output.y) +
               labs(title = paste0("Intervals = ", paste0(paste0(rev(quant[, 2]) * 100, "%"), collapse = ", ")))
           
           ## add observed data if specified
           if(!is.null(data[1])){
               data <- data %>%
                   rename(output.y = output, med = value)
               if(any(simNames %in% which)) {
                    p <- p + geom_line(aes(x = t, y = med), data = data, linetype = "dashed")
               }
           }  
        }
        ## add bootstrap end indicator if required
        if(!is.na(x$bootEnd)){
            p <- p + geom_vline(xintercept = x$bootEnd, linetype = "dashed", colour = "blue")
        }
    }
    p
}
    
