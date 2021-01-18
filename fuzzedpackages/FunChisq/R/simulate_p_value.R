# simulate_p_value.R
#
# Simulate p value for FunChisq
# Created by: Hua Zhong, Ruby Sharma and Dr. Joe Song
# Date : May 24 2017

# x: observed contingency table
simulate.p.value <- function(x, simulate.nruns)
{
  if(is.matrix(x) || is.data.frame(x)){
    sample.size <- sum(x)
    nr <- nrow(x)
    nc <- ncol(x)
    #simulate.nruns <- 500 # number of simulated tables
    row.marginal <- rowSums(x) / sample.size
    col.marginal <- colSums(x) / sample.size
    if(sum(row.marginal) != 1) row.marginal[nr] <- 1 - sum(row.marginal[1:(nr-1)])
    if(sum(col.marginal) != 1) col.marginal[nc] <- 1 - sum(col.marginal[1:(nc-1)])

    #noise.multinomial <- sample(c(0,0.1,0.2,0.3,0.4,0.5), simulate.nruns, replace = TRUE)
    noise.multinomial <- runif(simulate.nruns, min = 0, max = 1)
    # tables <- lapply(c(1:simulate.nruns), function(y){
    #   set.seed((nr+nc+sample.size)*y)
    #   sim.table <- simulate_tables(n = sample.size, type = "independent",
    #                          nrow = nr, ncol = nc, noise = noise.multinomial[y],
    #                          row.marginal = row.marginal,
    #                          col.marginal = col.marginal)$noise.list[[1]]
    #   return(sim.table)
    # })
    tables <- simulate.null.tables.with.marginal.prob(sample.size, simulate.nruns,
                                                      row.marginal, col.marginal,
                                                      noise = noise.multinomial)
    #tables <- simulate.null.tables.with.cell.prob(x, simulate.nruns, noise = noise.multinomial)

    tables.FunChisq <- as.numeric(unlist(sapply(c(1:simulate.nruns), function(y){
      fun.chisq.test(tables[[y]])$statistic
    }), use.names = FALSE))

    observed.FunChisq <- fun.chisq.test(x)$statistic
    p.value <- (1 + sum(tables.FunChisq >= observed.FunChisq)) / (1 + simulate.nruns)

    return(p.value)
  }else{
    stop("Wrong input contingency table \"x\" in simulate.p.value.R!")
  }
}


# simulate.null.tables.with.cell.prob <- function(x, simulate.nruns, noise){
#   sample.size <- sum(x)
#   nr <- nrow(x)
#   nc <- ncol(x)
#   noise.null <- noise
#   if(length(noise.null) == 1){
#     noise.null <- rep(noise.null, simulate.nruns)
#   }else if(length(noise.null) == simulate.nruns){
#
#   }else{
#     stop("Wrong length of noise values in \"simulate.null.tables\"!")
#   }
#   prob.matrix <- x / sum(x)
#   prob.matrix.vector <- as.vector(prob.matrix)
#   output.tables <- lapply(c(1:simulate.nruns), function(m){
#     sample.matrix.vector <- sample(c(1:length(prob.matrix.vector)), size = sample.size, replace = TRUE, prob = prob.matrix.vector)
#     sample.matrix.vector <- as.numeric(table(factor(sample.matrix.vector, levels = c(1:length(prob.matrix.vector)))))
#     sample.matrix <- matrix(sample.matrix.vector, nrow = nr, ncol = nc)
#     sample.matrix.noised <- add.noise(tables = sample.matrix, noise.model = "house", margin = 0, u = noise.null[m])
#     return(sample.matrix.noised)
#   })
#   return(output.tables)
# }

simulate.null.tables.with.marginal.prob <-
  function(sample.size, simulate.nruns,
           row.marginal, col.marginal, noise)
{
  noise.null <- noise
  if(length(noise.null) == 1){
    noise.null <- rep(noise.null, simulate.nruns)
  }else if(length(noise.null) == simulate.nruns){

  }else{
    stop("Wrong length of noise values in \"simulate.null.tables\"!")
  }
  prob.matrix <- row.marginal %o% col.marginal
  prob.matrix.vector <- as.vector(prob.matrix)
  output.tables <- lapply(c(1:simulate.nruns), function(m){
    sample.matrix.vector <- sample(c(1:length(prob.matrix.vector)),
                                   size = sample.size,
                                   replace = TRUE,
                                   prob = prob.matrix.vector)
    sample.matrix.vector <- as.numeric(table(factor(sample.matrix.vector,
                                                    levels = c(1:length(prob.matrix.vector)))))
    sample.matrix <- matrix(sample.matrix.vector,
                            nrow = length(row.marginal),
                            ncol = length(col.marginal))
    sample.matrix.noised <- add.noise(tables = sample.matrix,
                                      noise.model = "house",
                                      margin = 0,
                                      u = noise.null[m])
    return(sample.matrix.noised)
  })
  return(output.tables)
}
