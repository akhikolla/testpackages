#' Simulation function for strategy estimation.
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param data a \code{stratEst.data} object. Alternatively, the arguments \code{num.ids}, \code{num.games}, and \code{num.periods} can be used if no data is available.
#' @param strategies a list of strategies. Each element if the list must be an object of class \code{stratEst.strategy}.
#' @param shares a numeric vector of strategy shares. The order of the elements corresponds to the order in \code{strategies}. \code{NA} values are not allowed. Use a list of numeric vectors if data has more than one sample and shares are sample specific.
#' @param coefficients a matrix of regression coefficients. Column names correspond to the names of the strategies, row names to the names of the covariates.
#' @param covariate.mat a matrix with the covariates in columns. The column names of the matrix indicate the names of the covariates. The matrix must have as many rows as there are individuals.
#' @param num.ids an integer. The number of individuals. Default is 100.
#' @param num.games an integer. The number of games. Default is 5.
#' @param num.periods a vector of integers with as many elements \code{num.games}. The elements specify the number of periods in each game. Default (\code{NULL}) means 5 periods in each game.
#' @param fixed.assignment a logical value. If \code{FALSE} individuals use potentially different strategies in each each game. If \code{TRUE}, individuals use the same strategy in each game. Default is \code{FALSE}.
#' @param input.na a logical value. If \code{FALSE} an input value is randomly selected for the first period. Default is \code{FALSE}.
#' @param sample.id a character string indicating the name of the variable which identifies the samples in data. Individual observations must be nested in samples. Default is \code{NULL}.
#' @return A \code{stratEst.data} object. A data frame in the long format with the following variables:
#' \item{id}{the variable that identifies observations of the same individual.}
#' \item{game}{the variable that identifies observations of the same game.}
#' \item{period}{the period of the game.}
#' \item{choice}{the discrete choices.}
#' \item{input}{the inputs.}
#' \item{sample}{the sample of the individual.}
#' \item{strategy}{the strategy of the individual.}
#' @description The simulation function of the package.
#' @examples
#' ## Simulate data of two strategies for choices "left" and "right".
#' lr <- c("left","right")
#' pi <- runif(1)
#' pr <- c(1,0,0,1)
#' tr <- c(1,2,1,2)
#' mixed <- stratEst.strategy(choices = lr, inputs = lr, prob.choices = c(pi, 1 - pi))
#' pure <- stratEst.strategy(choices = lr, inputs = lr, prob.choices = pr, tr.inputs = tr)
#' gamma <- runif(1)/4
#' pure$tremble <- gamma
#' beta <- rnorm(1)
#' p <- 1/sum(1 + exp(beta))
#' sim.shares <- c(p, 1-p)
#' sim.strategies <- list("mixed" = mixed, "pure" = pure)
#' sim.data <- stratEst.simulate(strategies = sim.strategies, shares = sim.shares)
#'
#' @export
stratEst.simulate <- function( data = NULL, strategies, shares = NULL, coefficients = NULL, covariate.mat = NULL, num.ids = 100, num.games = 5, num.periods = NULL, fixed.assignment = TRUE, input.na = FALSE, sample.id = NULL ){

  ##################################################################################################
  # CHECK INPUTS
  ##################################################################################################

  # check strategies
  if( missing(strategies) ) {
    stop("stratEst error: Mandatory input object 'strategies' is missing. Specify an integer or create a data.frame object.")
  }else{
    stratEst.simulate.check.strategies.return <- stratEst.simulate.check.strategies( strategies )
    strategies_matrix_list <- stratEst.simulate.check.strategies.return$strategies.matrix.list
    num_samples_strategies <- length(strategies_matrix_list)
    trembles_list <- stratEst.simulate.check.strategies.return$trembles.list
    num_strats <- stratEst.simulate.check.strategies.return$num_strats
    unique_inputs <- stratEst.simulate.check.strategies.return$unique.inputs
    unique_outputs <- stratEst.simulate.check.strategies.return$unique.outputs
    num_inputs <- stratEst.simulate.check.strategies.return$num.unique.inputs
    num_outputs <- stratEst.simulate.check.strategies.return$num.unique.outputs
    sid <- stratEst.simulate.check.strategies.return$sid
    names_strategies <- stratEst.simulate.check.strategies.return$names.strategies
    names_strategies_lcr <- stratEst.simulate.check.strategies.return$names.strategies
  }

  if( is.null( shares ) & ( is.null( coefficients ) | is.null( covariate.mat ) ) ){
    stop("stratEst error: Either the input object 'shares' or the the inputs objects 'coefficients' and 'covariates' must be supplied.")
  }

  # check covariates
  if( is.null(covariate.mat) ){
    LCR = FALSE
  }
  else{
    stratEst.simulate.check.covariates.return <- stratEst.simulate.check.covariates(covariate.mat)
    covariate_mat <- stratEst.simulate.check.covariates.return$covariate.mat
    num.ids <- stratEst.simulate.check.covariates.return$num.ids
    num.covariates <- stratEst.simulate.check.covariates.return$num.covariates
    names.covariates <- stratEst.simulate.check.covariates.return$names.covariates
    LCR = TRUE
  }

  #check shares
  if( is.null(shares) == F ) {
    stratEst.simulate.check.shares.return <- stratEst.simulate.check.shares( shares , LCR , num_strats )
    shares = stratEst.simulate.check.shares.return$shares
    num_samples_shares = stratEst.simulate.check.shares.return$num_samples
    names_samples = stratEst.simulate.check.shares.return$names_samples
    if( num_samples_shares > 1 ){
      shares <- do.call(rbind,shares)
    }else{
      shares <- matrix(shares,1,length(shares))
    }
  }else{
    num_samples_shares = 1
    names_samples = "1"
  }

  #check coefficients
  if( is.null(coefficients) ) {
    if( LCR ){
      stop("stratEst error: Input object 'coefficients' must be supplied.")
    }
  }else{
    coefficient_mat <- stratEst.simulate.check.coefficients( coefficients , covariate_mat , num_strats , names_strategies_lcr )
    if( LCR == F ){
      warning("stratEst warning: No covariates specified. The input object 'coefficients' is ignored.");
    }
  }

  # check data
  if( is.null( data ) == F ){
    if( "stratEst.data" %in% class(data) | "data.frame" %in% class(data) ){
      if( is.null(covariate.mat) == F & num.ids != length(unique(data$id)) ){
        stop("stratEst error: The object 'covariate.mat' must have as many rows as there are unique ids in 'data'.")
      }
      num.ids <- length(unique(data$id))
      sample_vec <- rep(1,nrow(data))
      if( is.null(sample.id) == F ){
        if( sample.id %in% colnames(data)){
          sample_vec <- data[,sample.id]
        }
        else{
          stop(paste("stratEst error: The variable ",sample.id," cannot be found in the data.",sep=""))
        }
      }
      unique.samples.data <- unique(sample_vec)
      num.samples.data <- length(unique.samples.data)
      if( num_samples_shares > 1 & num_samples_shares != num.samples.data ){
        stop("stratEst error: Number of samples in 'data' and 'shares' must be the same.")
      }
      if( num_samples_strategies > 1 & num_samples_strategies != num.samples.data ){
        stop("stratEst error: Number of samples in 'data' and 'strategies' must be the same.")
      }
    }else{
      stop("stratEst error: Input object 'data' has to be of class 'stratEst.data'.")
    }
    stratEst.return.data <- stratEst.check.data( data )
    id <- stratEst.return.data$id
    game <- stratEst.return.data$game
    period <- stratEst.return.data$period
    input <- stratEst.return.data$input
    output <- stratEst.return.data$output
    input_factor <- stratEst.return.data$input.factor
    output_factor <- stratEst.return.data$output.factor
    levels_input <- stratEst.return.data$levels.input
    levels_output <- stratEst.return.data$levels.output
  }

  if( is.null( data ) ){
    # check num.ids
    num.ids <- as.integer(num.ids)
    if( "integer" %in% class(num.ids) == F ){
      stop("stratEst error: Input object 'num.ids' has to be of class 'integer'.")
    }else{
      if( num.ids < 0 ){
        stop("stratEst error: Input object 'num.ids' cannot be negative.")
      }
    }


    # check num.games
    if( "numeric" %in% class(num.games) == F ){
      stop("stratEst error: Input object 'num.games' has to be of class 'numeric'.")
    }else{
      num.games <- round(num.games)
      if( num.games < 0 ){
        stop("stratEst error: Input object 'num.games' cannot be negative.")
      }
    }

    # check num.periods
    if( is.null( num.periods ) ){
      num.periods <- rep(5,num.games)
    }else{
      if( "numeric" %in% class(num.periods) == F ){
        stop("stratEst error: Input object 'num.periods' has to be of class 'numeric'.")
      }else{
        if( length(num.periods) != num.games ){
          stop("stratEst error: Input object 'num.periods' has to be a vector of length 'num.games'.")
        }
        num.periods <- round(num.periods)
        if( any(num.periods < 0 ) ){
          stop("stratEst error: Elements of the input object 'num.periods' cannot be negative.")
        }
      }
    }

    # check fixed.assignment
    if( "logical" %in% class(fixed.assignment) == F ){
      stop("stratEst error: Input object 'fixed.assignment' has to be of class 'logical'.")
    }
  }

  if( is.null( data ) ){
    # balance num sample information
    if( num_samples_strategies == num_samples_shares ){
      num_samples = num_samples_strategies
    }
    else if( num_samples_strategies == 1 & num_samples_shares > 1 ){
      num_samples = num_samples_shares
      for( sam in 2:num_samples ){
        strategies_matrix_list[[sam]] <- strategies_matrix_list[[1]]
        trembles_list[[sam]] <- trembles_list[[1]]
      }
    }
    else if( num_samples_strategies > 1 & num_samples_shares == 1 ){
      num_samples = num_samples_strategies
      shares <- t(matrix( shares , length(shares) , num_samples_strategies ))
    }
    else{
      stop("stratEst error: The number of samples indicated in input object 'strategies' and 'shares' does not match.")
    }
  }else{

    num_samples = length(unique(sample_vec))
    # balance num sample information
    if( num_samples_strategies == 1 & num_samples > 1 ){
      for( sam in 2:num_samples ){
        strategies_matrix_list[[sam]] <- strategies_matrix_list[[1]]
        trembles_list[[sam]] <- trembles_list[[1]]
      }
    }
    if( num_samples_shares == 1 & num_samples > 1 ){
      shares <- t(matrix( shares , length(shares) , num_samples ))
    }
    if( ( num_samples_shares > 1 & num_samples_shares != num_samples ) ){
      stop("stratEst error: The number of samples indicated in input object 'shares' does not match the number of samples in the data.")
    }
    if( ( num_samples_strategies > 1 & num_samples_strategies != num_samples ) ){
      stop("stratEst error: The number of samples indicated in input object 'strategies' does not match the number of samples in the data.")
    }
  }

  ##################################################################################################
  # GENERATE sim.data
  ##################################################################################################

  # check fixed.assignment
  if( is.null( data ) ){

    sim.data <- NULL
    max_id_sample <- 0

    for( sam in 1:num_samples ){

      strategies_matrix <- strategies_matrix_list[[sam]]
      trembles <- trembles_list[[sam]]
      trembles[is.na(trembles)] <- 0
      trembles <- trembles
      response_mat <- matrix(strategies_matrix[,2:(num_outputs+1)],nrow(strategies_matrix),num_outputs)
      tremble_mat <- matrix(rep(trembles,ncol(response_mat)),nrow(response_mat),ncol(response_mat))
      response_mat <- response_mat*(1-tremble_mat) + (1-response_mat)*tremble_mat/(ncol(response_mat)-1)
      transition_mat <- matrix(strategies_matrix[,(num_outputs+2):(num_inputs+1+num_outputs)],nrow(strategies_matrix),num_inputs)

      if( LCR ){
        priors_individuals = exp( covariate_mat %*% coefficient_mat)/( apply(exp(covariate_mat %*% coefficient_mat),1,sum) )
      }else{
        priors_individuals = t( matrix( shares[sam,] , ncol(shares) , num.ids ) )
      }

      # assign ids to strategies (strategies assignmen mat)
      strategy_assignment_mat <- matrix(NA,num.ids,num.games)
      if( fixed.assignment ){
        rand_vec = stats::runif(num.ids)
        for( i in 1:num.ids){
          assigned = F
          for( s in 1:num_strats ){
            if( assigned == F ){
              if( rand_vec[i] <= sum( priors_individuals[i,(1:s)]) ){
                strategy_assignment_mat[i,] <- rep( s , num.games )
                assigned = T
              }
            }
          }
        }
      }else{
        rand_mat = matrix(stats::runif(num.ids*num.games),num.ids,num.games)
        for( i in 1:num.ids){
          for( j in 1:num.games ){
            assigned = F
            for( s in 1:num_strats ){
              if( assigned == F ){
                if( rand_mat[i,j] <= sum( priors_individuals[i,(1:s)]) ){
                  strategy_assignment_mat[i,j] <- s
                  assigned = T
                }
              }
            }
          }
        }
      }

      # generate variables
      if( fixed.assignment ){
        id = max_id_sample + rep(c(1:num.ids),each=sum(num.periods))
      }else{
        id <- NULL
        for( i in 1:num.games ){
          id = max_id_sample + c(id,rep(c(((i-1)*num.ids + 1):((i-1)*num.ids + num.ids)),each=num.periods[i]))
        }
      }
      game = rep( NA,num.ids*sum(num.periods) )
      period = rep(NA,num.ids*sum(num.periods))
      input = rep(0,num.ids*sum(num.periods))
      output = rep(NA,num.ids*sum(num.periods))
      sample = rep(names_samples[sam],num.ids*sum(num.periods))
      strategy.id = rep(NA,num.ids*sum(num.periods))
      if( LCR ){
        covariate_vars = matrix(NA,num.ids*sum(num.periods),num.covariates)
      }

      # generate output
      counter <- 1
      for( i in 1:num.ids ){
        for( g in 1:num.games ){
          sid_sbj = strategy_assignment_mat[i,g]
          indices_strat <- sid == sid_sbj
          response_mat_sbj <- matrix(response_mat[ indices_strat , ],sum(indices_strat),ncol(response_mat))
          response_mat_sbj <- abs(response_mat_sbj)
          transition_mat_sbj <- matrix(transition_mat[ indices_strat , ],sum(indices_strat),ncol(transition_mat))
          state <- 1
          for( p in 1:num.periods[g] ){
            strategy.id[counter] <- names_strategies[sid_sbj]
            game[counter] <- g
            period[counter] <- p
            if( p > 1 | input.na == FALSE ){
              input[counter] <- t(stats::rmultinom( 1 , 1 , rep(1,num_inputs)/num_inputs)) %*% c(1:num_inputs)
              state <- transition_mat_sbj[state,input[counter]]
            }
            output[counter] = t(stats::rmultinom( 1 , 1 , response_mat_sbj[state,] )) %*% c(1:num_outputs)
            if( LCR ){
                covariate_vars[counter,] <- covariate_mat[i,]
            }
            counter <- counter + 1
          }
        }
      }

      strategy <- strategy.id

      unique_output_values <- gsub("prob.", "", unique_outputs)
      unique_input_values <- gsub("tr\\(", "", unique_inputs)
      unique_input_values <- gsub("\\)", "", unique_input_values)

      choice <- unique_output_values[output]
      unique_input_values <- c(NA,unique_input_values)
      input <- unique_input_values[input+1]

      max_id_sample <- max(id)

      if( LCR ){
        sim.data_sample <- data.frame(id,game,period,choice,input,sample,strategy,covariate_vars)
      }else{
        sim.data_sample <- data.frame(id,game,period,choice,input,sample,strategy)
      }

      sim.data <- rbind(sim.data,sim.data_sample)

    }

    # define classes
    sim.data$input <- as.factor(sim.data$input)
    sim.data$choice <- as.factor(sim.data$choice)
    sim.data$sample <- as.factor(sim.data$sample)
    sim.data$id <- as.integer(sim.data$id)
    sim.data$strategy <- as.factor(sim.data$strategy)


    if( LCR ){
      colnames(sim.data)[8:(8+num.covariates-1)] <- colnames(covariate.mat)
    }

    # make object of class stratEst.sim.data
    attr(sim.data, "class") <- c("stratEst.data","data.frame")

  }else{ # if data is not NULL
    sim.data <- NULL

    # generate variables
    output = rep(NA,nrow(data))
    sample = sample_vec
    strategy.id = rep(NA,nrow(data))
    if( LCR ){
      covariate_vars = matrix(NA,nrow(data),num.covariates)
    }

    for( sam in 1:num_samples ){

      unique.ids.sample <- unique(data$id[sample_vec==unique.samples.data[sam]])
      num.ids.sample <- length(unique.ids.sample)

      strategies_matrix <- strategies_matrix_list[[sam]]
      trembles <- trembles_list[[sam]]
      trembles[is.na(trembles)] <- 0
      trembles <- trembles
      response_mat <- matrix(strategies_matrix[,2:(num_outputs+1)],nrow(strategies_matrix),num_outputs)
      tremble_mat <- matrix(rep(trembles,ncol(response_mat)),nrow(response_mat),ncol(response_mat))
      response_mat <- response_mat*(1-tremble_mat) + (1-response_mat)*tremble_mat/(ncol(response_mat)-1)
      transition_mat <- matrix(strategies_matrix[,(num_outputs+2):(num_inputs+1+num_outputs)],nrow(strategies_matrix),num_inputs)

      if( LCR ){
        priors_individuals = exp( covariate_mat %*% coefficient_mat)/( apply(exp(covariate_mat %*% coefficient_mat),1,sum) )
      }else{
        priors_individuals = t( matrix( shares[sam,] , ncol(shares) , num.ids.sample ) )
      }

      # assign ids to strategies (strategies assignmen mat)
      strategy_assignment_vec <- rep(NA,num.ids.sample)
      rand_vec = stats::runif(num.ids.sample)
      for( i in 1:num.ids.sample){
        assigned = F
        for( s in 1:num_strats ){
          if( assigned == F ){
            if( rand_vec[i] <= sum( priors_individuals[i,(1:s)]) ){
              strategy_assignment_vec[i] <- s
              assigned = T
            }
          }
        }
      }

      # generate output
      for( i in 1:num.ids.sample ){
        unique.games <- unique(game[id==unique.ids.sample[i]])
        num.games <- length(unique.games)
        sid_sbj = strategy_assignment_vec[i]
        indices_strat <- sid == sid_sbj
        response_mat_sbj <- matrix(response_mat[ indices_strat , ],sum(indices_strat),ncol(response_mat))
        response_mat_sbj <- abs(response_mat_sbj)
        transition_mat_sbj <- matrix(transition_mat[ indices_strat , ],sum(indices_strat),ncol(transition_mat))
        for( g in 1:num.games ){
          unique.periods <- unique(period[id==unique.ids.sample[i] & game == unique.games[g]])
          num.periods <- length(unique.periods)
          state <- 1
          for( p in 1:num.periods ){
            index <- c(1:nrow(data))[id == unique.ids.sample[i] & game == unique.games[g] & period == unique.periods[p]]
            strategy.id[index] <- names_strategies[sid_sbj]
            if( input[index] > 0 ){
              state <- transition_mat_sbj[state,input[index]]
            }
            output[index] = t(stats::rmultinom( 1 , 1 , response_mat_sbj[state,] )) %*% c(1:num_outputs)
            if( LCR ){
              covariate_vars[index,] <- covariate_mat[i,]
            }
          }
        }
      }

    }

    strategy <- strategy.id

    unique_output_values <- gsub("prob.", "", unique_outputs)
    unique_input_values <- gsub("tr\\(", "", unique_inputs)
    unique_input_values <- gsub("\\)", "", unique_input_values)

    choice <- unique_output_values[output]
    sample <- unique.samples.data[sample]
    unique_input_values <- c(NA,unique_input_values)
    suppressWarnings(unique_input_values <- as.numeric(unique_input_values))
    input <- unique_input_values[input+1]

    if( LCR ){
      sim.data <- data.frame(id,game,period,choice,input,sample,strategy,covariate_vars)
    }else{
      sim.data <- data.frame(id,game,period,choice,input,sample,strategy)
    }

    # define classes
    sim.data$choice <- as.factor(sim.data$choice)
    sim.data$sample <- as.factor(sim.data$sample)
    sim.data$strategy <- as.factor(sim.data$strategy)

    data$choice <- sim.data$choice
    data$sample <- sim.data$sample
    data$strategy <- sim.data$strategy

    if(LCR){
      num.cols.data <- ncol(data)
      data <- cbind(data,sim.data[8:(8+num.covariates-1)])
      colnames(data)[(num.cols.data+1):(num.cols.data+num.covariates)] <- colnames(covariate.mat)
    }

    sim.data <- data

    # make object of class stratEst.sim.data
    attr(sim.data, "class") <- c("stratEst.data","data.frame")


  }

  return(sim.data)
}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
