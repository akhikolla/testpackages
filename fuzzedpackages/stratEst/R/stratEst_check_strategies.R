# checks the input object strategies

stratEst.check.strategies <- function( strategies , input_factor , output_factor , input , output , select_strategies ){

  #check strategies
  unique_outputs = sort( unique( output_factor ) )
  unique_inputs = sort( unique( input_factor ) )
  num_unique_inputs = length( unique_inputs )
  num_unique_outputs = length( unique_outputs )
  response_mat_col_index <- NULL
  if( class(strategies) == "numeric" ){
    strategies = as.integer(strategies)
    integer_strategies = T
    num_strats = strategies
    input_has_na <- as.numeric(any( input == 0 ) )
    strategy_states = rep(c(1:(num_unique_inputs + input_has_na)),num_strats)
    response_par = matrix(NA,(num_unique_inputs + input_has_na)*num_strats,num_unique_outputs)
    transition_mat = rep(1,(num_unique_inputs + input_has_na)*num_strats) %*% t.default(c(( 1 + input_has_na ):( num_unique_inputs + input_has_na )))
    strategies_matrix <- cbind(strategy_states,response_par,transition_mat)
    trembles <- rep( NA , nrow(strategies_matrix) )
    sid <- rep(c(1:num_strats), each = num_unique_inputs + input_has_na )
    names_strategies <- NULL
    for( s in 1:num_strats){
        names_strategies <-  c(names_strategies,paste("strategy.",s,sep = ""))
    }
  }
  else if( class(strategies) == "list" ){
    integer_strategies = F
    num_strats = length(strategies)
    state <- NULL
    transition_mat <- NULL
    response_mat <- NULL
    trembles <- NULL
    sid <- NULL
    names_strategies <- names(strategies)
    if( is.null(names_strategies) ){
      for( s in 1:num_strats){
        names_strategies <-  c(names_strategies,paste("strategy.",s,sep = ""))
      }
    }
    response_mat_col_index <- matrix(NA,num_strats,num_unique_outputs)
    transition_mat_col_index <- matrix(NA,num_strats,num_unique_inputs)
    for( le in 1:num_strats ){
      strategy <- strategies[[le]]
      if( "stratEst.strategy" %in%  class(strategy) == F & "data.frame" %in% class(strategy) == F ){
        stop(paste("stratEst error: Strategy ",le," in the list 'strategies' is not an object of class stratEst.strategy.",sep=""))
      }
      state_strategy <- c(1:nrow(strategy))
      state <- c(state,state_strategy)
      # check and generate responses
      response_mat_strategy <- matrix(NA,nrow(strategy),num_unique_outputs)
      for( out in 1:num_unique_outputs ){
        r_string <- paste("prob.",as.character(unique_outputs[out]),sep="")
        if( r_string %in% colnames(strategy) ){
          response_mat_strategy[,out] <- strategy[,r_string]
          response_mat_col_index[le,out] <- grep(paste("^",r_string,"$",sep=""), colnames(strategy))
        }
        else{
          message <- paste("stratEst error: There is an output with value ", as.character(unique_outputs[out]) , " in the data but there is no column named '", r_string , "' in strategy",le,".",sep="")
          stop(message)
        }
      }
      response_mat <- rbind(response_mat,response_mat_strategy)
      # check and generate transitions
      transition_mat_strategy <- matrix(NA,nrow(strategy),num_unique_inputs)
      for( ins in 1:num_unique_inputs ){
        t_string <- paste("tr(",as.character(unique_inputs[ins]), ")" , sep="")
        if( t_string %in% colnames(strategy) ){
          transition_mat_strategy[,ins] <- strategy[,t_string]
          transition_mat_col_index[le,ins] <- c(1:ncol(strategy))[t_string == colnames(strategy)]
        }
        else{
          transition_mat_strategy[,ins] <- rep(1,nrow(strategy))
          #transition_mat_col_index[le,ins] <- grep(paste("^",t_string,"$",sep=""), colnames(strategy))
          # message <- paste("stratEst error: There is an input with value ", as.character(unique_inputs[ins]) , " in the data but there is no column named '", t_string , "' in strategy",le,".",sep="")
          # stop(message)
        }
        if( sum( transition_mat_strategy[,ins]%%1==0 ) < nrow(strategy) ){
          message <- paste("stratEst error: The transition columns in 'strategies' must be integers. Check the column named '", t_string , "'.",sep="")
          stop(message)
        }
      }
      if( paste("tr(",as.character(unique_inputs[1]), ")" , sep="") %in% colnames(strategy) == F ){
        strategy_with_transitions <- cbind(strategies[[le]],transition_mat_strategy)
        colnames( strategy_with_transitions) <- c(colnames(strategy),paste("tr(",as.character(unique_inputs),")",sep=""))
        strategies[[le]] <- strategy_with_transitions
      }
      transition_mat <- rbind(transition_mat,transition_mat_strategy)
      if( any( transition_mat <= 0 ) ){
        stop("stratEst error: No negative values allowed in the input columns of 'strategies'.");
      }

      # check for tremble column
      if( "tremble" %in% colnames(strategy) ){
        trembles_strategy <- strategy[,"tremble"]
      }
      else{
        trembles_strategy <- rep( NA , nrow(strategy) )
      }
      trembles <- c(trembles,trembles_strategy)

      # generate strategy id if missing
      sid_strategy <- rep(le,nrow(strategy))
      sid <- c(sid,sid_strategy)

    }
    strategies_matrix = cbind(state,response_mat,transition_mat)

  }else{
    stop("stratEst error: The input object 'strategies' must be an integer or a list of data frames.");
  }

  if ( select_strategies && num_strats == 1 ){
    stop("stratEst error: Strategies cannot be selected if there is only one strategy.");
  }

  stratEst.check.strategies.return <- list( "strategies" = strategies , "strategies.matrix" = strategies_matrix , "trembles" = trembles , "num_strats" = num_strats , "unique.inputs" = unique_inputs , "unique.outputs" = unique_outputs , "num.unique.inputs" = num_unique_inputs , "num.unique.outputs" = num_unique_outputs , "sid" = sid , "integer.strategies" = integer_strategies , "response.mat.col.index" = response_mat_col_index , "names.strategies" = names_strategies )

  return( stratEst.check.strategies.return )

}
