# checks the input object strategies

stratEst.simulate.check.strategies <- function( strategies  ){

  strategies_matrix_list <- list(NULL)
  trembles_list <- list(NULL)

  #check strategies
  if( "list" %in% class(strategies) ){
    if( "stratEst.strategy" %in%  class(strategies[[1]]) | "data.frame" %in% class(strategies[[1]]) ){
      num_samples_strategies = 1
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
      unique_outputs <- sort(grep("prob.", names(strategies[[1]]), value=TRUE))
      num_unique_outputs <- length(unique_outputs)
      unique_inputs <- sort(grep("tr\\(", names(strategies[[1]]), value=TRUE))
      num_unique_inputs <- length(unique_inputs)
      response_mat_col_index <- matrix(NA,num_strats,num_unique_outputs)
      transition_mat_col_index <- matrix(NA,num_strats,num_unique_inputs)
      for( le in 1:num_strats ){
        strategy <- strategies[[le]]
        if( "stratEst.strategy" %in%  class(strategy) == F & "data.frame" %in% class(strategy) == F ){
          stop(paste("stratEst error: Strategy ",le," in the list 'strategies' is not an object of class data.frame.",sep=""))
        }
        state_strategy <- c(1:nrow(strategy))
        state <- c(state,state_strategy)
        # check and generate responses
        response_mat_strategy <- matrix(NA,nrow(strategy),num_unique_outputs)
        for( out in 1:num_unique_outputs ){
          r_string <- unique_outputs[out]
          response_mat_strategy[,out] <- strategy[,r_string]
        }
        if( any(is.na(response_mat_strategy) ) ){
          stop("stratEst error: Responses in input object 'strategies' cannot contain NA values.")
        }else{
          response_mat <- rbind(response_mat,response_mat_strategy)
        }
        # check and generate transitions
        transition_mat_strategy <- matrix(NA,nrow(strategy),num_unique_inputs)
        for( ins in 1:num_unique_inputs ){
          t_string <- unique_inputs[ins]
          transition_mat_strategy[,ins] <- strategy[,t_string]
        }
        transition_mat <- rbind(transition_mat,transition_mat_strategy)

        # check for tremble column
        if( "tremble" %in% colnames(strategy) ){
          trembles_strategy <- strategy[,"tremble"]
        }
        else{
          trembles_strategy <- rep( 0 , nrow(strategy) )
        }
        trembles <- c(trembles,trembles_strategy)
        if( any(is.na(trembles) ) ){
          trembles[is.na(trembles)] <- 0
          #warning("stratEst warning: NA values in trembles in the input object 'strategies' are converted to zero.")
        }


        # generate strategy id if missing
        sid_strategy <- rep(le,nrow(strategy))
        sid <- c(sid,sid_strategy)

      }
      strategies_matrix = cbind(state,response_mat,transition_mat)
      strategies_matrix_list[[1]] <- strategies_matrix
      trembles_list[[1]] <- trembles
    }
    else if( "list" %in% class( strategies[[1]] ) ){
      strategies_list <- strategies
      num_samples_strategies = length(strategies_list)
      names_strategies <- names(strategies_list[[1]])
      if( is.null(names_strategies) ){
        for( s in 1:strategies_list[[1]]){
          names_strategies <-  c(names_strategies,paste("strategy.",s,sep = ""))
        }
      }
      for( sam in 1:num_samples_strategies){
        if( "list" %in% class( strategies_list[[sam]] ) == F ){
          stop(paste("stratEst error: Element ",sam, "of list 'strategies' is not an object of class list.",sep=""))
        }
        strategies <- strategies_list[[sam]]
        num_strats <- length(strategies)
        state <- NULL
        transition_mat <- NULL
        response_mat <- NULL
        trembles <- NULL
        sid <- NULL
        unique_outputs <- sort(grep("prob.", names(strategies[[1]]), value=TRUE))
        num_unique_outputs <- length(unique_outputs)
        unique_inputs <- sort(grep("tr\\(", names(strategies[[1]]), value=TRUE))
        num_unique_inputs <- length(unique_inputs)
        response_mat_col_index <- matrix(NA,num_strats,num_unique_outputs)
        transition_mat_col_index <- matrix(NA,num_strats,num_unique_inputs)
        for( le in 1:num_strats ){
          strategy <- strategies[[le]]
          if( "data.frame" %in% class(strategy) == F ){
            stop(paste("stratEst error: Strategy ",le," in list", sam ,"of list 'strategies' is not an object of class data.frame.",sep=""))
          }
          state_strategy <- c(1:nrow(strategy))
          state <- c(state,state_strategy)
          # check and generate responses
          response_mat_strategy <- matrix(NA,nrow(strategy),num_unique_outputs)
          for( out in 1:num_unique_outputs ){
            r_string <- unique_outputs[out]
            if( r_string %in% colnames(strategy) ){
              response_mat_strategy[,out] <- strategy[,r_string]
              #response_mat_col_index[le,out] <- grep(paste("^",r_string,"$",sep=""), colnames(strategy))
            }
            else{
              message <- paste("stratEst error: There is an output with value ", as.character(unique_outputs[out]) , " in the data but there is no column named '", r_string , "' in strategy",le,".",sep="")
              stop(message)
            }
          }
          if( any(is.na(response_mat_strategy)) ){
            stop("stratEst error: Response parameters in input object 'strategies' cannot be NA.");
          }else{
            response_mat <- rbind(response_mat,response_mat_strategy)
          }

          # check and generate transitions
          transition_mat_strategy <- matrix(NA,nrow(strategy),num_unique_inputs)
          for( ins in 1:num_unique_inputs ){
            t_string <- unique_inputs[ins]
            if( t_string %in% colnames(strategy) ){
              transition_mat_strategy[,ins] <- strategy[,t_string]
              #transition_mat_col_index[le,ins] <- grep(paste("^",t_string,"$",sep=""), colnames(strategy))
            }
            else{
              message <- paste("stratEst error: There is an input with value ", as.character(unique_inputs[ins]) , " in the data but there is no column named '", t_string , "' in strategy",le,".",sep="")
              stop(message)
            }
            if( sum( transition_mat_strategy[,ins]%%1==0 ) < nrow(strategy) ){
              message <- paste("stratEst error: The transition columns in 'strategies' must be integers. Check the column named '", t_string , "'.",sep="")
              stop(message)
            }
          }
          transition_mat <- rbind(transition_mat,transition_mat_strategy)

          # check for tremble column
          if( "tremble" %in% colnames(strategy) ){
            trembles_strategy <- strategy[,"tremble"]
          }
          else{
            trembles_strategy <- rep( 0 , nrow(strategy) )
          }
          trembles_strategy[is.na(trembles_strategy)] = 0
          trembles <- c(trembles,trembles_strategy)

          # generate strategy id if missing
          sid_strategy <- rep(le,nrow(strategy))
          sid <- c(sid,sid_strategy)

        }
        strategies_matrix = cbind(state,response_mat,transition_mat)
        strategies_matrix_list[[sam]] <- strategies_matrix
        trembles_list[[sam]] <- trembles
      }
    }else{
      stop("stratEst error: The elements of the list 'strategies' must be data frames or lists of data frames.");
    }

  }else{
    stop("stratEst error: The input object 'strategies' must be a list of data frames for one sample or a list of a lists of data frames for more than one sample.");
  }

  stratEst.simulate.check.strategies.return <- list( "strategies.matrix.list" = strategies_matrix_list , "trembles.list" = trembles_list , "num_strats" = num_strats , "unique.inputs" = unique_inputs , "unique.outputs" = unique_outputs , "num.unique.inputs" = num_unique_inputs , "num.unique.outputs" = num_unique_outputs , "sid" = sid , "response.mat.col.index" = response_mat_col_index , "names.strategies" = names_strategies )

  return( stratEst.simulate.check.strategies.return )

}
