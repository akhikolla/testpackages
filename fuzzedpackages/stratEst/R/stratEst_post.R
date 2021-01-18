# This function executes the post-processing after the estimation

stratEst.post <- function( data , cpp.output , stratEst.return , strategies , covariates , response , unique_ids , num_unique_ids , input , output , unique_inputs , unique_outputs , num_unique_inputs , num_unique_outputs , sample , sample.id , sample_factor , num_samples , specific_shares , specific_probs , specific_trembles , sample_is_factor , integer_strategies , LCR , probs_mat_col_index , crit , num_obs , se , quantiles, output_factor, input.is.null ){

  num_probs_to_est = length(stratEst.return$probs.par)
  num_trembles_to_est = length(stratEst.return$trembles.par)

  if( num_samples == 1 | ( (specific_probs == F | num_probs_to_est == 0) & (specific_trembles == F | num_trembles_to_est == 0) ) ){
    # strategies post-processing
    stratEst.return$strategies[ stratEst.return$strategies == -1 ] = NA
    num_strategies = sum(stratEst.return$strategies[,1] == 1)
    num_strategies_sample = num_strategies
    num_unique_inputs <- length(unique(input[input!=0]))
    post_sid <- cpp.output$sid
    unique_post_sids <- unique(post_sid)
    if( integer_strategies ){
      state <- stratEst.return$strategies[,1]
      prob_mat <- matrix(stratEst.return$strategies[,2:(1+num_unique_outputs)],nrow(stratEst.return$strategies), num_unique_outputs)
      r_names <- rep(NA,num_unique_outputs)
      for( outs in 1:num_unique_outputs ){
        r_names[outs] <- paste("prob.",as.character(unique_outputs[outs]),sep="")
      }
      colnames(prob_mat) <- r_names

      transition_mat <- matrix(stratEst.return$strategies[, ((1+num_unique_outputs+1):(1+num_unique_outputs+num_unique_inputs))],nrow(stratEst.return$strategies), num_unique_inputs)
      t_names <- rep(NA,num_unique_inputs)
      for( ins in 1:num_unique_inputs ){
        t_names[ins] <- paste( "tr(" , as.character(unique_inputs[ins]) , ")" , sep ="" )
      }
      colnames(transition_mat) <- t_names
      tremble_vec <- cpp.output$trembles
      colnames(tremble_vec) <- "tremble"
      if( input.is.null == F ){
        strategies_mat <- as.data.frame(cbind(prob_mat,tremble_vec,transition_mat))
      }else{
        strategies_mat <- as.data.frame(cbind(prob_mat,tremble_vec))
      }

      names_of_strategies <- NULL
      for( s in 1:num_strategies ){
        names_of_strategies <- c(names_of_strategies,paste("strategy.",as.character(s),sep=""))
      }
      stratEst.return$strategies <- strategies_mat
    }
    else if( "list" %in% class(strategies) ){
      strategies_list <- list(NULL)
      names_of_strategies <- names(strategies)
      if( is.null(names_of_strategies) ){
        for( s in 1:length(strategies) ){
          names_of_strategies <- c(names_of_strategies,paste("strategy.",as.character(s),sep=""))
        }
      }
      for( strs in 1:num_strategies){
        strategy <- strategies[[unique_post_sids[strs]]]
        for( out in 1:num_unique_outputs ){
          strategy[,probs_mat_col_index[unique_post_sids[strs],out]] <- stratEst.return$strategies[post_sid == unique_post_sids[strs],(1+out)]
        }
        prob_mat <- strategy[,probs_mat_col_index[unique_post_sids[strs],1:num_unique_outputs]]
        trembles_strategy <- cpp.output$trembles[post_sid == unique_post_sids[strs]]
        if( any(trembles_strategy != 0) ){
          strategy$tremble <- trembles_strategy
        }else{
          strategy$tremble <- rep(NA,length(trembles_strategy))
        }
        for( line in 1:nrow(prob_mat)){
          if( any( is.na(prob_mat[line,] )) ){
            strategy$tremble[line] = NA
          }
          else{
            if( is.na(strategy$tremble[line]) == F  ){
              if( any(prob_mat[line,] != 0  &  prob_mat[line,] != 1)  & strategy$tremble[line] == 0 ){
                strategy$tremble[line] = NA
              }
            }
          }
        }
        if( input.is.null ){
          strategy <- strategy[, -grep(")", colnames(strategy))]
        }
        strategies_list[[strs]] <- strategy
      }
      stratEst.return$strategies <- strategies_list
    }

    # create list of strategies
    if( "list" %in% class(strategies) == F ){
      strategies_list <- list(NULL)
      for( i in 1:num_strategies ){
        strategy <- stratEst.return$strategies[ post_sid == unique_post_sids[i] , ]
        rownames(strategy) <- c(1:nrow(strategy))
        strategies_list[[i]] <- strategy
      }
    }

  }else{  # num_samples > 1
    # strategies post-processing
    stratEst.return$strategies[ stratEst.return$strategies == -1 ] = NA
    num_strategies = sum(stratEst.return$strategies[,1] == 1)
    num_strategies_sample <- (num_strategies/num_samples)
    num_unique_inputs <- length(unique(input[input!=0]))
    original_sid <- cpp.output$sid
    post_sid <- original_sid
    post_sid <- rep( post_sid , num_samples ) + rep( c(0:(num_samples-1)) , each = length(post_sid))*max(post_sid)
    unique_post_sids <- unique(post_sid)
    unique_original_sids <- unique(original_sid)
    if( integer_strategies ){
      prob_mat <- cpp.output$responses
      r_names <- rep(NA,num_unique_outputs)
      for( outs in 1:num_unique_outputs ){
        r_names[outs] <- paste("prob.",as.character(unique_outputs[outs]),sep="")
      }
      colnames(prob_mat) <- r_names

      transition_mat <- matrix(stratEst.return$strategies[, ((1+num_unique_outputs+1):(1+num_unique_outputs+num_unique_inputs))],nrow(stratEst.return$strategies), num_unique_inputs)
      t_names <- rep(NA,num_unique_inputs)
      for( ins in 1:num_unique_inputs ){
        t_names[ins] <- paste( "tr(" , as.character(unique_inputs[ins]) , ")" , sep ="" )
      }
      colnames(transition_mat) <- t_names
      tremble_vec <- cpp.output$trembles
      colnames(tremble_vec) <- "tremble"
      if( input.is.null == F ){
        strategies_mat <- as.data.frame(cbind(prob_mat,tremble_vec,transition_mat))
      }else{
        strategies_mat <- as.data.frame(cbind(prob_mat,tremble_vec))
      }

      names_of_strategies <- NULL

      for( s in 1:num_strategies_sample ){
        names_of_strategies <- c(names_of_strategies,paste("strategy.",as.character(s),sep=""))
      }
      stratEst.return$strategies <- strategies_mat

    }
    else if( "list" %in% class(strategies) ){
      names_of_strategies <- names(strategies)
      if( is.null(names_of_strategies) ){
        for( s in 1:length(strategies) ){
          names_of_strategies <- c(names_of_strategies,paste("strategy.",as.character(unique_original_sids[s]),sep=""))
        }
      }
      strategies_list <- list(NULL)
      strategies_list_sample <- list(NULL)
      for( sam in 1:num_samples ){
        for( str in 1:num_strategies_sample){
          strategy <- strategies[[unique_post_sids[str]]]
          strs <- (sam-1)*num_strategies_sample + str
          for( out in 1:num_unique_outputs ){
            strategy[,probs_mat_col_index[unique_post_sids[str],out]] <- stratEst.return$strategies[post_sid == unique_post_sids[strs],(1+out)]
          }
          prob_mat <- strategy[,probs_mat_col_index[unique_post_sids[str],1:num_unique_outputs]]
          trembles_strategy <- cpp.output$trembles[post_sid == unique_post_sids[strs]]
          if( any(trembles_strategy != 0) ){
            strategy$tremble <- trembles_strategy
          }else{
            strategy$tremble <- rep(NA,length(trembles_strategy))
          }
          for( line in 1:nrow(prob_mat)){
            if( any( is.na(prob_mat[line,] )) ){
              strategy$tremble[line] = NA
            }
            else{
              if( is.na(strategy$tremble[line]) == F ){
                if( any(prob_mat[line,] != 0  &  prob_mat[line,] != 1)  & strategy$tremble[line] == 0 ){
                  strategy$tremble[line] = NA
                }
              }
            }
          }
          if( input.is.null ){
            strategy <- strategy[, -grep(")", colnames(strategy))]
          }
          strategies_list_sample[[str]] <- strategy
        }
        names(strategies_list_sample) <- names_of_strategies[unique_original_sids]
        strategies_list[[sam]] <- strategies_list_sample
      }
      stratEst.return$strategies <- strategies_list
    }
  }

  if( "list" %in% class(strategies)){
    names_of_strategies <- names_of_strategies[ cpp.output$selected.strats ]
  }

  # kill tremble iff all trembles are NA
  all.trembles.na = T
  strategies.clean <- stratEst.return$strategies
  if( "list" %in% class(stratEst.return$strategies[[1]]) ){
    for( i in 1:length(stratEst.return$strategies) ){
      for( j in 1:length(stratEst.return$strategies[[i]])){
        if( any( is.na(stratEst.return$strategies[[i]][[j]]$tremble ) == F ) ){
          all.trembles.na = F
        }
        strategies.clean[[i]][[j]] <- stratEst.return$strategies[[i]][[j]][,!colnames(stratEst.return$strategies[[i]][[j]])=="tremble"]
      }
    }
  }else{
    for( i in 1:length(stratEst.return$strategies) ){
      if( any( is.na(stratEst.return$strategies[[i]]$tremble ) == F ) ){
        all.trembles.na = F
      }
      strategies.clean[[i]] <- stratEst.return$strategies[[i]][,!colnames(stratEst.return$strategies[[i]])=="tremble"]
    }
  }
  if( all.trembles.na ){
    stratEst.return$strategies <- strategies.clean
  }


  # shares post-processing
  shares_mat <- stratEst.return$shares
  names_of_samples <- rep(NA,num_samples)
  if( num_samples > 1 ){
    if( sample_is_factor ){
      unique_factors <- sort(unique(sample_factor))
      for( smps in 1:num_samples ){
        names_of_samples[smps] <- paste(sample.id,as.character(unique_factors[smps]),sep=".")
      }
    }else{
      unique_samples <- sort(unique(sample))
      for( smps in 1:num_samples ){
        names_of_samples[smps] <-paste(sample.id,as.character(unique_samples[smps]),sep=".")
      }
    }
  }
  if( num_samples == 1 | specific_shares == F ){
    frame_shares_sample <- as.data.frame(t(shares_mat[,1]))
    colnames(frame_shares_sample) <- names_of_strategies
    rownames(frame_shares_sample) <- "share"
    stratEst.return$shares <- frame_shares_sample
    stratEst.return$shares = as.matrix(stratEst.return$shares)
  }else{
    shares <- list(NULL)
    for( i in 1:num_samples ){
      frame_shares_sample <- t(shares_mat[,i])
      rownames(frame_shares_sample) <- names_of_samples[i]

      colnames(frame_shares_sample) <- names_of_strategies
      shares[[i]] <- frame_shares_sample
    }
    names(shares) <- names_of_samples
    stratEst.return$shares <- shares
  }

  if( length(stratEst.return$shares.par) > 0 ){
    rownames(stratEst.return$shares.indices) = names_of_strategies
    if( num_samples == 1 | specific_shares == F ){
      colnames(stratEst.return$shares.indices) = "shares.par"
    }
    else{
      colnames(stratEst.return$shares.indices) = names_of_samples
    }
  }

  if( num_samples > 1 & ( (specific_probs & num_probs_to_est > 0) | (specific_trembles & num_trembles_to_est > 0) ) ){
    # create list of strategies
    if( "list" %in% class(strategies) == F ){
      strategies_list <- list(NULL)
      for( j in 1:num_samples ){
        strategies_list_sample <- list(NULL)
        counter_names <- 0
        for( i in 1:num_strategies_sample ){
          strategy <- stratEst.return$strategies[ post_sid == unique_post_sids[(num_strategies_sample*(j-1)+i)] , ]
          rownames(strategy) <- c(1:nrow(strategy))
          strats_sample_list <- strategy
          counter_names <- counter_names + nrow(strats_sample_list)
          strategies_list_sample[[i]] <- strats_sample_list
        }
        names(strategies_list_sample) <- names_of_strategies
        strategies_list[[j]] <- strategies_list_sample
      }
    }else{
      for( j in 1:num_samples ){
        strat_list_sample <- stratEst.return$strategies[[j]]
        names(strat_list_sample)  <- names_of_strategies
        stratEst.return$strategies[[j]] <- strat_list_sample
      }
    }
  }

  # assignments post-processing
  assignment_row_names <- NULL
  for( i in 1:num_unique_ids ){
    assignment_row_names <- c(assignment_row_names,paste( "id" , as.character(unique_ids[i]) , sep = "" ) )
  }
  rownames(stratEst.return$post.assignment) <- assignment_row_names
  colnames(stratEst.return$post.assignment) <- names_of_strategies

  # coefficients, covariate.mat and priors post-processing
  stratEst.return$coefficients.par <- cpp.output$coefficients.list$coefficients.par
  if( LCR & num_strategies > 1 ){
    colnames(stratEst.return$coefficients) <- names_of_strategies
    rownames(stratEst.return$coefficients) <- covariates

    colnames(stratEst.return$prior.assignment) <- names_of_strategies
    rownames(stratEst.return$prior.assignment) <- assignment_row_names
  }
  if( LCR ){
    rownames(stratEst.return$covariate.mat) <- assignment_row_names
    colnames(stratEst.return$covariate.mat) <- covariates
  }else{
    stratEst.return$covariate.mat <- NULL
  }

  # numbers
  stratEst.return$num.obs <- num_obs
  stratEst.return$num.ids <- num_unique_ids
  stratEst.return$num.par = length(stratEst.return$shares.par) + length(stratEst.return$probs.par) + length(stratEst.return$trembles.par) + length(stratEst.return$coefficients.par)
  stratEst.return$free.par = cpp.output$n.par
  stratEst.return$res.degrees = stratEst.return$num.ids - stratEst.return$free.par
  if( stratEst.return$res.degrees < 0 ){
    warning("stratEst warning: Residual degrees of freedom are negative. Parameter inference invalid.")
  }

  if( length(stratEst.return$shares.par) > 0 ){
    if( se == "bootstrap" ){
      stratEst.return$shares.quantiles = cpp.output$shares.list$shares.quantiles
    }
    else{
      shares.quantiles = matrix(NA,length(stratEst.return$shares.par),length(quantiles))
      for( i in 1:length(quantiles)){
        shares.quantiles[,i] =  stratEst.return$shares.par + stats::qt( quantiles[i] , df = stratEst.return$res.degrees )*stratEst.return$shares.se
      }
      stratEst.return$shares.quantiles = shares.quantiles
    }
    rownames(stratEst.return$shares.quantiles) <- paste("par.",as.character(seq(1,nrow(stratEst.return$shares.quantiles),by = 1)),sep="")
    colnames(stratEst.return$shares.quantiles) <- paste(as.character(quantiles*100),"%",sep="")
  }
  if( length(stratEst.return$probs.par) > 0 ){
    if( se == "bootstrap" ){
      stratEst.return$probs.quantiles = cpp.output$responses.list$responses.quantiles
    }
    else{
      probs.quantiles = matrix(NA,length(stratEst.return$probs.par),length(quantiles))
      for( i in 1:length(quantiles)){
        probs.quantiles[,i] =  stratEst.return$probs.par + stats::qt( quantiles[i] , df = stratEst.return$res.degrees )*stratEst.return$probs.se
      }
      stratEst.return$probs.quantiles = probs.quantiles
    }
    rownames(stratEst.return$probs.quantiles) <- paste("par.",as.character(seq(1,nrow(stratEst.return$probs.quantiles),by = 1)),sep="")
    colnames(stratEst.return$probs.quantiles) <- paste(as.character(quantiles*100),"%",sep="")
  }
  if( length(stratEst.return$trembles.par) > 0 ){
    if( se == "bootstrap" ){
      stratEst.return$trembles.quantiles = cpp.output$trembles.list$trembles.quantiles
    }
    else{
      trembles.quantiles = matrix(NA,length(stratEst.return$trembles.par),length(quantiles))
      for( i in 1:length(quantiles)){
        trembles.quantiles[,i] =  stratEst.return$trembles.par + stats::qt( quantiles[i] , df = stratEst.return$res.degrees )*stratEst.return$trembles.se
      }
      stratEst.return$trembles.quantiles = trembles.quantiles
    }
    rownames(stratEst.return$trembles.quantiles) <- paste("par.",as.character(seq(1,nrow(stratEst.return$trembles.quantiles),by = 1)),sep="")
    colnames(stratEst.return$trembles.quantiles) <- paste(as.character(quantiles*100),"%",sep="")
  }
  if( length(stratEst.return$coefficients.par) > 0 ){
    if( se == "bootstrap" ){
      stratEst.return$coefficients.quantiles = cpp.output$coefficients.list$coefficients.quantiles
    }
    else{
      coefficients.quantiles = matrix(NA,length(stratEst.return$coefficients.par),length(quantiles))
      for( i in 1:length(quantiles)){
        coefficients.quantiles[,i] =  stratEst.return$coefficients.par + stats::qt( quantiles[i] , df = stratEst.return$res.degrees )*stratEst.return$coefficients.se
      }
      stratEst.return$coefficients.quantiles = coefficients.quantiles
    }
    rownames(stratEst.return$coefficients.quantiles) <- paste("par.",as.character(seq(1,nrow(stratEst.return$coefficients.quantiles),by = 1)),sep="")
    colnames(stratEst.return$coefficients.quantiles) <- paste(as.character(quantiles*100),"%",sep="")
  }

  # gammas
  gammas = rep(NA,length(stratEst.return$trembles))
  indices.ok = (is.na(stratEst.return$trembles) == F & stratEst.return$trembles < 1 & stratEst.return$trembles != -1 )
  gammas[indices.ok] = -1/log(stratEst.return$trembles[indices.ok]/(1-stratEst.return$trembles[indices.ok]))
  gammas[indices.ok==F] = NA
  stratEst.return$gammas = matrix(gammas,length(gammas),1)
  #gamms.par
  gammas.par = rep(NA,length(stratEst.return$trembles.par))
  indices.ok = (is.na(stratEst.return$trembles.par) == F & stratEst.return$trembles.par < 1 & stratEst.return$trembles.par != -1)
  gammas.par[indices.ok] = -1/log(stratEst.return$trembles.par[indices.ok]/(1-stratEst.return$trembles.par[indices.ok]))
  gammas.par[indices.ok==F] = NA
  stratEst.return$gammas.par = matrix(gammas.par,length(gammas.par),1)
  # gammas.se
  gammas.se = rep(NA,length(stratEst.return$trembles.se))
  indices.ok = (is.na(stratEst.return$trembles.se) == F & stratEst.return$trembles.se < 1 & stratEst.return$trembles.se != -1)
  gammas.se[indices.ok] = -1/log(stratEst.return$trembles.se[indices.ok]/(1-stratEst.return$trembles.se[indices.ok]))
  gammas.se[indices.ok==F] = NA
  stratEst.return$gammas.se = matrix(gammas.se,length(gammas.se),1)

  #convergence post-processing
  stratEst.return$convergence[stratEst.return$convergence == -1] = NA
  if( num_strategies == 1 | LCR ){
    stratEst.return$convergence[1] = NA
  }
  if( response == "pure"){
    stratEst.return$convergence[2] = NA
  }
  if( sum(is.na(stratEst.return$convergence)) < length(stratEst.return$convergence) ){
    conv_shares <- stratEst.return$convergence[1]
    conv_probs <- stratEst.return$convergence[2]
    conv_trembles <- stratEst.return$convergence[3]
    conv_coefficients <- stratEst.return$convergence[4]
    convergence_vec <- cbind(conv_shares,conv_probs,conv_trembles,conv_coefficients)
    stratEst.return$convergence <- matrix(convergence_vec,1,length(convergence_vec))
    rownames(stratEst.return$convergence) <- "max.abs.score"
    colnames(stratEst.return$convergence) <- c("shares","probs","trembles","coefficients")
    #stratEst.return$convergence <- stratEst.return$convergence[, colSums(is.na(stratEst.return$convergence)) != nrow(stratEst.return$convergence)]
  }else{
    stratEst.return$convergence = NULL
  }


  convergence_string <- "no parameters estimated"

  # add names to strategies list
  if( "list" %in% class(strategies) == F ){
    stratEst.return$strategies <- strategies_list
    if( num_samples == 1 | ( (specific_probs == F | num_probs_to_est == 0) & (specific_trembles == F | num_trembles_to_est == 0) ) ){
      names(stratEst.return$strategies) <- names_of_strategies
    }else{
      names(stratEst.return$strategies) <- names_of_samples

    }
  }
  else{
    if( num_samples == 1 | ( (specific_probs == F | num_probs_to_est == 0) & (specific_trembles == F | num_trembles_to_est == 0) ) ){
      names(stratEst.return$strategies) <- names_of_strategies
    }
    else{
      names(stratEst.return$strategies) <- names_of_samples

    }
  }

  # strategy matrix row names
  if( "list" %in% class(stratEst.return$strategies[[1]])  ){
    strategies_sample_list <- NULL
    strategies_print <- stratEst.return$strategies
    names_samples <- names(stratEst.return$strategies)
    for( i in 1:length(stratEst.return$strategies) ){
      strategies_sample <- do.call(rbind,strategies_print[[i]])
      row_names_strategies_sample <- rownames(strategies_sample)
      row.names(strategies_sample) =  paste( names_samples[i], "." , row_names_strategies_sample, sep="")
      strategies_sample_list <- rbind( strategies_sample_list , strategies_sample )
    }
    strategies_matrix_names = rownames(strategies_sample_list)
  }else{
    strategies_matrix_names = rownames(do.call(rbind,stratEst.return$strategies))
  }

  # state.obs post-processing
  obs_names <- rep(NA,num_unique_outputs)
  for( outs in 1:num_unique_outputs ){
    obs_names[outs] <- paste("obs.",as.character(unique_outputs[outs]),sep="")
  }

  rownames(stratEst.return$state.obs) <- strategies_matrix_names
  colnames(stratEst.return$state.obs) <- obs_names

  r_names <- rep(NA,num_unique_outputs)
  for( outs in 1:num_unique_outputs ){
    r_names[outs] <- paste("prob.",as.character(unique_outputs[outs]),sep="")
  }
  rownames(stratEst.return$probs) = strategies_matrix_names
  colnames(stratEst.return$probs) = r_names

  rownames(stratEst.return$trembles) = strategies_matrix_names
  colnames(stratEst.return$trembles) = "tremble"

  rownames(stratEst.return$gammas) = strategies_matrix_names
  colnames(stratEst.return$gammas) = "gamma"

  if( length( stratEst.return$shares.par ) > 0 ){
    row_names_shares = paste("par.",as.character(seq(1,length(stratEst.return$shares.par),by=1)),sep="")
    rownames(stratEst.return$shares.par) = row_names_shares
    rownames(stratEst.return$shares.se) = row_names_shares
    rownames(stratEst.return$shares.score) = row_names_shares
    rownames(stratEst.return$shares.covar) = row_names_shares
    rownames(stratEst.return$shares.fisher) = row_names_shares

    colnames(stratEst.return$shares.par) = "probability"
    colnames(stratEst.return$shares.se) = "standard error"
    colnames(stratEst.return$shares.score) = "score"
    colnames(stratEst.return$shares.covar) = row_names_shares
    colnames(stratEst.return$shares.fisher) = row_names_shares
  }
  if( length( stratEst.return$probs.par ) > 0 ){
    rownames(stratEst.return$probs.indices) = strategies_matrix_names
    colnames(stratEst.return$probs.indices) = r_names

    row_names_probs = paste("par.",as.character(seq(1,length(stratEst.return$probs.par),by=1)),sep="")
    rownames(stratEst.return$probs.par) = row_names_probs
    rownames(stratEst.return$probs.se) = row_names_probs
    rownames(stratEst.return$probs.score) = row_names_probs
    rownames(stratEst.return$probs.covar) = row_names_probs
    rownames(stratEst.return$probs.fisher) = row_names_probs

    colnames(stratEst.return$probs.par) = "probability"
    colnames(stratEst.return$probs.se) = "standard error"
    colnames(stratEst.return$probs.score) = "score"
    colnames(stratEst.return$probs.covar) = row_names_probs
    colnames(stratEst.return$probs.fisher) = row_names_probs
  }
  if( length( stratEst.return$trembles.par ) > 0 ){
    rownames(stratEst.return$trembles.indices) = strategies_matrix_names
    colnames(stratEst.return$trembles.indices) = "tremble"

    row_names_trembles = paste("par.",as.character(seq(1,length(stratEst.return$trembles.par),by=1)),sep="")
    rownames(stratEst.return$trembles.par) = row_names_trembles
    rownames(stratEst.return$trembles.se) = row_names_trembles
    rownames(stratEst.return$trembles.score) = row_names_trembles
    rownames(stratEst.return$trembles.covar) = row_names_trembles
    rownames(stratEst.return$trembles.fisher) = row_names_trembles

    rownames(stratEst.return$gammas.par) = row_names_trembles
    rownames(stratEst.return$gammas.se) = row_names_trembles

    colnames(stratEst.return$trembles.par) = "probability"
    colnames(stratEst.return$trembles.se) = "standard error"
    colnames(stratEst.return$trembles.score) = "score"
    colnames(stratEst.return$trembles.covar) = row_names_trembles
    colnames(stratEst.return$trembles.fisher) = row_names_trembles

    colnames(stratEst.return$gammas.par) = "gamma"
    colnames(stratEst.return$gammas.se) = "standard error"
  }
  if( length( stratEst.return$coefficients.par ) > 0 ){
    row_names_coefficients = paste("par.",as.character(seq(1,length(stratEst.return$coefficients.par),by=1)),sep="")
    rownames(stratEst.return$coefficients.par) = row_names_coefficients
    rownames(stratEst.return$coefficients.se) = row_names_coefficients
    rownames(stratEst.return$coefficients.score) = row_names_coefficients
    rownames(stratEst.return$coefficients.covar) = row_names_coefficients
    rownames(stratEst.return$coefficients.fisher) = row_names_coefficients

    colnames(stratEst.return$coefficients.par) = "coefficient"
    colnames(stratEst.return$coefficients.se) = "standard error"
    colnames(stratEst.return$coefficients.score) = "score"
    colnames(stratEst.return$coefficients.covar) = row_names_coefficients
    colnames(stratEst.return$coefficients.fisher) = row_names_coefficients
  }

  # crit values
  stratEst.return$aic = cpp.output$fit[1,4]
  stratEst.return$bic = cpp.output$fit[1,5]
  stratEst.return$icl = cpp.output$fit[1,6]
  stratEst.return$chi.global = cpp.output$fit[1,7]

  # model.entropy
  if( "list" %in% class(stratEst.return$shares) ){
    num.strats <- length(stratEst.return$shares[[1]])
  }
  else{
    num.strats <- length(stratEst.return$shares)
  }

  stratEst.return$entropy.model <- cpp.output$fit[1,3]/(stratEst.return$num.ids*log(num.strats))
  stratEst.return$entropy.assignments <- cpp.output$fit[1,3]
  stratEst.return$chi.local = matrix(cpp.output$fit[1,8:(8+num.strats-1)],1,num.strats)
  colnames(stratEst.return$chi.local) = names_of_strategies
  rownames(stratEst.return$chi.local) = "chi^2"

  # fit vector
  stratEst.return$fit <- matrix(c(stratEst.return$loglike,stratEst.return$free.par,stratEst.return$aic,stratEst.return$bic,stratEst.return$icl),1,5)
  colnames(stratEst.return$fit) <- c("loglike","free.par","aic","bic","icl")

  # delete empty list entries
  if( length(stratEst.return$shares.par) == 0 ){
    stratEst.return$shares.indices = NULL
  }
  if( length(stratEst.return$probs.par) == 0 ){
    stratEst.return$probs.indices = NULL
  }
  else{
    # if( all( stratEst.return$probs.par == 0 | stratEst.return$probs.par == 1 ) ){
    #   stratEst.return$probs.se = NULL
    #   stratEst.return$probs.score = NULL
    #   stratEst.return$probs.covar = NULL
    #   stratEst.return$probs.fisher = NULL
    # }
  }
  if( length(stratEst.return$trembles.par) == 0 ){
    stratEst.return$trembles.indices = NULL
  }

  if( LCR == F ){
    stratEst.return$coefficients = NULL
    stratEst.return$coefficients.quantiles = NULL
  }

  stratEst.return <- stratEst.return[lapply(stratEst.return,length)>0]

  # return result
  return(stratEst.return)
}

