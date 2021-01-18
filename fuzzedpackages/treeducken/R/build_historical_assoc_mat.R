#' Reconstruct historical association matrix
#'
#' @details Given a time and a tree pair object produced by the `sim_cophylo_bdp`
#'     object will produce the association matrix at that time point for the
#'     tree object.
#' USER WARNING: this is still in development, and likely will not work all the time.
#'
#' @param t The time of interest
#' @param tr_pair_obj The tree pair object from `sim_cophylo_bdp`
#' @return Matrix of the associations at given time
#' @examples
#' host_mu <- 1.0 # death rate
#' host_lambda <- 2.0 # birth rate
#' numb_replicates <- 1
#' time <- 1.0
#' symb_mu <- 0.2
#' symb_lambda <- 0.4
#' host_shift_rate <- 0.0
#' cosp_rate <- 2.0
#'
#' cophylo_pair <- sim_cophylo_bdp(hbr = host_lambda,
#'                            hdr = host_mu,
#'                            cosp_rate = cosp_rate,
#'                            host_exp_rate = host_shift_rate,
#'                            sdr = symb_mu,
#'                            sbr = symb_lambda,
#'                            numbsim = numb_replicates,
#'                            time_to_sim = time)
#' time <- 1.0
#' assoc_mat_at_t <- build_historical_association_matrix(t=time, tr_pair_obj = cophylo_pair[[1]])
#'
#'
build_historical_association_matrix <- function(t, tr_pair_obj){
    times <- unique(tr_pair_obj$event_history$Event_Time)

    ##Error checking with regards to 't'
    if(!is.numeric(t)){
      stop("'t' needs to be a number")
    }
    if(t<=0){
      stop("'t' needs to be positive")
    }
    if( t > (max(ape::branching.times(tr_pair_obj$host_tree)) + tr_pair_obj$host_tree$root.edge )){
      stop("The chosen 't' is beyond the duration of the cophylogeny. We don't know the future.")
    }
    if(t<times[2]){
      warning("you chose a 't' before the first event")
      return(1) ##this seems like a valid choice to me and I don't think it should error out but I don't know what it should return.
      ##it maybe should come with a wanring or message that this is a sort of dubious time that they picked
    }


    times <- times[times <= t] ###only consider times that are before or at 't'
    events <- tr_pair_obj$event_history
    curr_indx <- which.max(times)
    init_mat <- matrix(1, nrow = 1, ncol = 1)
    colnames(init_mat) <- as.character(
                            length(tr_pair_obj$host_tree$tip.label) + 1)
    rownames(init_mat) <- as.character(
                            length(tr_pair_obj$symb_tree$tip.label) + 1)
    prev_mat <- init_mat
    for(i in 2:curr_indx) {
        curr_events <- subset(events, events[,4] == times[i])
        # what is first row of curr_events
        curr_events
        main_event <- curr_events[1,]$Event_Type
        if(main_event == "HG") {
            # make new matrix with dimensions of prev_mat same matrix
            # remove column indicated by curr_events[1,]$Host.Index
            # make 2 new column vectors name them correctly based on
            # curr_events[2,]$Host.Index whichever rows correspond to 4 in matrx are 1
            # curr_events[3,]$Host.Index
            curr_mat <- matrix(0,
                               nrow = nrow(prev_mat),
                               ncol = ncol(prev_mat) + 1)
            rownames(curr_mat) <- as.character(rownames(prev_mat))
            prev_col_names <- colnames(prev_mat)
            prev_mat <- as.matrix(prev_mat[,-which(colnames(prev_mat) == as.character(curr_events[1,]$Host_Index))])
            if(length(prev_mat) > 0)
                curr_mat[seq_len(nrow(prev_mat)),
                         seq_len(ncol(prev_mat))] <- prev_mat

            prev_col_names_cut <- prev_col_names[which(prev_col_names != as.character(curr_events[1,]$Host_Index))]
            # if(length(prev_col_names_cut) > 0)
            #     colnames(prev_mat) <- prev_col_names_cut
            col_names <- c(prev_col_names_cut, curr_events[2:3,]$Host_Index)
            colnames(curr_mat) <- as.character(col_names)
            for(j in 1:nrow(curr_events)){
                type <- curr_events[j,]$Event_Type
                if(type == "AG"){
                    curr_mat[as.character(curr_events[j,]$Symbiont_Index),
                             as.character(curr_events[j,]$Host_Index)] <- 1
                }
            }
        }
        else if(main_event == "HL"){
            curr_mat <- matrix(0, nrow = nrow(prev_mat), ncol = ncol(prev_mat) - 1)
            rownames(curr_mat) <- rownames(prev_mat)
            prev_col_names <- colnames(prev_mat)
            prev_mat <- as.matrix(prev_mat[,-which(colnames(prev_mat) == as.character(curr_events[1,]$Host_Index))])
            prev_col_names_cut <- prev_col_names[which(prev_col_names != as.character(curr_events[1,]$Host_Index))]
            colnames(prev_mat) <- prev_col_names_cut
            curr_mat <- prev_mat
            colnames(curr_mat) <- as.character(prev_col_names_cut)
        }
        else if(main_event == "SG"){
            curr_mat <- matrix(0, nrow = nrow(prev_mat) + 1, ncol = ncol(prev_mat))
            colnames(curr_mat) <- as.character(colnames(prev_mat))
            prev_row_names <- rownames(prev_mat)
            prev_mat <- as.matrix(prev_mat[-which(rownames(prev_mat) == as.character(curr_events[1,]$Symbiont_Index)),], )
            if(length(prev_mat) > 0)
                curr_mat[seq_len(nrow(prev_mat)), seq_len(ncol(prev_mat))] <- prev_mat

            prev_row_names_cut <- prev_row_names[which(prev_row_names != as.character(curr_events[1,]$Symbiont_Index))]
            # if(length(prev_row_names_cut) > 0)
            #     rownames(prev_mat) <- prev_row_names_cut
            row_names <- c(prev_row_names_cut, curr_events[2:3,]$Symbiont_Index)
            rownames(curr_mat) <- as.character(row_names)
            for(j in 2:nrow(curr_events)){
                type <- curr_events[j,]$Event_Type
                if(type == "AG"){
                    curr_mat[as.character(curr_events[j,]$Symbiont_Index),
                             as.character(curr_events[j,]$Host_Index)] <- 1
                }
            }
        }
        else if(main_event == "SL"){
            curr_mat <- matrix(0, nrow = nrow(prev_mat) - 1, ncol = ncol(prev_mat))
            colnames(curr_mat) <- colnames(prev_mat)
            prev_row_names <- rownames(prev_mat)
            prev_mat <- as.matrix(prev_mat[-which(rownames(prev_mat) == as.character(curr_events[1,]$Symbiont_Index)),])
            prev_row_names_cut <- prev_row_names[which(prev_row_names != as.character(curr_events[1,]$Symbiont_Index))]
            rownames(prev_mat) <- rownames(prev_row_names_cut)
            curr_mat <- prev_mat
            rownames(curr_mat) <- as.character(prev_row_names_cut)
        }


        else if(main_event %in% c("AL","AG")){ ##We have a series of association losses and gains not tied to a specific event
          curr_mat<-prev_mat ##we don't change the size of the matrix at all


          for(j in 1:nrow(curr_events)){ ##Go thru all loseses and gains
            ev<-curr_events[j,]
            if(ev$Event_Type== "AG"){
              ## If the event is a gain we puta 1 at their spot in the matrix
              curr_mat[ev$Symbiont_Index,ev$Host_Index]<-1
            }else if(ev$Event_Type=="AL"){
              ##If the event is a loss we put a 0 at their spot in the matrix
              curr_mat[ev$Symbiont_Index,ev$Host_Index]<-0
            }
          }
        }

        else{
            curr_mat <- matrix(0, nrow = nrow(prev_mat) + 1, ncol = ncol(prev_mat) + 1)
            curr_mat[nrow(prev_mat), ncol(prev_mat)] <- 1
            curr_mat[nrow(prev_mat) + 1, ncol(prev_mat) + 1] <- 1
            prev_col_names <- colnames(prev_mat)
            prev_row_names <- rownames(prev_mat)

            prev_mat <- as.matrix(prev_mat[-which(rownames(prev_mat) == as.character(curr_events[1,]$Symbiont_Index)),
                                           -which(colnames(prev_mat) == as.character(curr_events[1,]$Host_Index))])

            prev_col_names_cut <- prev_col_names[which(prev_col_names != as.character(curr_events[1,]$Host_Index))]
            prev_row_names_cut <- prev_row_names[which(prev_row_names != as.character(curr_events[1,]$Symbiont_Index))]

            if(length(prev_mat) > 0)
                curr_mat[1:nrow(prev_mat), 1:ncol(prev_mat)] <- prev_mat

            # if(length(prev_col_names_cut) > 0)
            #     colnames(prev_mat) <- prev_col_names_cut
            # if(length(prev_row_names_cut) > 0)
            #     rownames(prev_mat) <- prev_row_names_cut
            col_names <- c(prev_col_names_cut, curr_events[2:3,]$Host_Index)
            row_names <- c(prev_row_names_cut, curr_events[2:3,]$Symbiont_Index)
            rownames(curr_mat) <- as.character(row_names)
            colnames(curr_mat) <- as.character(col_names)
            for(j in 2:nrow(curr_events)){
                type <- curr_events[j,]$Event_Type
                if(type == "AG"){
                    curr_mat[as.character(curr_events[j,]$Symbiont_Index),
                             as.character(curr_events[j,]$Host_Index)] <- 1
                }
            }
        }
        # if C delete both, add 2 new columns and rows
        # else if HG delete col, add 2 new columns
        # else if SG delete col, add 2 new columns
        # else if SL delete col
        # else if HL delete col
        if(nrow(curr_mat) || ncol(curr_mat))
          return(curr_mat)
        prev_mat <- curr_mat
    }

    ##if we chose a time past all our events then the built matrix should be the same as what we have stored in our $association_mat
    if(t > max(times)){
      # if(!all(tr_pair_obj$association_mat==prev_mat)){
      #   warning(paste("with a chosen time of ", t, "the built matrix should be equal to tr_pair_ob$association_mat but it is not"))
      # }
      return(tr_pair_obj$association_mat)
    }

    prev_mat
}
