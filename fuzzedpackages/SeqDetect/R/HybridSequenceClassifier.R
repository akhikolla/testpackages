HybridSequenceClassifier <- setRefClass("HybridSequenceClassifier",
                                        fields = list(
                                          pp = "ANY",
                                          pc = "ANY",
                                          fds = "vector",
                                          ts_f = "character",
                                          tf_f = "character",
                                          ctx_f = "ANY",
                                          ew = "ANY",
                                          pattern_f = "ANY",
                                          cache = "ANY",
                                          tsss = "logical"
                                        ),
                                        methods = list(
                                          initialize = function(fields, timestamp_start_field, timestamp_finish_field,
                                                                context_field=NULL, preclassifier=NULL, preprocessor=NULL,
                                                                decay_descriptors=NULL, pattern_field=NULL,
                                                                time_series_sequence_stats=FALSE, reuse_states=TRUE,
                                                                parallel_execution=FALSE) {
                                            setInputDefinitions(fields, timestamp_start_field, timestamp_finish_field,
                                                                context_field, preclassifier, preprocessor, pattern_field)
                                            ew <<- create_ETT_wrapper(decay_descriptors,reuse_states,parallel_execution,time_series_sequence_stats);
                                            cache <<- NULL
                                            tsss <<- time_series_sequence_stats
                                          }
                                        )
)

HybridSequenceClassifier$methods(list(
    process = function(streams,learn=TRUE,give_explain=TRUE,threshold=NULL,debug=FALSE,out_filename=NULL) {
      deserialize()
      pp_res <- preprocess(pp,streams)
      pp <<- pp_res$obj
      event_stream <- pp_res$res
      cons_stream <- classify(pc,event_stream)
      explDataList <- list()
      sequences <- list()
      if(is.null(threshold)) threshold=-1
      key_f <- ctx_f
      if(".key" %in% names(cons_stream)) key_f <- ".key"
      if(is.null(key_f)) stop("Key field is not defined");
      res <- ew$process(cons_stream,".clazz",key_f,ts_f,tf_f,pattern_f,out_filename,NULL,debug,threshold,learn)
      for(j in 1:length(res)) {
        if(tsss && res[[j]]$from>=0 && res[[j]]$to>=0) res[[j]]$seq_stats <- TRUE
        else res[[j]]$seq_stats <- FALSE
        if(!"actual" %in% names(res[[j]])) res[[j]]$actual <- c()
        if(!"potential" %in% names(res[[j]])) res[[j]]$potential <- c()
        res[[j]]$.clazz <- as.character(cons_stream[j,".clazz"])
        res[[j]]$.key <- as.character(cons_stream[j,key_f])
      }
      if(tsss) sequences <- .updateSequences(cons_stream,res,sequences,learn)$sequences
      resOut <- list(stream=cons_stream)
      if(give_explain) {
        resOut$explanation <- res
        if(tsss) resOut$sequences <- sequences
      }
      return(resOut)
    },
    printMachines = function(machine_id=NULL,state=NULL,print_cache=TRUE,print_keys=TRUE) {
      deserialize()
      ew$printMachines(machine_id,state,print_cache,print_keys)
    },
    getMachineIdentifiers = function() {
      return(ew$getMachineIdentifiers())
    },
    setOutputPattern = function(states=c(),transitions=c(),pattern,machine_id=NULL) {
      deserialize()
      if(!is.null(machine_id)) {
        for(state in states) ew$setStatePattern(machine_id,state,pattern)
        for(trans in transitions) ew$setTransitionPattern(machine_id,trans,pattern)
      } else {
        message(getMachineIdentifiers())
        for(machine_id in getMachineIdentifiers()) {
          for(state in states) ew$setStatePattern(machine_id,state,pattern)
          for(trans in transitions) ew$setTransitionPattern(machine_id,trans,pattern)
        }
      }
    },
    setPreprocessor = function(preprocessor) {
      pp <<- preprocessor
    },
    setPreclassifier = function(preclassifier) {
      pc <<- preclassifier
    },
    setInputDefinitions = function(fields, timestamp_start_field, timestamp_finish_field, context_field=NULL, 
                                   preclassifier=NULL, preprocessor=NULL, pattern_field=NULL) {
        if(is.null(preclassifier)) preclassifier <- HSC_PC_None()
        if(!is(preclassifier,"HSC_PC"))
          stop("HSC: preclassifier must be an object of class that inherits HSC_PC")
        pc <<- preclassifier
        if(!is.vector(fields) || length(fields)==0)
          stop("HSC: fields list must be of type vector, and contain at least one element")
        fds <<- fields
        if(!timestamp_start_field %in% fields)
          stop("HSC: starting timestamp field is not found in the fields vector")
        ts_f <<- timestamp_start_field
        if(!timestamp_finish_field %in% fields)
          stop("HSC: finish timestamp field is not found in the fields list")
        tf_f <<- timestamp_finish_field
        if(!is.null(context_field) && !context_field %in% fields)
          stop("HSC: context field not found in the fields vector")
        ctx_f <<- context_field
        if(!is.null(preprocessor)) {
          if(!is(preprocessor,"HSC_PP"))
            stop("HSC: preprocessor must be an object of class that inherits HSC_PP")
        } else preprocessor <- HSC_PP(fds, ts_f)
        pp <<- preprocessor
        if (!is.null(pattern_field) && !is.character(pattern_field))
          stop("HSC: pattern field must of string type")
        pattern_f <<- pattern_field
    },
    induceSubmachine = function(threshold=10,isolate=FALSE) {
      deserialize()
      return(ew$induceSubmachine(threshold,isolate))
    },
    plotMachines = function(machine_id=NULL,label_aspect=1.0) {
      deserialize()
      if(is.null(machine_id)) {
        for(mid in ew$getMachineIdentifiers()) {
          x <- ew$getCoincidenceValues(mid,TRUE)
          if(length(x$coincidence)>0) {
            g <- graph_from_data_frame(x$coincidence)
            plot.igraph(g,edge.label=x$coincidence[,"weight"],vertex.color=NA,
                        vertex.frame.color="gray",vertex.label.font=2,vertex.label.cex=label_aspect,edge.label.cex=label_aspect)
          } else warning(paste0("No plot for the machine ",mid))
        }
      } else {
        x <- ew$getCoincidenceValues(machine_id,TRUE)
        if(length(c$coincidence)>0) {
          g <- graph_from_data_frame(x$coincidence)
          plot.igraph(g,edge.label=x$coincidence[,"weight"],vertex.color=NA,
                      vertex.frame.color="gray",vertex.label.font=2,vertex.label.cex=label_aspect,edge.label.cex=label_aspect)
        } else warning(paste0("No plot for the machine ",machine_id))
      }
    },
    cleanKeys = function(machine_id=NULL) {
      deserialize()
      ew$cleanMachineKeys(machine_id=machine_id)
    },
    mergeMachines = function() {
      deserialize()
      ew$mergeAllMachines()
    },
    compressMachines = function(ratio=0.5) {
      deserialize()
      ew$compressMachines(ratio)
    },
    clone = function() {
      deserialize()
      hsc <- HybridSequenceClassifier(fds,ts_f,tf_f,ctx_f,pc,preprocessor=pp,pattern_field=pattern_f,
                                      time_series_sequence_stats=tsss)
      hsc$ew <- ew$clone()
      return(hsc)
    },
    .updateSequences = function(consDataFrame,explDataList,sequences,learn) {
      #.GlobalEnv$cdf <- consDataFrame
      #.GlobalEnv$edl <- explDataList
      for(i in 1:nrow(consDataFrame)) {
        element <- consDataFrame[i,]
        ett <- explDataList[[i]]
        if(!learn) {
          clazz <- element[,".clazz",drop=TRUE]
          key <- ett$.key
          pattern <- element[,pattern_f,drop=TRUE]
          if(!key %in% names(sequences)) sequences[[key]] <- list(full_chain=c(),chains=0,total_chain_length=0)
          seq <- sequences[[key]]
          seq$full_chain[length(seq$full_chain)+1] <- pattern
          active_ttok <- c()
          if(ett$seq_stats) {
            stok <- as.character(ett$from)
            ttok <- as.character(ett$to)
            if(stok %in% names(seq) && seq[[stok]]$active) {
              chain <- seq[[stok]]
              chain$sequence[[length(chain$sequence)+1]] <- clazz
              chain$end <- element[,ts_f,drop=TRUE]
              seq[[ttok]] <- chain
              seq[[stok]] <- NULL
              seq$total_chain_length <- seq$total_chain_length+1
              active_ttok <- c(active_ttok,ttok)
            } else if (!ttok %in% names(seq)) {
              chain <- list(active=TRUE,sequence=c(),start=seq$temp_start,end=element[,tf_f,drop=TRUE])
              for(sp in ett$from_patterns) chain$sequence[[length(chain$sequence)+1]] <- sp
              for(tp in ett$to_patterns) chain$sequence[[length(chain$sequence)+1]] <- tp
              seq[[ttok]] <- chain
              seq$chains <- seq$chains+1
              seq$total_chain_length <- seq$total_chain_length+1
              seq$temp_start <- NULL
              active_ttok <- c(active_ttok,ttok)
            }
          } else seq$temp_start <- element[,ts_f,drop=TRUE]
          for(ttok in base::setdiff(names(seq),c(active_ttok,"full_chain","total_chain_length","chains","temp_start"))) seq[[ttok]]$active <- FALSE
          sequences[[key]] <- seq
        }
      }
      return(list(sequences=sequences))
    },
    
    serialize = function() {
      if(is.null(cache) && !is.null(ew)) {
        cache <<- ew$serialize_ETT()
        ew <<- NULL
      }
    },
    serializeToList = function() {
      if(is.null(cache) && !is.null(ew)) {
        res <- ew$serialize_ETT()
        res$hsc_options <- list(preprocessor=pp,preclassifier=pc,fields=fds,time_start_field=ts_f,time_end_field=tf_f,
                                context_field=ctx_f,pattern_field=pattern_f,time_series_sequence_stats=tsss)
        return(res)
      } else stop("Cannot serialize this HSC");
    },
    
    deserialize = function() {
      if(is.null(ew) && !is.null(cache)) {
        ew <<- deserialize_ETT(cache)
        cache <<- NULL
      }
    }
))


deserializeFromList <- function(l) {
  if(!"hsc_options" %in% names(l)) stop("List does not contain a seralized HSC")
  hsco <- l$hsc_options;
  l$hsc_options <- NULL
  hsc <- HybridSequenceClassifier(hsco$fields,hsco$time_start_field,hsco$time_end_field,hsco$context_field,
                                  hsco$preclassifier,preprocessor=hsco$preprocessor,pattern_field=hsco$pattern_field,
                                  time_series_sequence_stats=hsco$time_series_sequence_stats)
  hsc$ew <- deserialize_ETT(l)
  hsc$cache <- NULL
  return(hsc)
}