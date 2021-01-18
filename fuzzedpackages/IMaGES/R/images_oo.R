#

IMaGESclass <- setRefClass("IMaGESclass",
                      fields = list(
                        matrices="list",
                        penalty="numeric",
                        .rawscores="list",
                        .graphs="list",
                        imscore = "numeric",
                        results="list",
                        scores = "list",
                        num.markovs = "numeric",
                        use.verbose = "logical"),
                   
                      methods = list(
                        
                        ## Purpose: cast input data into graph object for use with IMaGES algorithm
                        ## ----------------------------------------------------------------------
                        ##' @param score 	scoring object to be used
                        ##' @param labels 	node labels
                        ##' @param fixedGaps 	logical matrix indicating forbidden edges
                        ##' @param adaptive sets the behaviour for adaptiveness in the forward phase (cf. "ARGES")
                        ##' @param phase  lists the phases that should be executed
                        ##' @param iterate  indicates whether the phases should be iterated. iterated = FALSE
                        ##'   means that the required phases are run just once
                        ##' @param turning	indicates whether the turning step should be included (DEPRECATED).
                        ##' @param maxDegree 	maximum vertex degree allowed
                        ##' @param verbose 	indicates whether debug output should be printed
                        ##' @param ... 		additional parameters (currently none)
                        ##' @param targets 	unique list of targets. Normally determined from the scoring object
                        ## ----------------------------------------------------------------------
                        ## return: IMGraph object
                        ## Author: Markus Kalisch, Date: 26 Jan 2006;  Martin Maechler; Modified by Noah Frazier-Logue
                        create.graph = function(
                          score, 
                          labels = score$getNodes(), 
                          targets = score$getTargets(),
                          fixedGaps = NULL, 
                          #adaptive = c("none", "vstructures", "triples"), 
                          #phase = c("forward", "backward", "turning"),
                          #iterate = length(phase) > 1,
                          turning = NULL, 
                          maxDegree = integer(0),
                          verbose = FALSE, 
                          ...)
                        {
                          # Catch calling convention of previous package versions:
                          # ges(p, targets, score, fixedGaps = NULL, ...)
                          # If this calling convention is used, issue a warning, but adjust the 
                          # arguments
                          if (is.numeric(score) && is.list(labels) && inherits(targets, "Score")) {
                            score <- targets
                            targets <- labels
                            labels <- as.character(1:length(score$getNodes()))
                            warning(paste("You are using a deprecated calling convention for gies()",
                                          "which will be disabled in future versions of the package;",
                                          "cf. ?gies.", sep = " "))
                          }

                          if (!missing(turning)) {
                            stopifnot(is.logical(turning))
                            warning(paste0("The argument 'turning' is deprecated; please use 'phase'",
                                           "instead (cf. ?ges)"))
                            
                          }
                          
                          # Error checks
                          if (!inherits(score, "Score")) {
                            #print(score)
                            stop("Argument 'score' must be an instance of a class inherited from 'Score'.")
                          }

                          if (is.numeric(score)) {
                            # This happens when the old calling convention is used with all 
                            # mandatory arguments unnamed
                            p <- score
                            if (is.list(labels) && is(targets, "Score")) {
                              score <- targets
                              targets <- labels
                              labels <- as.character(1:p)
                              warning(paste("You are using a DEPRECATED calling convention for",
                                            "gies(), gds() or simy(); please refer to the documentation",
                                            "of these functions to adapt to the new calling conventions."))
                            } else if (is(labels, "Score")) {
                              score <- labels
                              labels <- as.character(1:p)
                              warning(paste("You are using a DEPRECATED calling convention for",
                                            "ges(); please refer to the documentation",
                                            "to adapt to the new calling convention."))
                            }
                          } else if (is.numeric(labels) && length(labels) == 1) {
                            # This happens when the old calling convention is used with only the
                            # 'score' argument named
                            labels <- as.character(1:labels)
                            warning(paste("You are using a DEPRECATED calling convention for",
                                          "gies(), ges(), gds() or simy(); please refer to the documentation",
                                          "of these functions to adapt to the new calling conventions."))
                          }
                          
                          if (!is(score, "Score")) {
                            stop("'score' must be of a class inherited from the class 'Score'.")
                          }
                          if (!is.character(labels)) {
                            stop("'labels' must be a character vector.")
                          }
                          if (!is.list(targets) || !all(sapply(targets, is.numeric))) {
                            stop("'targets' must be a list of integer vectors.")
                          }
                          
                          
                          #print("You made it this far")
                          imgraph <- new("IMGraph", nodes = labels, targets = targets, score = score)
                          return(imgraph)
                        },
                        
                        
                        ## Purpose: run a particular GES phase (forward, backward, turning) on a 
                        ##          given graph
                        ## ----------------------------------------------------------------------
                        ##' @param phase GES phase to run on a given graph
                        ##' @param j index of graph to run supplied GES phase on
                        ## ----------------------------------------------------------------------
                        ##
                        ## Author: Noah Frazier-Logue
                        run_phase = function(phase="GIES-F", j) {
                          .graphs[[j]]$greedy.step(alg.name=phase, direction = phase, verbose = FALSE)
                        },
                        
                        ## Purpose: find optimal IMaGES phase to run by calculating the mode of
                        ##          the phases list found
                        ## NOTE: this function is no longer used
                        ## ----------------------------------------------------------------------
                        ##' @param opt_phases list of optimal phases selected for each individual
                        ##'                   graph
                        ## ----------------------------------------------------------------------
                        ##' @return: value of mode
                        ## Author: Noah Frazier-Logue
                        find_opt = function(opt_phases) {
                          phases <- list()
                          phase_counts = list()
                          for (i in 1:length(opt_phases)) {
                            if ((length(phases) == 0) || !(is.element(opt_phases[[i]], phases))) {
                              phases[[length(phases) + 1]] <- opt_phases[[i]]
                              phase_counts[[length(phase_counts) + 1]] = 1
                            }
                            else {
                              phase_counts[[match(opt_phases[[i]], phases)]] = phase_counts[[match(opt_phases[[i]], phases)]] + 1
                            }
                            
                          }
                          max_val = max(unlist(phase_counts))
                          index = match(max_val, phase_counts)

                          
                          return(phases[[index]])
                        },
                        
                        ## Purpose: update global IMScore
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        update_score = function() {
                          
                          if (use.verbose) {
                            print("Updating IMScore...")
                          }
                          
                          im.score <- IMScore()
                          #assign("score", imscore, envir=trueIM)
                          
                          if (use.verbose) {
                            print(paste("Updated IMScore: ", im.score))
                          }
                          trueIM$score <- im.score
                        },
                        
                        ## Purpose: run a particular GES phase (forward, backward, turning) on a 
                        ##          given graph
                        ## ----------------------------------------------------------------------
                        ##' @param x list to find mode
                        ## ----------------------------------------------------------------------
                        ##' @return: index of mode
                        ## Author: Noah Frazier-Logue
                        mode = function(x) {
                          ux <- unique(x)
                          #ux[which.max(tabulate(match(x, ux)))]
                          tab <- tabulate(match(x, ux)); ux[tab == max(tab)]
                        },
                        
                        ## Purpose: runs one phase of IMaGES on all the graphs. Works by 
                        ##          determining the optimal step (forward, backward, turning) for
                        ##          each graph and picking the step that most increases the IMscore
                        ##          across all graphs
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        run = function() {
                          
                          update_score()
                          
                          opt_phases = list()
                          
                          if (use.verbose) {
                            print("Selecting optimal step...")
                          }
                          
                          #call C++ function that determines best step for each graph
                          for (j in 1:length(.graphs)) {
                            opt_phases[[j]] <- .Call("greedyStepRFunc",
                                                     .graphs[[j]]$.in.edges,
                                                     .graphs[[j]]$.score$pp.dat,
                                                     .graphs[[j]]$.score$c.fcn,
                                                     .graphs[[j]]$causal.inf.options(caching = FALSE, maxSteps = 1),
                                                     PACKAGE = "IMaGES")
                          }

                          #find mode of optimal phase list
                          opt <- mode(opt_phases)[[1]]
                          
                          str_opt <- ''
                          
                          #map numeric value to text value for use in
                          #run_phase
                          if (opt == 1) {
                            str_opt <- 'GIES-F'
                            if (use.verbose) {
                              print("Forward pass selected")
                            }
                          }
                          else if (opt == 2) {
                            str_opt <- 'GIES-B'
                            if (use.verbose) {
                              print("Backward pass selected")
                            }
                          }
                          else if (opt == 3) {
                            str_opt <- 'GIES-T'
                            if (use.verbose) {
                              print("Turning pass selected")
                            }
                          }
                          else if (opt == 0) {
                            str_opt <- 'none'
                            if (use.verbose) {
                              print("No step selected")
                            }
                          }
                          
                          temp.scores <- vector()
                          if (!(str_opt == "none")) {
                            for (j in 1:length(.graphs)) {
                              #run single phase on graph j
                              run_phase(phase=str_opt, j)
                              #save score value generated by individual graph change
                              temp.scores[[j]] <- IMScore()
                              #undo step so next iteration has clean comparison
                              .graphs[[j]]$undo.step()
                            }
                            
                            #re-enable all steps after seeing which best impacts IMScore
                            for (j in 1:length(.graphs)) {
                              .graphs[[j]]$redo.step()
                              #see if graph j is in MEC
                              if (length(.graphs[[j]]$.in.edges) > 0) {
                                
                                if (use.verbose) {
                                  print("Updating markov equivalence class")
                                }
                                update.markovs(.graphs[[j]], temp.scores[[j]])
                              }
                            }
                            #find.best.step(temp.scores)
                          }
                        },
                        
                        
                        ## Purpose: finds ideal step/edge to insert into global structure
                        ## Note: this funciton is no longer used
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        find.best.step = function(score.list) {
                          # for (i in 1:length(score.list)) {
                          #   .graphs[[i]]$undo.step()
                          # }
                          inf <- 0
                          for (i in 1:length(score.list)) {
                            if (is.infinite(score.list[[i]])) {
                              inf <- inf + 1
                            }
                          }
                          if (inf == length(score.list)) {
                            return()
                          }
                          
                          #find index of the graph that had greatest positive effect
                          #on IMScore
                          #Note: can return two values if two edge changes have an identical
                          #      impact on the IMScore
                          best.graph.index <- which(score.list == max(score.list))


                          #update IMScore with all changes
                          update_score()
                          
                          #just use first index of best.graph.index for now
                          #TODO: find smarter way to do this
                          if (!update.global(.graphs[[best.graph.index[[1]]]]$.edge.change)) {
                            #make the score lower than the rest so it considers the next highest score
                            #better way to do this? probably
                            score.list[[best.graph.index[[1]]]] <- -Inf
                            find.best.step(score.list)
                          }
                          
                        },
                        
                        
                        ## Purpose: insert edge into global graph, 
                        ##          where dst contains the list of edges going towards it
                        ## ----------------------------------------------------------------------
                        ##' @param src source vertex, represented as an int
                        ##' @param dst destination vertex, represented as an int
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        insert.global = function(src, dst) {
                          insert <- src
                          insert.point <- which(order(c(insert, trueIM$global.edges[[dst]]))==1)
                          #insert the edge into the edgelist for that vertex
                          trueIM$global.edges[[dst]] <- append(trueIM$global.edges[[dst]], insert, insert.point - 1)
                        },
                        
                        ## Purpose: remove edge from global graph, where *dst* contains the list of edges going towards it
                        ##          where dst contains the list of edges going towards it
                        ## ----------------------------------------------------------------------
                        ##' @param src source vertex, represented as an int
                        ##' @param dst destination vertex, represented as an int
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        remove.global = function(src, dst) {
                          #remove edge by reassigning edge list to itself, where none of the values are *src*
                          trueIM$global.edges[[dst]] <- trueIM$global.edges[[dst]][trueIM$global.edges[[dst]] != src]
                        },
                        
                        ## Purpose: perform deletion of (src,dst) edge and insertion of (dst,src)
                        ##          edge to simulate "turning" mode for global graph
                        ## ----------------------------------------------------------------------
                        ##' @param src source vertex, represented as an int
                        ##' @param dst destination vertex, represented as an int
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        turn.global = function(src, dst) {
                          remove.global(src, dst)
                          insert.global(dst,src)
                        },
                        
                        ## Purpose: determine whether or not an edge exists in the global graph
                        ##          given the src and dst vertices
                        ## ----------------------------------------------------------------------
                        ##' @param src source vertex, represented as an int
                        ##' @param dst destination vertex, represented as an int
                        ##' @return TRUE if edge exists in global and FALSE otherwise
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        edge.exists = function(src, dst) {
                          return(src %in% trueIM$global.edges[[dst]])
                        },
                        
                        ## Purpose: determine whether or not an edge is legal to insert
                        ##          given the src and dst vertices
                        ## ----------------------------------------------------------------------
                        ##' @param src source vertex, represented as an int
                        ##' @param dst destination vertex, represented as an int
                        ##' @return TRUE if edge is able to be inserted and FALSE otherwise
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        is.legal.edge = function(src, dst) {
                          if (src > 0 && src <= ncol(.graphs[[1]]$.score$pp.dat$data)) {
                            if (dst > 0 && dst <= ncol(.graphs[[1]]$.score$pp.dat$data)) {
                              return(TRUE)
                            }
                          }
                          return(FALSE)
                        },
                        
                        
                        ## Purpose: handles updating of global graph. calls insert.global, remove.global,
                        ##          or turn.global depending on what the edge change specifies
                        ##
                        ## Note: this function is no longer used
                        ## ----------------------------------------------------------------------
                        ##' @param edge.change list containing src vertex, dst vertex, and step
                        ##'                    to run on that edge
                        ##' @return TRUE if global structure is able to be updated and FALSE otherwise
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        update.global = function(edge.change) {
                          src <- edge.change[[1]]
                          dst <- edge.change[[2]]
                          dir <- edge.change[[3]]
                          if (dir == 'GIES-F') {
                            #insert
                            if (is.legal.edge(src,dst) && !(edge.exists(src, dst))) {
                              #print("Inserting edge")
                              insert.global(src, dst)
                              return(TRUE)
                            }
                            else {
                              return(FALSE)
                            }
                          }
                          else if (dir == 'GIES-B') {
                            #remove
                            if (is.legal.edge(src,dst) && edge.exists(src, dst)) {
                              #print("Removing edge")
                              remove.global(src, dst)
                              return(TRUE)
                            }
                            else {
                              return(FALSE)
                            }
                          }
                          else if (dir == 'GIES-T') {
                            #turn
                            #print("Made it to TURN")
                            #invisible(readline(prompt="Press [enter] to continue"))
                            if (is.legal.edge(src,dst) && edge.exists(src, dst)) {
                              #print("Turning edge")
                              turn.global(src, dst)
                              return(TRUE)
                            }
                            else {
                              #print("Something messed up here")
                              #invisible(readline(prompt="Press [enter] to continue"))
                              return(FALSE)
                              
                            }
                          }
                          else {
                            return(TRUE)
                          }
                        },
                        
                        ## Purpose: initializes global graph prior to IMaGES run
                        ##
                        ## Note: this function is no longer used
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        initialize.global = function() {
                          edges <- list()
                          
                          for (i in 1:ncol(.graphs[[1]]$.score$pp.dat$data)) {
                            edges[[i]] <- vector()
                          }
                          
                          #assign("global.edges", edges, env=trueIM)
                          trueIM$global.edges <- edges
                        },
                        
                        
                        ## Purpose: calculates IMScore, the global replacement for the local score.
                        ##          Equation from Six Problems paper used.
                        ##
                        ## ----------------------------------------------------------------------
                        ##' @return IMScore
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        IMScore = function() {
                          penalty <<- penalty
                          m <- length(.graphs)
                          n <- ncol(.graphs[[1]]$.score$pp.dat$data) * nrow(.graphs[[1]]$.score$pp.dat$data)
                          sum <- 0
                          k <- nrow(.graphs[[1]]$.score$pp.dat$data)
                          trueIM$isLocalIM <- FALSE

                          for (i in 1:length(.graphs)) {
                            sum <- sum + .graphs[[i]]$.score$global.score.int(.graphs[[i]]$.in.edges)
                          }
                          imscore <- ((2/m) *  sum) + ((penalty * k) * log(n))
                          return(imscore)
                        },
                        
                        ## Purpose: fixes edge directions so there are no double directed edges.
                        ##          
                        ## TODO: improve this approach
                        ## ----------------------------------------------------------------------
                        ##' @param graph graph object to be modified
                        ##' @return modified graph object
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        fix.edges = function(graph) {
                          new.graph <- rep(list(c()), length(graph))
                          for (i in 1:length(graph)) {
                            #if there are edges going to this vertex
                            if (length(graph[[i]]) > 0) {
                              #for each edge going to vertex i
                              for (j in (1:length(graph[[i]]))) {
                                if (!((i %in% graph[[graph[[i]][[j]]]]) && (i < graph[[i]][[j]]))) {
                                  #adding edges that only point forward -- hacky fix for now
                                  new.graph[[i]] <- append(new.graph[[i]], graph[[i]][[j]])
                                }
                              }
                            }
                          }
                          new.graph
                        },
                        
                        ## Purpose: compares two graphs and returns TRUE if the graphs are 
                        ##          identical and FALSE otherwise
                        ## ----------------------------------------------------------------------
                        ##' @param graph1 first graph to be compared
                        ##' @param graph2 second graph to be compared
                        ##' @return TRUE if the graphs are identical and FALSE otherwise
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        is.identical = function(graph1, graph2) {
                          
                          ####
                          #this implementation, while very compact, was extremely slow. my 
                          #spaghetti implementation, while spaghetti, ends up making the program
                          #~3-12x faster (scales to the size of the input)
                          ####
                          # graph1 <- igraph::as_adjacency_matrix(igraph::igraph.from.graphNEL(convert(graph1)), names=FALSE)
                          # graph2 <- igraph::as_adjacency_matrix(igraph::igraph.from.graphNEL(convert(graph2)), names=FALSE)
                          # 
                          # truth.matrix <- graph1 == graph2
                          # 
                          # if (length(truth.matrix[truth.matrix == FALSE]) > 0) {
                          #   return(FALSE)
                          # }
                          # return(TRUE)
                          

                          for (i in 1:length(graph1)) {
                            #vertex i in graph1 has edges but vertex i in 
                            #graph 2 doesn't
                            if ((length(graph1[[i]]) > 0)) {
                              if (length(graph2[[i]]) == 0) {
                                return(FALSE)
                              }
                            }
                            
                            #vertex i in graph1 has no edges but vertex i in 
                            #graph 2 does
                            if ((length(graph2[[i]]) > 0)) {
                              if (length(graph1[[i]]) == 0) {
                                return(FALSE)
                              }
                            }
                            
                            #if they actually both have edges
                            if ((length(graph1[[i]]) > 0) && (length(graph2[[i]]) > 0)) {
                              #if they have differing numbers of edges
                              if (length(graph1[[i]]) != length(graph2[[i]])) {
                                return(FALSE)
                              }
                              
                              #list containing bools from element-wise comparison between
                              #graph1[[i]] and graph2[[i]]
                              comp.bools <- graph1[[i]] == graph2[[i]]
                              if (!all(comp.bools)) {
                                return(FALSE)
                              }
                            }
                          }
                          return(TRUE)
                        },
                        
                        ## Purpose: update the Markov Equivalence class throughout the IMaGES
                        ##          process with the top n graphs, where n is the user-specified
                        ##          max size of the MEC
                        ## ----------------------------------------------------------------------
                        ##' @param graph graph that might be added to the MEC
                        ##' @param score IMScore of the graph
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        update.markovs = function(graph, score) {
                          #iterate backwards to find lowest score that it's higher than
                          sum <- 0
                          graph1 <- list(.in.edges=graph$.in.edges, .nodes=graph$.nodes)
                          for (i in 1:length(trueIM$markovs)) {
                            if (!is.null(trueIM$markovs[[i]]$.graph)) {
                              true <- list(.in.edges=trueIM$markovs[[i]]$.graph, .nodes=graph$.nodes)
                              #don't include duplicate graphs 
                              if (is.identical(graph1$.in.edges, true$.in.edges)) {
                                rm(true)
                                if (use.verbose) {
                                  print("No new MEC found") 
                                }
                                return()
                              }
                            }
                          }
                          
                          for (i in 1:length(trueIM$markovs)) {
                            # if it is higher, bump lowest markov and sort in the new one
                            #exclude graphs that are identical since that's not too helpful
                            res <- (score > trueIM$markovs[[i]]$.score)
                            if (score > trueIM$markovs[[i]]$.score && !is.null(score) && !is.null(trueIM$markovs[[i]]$.score)) {
                              
                              if (use.verbose) {
                                print("New graph being inserted into MEC") 
                              }
                              markov <- list(.graph=graph$.in.edges, .score=score, .data=graph$.score$pp.dat$data)#,
                                                          #.params=apply.sem(converted, graph$.score$pp.dat$data))
                              num.markovs <<- num.markovs
                              
                              for (k in num.markovs - 1:i + 1) {
                                if (i == num.markovs) {
                                  trueIM$markovs[[i]] <- markov
                                  break
                                }
                                trueIM$markovs[[k]] <- trueIM$markovs[[k - 1]]
                              }
                              trueIM$markovs[[i]] <- markov
                              break
                            }
                          }
                          
                        },
                        
                        
                        ## Purpose: convert graph objects from edge lists to graphNEL objects
                        ## ----------------------------------------------------------------------
                        ##' @param graph graph object (containing .in.edges and .nodes)
                        ##' @return graphNEL object
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        convert = function(from) {
                          edgeList <- lapply(from$.in.edges, function(v) from$.nodes[v])
                          names(edgeList) <- from$.nodes
                          result <- new("graphNEL",
                                        nodes = from$.nodes,
                                        edgeL = edgeList,
                                        edgemode = "directed")
                          return(graph::reverseEdgeDirections(result))
                        },
                        
                        ## Purpose: apply structural equation modeling to a given graph
                        ## ----------------------------------------------------------------------
                        ##' @param converted graphNEL object to receive SEM modeling
                        ##' @param dataset original data off of which the graph is based
                        ##' @return SEM data
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        apply.sem = function(converted, dataset) {
                          graph.nodes <- igraph::igraph.from.graphNEL(converted)
                          edge.list <- igraph::get.edgelist(graph.nodes)
                          model <- paste(edge.list[,1], "~", edge.list[,2])

                          fit <- lavaan::sem(model, data=data.frame(dataset))
                          estimate <- lavaan::partable(fit)$est
                          estimate <- round(estimate,2)

                          names(estimate) <- graph::edgeNames(converted)
                          return(estimate)
                        },
                        
                        ## Purpose: averages the SEM data from all of the individual graphs and
                        ##          applies the results to the global graph
                        ## ----------------------------------------------------------------------
                        ##' @param params.list list of SEM data/parameters for all the graphs
                        ##' @return averaged SEM data for global graph
                        ## ----------------------------------------------------------------------
                        ## Author: Noah Frazier-Logue
                        average.sem = function(params.list) {
                          mat <- matrix(unlist(params.list), length(params.list), 
                                        length(params.list[[1]]))
                          means <- apply(rbind(mat), 2, "mean")
                          names(means) <- names(params.list[[1]])
                          
                          return(means)
                        }
                      )
)

IMaGESclass$methods(
  initialize = function(matrices = NULL,
                        scores = NULL,
                        penalty = 3,
                        imscore = NULL,
                        num.markovs=5, 
                        use.verbose=FALSE) {
    
    rawscores <- list()
    #handling raw matrices
    if (!is.null(matrices)) {
      trueIM$numDatasets <- length(matrices)
      #create score objects from GES
      for (i in 1:length(matrices)) {
        rawscores[[i]] <- new("GaussL0penObsScore", matrices[[i]], lambda = penalty)
      }  
    }
    #handle score objects
    else {
      assign("numDatasets", length(scores), envir=trueIM)
      for (i in 1:length(scores)) {
        scores[[i]]$pp.dat$lambda <- penalty
        rawscores[[i]] <- scores[[i]]
      }
    }
    
    use.verbose <<- use.verbose
    penalty <<- penalty
    .rawscores <<- rawscores
    graphs <- list()
    
    for (i in 1:length(.rawscores)) {
      #print("creating graph")
      graphs[[i]] <- create.graph(rawscores[[i]])
    }
    
    .graphs <<- graphs
    
    initialize.global()
    
    if (use.verbose) {
      print("Verbose mode selected.")
    }
    
    print("Running...")
    
    #create list of size num.markovs
    #initialize to lowest signed int value because you can't compare
    #ints and NULL
    num.markovs <<- num.markovs
    
    trueIM$markovs <- rep(list(list(.graph=NULL, .score=-2147483648)), num.markovs)
    num.runs <- ncol(.graphs[[1]]$.score$pp.dat$data) * ncol(.graphs[[1]]$.score$pp.dat$data)
    num.equal <- 0
    prev.score <- -2147483648
    for (i in 1:num.runs) {
      #run IMaGES
      run()
      if(use.verbose) {
        print(paste("Run ", i, " of ", num.runs))
      }
      #determine if early stopping should be used
      if ( !is.nan(IMScore()) && !is.nan(prev.score) && IMScore() == prev.score) {
        num.equal = num.equal + 1
      }
      else {
        num.equal = 0
        prev.score = IMScore()
      }
      if (num.equal == 5) {
        print("Stopping early. IMaGES run has converged on a representative graph.")
        break
      }
    }
    
    #de-double the edge directions
    .graphs[[1]]$.in.edges <- fix.edges(.graphs[[1]]$.in.edges)
    
    
    #apply double-edge fix to the rest of the graphs
    # if (length(.graphs) > 1) {
    #   for (i in 2:length(.graphs)) {
    #     .graphs[[i]]$.in.edges <- .graphs[[1]]$.in.edges
    #   }
    # }
    
    if (use.verbose) {
      print("Applying structural equation modeling to graphs.")
    }
    
    #apply SEM and structure the results 
    single.graphs <- list()
    params.list <- list()
    converted <- convert(list(.in.edges = graphs[[1]]$.in.edges, .nodes = .graphs[[1]]$.nodes))
    for (i in 1:length(.graphs)) {
      if (use.verbose) {
        print(paste("Applied to graph ", i))
      }
      params <- apply.sem(converted, .graphs[[i]]$.score$pp.dat$raw.data)
      #removing NA data
      for (k in 1:length(params)) {
        if (is.na(names(params[k]))) {
          params <- params[1:k - 1]
          break
        }
      }
      params.list[[i]] <- params
      single.converted <- convert(list(.in.edges = .graphs[[i]]$.in.edges, .nodes = .graphs[[i]]$.nodes))
      #single.converted <- converted
      single.graphs[[i]] <-list(.graph = single.converted, .params = params)
    }
    
    print("Done with IMaGES run.")
    
    means <- single.graphs[[1]]$.params
    
    
    if (length(single.graphs) > 1) {
      for (i in 2:length(single.graphs)) {
        for (k in 1:length(single.graphs[[i]]$.params)) {
          means[k] <- as.numeric(means[k]) + as.numeric(single.graphs[[i]]$.params[k])
        }
      }
    }
    
    #means.fixed <- means
    for (i in 1:length(means)) {
      means[i] <- as.numeric(format(round(as.numeric(means[i]) / length(.graphs), 2), nsmall=2))
    }
    
    
    

    print(paste("Final IMScore: ", trueIM$score))
    #global <- list(.graph = converted, .params = average.sem(params.list))
    global <- list(.graph = converted, .params = means)
    markovs <- list()

    if (use.verbose) {
      print("Applying structural equation modeling to markov equivalence class.")
    }
    
    #only add markovs from list that aren't NULL
    for (i in 1:num.markovs) {
      attempt <- tryCatch(
        {
          converted.markov <- convert(list(.in.edges = fix.edges(trueIM$markovs[[i]]$.graph), .nodes = .graphs[[1]]$.nodes))
          markovs[[i]] <- list(.graph=converted.markov, .params = apply.sem(converted.markov, trueIM$markovs[[i]]$.data))
          if (use.verbose) {
            print(paste("Applied to MEC ", i))
          }
        },
        error = function(e) {
          print(paste("MEC ", i, " encountered an error while calculating SEM data. Set to NULL so graph can still be displayed."))
          markovs[[i]] <- NULL
        }
      )

    }
    
    if (length(.graphs) == 1) {
      results <<- list(.global = global, .single.graphs = single.graphs, .markovs = markovs, .means = NULL, .std.errs = NULL)
      return(results)
    }
    
    
    #means <- average.sem(params)
    
    # means <- single.graphs[[1]]$.params
    # 
    # 
    # if (length(single.graphs) > 1) {
    #   for (i in 2:length(single.graphs)) {
    #     for (k in 1:length(single.graphs[[i]]$.params)) {
    #       means[k] <- means[k] + single.graphs[[i]]$.params[k]
    #     }
    #   }
    # }
    # 
    # #means.fixed <- means
    # for (i in 1:length(means)) {
    #   means[i] <- format(round(means[i] / length(.graphs), 2), nsmall=2)
    # }
    
    
    if (use.verbose) {
      print("Calculating means and standard errors")
    }
    
    std.errs <- means
    
    for (i in 1:length(means)) {
      sum <- 0
      #sum from 1 to k of (c_n - c_i)^2
      #for each edge
      for (k in 1:length(single.graphs)) {
        #(c_n - c_i)^2
        difference <- (means[i] - single.graphs[[k]]$.params[i])^2
        sum <- sum + difference
      }
      sigma <- sqrt(sum / (length(.graphs) - 1))
      std.err <- sigma / sqrt(length(.graphs))
      std.errs[i] <- std.err
    }
    if (use.verbose) {
      print("Finished.")
    }
    results <<- list(.global = global, .single.graphs = single.graphs, .markovs = markovs, .means = means, .std.errs = std.errs)
    return(results)
  }
)

#' Interventional essential graph for IMaGES
setRefClass("IMGraph",
            fields = list(
              .nodes = "vector",
              .in.edges = "list",
              .targets = "list",
              .score = "Score",
              #.current_repr = "list",
              .old.edges = "list",
              .edge.change = "list"
              
            ),
            
            validity = function(object) {
              ## Check in-edges
              if (!all(sapply(object$.in.edges, is.numeric))) {
                return("The vectors in 'in.edges' must contain numbers.")
              }
              if (!all(unique(unlist(object$.in.edges)) %in% 1:object$node.count())) {
                return(sprintf("Invalid edge source(s): edge sources must be in the range 1:%d.",
                               object$node.count()))
              }
              
              ## Check targets
              if (anyDuplicated(object$.targets)) {
                return("Targets are not unique.")
              }
              if (!all(unique(unlist(object$.targets)) %in% 1:object$node.count())) {
                return(sprintf("Invalid target(s): targets must be in the range 1:%d.",
                               object$node.count()))
              }
              
              return(TRUE)
            },
            
            methods = list(
              #' Constructor
              initialize = function(nodes,
                                    in.edges = replicate(length(nodes), integer(0), simplify = FALSE),
                                    targets = list(integer(0)),
                                    score = NULL) {
                ## Store nodes names
                if (missing(nodes)) {
                  stop("Argument 'nodes' must be specified.")
                }
                .nodes <<- as.character(nodes)
                
                ## Store in-edges
                stopifnot(is.list(in.edges) && length(in.edges) == length(nodes))
                # More error checking is done in validity check
                .in.edges <<- in.edges
                .current.repr <<- new("GaussParDAG", nodes=.nodes, in.edges=.in.edges)
                #names(.in.edges) <<- NULL
                
                ## Store targets
                setTargets(targets)
                
                ## Store score
                setScore(score)
              },
              
              #' Yields the number of nodes
              node.count = function() {
                length(.nodes)
              },
              
              #' Yields the total number of edges in the graph
              edge.count = function() {
                sum(vapply(.in.edges, length, 1L))
              },
              
              #' Getter and setter functions for score object
              getScore = function() {
                .score
              },
              
              setScore = function(score) {
                if (!is.null(score)) {
                  .score <<- score
                }
              },
              
              #' Getter and setter functions for targets list
              getTargets = function() {
                .targets
              },
              
              setTargets = function(targets) {
                .targets <<- lapply(targets, sort)
              },
              
              #' Creates a list of options for the C++ function "causalInference";
              #' internal function
              causal.inf.options = function(
                caching = TRUE,
                phase = c("forward"),
                iterate = length(phase) > 1,
                maxDegree = integer(0),
                maxSteps = 0,
                childrenOnly = integer(0),
                fixedGaps = NULL,
                adaptive = c("none", "vstructures", "triples"),
                verbose = 0,
                p = 0) {
                # Check for deprecated calling convention and issue a warning
                if (p > 0) {
                  warning(paste("Argument 'p' is deprecated in calls of ges() or gies",
                                "and will be disabled in future package versions;",
                                "please refer to the corresponding help page.", sep = " "))
                }
                
                # Error checks for supplied arguments
                # TODO extend!
                if (is.logical(adaptive)) {
                  adaptive <- ifelse(adaptive, "vstructures", "none")
                  warning(paste("The parameter 'adaptive' should not be provided as logical anymore;",
                                "cf. ?ges or gies", sep = " "))
                }
                #phase <- match.arg(phase, several.ok = TRUE)
                stopifnot(is.logical(iterate))
                adaptive <- match.arg(adaptive)
                if (is.null(fixedGaps)) {
                  adaptive <- "none"
                }
                list(caching = caching,
                     phase = phase,
                     iterate = iterate,
                     maxDegree = maxDegree,
                     maxSteps = maxSteps,
                     childrenOnly = childrenOnly,
                     fixedGaps = fixedGaps,
                     adaptive = adaptive,
                     DEBUG.LEVEL = as.integer(verbose))
              },
              
              #' Performs one greedy step
              greedy.step = function(alg.name="GIES-F", direction = c("forward"), verbose = FALSE, ...) {
                stopifnot(!is.null(score <- getScore()))
                
                ## Cast direction
                #direction <- match.arg(direction)
                #alg.name <- switch(direction,
                #                   forward = "GIES-F",
                #                   backward = "GIES-B",
                #                   turning = "GIES-T")
                #alg.name <- match.arg(alg.name)
                #print("calling GES")
                new.graph <- .Call("causalInferenceEdge",
                                   #.current_repr$.in.edges,
                                   .in.edges,
                                   score$pp.dat,
                                   alg.name,
                                   score$c.fcn,
                                   causal.inf.options(caching = FALSE, maxSteps = 1, verbose = verbose, phase=direction, ...),
                                   PACKAGE = "IMaGES")
                
                .old.edges <<- .in.edges
                .in.edges <<- new.graph$in.edges
                .saved.edges <<- new.graph$in.edges
                names(.in.edges) <<- .nodes
                .edge.change <<- list(new.graph$edge.change$src, new.graph$edge.change$dst, new.graph$edge.change$alg)
                return(FALSE)
              },
              
              #' undo individual step. used for calculating best step in IMaGES run
              undo.step = function() {
                .in.edges <<- .old.edges
              },
              
              #' redo individual step. used for calculating best step in IMaGES run
              redo.step = function() {
                .in.edges <<- .saved.edges
              },
              
              
              #' Yields a representative (estimating parameters via MLE)
              repr = function() {
                stopifnot(!is.null(score <- getScore()))
                
                result <- new("GaussParDAG", nodes =.nodes)
                # result$.in.edges <- .current_repr$.in.edges
                result$.in.edges <- .in.edges
                result$.params <- score$global.fit(result)
                
                return(result)
              },
              
              #' Calculates an optimal intervention target
              #'
              #' @param   max.size    maximum target size; allowed values: 1, p (= # nodes)
              ## TODO document that function... or better: provide a documented wrapper function
              opt.target = function(max.size) {
                .Call("optimalTarget", .in.edges, max.size, PACKAGE = "IMaGES")
              }
            ))


## Purpose: finds lowest two values that, when multiplied, create square-ish grid for
##          plotting graphs
## ----------------------------------------------------------------------
##' @param number number of graphs that need to be plotted
##' @return vector containing the two (very similar) values that, when multiplied together,
##'         result in a value that is >= to the input
## ----------------------------------------------------------------------
## Author: Noah Frazier-Logue
find_dimensions = function(number) {
  val1 <- floor(sqrt(number))
  val2 <- ceiling(sqrt(number))
  
  
  while (val1*val2 < number) {
    val2 <- val2 + 1
    if (val1 * val2 < number) {
      val1 <- val1 - 1
    }
  }
  return(c(val1, val2))
}



## Purpose: plot function for graphs generated by IMaGES
## ----------------------------------------------------------------------
##' @param graph.list graph object (containing .graph and .params) to plot
##' @param title plot title. defaults to "Global"
## ----------------------------------------------------------------------
## Author: Noah Frazier-Logue
plotIMGraph = function(graph.object, title="Global") {#, layout=NULL) {
  
  #TODO: convert to igraph calls
  #figure out how to make graphs scale to fill individual cells
  graph <- graph.object$.graph
  params <- graph.object$.params
  finalgraph <- Rgraphviz::agopen(graph, "", attrs=list(node=list(shape="ellipse")), edgeAttrs=list(label=params))
  Rgraphviz::plot(finalgraph, main=title)
  # graph <- igraph::igraph.from.graphNEL(graph.object$.graph)
  # params <- graph.object$.params
  # 
  # if (is.null(layout)) {
  #   layout <- igraph::layout_(graph, with_fr())
  # }
  # 
  # igraph::plot.igraph(graph, edge.label=params, layout=layout, main=title, edge.arrow.size=0.1, vertex.size=30)
}

## Purpose: plot global and individual graphs in a grid using the object
##          returned by IMaGES. 
## ----------------------------------------------------------------------
##' @param im.fits object returned by IMaGES
## ----------------------------------------------------------------------
## Author: Noah Frazier-Logue
plotAll = function(im.fits) {
  single.length <- length(im.fits$.single.graphs)
  plot.vals <- find_dimensions(single.length + 1)
  graphics::par(mfrow=plot.vals)
  
  # graph <- igraph::igraph.from.graphNEL(im.fits$.global$.graph)
  # layout <- igraph::layout_(graph, igraph::with_fr())
  #plotIMGraph(im.fits$.global, layout=layout)
  plotIMGraph(im.fits$.global)
  

  
  
  for (i in 1:single.length) {
    plotIMGraph(im.fits$.single.graphs[[i]], title=paste("Graph ", i))#, layout=layout)
  }
}

## Purpose: plot global and Markov Equivalence Class graphs in a grid using the object
##          returned by IMaGES. 
## ----------------------------------------------------------------------
##' @param im.fits object returned by IMaGES
## ----------------------------------------------------------------------
## Author: Noah Frazier-Logue
plotMarkovs = function(im.fits) {
  single.length <- length(im.fits$.markovs)
  #print(single.length)
  plot.vals <- find_dimensions(single.length + 1)
  graphics::par(mfrow=plot.vals)
  plotIMGraph(im.fits$.global)
  
  for (i in 1:single.length) {
    plotIMGraph(im.fits$.markovs[[i]], title=paste("MEC Graph ", i))
  }
  
}


## Purpose: wrapper function for IMaGESclass. used to control output
## so intermediary fields used in IMaGESclass are not accessible
## ----------------------------------------------------------------------
##' @param matrices list of datasets in matrix form to be used in IMaGES
##' @param penalty penalty value to be applied to BIC calculations
##' @param num.markovs number of Markov Equivalence Classes to return
##' @param use.verbose boolean indicating whether or not the user wants
##'                    verbose mode enabled
##' @return results field from IMaGESclass
## ----------------------------------------------------------------------
## Author: Noah Frazier-Logue
IMaGES = function(matrices, penalty=3, num.markovs=5, use.verbose=FALSE) {
  result <- IMaGESclass(matrices, penalty, num.markovs, use.verbose)
  return(result$results)
}

### GLOBAL VAR

#trueIM <- 100
trueIM <- new.env()
assign("score", 100, env=trueIM)
assign("isLocalIM", FALSE, env=trueIM)
assign("numDatasets", 1, env=trueIM)
assign("global.edges", list(), env=trueIM)
assign("markovs", list(), env=trueIM)