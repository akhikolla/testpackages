.check_nb_cores <- function(nb_cores=NULL) {
  if (is.null(nb_cores)) nb_cores <- blackbox.getOption("coreNbr") ## may be NULL
  machine_cores <- parallel::detectCores()
  if (is.null(nb_cores)) {
    nb_cores <- 1L ## default
    if (machine_cores>1L && interactive()) {
      if (! identical(blackbox.getOption("cores_avail_warned"),TRUE)) {
        message(paste(machine_cores,"cores are available for parallel computation\n(you may be allowed to fewer of them on a shared cluster).\n"))
        if ("Migraine" %in% blackbox.getOption("usedBy")) {
          message("Use 'CoreNumberForR' keyword in Migraine settings file to set the number of cores to be used by Migraine during the R blackbox analysis.\n Or use blackbox.options(coreNbr=<n>) to control coreNbr during the execution of the R script.")
        } else {
          message("Change 'coreNbr' argument to use some of them.\nUse blackbox.options(coreNbr=<n>) to control coreNbr globally.")
        }
        blackbox.options(cores_avail_warned=TRUE)
      }
    }
  } else if (nb_cores>machine_cores) {
    if (! identical(blackbox.getOption("nb_cores_warned"),TRUE)) {
      if ("Migraine" %in% blackbox.getOption("usedBy")) {
        warning(paste("More cores were requested than found by parallel::detectCores(). Check 'CoreNumberForR' keyword in Migraine settings file.\n",
                      "Number of available cores automatically reduced to the number of cores found by parallel::detectCores(). I continue."))
      } else {
        warning(paste("More cores were requested than found by parallel::detectCores(). Check blackbox.getOption(\"coreNbr\") argument.\n",
                      "Number of available cores automatically reduced to the number of cores found by parallel::detectCores(). I continue."))
      } 
      blackbox.options(nb_cores_warned=TRUE)
    }
    nb_cores <- machine_cores
  }
  return(nb_cores)
}

.init_cores <- local({
  doSNOW_warned <- FALSE
  function(nb_cores=NULL, ## passing explicit value from user
           ...) {  ## ... are arguments used by functions called by the loc_calc_logL function
    nb_cores <- .check_nb_cores(nb_cores=nb_cores)
    cores_info <- list(nb_cores=nb_cores)
    if (nb_cores > 1L) {
      cores_info$cl <- parallel::makeCluster(nb_cores) 
      dotenv <- list2env(list(...))
      parallel::clusterExport(cl=cores_info$cl, as.list(ls(dotenv)),envir=dotenv) 
      ## foreach is NOT a parallel backend so there is no point using it if doSNOW is not available
      if (cores_info$has_doSNOW <- ("package:doSNOW" %in% search())) {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        ## allows progressbar but then requires foreach
        assign(".Random.seed", R.seed, envir = .GlobalEnv) # loading (?) the namespace of 'snow' changes the global RNG state!
        fn <- get("registerDoSNOW", asNamespace("doSNOW"))
        do.call(fn,list(cl=cores_info$cl)) 
      } else {
        if ( ! doSNOW_warned) {
          message("If the 'doSNOW' package were attached, better load-balancing might be possible.")
          doSNOW_warned <<- TRUE
        } 
      }
    }
    return(cores_info)
  }
})

.run_cores <- function(loc_calc_logL, ## called for loc_calc_logL=indepCalcProfileLRforeachPoint which alway returns a list
                       pargrid, profileMethod, cores_info) {
  #prevmsglength <- 0 ## no longer used ?
  if (cores_info$nb_cores > 1L) {
    blackboxOptions <- blackbox.options()
    if (cores_info$has_doSNOW) {
      pb <- txtProgressBar(max = nrow(pargrid), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      parallel::clusterExport(cl=cores_info$cl, list("progress"),envir=environment()) ## slow! why?
      ii <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i'
      foreach_args <- list(
        ii = seq_len(nrow(pargrid)), 
        .packages= "blackbox",
        .options.snow = list(progress = progress),
        .inorder = TRUE, .errorhandling = "remove"
        #                                 "pass"## "pass" to see error messages
      )
      foreach_blob <- do.call(foreach::foreach,foreach_args)
      grid_obj <- foreach::`%dopar%`(foreach_blob,
                                     loc_calc_logL(pargrid[ii,], profileMethod, blackboxOptions))
      if (cores_info$has_doSNOW) close(pb)
    } else {
      parallel::clusterEvalQ(cores_info$cl, {library("blackbox")}) 
      pbopt <- pboptions(nout=min(100,2*nrow(pargrid)),type="timer")
      grid_obj <- pbapply(X=pargrid,MARGIN = 1L,FUN = loc_calc_logL, cl=cores_info$cl, profileMethod=profileMethod, blackboxOptions=blackboxOptions)
      # pblapply on a function that returns a list should (?) always retrun a list of lists
      pboptions(pbopt)
      # for ( i in 1:nrow(pargrid) ) {
      #   vec_obj <- vector(mode = "list",length(grid_obj_list[[1]]))
      #   names(vec_obj) <- names(grid_obj_list[[1]])
      #   for ( j in 1:length(grid_obj_list[[1]]) ) vec_obj[j] <- list(grid_obj_list[[i]][[j]])
      #   if ( i==1) {grid_obj <- vec_obj} else {grid_obj <- rbind(grid_obj,vec_obj)}
      # }
      #browser() # tres utile pour voir les sorties des calculs en parallele, débugger avec des return a chaque point de blocage, et pour voir l'avancement
    }
  } else { ## F I X M E: but is this ever run ? if nb_cores=1 .run_cores is not called and grid_list has another format ?
    grid_obj <- pbapply(pargrid,MARGIN = 1L,FUN = loc_calc_logL,cl=NULL, profileMethod=profileMethod, blackboxOptions=NULL)
    #grid_obj <- t(grid_obj) 
  }
  return(grid_obj) ## a list of lists
}

