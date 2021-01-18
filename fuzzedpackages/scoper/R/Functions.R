#### Classes ####

#' S4 class containing clonal assignments and summary data
#' 
#' \code{ScoperClones} stores output from \link{identicalClones}, \link{hierarchicalClones} and
#' \link{spectralClones} functions.
#'
#' @slot   db              \code{data.frame} of repertoire data including with clonal identifiers in 
#'                         the column specified during processing.
#' @slot   vjl_groups      \code{data.frame} of clonal summary, including sequence count, V gene, 
#'                         J gene, junction length, and clone counts.
#' @slot   inter_intra     \code{data.frame} containing minimum inter (between) and maximum intra 
#'                         (within) clonal distances.
#' @slot   eff_threshold   effective cut-off separating the inter (between) and intra (within) clonal 
#'                         distances.
#'
#' @seealso      \link{identicalClones}, \link{hierarchicalClones} and \link{spectralClones}
#'
#' @name         ScoperClones-class
#' @rdname       ScoperClones-class
#' @aliases      ScoperClones
#' @exportClass  ScoperClones
setClass("ScoperClones",
         slots=c(db="data.frame",
                 vjl_groups="data.frame",
                 inter_intra="data.frame",
                 eff_threshold="numeric"))

#### Methods ####

#' @param    x      ScoperClones object
#' 
#' @rdname   ScoperClones-class
#' @aliases  ScoperClones-method
#' @export
setMethod("print", c(x="ScoperClones"), function(x) { print(x@eff_threshold) })

#' @param    object  ScoperClones object
#' 
#' @rdname   ScoperClones-class
#' @aliases  ScoperClones-method
#' @export
setMethod("summary", c(object="ScoperClones"),
          function(object) { object@vjl_groups })

#' @param    y      ignored.
#' @param    ...    arguments to pass to \link{plotCloneSummary}.
#' 
#' @rdname   ScoperClones-class
#' @aliases  ScoperClones-method
#' @export
setMethod("plot", c(x="ScoperClones", y="missing"),
          function(x, y, ...) { plotCloneSummary(x, ...) })

#' @rdname   ScoperClones-class
#' @aliases  ScoperClones-method
#' @export
setMethod("as.data.frame", c(x="ScoperClones"),
          function(x) { as.data.frame(x@db) })

#### Internal functions ####

# find density gap
findGapSmooth <- function(vec) {
    # bandwidth <- kedd::h.ucv(vec, 4)$h
    # bandwidth <- density(vec)$bw
    # dens <- KernSmooth::bkde(vec, canonical=TRUE) #, bandwidth=bandwidth
    # suppressWarnings(dens <- density(vec, kernel="gaussian", adjust=1, bw="ucv"))  #"nrd0"
    dens <- density(vec)
    tryCatch({
        idy <- which(diff(sign(diff(dens$y))) == 2) + 1
        idx <- idy[which.min(dens$y[idy])]
        d <- ifelse(length(idx) != 0, dens$x[idx], NA)
    },
    error = function(e) {
        warning('No minimum was found between two modes.')
        return(NULL) })
    return(d)
}
# *****************************************************************************

# *****************************************************************************
# epsilon calculator by "infer"
infer <- function(vec) {
    vec <- sort(vec)
    # vec[1] <- vec[2]/2
    n <- length(vec)
    d <- NA
    # upper level search
    if (n > 2) {
        d <- findGapSmooth(vec=vec)    
        if (!is.na(d)) { 
            # d <- max(vec[vec <= d]) 
            d <- ifelse(max(vec[vec <= d]) == 0, d, max(vec[vec <= d]))
        }
    }
    # lower level search
    if (is.na(d)) {
        diffVec <- diff(vec)
        if (length(unique(diffVec[diffVec > 0])) == 1) {
            d <- ceiling(mean(vec))
        } else {
            x <- which.max(diffVec)
            d <- ifelse(vec[x] == 0, mean(c(vec[x], vec[vec>0][1])), vec[x])
        }   
    }
    return(d)
}
# *****************************************************************************

# *****************************************************************************
# kernel matrix calculator
krnlMtxGenerator <- function(mtx) {
    # Radial basis function kernel: In the Gaussian Kernel if two points are
    # close then K_ij≈1 and when two points are far apart then Kij≈0
    n <- nrow(mtx)
    # calculate epsilons
    epsilon <- rep(0, length=n)
    for (i in 1:n) {
        epsilon[i] <- infer(vec=mtx[i,])
    }    
    # calculate kernel matrix
    krnl_mtx <- matrix(data=1, nrow=n, ncol=n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            krnl_mtx[i,j] <- exp(-mtx[i,j]^2/(epsilon[i]*epsilon[j]))  #2*
            krnl_mtx[j,i] <- krnl_mtx[i,j]
        }
    }
    krnl_mtx[is.nan(krnl_mtx)] <- 1  # if mtx[i,j] and epsilon == 0
    krnl_mtx[krnl_mtx < 0.05] <- 0
    # krnl_mtx <- round(krnl_mtx, 3)
    return(krnl_mtx)
}
# *****************************************************************************

# *****************************************************************************
# affinity matrix calculator
# Disconnect those edges with distance larger than threshold
makeAffinity <- function(mtx_o, mtx_k, thd) {
    mtx_k[mtx_o > thd] <- 0
    return(mtx_k)
}
# *****************************************************************************

# *****************************************************************************
# laplacian matrix calculator
laplacianMtx <- function(entry) {
    # Calculate unnormalised Laplacian matrix and its eigenfunctions
    D <- diag(apply(entry, 1, sum))
    L <- D - entry
    return(L)
}
# *****************************************************************************

# *****************************************************************************
# range a vector from a to b
rangeAtoB <- function(x, a, b){
    return((b-a)*(x-min(x))/(max(x)-min(x)) + a)
}
# *****************************************************************************

# *****************************************************************************
# calculate likelihoods
likelihoods <- function(tot_mtx, sh_mtx, mutab_mtx) {
    sigma_sh <- sd(sh_mtx[upper.tri(sh_mtx)])
    sigma_tot <- sd(tot_mtx[upper.tri(tot_mtx)])
    if (sigma_sh %in% c(NA, 0)) { # there is no any prefences among pairs of sequences, therefore all likelihoods are zero
        z_mtx <- matrix(0, nrow=nrow(sh_mtx), ncol=nrow(sh_mtx))
    } else if (sigma_tot %in% c(NA, 0)) { # pairs with more shared mutations are more likly to belong to the same clone
        x_mtx <- 1.0 - exp(-sh_mtx^2/(2.0*sigma_sh^2))
        z_mtx <- mutab_mtx*x_mtx
    } else {
        x_mtx <- 1.0 - exp(-sh_mtx^2/(2.0*sigma_sh^2))
        y_mtx <- exp(-(tot_mtx-sh_mtx)^2/(2.0*sigma_tot^2))
        z_mtx <- mutab_mtx*x_mtx*y_mtx
    }
    diag(z_mtx) <- 1 # diagonals should have liklihoods equal to one
    return(z_mtx)
}
# *****************************************************************************

# *****************************************************************************
pairwiseMutions <- function(germ_imgt, 
                            seq_imgt, 
                            junc_length, 
                            len_limit = NULL,
                            cdr3 = FALSE,
                            mutabs = NULL,
                            norm_fact = TRUE) {
    
    ##### get number of seqs
    n <- unique(c(length(seq_imgt), length(germ_imgt)))
    ##### check number of sequences
    if (length(n) > 1) stop("germ_imgt and seq_imgt number should be the same")
    if (n == 1) stop("there should be at least two seqs")
    # check consensus length
    if (!is.null(len_limit)) {
        lenConsensus <- len_limit@seqLength
        seq_imgt  <- substr(seq_imgt, start = 1, stop = lenConsensus)
        germ_imgt <- substr(germ_imgt, start = 1, stop = lenConsensus)
        eff_germ <- ifelse(length(unique(germ_imgt)) == 1, 
                           unique(germ_imgt), 
                           consensusSequence(sequences = unique(germ_imgt),
                                             muFreqColumn = NULL, 
                                             lenLimit = lenConsensus,
                                             method = "catchAll",
                                             minFreq = NULL,
                                             includeAmbiguous = FALSE,
                                             breakTiesStochastic = FALSE,
                                             breakTiesByColumns = NULL, 
                                             db = NULL)$cons)
    } else {
        ##### constants
        lv <- ifelse(cdr3, shazam::IMGT_V@seqLength, shazam::IMGT_V@seqLength - 3)
        trim_l <- junc_length
        ##### trim out junction/cdr3 segments from seq_imgt
        seq_imgt <- sapply(1:length(seq_imgt), function(i){
            x <- strsplit(seq_imgt[i], split="")[[1]]
            x[(lv+1):(lv+trim_l)] <- ""   # x[(lv+1):(lv+trim_l[i])] <- ""
            return(paste(x, collapse=""))
        })
        ##### Pads ragged ends
        l <- unique(nchar(seq_imgt))
        if (length(l) > 1) {
            seq_imgt <- padSeqEnds(seq = seq_imgt, len = NULL, start = FALSE, pad_char = "N")
        }
        ##### trim out junction/cdr3 segments from germ_imgt
        germ_imgt <- sapply(1:length(germ_imgt), function(i){
            x <- strsplit(germ_imgt[i], split="")[[1]]
            x[(lv+1):(lv+trim_l)] <- ""  # x[(lv+1):(lv+trim_l[i])] <- ""
            return(paste(x, collapse=""))
        })
        ##### Pads ragged ends
        l <- unique(nchar(germ_imgt))
        if (length(l) > 1) {
            germ_imgt <- padSeqEnds(seq = germ_imgt, len = NULL, start = FALSE, pad_char = "N")
        }
        ##### find consensus germline (allel level grouping)
        # see arg "method" from shazam::collapseClones function
        eff_germ <- ifelse(length(unique(germ_imgt)) == 1, 
                           unique(germ_imgt), 
                           consensusSequence(sequences = unique(germ_imgt),
                                             muFreqColumn = NULL, 
                                             lenLimit = NULL,
                                             method = "catchAll",
                                             minFreq = NULL,
                                             includeAmbiguous = FALSE,
                                             breakTiesStochastic = FALSE,
                                             breakTiesByColumns = NULL, 
                                             db = NULL)$cons)
        ##### check germ and seqs lengths
        seq_imgt_lent <- unique(nchar(seq_imgt))
        germ_imgt_lent <- unique(nchar(germ_imgt))
        eff_germ_lent <- nchar(eff_germ)
        lenConsensus <- min(seq_imgt_lent, germ_imgt_lent, eff_germ_lent) 
        ##### trim extra characters
        if ( seq_imgt_lent > lenConsensus)  { seq_imgt <- substr( seq_imgt, start = 1, stop = lenConsensus) }
        if (germ_imgt_lent > lenConsensus) { germ_imgt <- substr(germ_imgt, start = 1, stop = lenConsensus) }
        if ( eff_germ_lent > lenConsensus)  { eff_germ <- substr( eff_germ, start = 1, stop = lenConsensus) }  
    }
    ##### count informative positions
    if (norm_fact) {
        informative_pos <- sapply(1:n, function(x){ sum(stri_count(seq_imgt[x], fixed = c("A","C","G","T"))) })   
    } else {
        informative_pos <- rep(1, n)
    }
    ##### convert eff_germ and seq_imgt to matrices
    seqsMtx <- matrix(NA, nrow=n, ncol=lenConsensus)
    effMtx <- matrix(NA, nrow=n, ncol=lenConsensus)
    for (i in 1:n) {
        seqsMtx[i, ] <- strsplit(seq_imgt[i], split = "")[[1]][1:lenConsensus]
        effMtx[i, ] <- strsplit(eff_germ, split = "")[[1]][1:lenConsensus]
    }
    ##### make a distance matrix
    dnaMtx <- getDNAMatrix(gap = 0)
    mutMtx <- matrix(NA, nrow=n, ncol=lenConsensus)
    for (i in 1:n) {
        mutMtx[i, ] <- sapply(1:lenConsensus, function(j) {
            return(dnaMtx[effMtx[i,j], seqsMtx[i,j]]) 
        })
    }
    ##### make a mutation matrix
    mutMtx <- matrix(paste0(effMtx, mutMtx, seqsMtx), nrow=n, ncol=lenConsensus)
    ##### clean non-mutated elements
    mutMtx[grepl(pattern="0", mutMtx)] <- NA
    ##### check mutabilities
    ##### make a motif matrix
    motifMtx <- matrix(0, nrow=n, ncol=lenConsensus)
    if (!is.null(mutabs)) {
        for (i in 1:n) {
            for (j in 3:(lenConsensus-2)) {
                motifMtx[i, j] <- mutabs[substr(germ_imgt[i], start = j-2, stop = j+2)]
            }
        }
        motifMtx[is.na(motifMtx)] <- 0
    }
    ##### calculate mutation matrix
    results <- pairwiseMutMatrix(informative_pos = informative_pos, 
                                 mutMtx = mutMtx, 
                                 motifMtx = motifMtx)
    sh_mtx <- results$sh_mtx
    tot_mtx <- results$tot_mtx
    mutab_mtx <- results$mutab_mtx
    ##### make symmetric matrix
    sh_mtx[lower.tri(sh_mtx)] <- t(sh_mtx)[lower.tri(sh_mtx)]
    tot_mtx[lower.tri(tot_mtx)] <- t(tot_mtx)[lower.tri(tot_mtx)]
    mutab_mtx[lower.tri(mutab_mtx)] <- t(mutab_mtx)[lower.tri(mutab_mtx)]
    # return results
    return_list <- list("pairWiseSharedMut" = sh_mtx,
                        "pairWiseTotalMut" = tot_mtx,
                        "pairWiseMutability" = mutab_mtx)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
### make a dataframe of unique seqs in each clone
uniqueSeq <- function(seqs) {
    seqs_db <- data.frame(value = seqs, name = names(seqs), stringsAsFactors = FALSE) %>%
        dplyr::group_by(!!!rlang::syms(c("name", "value"))) %>% # alternatively: group_by(name) if name value pair is always unique
        dplyr::slice(1) %>%
        dplyr::ungroup()
    seqs <- seqs_db$value
    names(seqs) <- seqs_db$name
    return(seqs)
}
# *****************************************************************************

# *****************************************************************************
# inter-clone-distance vs intra-clone-distance
calculateInterVsIntra <- function(db,
                                  clone,
                                  vjl_gps,
                                  junction = "junction",
                                  cdr3 = FALSE,
                                  cdr3_col = NA,
                                  nproc = 1,
                                  verbose = FALSE) {
    ### Create cluster of nproc size and export namespaces
    if(nproc == 1) {
        # If needed to run on a single core/cpu then, register DoSEQ
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK", outfile = "")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    ### expoer function to clusters
    if (nproc > 1) { 
        export_functions <- list("pairwiseDist", "getDNAMatrix", "uniqueSeq")
        parallel::clusterExport(cluster, export_functions, envir=environment())
    }
    
    n_groups <- nrow(vjl_gps)  
    ### check the progressbar
    if (verbose) {
        pb <- progressBar(n_groups)
    }
    
    k <- NULL
    # open dataframes
    vec_ff <- foreach(k=1:n_groups,
                      .combine="c",
                      .errorhandling='stop') %dopar% {
                          
                          # *********************************************************************************
                          clones <- strsplit(vjl_gps$clone_id[k], split=",")[[1]]
                          l <- vjl_gps$junction_length[k]
                          n_clones <- length(clones)
                          seqs <- db[[ifelse(cdr3, cdr3_col, junction)]][db[[clone]] %in% clones]
                          names(seqs) <- db[[clone]][db[[clone]] %in% clones]
                          seqs <- uniqueSeq(seqs)
                          ### calculate distance matrix among all seqs
                          dist_mtx <- pairwiseDist(seqs, dist_mat=getDNAMatrix(gap = 0))
                          ### prealoocate a vector = no. of max-dist in each clone (intra) + no. of min-dist between clones (inter)
                          nrow_f <- n_clones + n_clones*(n_clones-1)/2
                          vec_f <- rep(NA, nrow_f)
                          ### calculate minimum and maximum distance in each clone
                          n <- 0
                          if (n_clones == 1) {
                              n <- n + 1
                              vec_f[n] <- max(dist_mtx)/l
                              names(vec_f)[n] <- paste(clones[1], "NA", "intra", sep="_")
                          } else {
                              for (i in 1:(n_clones-1)) {
                                  xx <- dist_mtx[rownames(dist_mtx) == clones[i], colnames(dist_mtx) == clones[i]]
                                  n <- n + 1
                                  vec_f[n] <- max(xx)/l
                                  names(vec_f)[n] <- paste(clones[i], "NA", "intra", sep="_")
                                  for (j in (i+1):n_clones) {
                                      xy <- dist_mtx[rownames(dist_mtx) == clones[i], colnames(dist_mtx) == clones[j]]
                                      n <- n + 1
                                      vec_f[n] <- min(xy)/l
                                      names(vec_f)[n] <- paste(clones[i], clones[j], "inter", sep="_")
                                  }
                              }
                              yy <- dist_mtx[rownames(dist_mtx) == clones[j], colnames(dist_mtx) == clones[j]]
                              n <- n + 1
                              vec_f[n] <- max(yy)/l
                              names(vec_f)[n] <- paste(clones[j], "NA", "intra", sep="_")
                          }
                          
                          # Update progress
                          if (verbose) { pb$tick() }
                          
                          return(vec_f)
                          # *********************************************************************************
                      }
    
    ### Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    # convert to a data.frame
    db_dff <- data.frame(keyName = names(vec_ff), 
                         distance = vec_ff, 
                         row.names=NULL)
    db_dff$label <- "intra"
    db_dff$label[grepl("inter", db_dff$keyName)] <- "inter"
    clones_xy <- data.frame(matrix(unlist(stri_split_fixed(db_dff$keyName, "_", n=3)),
                                   nrow=nrow(db_dff),
                                   byrow=T),
                            stringsAsFactors=FALSE)
    db_dff <- cbind(clones_xy, db_dff)
    db_dff$keyName <- NULL
    colnames(db_dff)[colnames(db_dff) == "X1"] <- "clone_id_x"
    colnames(db_dff)[colnames(db_dff) == "X2"] <- "clone_id_y"
    db_dff$X3 <- NULL
    
    # return results
    return(db_dff)
}
# *****************************************************************************
### Define verbose reporting function
printVerbose <- function(n_groups, vjl_gp, model, method, linkage, cdr3,
                         gp_vcall, gp_jcall, gp_lent, gp_size, n_cluster) {
    method <- ifelse(model == "hierarchical", paste(linkage, "linkage", method, sep="-"), method)
    cat("     TOTAL_GROUPS> ", n_groups,  "\n", sep=" ")
    cat("            GROUP> ", vjl_gp, "\n", sep=" ")
    cat("   SEQUENCE_COUNT> ", gp_size,   "\n", sep=" ")
    cat("           V_CALL> ", gp_vcall,  "\n", sep=" ")
    cat("           J_CALL> ", gp_jcall,  "\n", sep=" ")
    cat("  JUNCTION_LENGTH> ", gp_lent,   "\n", sep=" ") 
    cat("            MODEL> ", model,     "\n", sep=" ")
    cat("           METHOD> ", method,    "\n", sep=" ")
    cat("             CDR3> ", cdr3,      "\n", sep=" ")
    cat("            CLONE> ", n_cluster, "\n", sep=" ")
    cat("", "\n", sep=" ")
}
# *****************************************************************************

# *****************************************************************************
logVerbose <- function(out_dir, log_verbose_name,
                       n_groups, vjl_gp, model, method, linkage, cdr3,
                       gp_vcall, gp_jcall, gp_lent, gp_size, n_cluster) {
    method <- ifelse(model == "hierarchical", paste(linkage, "linkage", method, sep="-"), method)
    cat("     TOTAL_GROUPS> ", n_groups,  "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("            GROUP> ", vjl_gp, "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("   SEQUENCE_COUNT> ", gp_size,   "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("           V_CALL> ", gp_vcall,  "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("           J_CALL> ", gp_jcall,  "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("  JUNCTION_LENGTH> ", gp_lent,   "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)   
    cat("            MODEL> ", model,     "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)   
    cat("           METHOD> ", method,    "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("             CDR3> ", cdr3,      "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("            CLONE> ", n_cluster, "\n", sep=" ", file=file.path(out_dir, log_verbose_name), append=TRUE)
    cat("", "\n", sep=" ", file=file.path(out_dir, log_verbose_name), append=TRUE)
}
# *****************************************************************************

# *****************************************************************************
prepare_db <- function(db, 
                       junction = "junction", v_call = "v_call", j_call = "j_call",
                       first = FALSE, cdr3 = FALSE, 
                       cell_id = NULL, locus = NULL, only_heavy = TRUE,
                       mod3 = FALSE, max_n = 0) {
    # add junction length column
    db$junction_l <- stri_length(db[[junction]])
    junction_l <- "junction_l"
    
    ### check for mod3
    # filter mod 3 junction lengths
    if (mod3) {
        n_rmv_mod3 <- sum(db[[junction_l]]%%3 != 0)
        db <- db %>% 
            dplyr::filter(!!rlang::sym(junction_l)%%3 == 0)
    } else {
        n_rmv_mod3 <- 0
    }
    
    ### check for cdr3
    # filter junctions with length > 6
    if (cdr3) {
        n_rmv_cdr3 <- sum(db[[junction_l]] <= 6)
        db <- db %>% 
            dplyr::filter(!!rlang::sym(junction_l) > 6)
        # add cdr3 column
        db$cdr3_col <- substr(db[[junction]], 4, db[[junction_l]]-3)
        cdr3_col <- "cdr3_col"
    } else {
        n_rmv_cdr3 <- 0
        cdr3_col <- NA
    }
    
    ### check for N's
    # Count the number of 'N's in junction
    if (!is.null(max_n)) {
        n_rmv_N <- sum(stri_count(db[[junction]], regex = "N") > max_n)
        db <- db %>% 
            dplyr::filter(stri_count(!!rlang::sym(junction), regex = "N") <= max_n)
    } else {
        n_rmv_N <- 0
    }
    
    ### Parse V and J columns to get gene
    db <- groupGenes(db,
                     v_call = v_call,
                     j_call = j_call,
                     junc_len = NULL,
                     cell_id = cell_id,
                     locus = locus,
                     only_heavy = only_heavy,
                     first = first)
    
    ### groups to use
    groupBy <- c("vj_group", junction_l)
    
    ### assign group ids to db
    db$vjl_group <- db %>%
        dplyr::group_by(!!!rlang::syms(groupBy)) %>%
        dplyr::group_indices()
    
    ### retrun results
    return_list <- list("db" = db, 
                        "n_rmv_mod3" = n_rmv_mod3, 
                        "n_rmv_cdr3" = n_rmv_cdr3,
                        "n_rmv_N" = n_rmv_N,
                        "junction_l" = junction_l,
                        "cdr3_col" =  cdr3_col)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
# @export
pairwiseMutMatrix <- function(informative_pos, mutMtx, motifMtx) {
    pairwiseMutMatrixRcpp(informative_pos, mutMtx, motifMtx)
}
# *****************************************************************************

#### plotCloneSummary ####

#' Plot clonal clustering summary
#' 
#' \code{plotCloneSummary} plots the results in a \code{ScoperClones} object returned 
#' by \code{spectralClones}, \code{identicalClones} or \code{hierarchicalClones}.  
#' Includes the minimum inter (between) and maximum intra (within) clonal distances 
#' and the calculated efective threshold.
#'
#' @param    data      \link{ScoperClones} object output by the \link{spectralClones}, 
#'                     \link{identicalClones} or \link{hierarchicalClones}.
#' @param    xmin      minimum limit for plotting the x-axis. If \code{NULL} the limit will 
#'                     be set automatically.
#' @param    xmax      maximum limit for plotting the x-axis. If \code{NULL} the limit will 
#'                     be set automatically.
#' @param    breaks    number of breaks to show on the x-axis. If \code{NULL} the breaks will 
#'                     be set automatically.
#' @param    binwidth  binwidth for the histogram. If \code{NULL} the binwidth 
#'                     will be set automatically.
#' @param    title     string defining the plot title.
#' @param    size      numeric value for lines in the plot.
#' @param    silent    if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                     object; if \code{FALSE} draw the plot.
#' @param    ...       additional arguments to pass to ggplot2::theme.
#' 
#' @return   A ggplot object defining the plot.
#'
#' @seealso  See \link{ScoperClones} for the the input object definition.  
#'           See \link{spectralClones}, \link{identicalClones} and \link{hierarchicalClones} 
#'           for generating the input object.
#'
#' @examples
#' # Find clones
#' results <- hierarchicalClones(ExampleDb, threshold=0.15)
#' 
#' # Plot clonal summaries 
#' plot(results, binwidth=0.02)
#' 
#' @export
plotCloneSummary <- function(data, xmin=NULL, xmax=NULL, breaks=NULL, 
                             binwidth=NULL, title=NULL, size=0.75, silent=FALSE, 
                             ...) {
    
    eff_threshold <- data@eff_threshold
    xdf <- select(data@inter_intra, c("distance", "label"))
    data_intra <- xdf %>% 
        dplyr::filter(!!rlang::sym("label") == "intra", !!rlang::sym("distance") > 0)
    data_inter <- xdf %>%
        dplyr::filter(!!rlang::sym("label") == "inter", !!rlang::sym("distance") > 0)
    
    # fill color
    fill_manual <- c("intra"="grey30",
                     "inter"="grey60")
    
    # ggplot workaround
    if (is.null(xmin)) { xmin <- NA }
    if (is.null(xmax)) { xmax <- NA }
    
    ### plot
    p <- ggplot() +
        baseTheme() +
        theme(plot.title = element_text(size = 13, hjust = 0.5),
              axis.title.x = element_blank(),
              axis.title = element_text(size=13),
              axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              legend.position = "bottom",
              legend.text=element_text(size=12)) +
        xlab("Normalized hamming distance") +
        ylab("Density") +
        scale_y_continuous(labels = abs) +
        scale_fill_manual(name="",
                          values=fill_manual,
                          labels=c("intra"="maximum-distance within clones  ",
                                   "inter"="minimum-distance between clones  "))
    
    # Plot within clonal distances
    if (nrow(data_intra) > 0) {
        p <- p + 
            geom_histogram(data = data_intra,
                           aes_string(x = "distance", y = "..density..", fill = "label"),
                           binwidth = binwidth, color = "white", alpha = 0.85) +
            geom_density(data = data_intra, 
                         aes_string(x = "distance"),
                         size = size, color = "grey30")
    }
    
    # Plot between clonal distances
    if (nrow(data_inter) > 0) {
        p <- p + 
            geom_histogram(data = data_inter,
                           aes_string(x = "distance", y = "..density..", fill = "label"),
                           binwidth = binwidth, color = "white", alpha = 0.75) +
            geom_density(data = data_inter, 
                         aes_string(x = "distance"),
                         size = size, color = "grey60")
    } else {
        warning("No inter clonal distance is detected. Each group of sequences with same V-gene, J-gene, and junction length may contain only one clone.")
    }
    
    # Plot vertical treshold line
    if (!is.na(eff_threshold)) {
        p <- p + 
            ggtitle(paste("Effective threshold=", eff_threshold)) +
            geom_vline(xintercept=eff_threshold, color="grey30", linetype=2, size=size)
    } else {
        p <- p + 
            ggtitle(paste("Effective threshold not found"))
    }
    
    # Add x limits
    if (is.null(breaks) & (!is.na(xmin) | !is.na(xmax))) {
        p <- p + coord_cartesian(xlim = c(xmin, xmax))
    }
    # Set breaks
    if (!is.null(breaks)) {
        p <- p + scale_x_continuous(breaks=pretty_breaks(n=breaks),
                                    limits=c(xmin, xmax))
    }
    # Add Title
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    
    # Add additional theme elements
    p <- p + do.call(theme, list(...))
    
    # Plot
    if (!silent) {
        plot(p)
    } else {
        return(p)
    }
}

#### identicalClones ####

#' Sequence identity method for clonal partitioning
#'
#' \code{identicalClones} provides a simple sequence identity based partitioning 
#' approach for inferring clonal relationships in high-throughput Adaptive Immune Receptor 
#' Repertoire sequencing (AIRR-seq) data. This approach partitions B or T cell receptor 
#' sequences into clonal groups based on junction region sequence identity within 
#' partitions that share the same V gene, J gene, and junction length, allowing for 
#' ambiguous V or J gene annotations.
#'
#' @param    db                 data.frame containing sequence data.
#' @param    method             one of the \code{"nt"} for nucleotide based clustering or 
#'                              \code{"aa"} for amino acid based clustering.
#' @param    junction           character name of the column containing junction sequences.
#'                              Also used to determine sequence length for grouping.
#' @param    v_call             character name of the column containing the V-segment allele calls.
#' @param    j_call             character name of the column containing the J-segment allele calls.
#' @param    clone              the output column name containing the clonal clustering identifiers.
#' @param    cell_id            name of the column containing cell identifiers or barcodes. 
#'                              If specified, grouping will be performed in single-cell mode
#'                              with the behavior governed by the \code{locus} and 
#'                              \code{only_heavy} arguments. If set to \code{NULL} then the 
#'                              bulk sequencing data is assumed.
#' @param    locus              name of the column containing locus information. 
#'                              Only applicable to single-cell data.
#'                              Ignored if \code{cell_id=NULL}.
#' @param    only_heavy         use only the IGH (BCR) or TRB/TRD (TCR) sequences 
#'                              for grouping. Only applicable to single-cell data.
#'                              Ignored if \code{cell_id=NULL}.
#' @param    split_light        split clones by light chains. Ignored if \code{cell_id=NULL}.
#' @param    first              specifies how to handle multiple V(D)J assignments for initial grouping. 
#'                              If \code{TRUE} only the first call of the gene assignments is used. 
#'                              If \code{FALSE} the union of ambiguous gene assignments is used to 
#'                              group all sequences with any overlapping gene calls.
#' @param    cdr3               if \code{TRUE} removes 3 nucleotides from both ends of \code{"junction"} 
#'                              prior to clustering (converts IMGT junction to CDR3 region). 
#'                              If \code{TRUE} this will also remove records with a junction length 
#'                              less than 7 nucleotides.
#' @param    mod3               if \code{TRUE} removes records with a \code{junction} length that is not divisible by 
#'                              3 in nucleotide space.
#' @param    max_n              The maximum number of N's to permit in the junction sequence before excluding the 
#'                              record from clonal assignment. Default is set to be zero. Set it as \code{"NULL"} for no 
#'                              action.
#' @param    nproc              number of cores to distribute the function over.
#' @param    verbose            if \code{TRUE} prints out a summary of each step cloning process.
#'                              if \code{FALSE} (default) process cloning silently.
#' @param    log                output path and filename to save the \code{verbose} log. 
#'                              The input file directory is used if path is not specified.
#'                              The default is \code{NULL} for no action.
#' @param    summarize_clones   if \code{TRUE} performs a series of analysis to assess the clonal landscape
#'                              and returns a \link{ScoperClones} object. If \code{FALSE} then
#'                              a modified input \code{db} is returned.
#'
#' @return
#' If \code{summarize_clones=TRUE} (default) a \link{ScoperClones} object is returned that includes the 
#' clonal assignment summary information and a modified input \code{db} in the \code{db} slot that 
#' contains clonal identifiers in the specified \code{clone} column.
#' If \code{summarize_clones=FALSE} modified \code{data.frame} is returned with clone identifiers in the 
#' specified \code{clone} column.
#' 
#' @section Single-cell data:
#' To invoke single-cell mode the \code{cell_id} argument must be specified and the \code{locus} 
#' column must be correct. Otherwise, clustering will be performed with bulk sequencing assumptions, 
#' using all input sequences regardless of the values in the \code{locus} column.
#' 
#' Values in the \code{locus} column must be one of \code{c("IGH", "IGI", "IGK", "IGL")} for BCR 
#' or \code{c("TRA", "TRB", "TRD", "TRG")} for TCR sequences. Otherwise, the operation will exit and 
#' return and error message.
#' 
#' Under single-cell mode with paired-chain sequences, there is a choice of whether 
#' grouping should be done by (a) using IGH (BCR) or TRB/TRD (TCR) sequences only or
#' (b) using IGH plus IGK/IGL (BCR) or TRB/TRD plus TRA/TRG (TCR) sequences. 
#' This is governed by the \code{only_heavy} argument. There is also choice as to whether 
#' inferred clones should be split by the light/short chain (IGK, IGL, TRA, TRG) following 
#' heavy/long chain clustering, which is governed by the \code{split_light} argument.
#' 
#' In single-cell mode, clonal clustering will not be performed on data were cells are 
#' assigned multiple heavy/long chain sequences (IGH, TRB, TRD). If observed, the operation 
#' will exit and return an error message. Cells that lack a heavy/long chain sequence (i.e., cells with 
#' light/short chains only) will be assigned a \code{clone_id} of \code{NA}.
#'
#' @seealso See \link{plotCloneSummary} for plotting summary results. See \link{groupGenes} for 
#' more details about grouping requirements.
#'
#' @examples
#' # Find clonal groups
#' results <- identicalClones(ExampleDb)
#' 
#' # Retrieve modified input data with clonal clustering identifiers
#' df <- as.data.frame(results)
#' 
#' # Plot clonal summaries 
#' plot(results, binwidth=0.02)
#' 
#' @export
identicalClones <- function(db, method=c("nt", "aa"), junction="junction", 
                            v_call="v_call", j_call="j_call", clone="clone_id",
                            cell_id=NULL, locus="locus", only_heavy=TRUE, split_light=TRUE,
                            first=FALSE, cdr3=FALSE, mod3=FALSE, max_n=0, nproc=1,
                            verbose=FALSE, log=NULL, 
                            summarize_clones=TRUE) {
    
    results <- defineClonesScoper(db = db,
                                  method = match.arg(method), model = "identical", 
                                  junction = junction, v_call = v_call, j_call = j_call, clone = clone,
                                  cell_id = cell_id, locus = locus, only_heavy = only_heavy, split_light = split_light,
                                  first = first, cdr3 = cdr3, mod3 = mod3, max_n = max_n, nproc = nproc,        
                                  verbose = verbose, log = log, 
                                  summarize_clones = summarize_clones)
    
    ### return results
    if (summarize_clones) {
        return_list <- new("ScoperClones",
                           db = results$db,
                           vjl_groups = results$vjl_gps,
                           inter_intra = results$inter_intra,
                           eff_threshold = results$eff_threshold)    
        
        return(return_list)
    } else {
        return(results)
    }
}


#### hierarchicalClones ####

#' Hierarchical clustering method for clonal partitioning
#'
#' \code{hierarchicalClones} provides a hierarchical agglomerative clustering 
#' approach to infer clonal relationships in high-throughput Adaptive Immune Receptor 
#' Repertoire sequencing (AIRR-seq) data. This approach clusters B or T cell receptor 
#' sequences based on junction region sequence similarity within partitions that share the 
#' same V gene, J gene, and junction length, allowing for ambiguous V or J gene annotations.
#'
#' @param    db                 data.frame containing sequence data.
#' @param    threshold          a numeric scalar where the tree should be cut (the distance threshold for clonal grouping).
#' @param    method             one of the \code{"nt"} for nucleotide based clustering or 
#'                              \code{"aa"} for amino acid based clustering.
#' @param    linkage            available linkage are \code{"single"}, \code{"average"}, and \code{"complete"}.
#' @param    normalize	        method of normalization. The default is \code{"len"}, which divides the distance by the length 
#'                              of the sequence group. If \code{"none"} then no normalization if performed.
#' @param    junction           character name of the column containing junction sequences.
#'                              Also used to determine sequence length for grouping.
#' @param    v_call             character name of the column containing the V-segment allele calls.
#' @param    j_call             character name of the column containing the J-segment allele calls.
#' @param    clone              the output column name containing the clonal cluster identifiers.
#' @param    cell_id            name of the column containing cell identifiers or barcodes. 
#'                              If specified, grouping will be performed in single-cell mode
#'                              with the behavior governed by the \code{locus} and 
#'                              \code{only_heavy} arguments. If set to \code{NULL} then the 
#'                              bulk sequencing data is assumed.
#' @param    locus              name of the column containing locus information. 
#'                              Only applicable to single-cell data.
#'                              Ignored if \code{cell_id=NULL}.
#' @param    only_heavy         use only the IGH (BCR) or TRB/TRD (TCR) sequences 
#'                              for grouping. Only applicable to single-cell data.
#'                              Ignored if \code{cell_id=NULL}.
#' @param    split_light        split clones by light chains. Ignored if \code{cell_id=NULL}.
#' @param    first              specifies how to handle multiple V(D)J assignments for initial grouping. 
#'                              If \code{TRUE} only the first call of the gene assignments is used. 
#'                              If \code{FALSE} the union of ambiguous gene assignments is used to 
#'                              group all sequences with any overlapping gene calls.
#' @param    cdr3               if \code{TRUE} removes 3 nucleotides from both ends of \code{"junction"} 
#'                              prior to clustering (converts IMGT junction to CDR3 region). 
#'                              If \code{TRUE} this will also remove records with a junction length 
#'                              less than 7 nucleotides.
#' @param    mod3               if \code{TRUE} removes records with a \code{junction} length that is not divisible by 
#'                              3 in nucleotide space.
#' @param    max_n              The maximum number of \code{N} characters to permit in the junction sequence 
#'                              before excluding the record from clonal assignment. Note, with 
#'                              \code{linkage="single"} non-informative positions can create artifactual 
#'                              links between unrelated sequences. Use with caution. 
#'                              Default is set to be zero. Set it as \code{"NULL"} for no action.
#' @param    nproc              number of cores to distribute the function over.
#' @param    verbose            if \code{TRUE} prints out a summary of each step cloning process.
#'                              if \code{FALSE} (default) process cloning silently.
#' @param    log                output path and filename to save the \code{verbose} log. 
#'                              The input file directory is used if path is not specified.
#'                              The default is \code{NULL} for no action.
#' @param    summarize_clones   if \code{TRUE} performs a series of analysis to assess the clonal landscape
#'                              and returns a \link{ScoperClones} object. If \code{FALSE} then
#'                              a modified input \code{db} is returned.
#'
#' @return
#' If \code{summarize_clones=TRUE} (default) a \link{ScoperClones} object is returned that includes the 
#' clonal assignment summary information and a modified input \code{db} in the \code{db} slot that 
#' contains clonal identifiers in the specified \code{clone} column.
#' If \code{summarize_clones=FALSE} modified \code{data.frame} is returned with clone identifiers in the 
#' specified \code{clone} column.
#' 
#' @section Single-cell data:
#' To invoke single-cell mode the \code{cell_id} argument must be specified and the \code{locus} 
#' column must be correct. Otherwise, clustering will be performed with bulk sequencing assumptions, 
#' using all input sequences regardless of the values in the \code{locus} column.
#' 
#' Values in the \code{locus} column must be one of \code{c("IGH", "IGI", "IGK", "IGL")} for BCR 
#' or \code{c("TRA", "TRB", "TRD", "TRG")} for TCR sequences. Otherwise, the operation will exit and 
#' return and error message.
#' 
#' Under single-cell mode with paired-chain sequences, there is a choice of whether 
#' grouping should be done by (a) using IGH (BCR) or TRB/TRD (TCR) sequences only or
#' (b) using IGH plus IGK/IGL (BCR) or TRB/TRD plus TRA/TRG (TCR) sequences. 
#' This is governed by the \code{only_heavy} argument. There is also choice as to whether 
#' inferred clones should be split by the light/short chain (IGK, IGL, TRA, TRG) following 
#' heavy/long chain clustering, which is governed by the \code{split_light} argument.
#' 
#' In single-cell mode, clonal clustering will not be performed on data were cells are 
#' assigned multiple heavy/long chain sequences (IGH, TRB, TRD). If observed, the operation 
#' will exit and return an error message. Cells that lack a heavy/long chain sequence (i.e., cells with 
#' light/short chains only) will be assigned a \code{clone_id} of \code{NA}.
#' 
#' @seealso 
#' See \link{plotCloneSummary} for plotting summary results. See \link{groupGenes} for 
#' more details about grouping requirements.
#'
#' @examples
#' # Find clonal groups
#' results <- hierarchicalClones(ExampleDb, threshold=0.15)
#' 
#' # Retrieve modified input data with clonal clustering identifiers
#' df <- as.data.frame(results)
#' 
#' # Plot clonal summaries 
#' plot(results, binwidth=0.02)
#' 
#' @export
hierarchicalClones <- function(db, threshold, method=c("nt", "aa"), linkage=c("single", "average", "complete"), 
                               normalize=c("len", "none"), junction="junction", 
                               v_call="v_call", j_call="j_call", clone="clone_id",
                               cell_id=NULL, locus="locus", only_heavy=TRUE, split_light=TRUE,
                               first=FALSE, cdr3=FALSE, mod3=FALSE, max_n=0, nproc=1,
                               verbose=FALSE, log=NULL,
                               summarize_clones=TRUE) {
    
    results <- defineClonesScoper(db = db, threshold = threshold, model = "hierarchical", 
                                  method = match.arg(method), linkage = match.arg(linkage), normalize = match.arg(normalize),
                                  junction = junction, v_call = v_call, j_call = j_call, clone = clone,
                                  cell_id = cell_id, locus = locus, only_heavy = only_heavy, split_light = split_light,
                                  first = first, cdr3 = cdr3, mod3 = mod3, max_n = max_n, nproc = nproc,   
                                  verbose = verbose, log = log, 
                                  summarize_clones = summarize_clones)
    
    # return results
    if (summarize_clones) {
        return_list <- new("ScoperClones",
                           db = results$db,
                           vjl_groups = results$vjl_gps,
                           inter_intra = results$inter_intra,
                           eff_threshold = results$eff_threshold)  
        return(return_list)
    } else {
        return(results)
    }
}

#### spectralClones ####

#' Spectral clustering method for clonal partitioning
#'
#' \code{spectralClones} provides an unsupervised spectral clustering 
#' approach to infer clonal relationships in high-throughput Adaptive Immune Receptor 
#' Repertoire sequencing (AIRR-seq) data. This approach clusters B or T cell receptor 
#' sequences based on junction region sequence similarity and shared mutations within 
#' partitions that share the same V gene, J gene, and junction length, allowing for 
#' ambiguous V or J gene annotations.
#'
#' @param    db                 data.frame containing sequence data.
#' @param    method             one of the \code{"novj"} or \code{"vj"}. See Details for description.
#' @param    germline           character name of the column containing the germline or reference sequence.
#' @param    sequence           character name of the column containing input sequences. 
#' @param    junction           character name of the column containing junction sequences.
#'                              Also used to determine sequence length for grouping.
#' @param    v_call             character name of the column containing the V-segment allele calls.
#' @param    j_call             character name of the column containing the J-segment allele calls.
#' @param    clone              the output column name containing the clone ids.
#' @param    cell_id            name of the column containing cell identifiers or barcodes. 
#'                              If specified, grouping will be performed in single-cell mode
#'                              with the behavior governed by the \code{locus} and 
#'                              \code{only_heavy} arguments. If set to \code{NULL} then the 
#'                              bulk sequencing data is assumed.
#' @param    locus              name of the column containing locus information. 
#'                              Only applicable to single-cell data.
#'                              Ignored if \code{cell_id=NULL}.
#' @param    only_heavy         use only the IGH (BCR) or TRB/TRD (TCR) sequences 
#'                              for grouping. Only applicable to single-cell data.
#'                              Ignored if \code{cell_id=NULL}.
#' @param    split_light        split clones by light chains. Ignored if \code{cell_id=NULL}.
#' @param    targeting_model    \link[shazam]{TargetingModel} object. Only applicable if 
#'                              \code{method="vj"}. See Details for description.
#' @param    len_limit          \link{IMGT_V} object defining the regions and boundaries of the Ig 
#'                              sequences. If NULL, mutations are counted for entire sequence. Only 
#'                              applicable if \code{method} = \code{"vj"}.
#' @param    first              specifies how to handle multiple V(D)J assignments for initial grouping. 
#'                              If \code{TRUE} only the first call of the gene assignments is used. 
#'                              If \code{FALSE} the union of ambiguous gene assignments is used to 
#'                              group all sequences with any overlapping gene calls.
#' @param    cdr3               if \code{TRUE} removes 3 nucleotides from both ends of \code{"junction"} 
#'                              prior to clustering (converts IMGT junction to CDR3 region). 
#'                              If \code{TRUE} this will also remove records with a junction length 
#'                              less than 7 nucleotides.
#' @param    mod3               if \code{TRUE} removes records with a \code{junction} length that is not divisible by 
#'                              3 in nucleotide space.
#' @param    max_n              the maximum number of N's to permit in the junction sequence before excluding the 
#'                              record from clonal assignment. Default is set to be zero. Set it as \code{"NULL"} 
#'                              for no action.
#' @param    threshold          the supervising cut-off to enforce an upper-limit distance for clonal grouping.
#'                              A numeric value between (0,1).
#' @param    base_sim           required similarity cut-off for sequences in equal distances from each other.
#' @param    iter_max	        the maximum number of iterations allowed for kmean clustering step.
#' @param    nstart	            the number of random sets chosen for kmean clustering initialization.
#' @param    nproc              number of cores to distribute the function over.
#' @param    verbose            if \code{TRUE} prints out a summary of each step cloning process.
#'                              if \code{FALSE} (default) process cloning silently.
#' @param    log                output path and filename to save the \code{verbose} log. 
#'                              The input file directory is used if path is not specified.
#'                              The default is \code{NULL} for no action.
#' @param    summarize_clones   if \code{TRUE} performs a series of analysis to assess the clonal landscape
#'                              and returns a \link{ScoperClones} object. If \code{FALSE} then
#'                              a modified input \code{db} is returned.
#'
#' @return
#' If \code{summarize_clones=TRUE} (default) a \link{ScoperClones} object is returned that includes the 
#' clonal assignment summary information and a modified input \code{db} in the \code{db} slot that 
#' contains clonal identifiers in the specified \code{clone} column.
#' If \code{summarize_clones=FALSE} modified \code{data.frame} is returned with clone identifiers in the 
#' specified \code{clone} column.
#'
#' @details
#' If \code{method="novj"}, then clonal relationships are inferred using an adaptive 
#' threshold that indicates the level of similarity among junction sequences in a local neighborhood. 
#' 
#' If \code{method="vj"}, then clonal relationships are inferred not only on 
#' junction region homology, but also taking into account the mutation profiles in the V 
#' and J segments. Mutation counts are determined by comparing the input sequences (in the 
#' column specified by \code{sequence}) to the effective germline sequence (IUPAC representation 
#' of sequences in the column specified by \code{germline}). 
#' 
#' While not mandatory, the influence of SHM hot-/cold-spot biases in the clonal inference 
#' process will be noted if a SHM targeting model is provided through the \code{targeting_model} 
#' argument. See \link[shazam]{TargetingModel} for more technical details.
#' 
#' If the \code{threshold} argument is specified, then an upper limit for clonal grouping will 
#' be imposed to prevent sequences with dissimilarity above the threshold from grouping together. 
#' Any sequence with a distance greater than the \code{threshold} value from the other sequences, 
#' will be assigned to a singleton group.
#' 
#' @section Single-cell data:
#' To invoke single-cell mode the \code{cell_id} argument must be specified and the \code{locus} 
#' column must be correct. Otherwise, clustering will be performed with bulk sequencing assumptions, 
#' using all input sequences regardless of the values in the \code{locus} column.
#' 
#' Values in the \code{locus} column must be one of \code{c("IGH", "IGI", "IGK", "IGL")} for BCR 
#' or \code{c("TRA", "TRB", "TRD", "TRG")} for TCR sequences. Otherwise, the operation will exit and 
#' return and error message.
#' 
#' Under single-cell mode with paired-chain sequences, there is a choice of whether 
#' grouping should be done by (a) using IGH (BCR) or TRB/TRD (TCR) sequences only or
#' (b) using IGH plus IGK/IGL (BCR) or TRB/TRD plus TRA/TRG (TCR) sequences. 
#' This is governed by the \code{only_heavy} argument. There is also choice as to whether 
#' inferred clones should be split by the light/short chain (IGK, IGL, TRA, TRG) following 
#' heavy/long chain clustering, which is governed by the \code{split_light} argument.
#' 
#' In single-cell mode, clonal clustering will not be performed on data were cells are 
#' assigned multiple heavy/long chain sequences (IGH, TRB, TRD). If observed, the operation 
#' will exit and return an error message. Cells that lack a heavy/long chain sequence (i.e., cells with 
#' light/short chains only) will be assigned a \code{clone_id} of \code{NA}.
#'
#' @seealso
#' See \link{plotCloneSummary} for plotting summary results. See \link{groupGenes} for 
#' more details about grouping requirements.
#'
#' @examples
#' # Subset example data
#' db <- subset(ExampleDb, sample_id == "-1h")
#' 
#' # Find clonal groups
#' results <- spectralClones(db, method="novj", germline="germline_alignment_d_mask")
#' 
#' # Retrieve modified input data with clonal clustering identifiers
#' df <- as.data.frame(results)
#'   
#' # Plot clonal summaries 
#' plot(results, binwidth=0.02)
#' 
#' @export
spectralClones <- function(db, method=c("novj", "vj"), germline="germline_alignment", sequence="sequence_alignment",
                           junction="junction", v_call="v_call", j_call="j_call", clone="clone_id",
                           cell_id=NULL, locus="locus", only_heavy=TRUE, split_light=TRUE,
                           targeting_model=NULL, len_limit=NULL, first=FALSE, cdr3=FALSE, mod3=FALSE, max_n=0, 
                           threshold=NULL, base_sim=0.95, iter_max=1000,  nstart=1000, nproc=1,
                           verbose=FALSE, log=NULL,
                           summarize_clones=TRUE) {
    
    results <- defineClonesScoper(db = db, method = match.arg(method), model = "spectral", 
                                  germline = germline, sequence = sequence,
                                  junction = junction, v_call = v_call, j_call = j_call, clone = clone, 
                                  cell_id = cell_id, locus = locus, only_heavy = only_heavy, split_light = split_light,
                                  targeting_model = targeting_model, len_limit = len_limit,
                                  first = first, cdr3 = cdr3, mod3 = mod3, max_n = max_n,
                                  threshold = threshold, base_sim = base_sim,
                                  iter_max = iter_max, nstart = nstart, nproc = nproc,
                                  verbose = verbose, log = log,
                                  summarize_clones = summarize_clones)
    
    # return results
    if (summarize_clones) {
        return_list <- new("ScoperClones",
                           db = results$db,
                           vjl_groups = results$vjl_gps,
                           inter_intra = results$inter_intra,
                           eff_threshold = results$eff_threshold)    
        return(return_list)
    } else {
        return(results)
    }
}

#### defineClonesScoper ####

# *****************************************************************************
defineClonesScoper <- function(db,
                               model = c("identical", "hierarchical", "spectral"),
                               method = c("nt", "aa", "novj", "vj"),
                               linkage = c("single", "average", "complete"), normalize = c("len", "none"),
                               germline = "germline_alignment", sequence = "sequence_alignment",
                               junction = "junction", v_call = "v_call", j_call = "j_call", clone = "clone_id",
                               cell_id = NULL, locus = NULL, only_heavy = TRUE, split_light = TRUE,
                               targeting_model = NULL, len_limit = NULL,
                               first = FALSE, cdr3 = FALSE, mod3 = FALSE, max_n = 0, 
                               threshold = NULL, base_sim = 0.95,
                               iter_max = 1000, nstart = 1000, nproc = 1,
                               verbose = FALSE, log = NULL,
                               summarize_clones = TRUE) {

    ### get model
    model <- match.arg(model)
    
    ### get method
    method <- match.arg(method)

    # Initial checks
    if (!is.data.frame(db)) {
        stop("'db' must be a data frame")
    }
    
    ### check model andmethod
    if (model == "identical") {
        if (!(method %in% c("nt", "aa"))) {
            stop(paste0("'method' should be one of 'nt' or 'aa' for model '", model, "'.")) 
        }
    } else if (model == "hierarchical") {
        ### get normalize
        normalize <- match.arg(normalize)
        if (!normalize %in% c("len", "none")) { 
            stop(paste0("'normalize' should be one of 'len' or 'none for model '", model, "'.")) 
        }
        ### get linkage
        linkage <- match.arg(linkage)
        if (!linkage %in% c("single", "average", "complete")) { 
            stop(paste0("'linkage' should be one of 'single', 'average', or 'complete' for model '", model, "'.")) 
        }
        ### check threshold
        if (is.null(threshold) | threshold > 1) {
            stop(paste0("'threshold' should be a positive value less than 1 for model '", model, "'.")) 
        }
    } else if (model == "spectral") {
        if (!method %in% c("novj", "vj")) { 
            stop(paste0("'method' should be one of 'novj' or 'vj' for model '", model, "'.")) 
        }
    }  else {
        stop("model must be one of 'identical', 'hierarchical', or 'spectral'.")
    }
    
    ### Check for invalid characters
    valid_chars <- colnames(getDNAMatrix(gap = 0))
    .validateSeq <- function(x) { all(unique(strsplit(x, "")[[1]]) %in% valid_chars) }
    valid_seq <- sapply(db[[junction]], .validateSeq)
    not_valid_seq <- which(!valid_seq)
    if (length(not_valid_seq) > 0) {
        stop("invalid sequence characters in the ", junction, " column. ",
             length(not_valid_seq)," sequence(s) found.", "\n Valid characters are: '",  valid_chars, "'")
    }
    
    ### temp columns
    temp_cols <- c("vj_group", "vjl_group", "junction_l",  "cdr3_col", "clone_temp")
    
    ### check for invalid columns
    invalid_cols <- c(clone, temp_cols)
    if (any(invalid_cols %in% colnames(db))) {
        stop("Column(s) '", paste(invalid_cols[invalid_cols %in% colnames(db)], collapse = "', '"), "' already exist.",
             "\n Invalid column names are: '", paste(invalid_cols, collapse = "', '"), "'.")
    }
    
    ### Check general required columns
    columns <- c(junction, v_call, j_call) #, fields
    check <- checkColumns(db, columns)
    if (is.character(check)) { 
        stop(check)
    }
    
    ### Check required columns for method "vj"
    if (model == "spectral" & method == "vj") {
        columns <- c(germline, sequence) #, fields
        check <- checkColumns(db, columns)
        if (is.character(check)) { 
            stop(check)
        }
    } 
    
    ### Check single-cell mode
    single_cell <- FALSE
    if (!is.null(cell_id) & !is.null(locus)) {
        # Check required columns for single-cell mode
        columns <- c(cell_id, locus) #, fields
        check <- checkColumns(db, columns)
        if (check != TRUE) { stop(check) }
        
        # check locus column
        valid_loci <- c("IGH", "IGI", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
        check <- !all(unique(db[[locus]]) %in% valid_loci)
        if (check) {
            stop("The locus column contains invalid loci annotations.")
        }
        # check multiple heavy chains
        x <- sum(table(db[[cell_id]][db[[locus]] == "IGH"]) > 1)
        if (x > 0) {
            stop(paste(x, "cell(s) with multiple heavy chains found. One heavy chain per cell is expected."))
        }
        # check multiple beta chains
        x <- sum(table(db[[cell_id]][db[[locus]] == "TRB"]) > 1)
        if (x > 0) {
            stop(paste(x, "cell(s) with multiple beta chains found. One beta chain per cell is expected."))
        }
        # check multiple delta chains
        x <- sum(table(db[[cell_id]][db[[locus]] == "TRD"]) > 1)
        if (x > 0) {
            stop(paste(x, "cell(s) with multiple delta chains found. One delta chain per cell is expected."))
        }
        # if passed
        single_cell <- TRUE
    } 
    
    ### check verbose and log
    verbose <- ifelse(verbose, 1, 0)
    if (!is.null(log)) {
        out_dir <- dirname(log)
        if (!dir.exists(out_dir)) stop("out_dir '", out_dir, "' does not exist.")
        log_verbose <- 1
        log_verbose_name <- basename(log)
        cat(file=file.path(out_dir, log_verbose_name), append=FALSE)
    } else {
        log_verbose <- 0
    }

    ### Prepare db
    results_prep <- prepare_db(db = db, 
                               junction = junction, v_call = v_call, j_call = j_call,
                               first = first, cdr3 = cdr3, 
                               cell_id = cell_id, locus = locus, only_heavy = only_heavy,
                               mod3 = mod3, max_n = max_n)
    db <- results_prep$db
    n_rmv_mod3 <- results_prep$n_rmv_mod3
    n_rmv_cdr3 <- results_prep$n_rmv_cdr3
    n_rmv_N <- results_prep$n_rmv_N
    junction_l <- results_prep$junction_l
    cdr3_col <-  results_prep$cdr3_col
    
    ### for single-cell mode: separates heavy and light chain data frames
    ### performs cloning only on heavy chains
    if (single_cell) {
        db_l <- db[db[[locus]] %in% c("IGK", "IGL", "TRA", "TRG"), ]
        db_h <- db[db[[locus]] %in% c("IGH", "TRB", "TRD"), ]
        db <- db_h
    }
    
    ### summary of the groups
    vjl_gps <- db %>% 
        dplyr::group_by(!!rlang::sym("vjl_group")) %>% 
        dplyr::summarise(group_v_call = paste(unique(!!rlang::sym(v_call)), collapse=","),
                         group_j_call = paste(unique(!!rlang::sym(j_call)), collapse=","),
                         group_junction_length = unique(!!rlang::sym(junction_l)),
                         group_size = n())
    vjl_gps$group_v_call <- sapply(1:nrow(vjl_gps), 
                                      function(i){ paste(unique(strsplit(vjl_gps$group_v_call[i], split=",")[[1]]), collapse=",") })
    vjl_gps$group_j_call <- sapply(1:nrow(vjl_gps), 
                                      function(i){ paste(unique(strsplit(vjl_gps$group_j_call[i], split=",")[[1]]), collapse=",") })
    n_groups <- nrow(vjl_gps)
    
    ### Create cluster of nproc size and export namespaces
    if(nproc == 1) {
        # If needed to run on a single core/cpu then, register DoSEQ
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK", outfile = "")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    ### expoer function to clusters
    if (nproc > 1) { 
        export_functions <- list("passToClustering_lev1", "passToClustering_lev2", "passToClustering_lev3", "passToClustering_lev4",
                                 "findGapSmooth", "infer", "krnlMtxGenerator", "makeAffinity", "laplacianMtx", 
                                 "rangeAtoB", "likelihoods", "pairwiseMutions", "pairwiseMutMatrix",
                                 "printVerbose", "logVerbose")
        parallel::clusterExport(cluster, export_functions, envir=environment())
    }
    
    ### perform clustering for each group
    gp <- NULL
    db_cloned <- foreach(gp=1:n_groups,
                         .final=dplyr::bind_rows,
                         .inorder=TRUE,
                         .errorhandling='stop') %dopar% { 
                             # *********************************************************************************
                             # filter each group
                             vjl_gp <- vjl_gps$vjl_group[gp]
                             gp_vcall <- vjl_gps$group_v_call[gp]
                             gp_jcall <- vjl_gps$group_j_call[gp]
                             gp_lent <- vjl_gps$group_junction_length[gp]
                             gp_size <- vjl_gps$group_size[gp]
                             db_gp <- dplyr::filter(db, !!rlang::sym("vjl_group") == vjl_gp)
                             
                             # pass the group for clustering
                             # cat(paste(vjl_gp, "here"), sep="\n")  # for tests
                             results <- passToClustering_lev1(db_gp,
                                                              model = model,
                                                              method = method,
                                                              linkage = ifelse(model == "hierarchical", linkage, NA),
                                                              normalize = ifelse(model == "hierarchical", normalize, NA),
                                                              germline = germline,
                                                              sequence = sequence,
                                                              junction = junction,
                                                              mutabs = targeting_model,
                                                              len_limit = len_limit,
                                                              cdr3 = cdr3,
                                                              cdr3_col = cdr3_col,
                                                              threshold = threshold,
                                                              base_sim = base_sim,
                                                              iter_max = iter_max, 
                                                              nstart = nstart)
                             idCluster <- results$idCluster
                             n_cluster <- results$n_cluster
                             # cat(paste(vjl_gp), sep="\n")  # for tests
                             
                             if (length(idCluster) == 0 | any(is.na(idCluster))) {
                                 stop(printVerbose(n_groups, vjl_gp, model, method, linkage, cdr3,
                                                   gp_vcall, gp_jcall, gp_lent, gp_size, n_cluster) )  
                             } 
                             
                             # check verbose
                             if (verbose) { printVerbose(n_groups, vjl_gp, model, method, linkage, cdr3,
                                                         gp_vcall, gp_jcall, gp_lent, gp_size, n_cluster) 
                             }
                             
                             # check log verbose
                             if (log_verbose) { logVerbose(out_dir, log_verbose_name,
                                                           n_groups, vjl_gp, model, method, linkage, cdr3,
                                                           gp_vcall, gp_jcall, gp_lent, gp_size, n_cluster) }
                             
                             # attache clones
                             db_gp[[clone]] <- paste(vjl_gp, idCluster, sep="_")   
                             
                             # return result from each proc
                             return(db_gp)
                             # *********************************************************************************
                         }
    
    ### Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    ### sort clone ids
    db_cloned$clone_temp <- db_cloned %>%
        dplyr::group_by(!!rlang::sym(clone)) %>%
        dplyr::group_indices()
    db_cloned[[clone]] <- db_cloned$clone_temp
    db_cloned <- db_cloned[order(db_cloned[[clone]]), ]
    db_cloned[[clone]] <- as.character(db_cloned[[clone]])
    db_cloned$clone_temp <- NULL
    
    ### report removed sequences
    if (mod3) {
        if (verbose) {
            cat("      MOD3_FILTER> ", n_rmv_mod3, "invalid junction length(s) (not mod3) in the", junction, "column removed.", "\n", sep=" ")   
        }
        if (log_verbose)  { 
            cat("      MOD3_FILTER> ", n_rmv_mod3, "invalid junction length(s) (not mod3) in the", junction, "column removed.", "\n", sep=" ",
                file = file.path(out_dir, log_verbose_name), append=TRUE) 
        }
    }
    if (cdr3) {
        if (verbose) {
            cat("      CDR3_FILTER> ", n_rmv_cdr3, "invalid junction length(s) (< 7) in the", junction, "column removed.", "\n", sep=" ")   
        }
        if (log_verbose)  { 
            cat("      CDR3_FILTER> ", n_rmv_cdr3, "invalid junction length(s) (< 7) in the", junction, "column removed.", "\n", sep=" ",
                file = file.path(out_dir, log_verbose_name), append=TRUE) 
        }
    }
    if (!is.null(max_n)) {
        if (verbose) {
            cat("     MAX_N_FILTER> ", n_rmv_N, "invalid junction(s) ( # of N >", max_n, ") in the", junction, "column removed.", "\n", sep=" ")   
        }
        if (log_verbose)  { 
            cat("     MAX_N_FILTER> ", n_rmv_N, "invalid junction(s) ( # of N >", max_n, ") in the", junction, "column removed.", "\n", sep=" ",
                file = file.path(out_dir, log_verbose_name), append=TRUE) 
        }
    }
    
    ### make summary 
    if (summarize_clones) {
        ### vjl group summary
        vjl_gps <- db_cloned %>%
            dplyr::group_by(!!rlang::sym("vjl_group")) %>%
            dplyr::summarise(sequence_count = n(),
                             v_call = paste(unique(!!rlang::sym(v_call)), collapse=","),
                             j_call = paste(unique(!!rlang::sym(j_call)), collapse=","),
                             junction_length = unique(!!rlang::sym(junction_l)),
                             clone_count = length(unique(!!rlang::sym(clone))),
                             clone_id = paste(unique(!!rlang::sym(clone)), collapse = ","))
        vjl_gps$v_call <- sapply(1:nrow(vjl_gps),
                                    function(i){ paste(unique(strsplit(vjl_gps$v_call[i], split=",")[[1]]), collapse=",") })
        vjl_gps$j_call <- sapply(1:nrow(vjl_gps),
                                    function(i){ paste(unique(strsplit(vjl_gps$j_call[i], split=",")[[1]]), collapse=",") })
        
        ### calculate inter and intra distances
        df_inter_intra <- calculateInterVsIntra(db = db_cloned,
                                                clone = clone,
                                                vjl_gps = vjl_gps,
                                                junction = junction,
                                                cdr3 = cdr3,
                                                cdr3_col = cdr3_col,
                                                nproc = nproc,
                                                verbose = verbose)
        
        ### calculate effective threshold
        data_eff <- select(df_inter_intra, c("distance", "label"))
        data_intra <- data_eff %>% 
            dplyr::filter(!!rlang::sym("label") == "intra", !!rlang::sym("distance") > 0)
        data_inter <- data_eff %>%
            dplyr::filter(!!rlang::sym("label") == "inter", !!rlang::sym("distance") > 0)
        
        eff_threshold <- as.numeric(NA)
        if (nrow(data_intra) > 5 & nrow(data_inter) > 5) {
            a <- data_intra$distance
            b <- data_inter$distance
            xlim = c(min(c(a,b)), max(c(a,b)))
            df <- merge(
                as.data.frame(density(a, from = xlim[1], to = xlim[2])[c("x", "y")]),
                as.data.frame(density(b, from = xlim[1], to = xlim[2])[c("x", "y")]),
                by = "x", suffixes = c(".a", ".b")
            )
            df$comp <- as.numeric(df$y.a > df$y.b)
            df$cross <- c(NA, diff(df$comp))
            df <- df[which(df$cross != 0), c("x", "y.a")]
            if (nrow(df) > 0) {
                eff_th <- df$x
                eff_th <- eff_th[mean(a) - sd(a) < eff_th & eff_th < mean(b) + sd(b)]
                if (length(eff_th) > 0) {
                    eff_threshold <- round(mean(eff_th), 2)
                } 
            }    
        }
    }
    
    ### remove extra columns
    db_cloned <- db_cloned[, !(names(db_cloned) %in% temp_cols)]
    
    ### singl cell pipeline
    if (single_cell) {
        db_l <- db_l[, !(names(db_l) %in% temp_cols)]
        db_l[[clone]] <- NA
        # copy clone ids from heavy chains into light chains
        cell_ids_h <- unique(db_cloned[[cell_id]])
        cell_ids_l <- unique(db_l[[cell_id]])
        for (cellid in cell_ids_l) {
            if (cellid %in% cell_ids_h) {
                db_l[[clone]][db_l[[cell_id]] == cellid] <- db_cloned[[clone]][db_cloned[[cell_id]] == cellid]
            } 
        }
        # bind heavy and light chain data.frames
        stopifnot(all(names(db_cloned) == names(db_l)))
        db_cloned <- bind_rows(db_cloned, db_l)
        # split clones by light chains
        if (split_light) {
            clones <- unique(db_cloned[[clone]])
            clones <- clones[!is.na(clones)]
            for (cloneid in clones) {
                db_c <- dplyr::filter(db_cloned, !!rlang::sym(clone) == cloneid)
                if (length(unique(db_c[[cell_id]])) == 1) next()
                db_c <- groupGenes(data = db_c,
                                   v_call = v_call,
                                   j_call = j_call,
                                   junc_len = NULL,
                                   cell_id = cell_id,
                                   locus = locus,
                                   only_heavy = FALSE,
                                   first = FALSE)
                if (length(unique(db_c$vj_group)) == 1) next()
                db_c[[clone]] <- paste(db_c[[clone]], db_c$vj_group, sep="_") 
                for (cellid in unique(db_c[[cell_id]])) {
                    db_cloned[[clone]][db_cloned[[clone]] == cloneid & db_cloned[[cell_id]] == cellid] <- 
                        db_c[[clone]][db_c[[cell_id]] == cellid]
                }
            }
        }
        # sort clone ids
        na.count <- sum(is.na(db_cloned[[clone]]))
        if (na.count > 0) {
            db_na <- db_cloned[is.na(db_cloned[[clone]]), ]
            db_cloned <- db_cloned[!is.na(db_cloned[[clone]]), ]        
        }
        db_cloned$clone_temp <- db_cloned %>%
            dplyr::group_by(!!rlang::sym(clone)) %>%
            dplyr::group_indices()
        db_cloned[[clone]] <- db_cloned$clone_temp
        db_cloned <- db_cloned[order(db_cloned[[clone]]), ]
        db_cloned[[clone]] <- as.character(db_cloned[[clone]])
        db_cloned$clone_temp <- NULL
        if (na.count > 0) {
            db_cloned <- bind_rows(db_cloned, db_na)
        }
    }
    
    # return results
    if (summarize_clones) {
        return_list <- list("db" = db_cloned,
                            "vjl_gps" = vjl_gps,
                            "inter_intra" = df_inter_intra,
                            "eff_threshold" = eff_threshold)  
        
        return(return_list)
    } else {
        return(db_cloned)
    }
}
# *****************************************************************************

# *****************************************************************************
passToClustering_lev1 <- function (db_gp, 
                                   model = c("identical", "hierarchical", "spectral"),
                                   method = c("nt", "aa", "novj", "vj"),
                                   linkage = c("single", "average", "complete"),
                                   normalize = c("len", "none"),
                                   germline = "germline_alignment",
                                   sequence = "sequence_alignment",
                                   junction = "junction",
                                   mutabs = NULL,
                                   len_limit = NULL,
                                   cdr3 = FALSE,
                                   cdr3_col = NA,
                                   threshold = NULL,
                                   base_sim = 0.95,
                                   iter_max = 1000, 
                                   nstart = 1000) {
    ### get model
    model <- match.arg(model)
    
    ### begin clustering
    if (model == "identical") {
        clone_results <- identicalClones_helper(db_gp,
                                                method = method,
                                                junction = junction,
                                                cdr3 = cdr3,
                                                cdr3_col = cdr3_col)
    } else if (model == "hierarchical") {
        clone_results <- hierarchicalClones_helper(db_gp,
                                                   method = method,
                                                   linkage = linkage,
                                                   normalize = normalize,
                                                   junction = junction,
                                                   cdr3 = cdr3,
                                                   cdr3_col = cdr3_col,
                                                   threshold = threshold)
    } else if (model == "spectral") {
        clone_results <- spectralClones_helper(db_gp,
                                               method = method,
                                               germline = germline,
                                               sequence = sequence,
                                               junction = junction,
                                               mutabs = mutabs,
                                               len_limit = len_limit,
                                               cdr3 = cdr3,
                                               cdr3_col = cdr3_col,
                                               threshold = threshold,
                                               base_sim = base_sim,
                                               iter_max = iter_max, 
                                               nstart = nstart)
    }
    
    ### retrun results
    return_list <- list("idCluster" = clone_results$idCluster, 
                        "n_cluster" = clone_results$n_cluster, 
                        "eigen_vals" = clone_results$eigen_vals)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
identicalClones_helper <- function(db_gp,
                                   method = c("nt", "aa"),
                                   junction = "junction",
                                   cdr3 = FALSE,
                                   cdr3_col = NA) {
    ### get method
    method <- match.arg(method)
    
    ### number of sequences
    n <- nrow(db_gp)
    
    ### cloning
    seq_col <- ifelse(cdr3, cdr3_col, junction)
    if (method == "aa") {
        db_gp[[seq_col]] <- translateDNA(db_gp[[seq_col]])
    }
    idCluster <- db_gp %>% 
        dplyr::group_by(!!rlang::sym(seq_col)) %>% 
        group_indices()
    n_cluster <- length(unique(idCluster))
    eigen_vals <- rep(0, n)
    
    ### retrun results
    return_list <- list("idCluster" = idCluster, 
                        "n_cluster" = n_cluster, 
                        "eigen_vals" = eigen_vals)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
hierarchicalClones_helper <- function(db_gp,
                                      method = c("nt", "aa"),
                                      linkage = c("single", "average", "complete"),
                                      normalize = c("len", "none"),
                                      junction = "junction",
                                      cdr3 = FALSE,
                                      cdr3_col = NA,
                                      threshold = NULL) {
    ### get method
    method <- match.arg(method)

    # get linkage
    linkage <- match.arg(linkage)
    
    # get normalize
    normalize <- match.arg(normalize)
    
    ### number of sequences
    n <- nrow(db_gp)
    
    # get sequences
    if (method == "nt") {
        seqs <- db_gp[[ifelse(cdr3, cdr3_col, junction)]]   
    } else if (method == "aa") {
        # translate amino acid for method "aa"
        seqs <- translateDNA(db_gp[[ifelse(cdr3, cdr3_col, junction)]])
    }
    
    # find unique seqs
    df <- as.data.table(seqs)[, list(list(.I)), by = seqs]
    n_unq <- nrow(df)
    ind_unq <- df$V1
    seqs_unq <- df$seqs
    if (n_unq == 1) {
        return(list("idCluster" = rep(1, n), 
                    "n_cluster" = 1, 
                    "eigen_vals" = rep(0, n)))
    }
    
    # calculate distance matrix
    if (method == "nt") {
        dist_mtx <- pairwiseDist(seq = seqs_unq, 
                                 dist_mat = getDNAMatrix(gap = 0))
    } else if (method == "aa") {
        dist_mtx <- pairwiseDist(seq = seqs_unq, 
                                 dist_mat = getAAMatrix(gap = 0))
    }
    
    # perform hierarchical clustering
    if (normalize == "len") {
        # calculate normalization factor
        junc_length <- unique(stri_length(seqs_unq))
        hc <- hclust(as.dist(dist_mtx/junc_length), method = linkage)    
    } else if (normalize == "none") {
        hc <- hclust(as.dist(dist_mtx), method = linkage)    
    }
    
    # cut the tree
    idCluster_unq <- cutree(hc, h = threshold)
    
    # back to reality
    idCluster <- rep(NA, n)
    for (i in 1:n_unq) {
        idCluster[ind_unq[[i]]] <- idCluster_unq[i]
    }
    n_cluster <- length(unique(idCluster))
    eigen_vals <- rep(0, n)
    
    ### retrun results
    return_list <- list("idCluster" = idCluster, 
                        "n_cluster" = n_cluster, 
                        "eigen_vals" = eigen_vals)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
spectralClones_helper <- function(db_gp,
                                  method = c("novj", "vj"),
                                  germline = "germline_alignment",
                                  sequence = "sequence_alignment",
                                  junction = "junction",
                                  mutabs = NULL,
                                  len_limit = NULL,
                                  cdr3 = FALSE,
                                  cdr3_col = NA,
                                  threshold = NULL,
                                  base_sim = 0.95,
                                  iter_max = 1000, 
                                  nstart = 1000) {
    
    ### get method
    method <- match.arg(method)
    
    ### number of sequences
    n <- nrow(db_gp)
    
    ### cloning
    if (method == "vj") {
        ### check targeting model
        if (!is.null(mutabs)) { 
            mutabs <- mutabs@mutability 
        } else {
            mutabs <- NULL
        }
        # get required info based on the method
        germs <- db_gp[[germline]]
        seqs <- db_gp[[sequence]]
        juncs <- db_gp[[ifelse(cdr3, cdr3_col, junction)]]
        junc_length <- unique(stri_length(juncs))
        # find unique seqs
        seqs <- paste(seqs, juncs, germs, sep = "|")
        df <- as.data.table(seqs)[, list(list(.I)), by=seqs] %>%
            tidyr::separate(col = seqs, into = c("seqs_unq", "juncs_unq", "germs_unq"), sep = "\\|")
        n_unq <- nrow(df)
        ind_unq <- df$V1
        if (n_unq == 1) {
            return(list("idCluster" = rep(1, n), 
                        "n_cluster" = 1, 
                        "eigen_vals" = rep(0, n)))
        }
        # find corresponding unique germs and junctions
        seqs_unq <- df$seqs_unq
        germs_unq <- df$germs_unq
        juncs_unq <- df$juncs_unq
        # calculate unique junctions distance matrix
        dist_mtx <- pairwiseDist(seq = juncs_unq, 
                                 dist_mat = getDNAMatrix(gap = 0))
        # count mutations from unique sequence imgt
        results <- pairwiseMutions(germ_imgt = germs_unq, 
                                   seq_imgt = seqs_unq,
                                   junc_length = junc_length, 
                                   len_limit = len_limit,
                                   cdr3 = cdr3,
                                   mutabs = mutabs)
        tot_mtx <- results$pairWiseTotalMut
        sh_mtx <- results$pairWiseSharedMut
        mutab_mtx <- results$pairWiseMutability
        # calculate likelihhod matrix
        lkl_mtx <- likelihoods(tot_mtx = tot_mtx, 
                               sh_mtx = sh_mtx, 
                               mutab_mtx = mutab_mtx)
        # calculate weighted matrix  
        disim_mtx <- dist_mtx * (1.0 - lkl_mtx)
        # check if vj method made any changes, otherwise go back to method "novj"
        # check if one of the rows is all zeros, meaning vj mathod cannot decide (a highly rare case)
        if (all(disim_mtx == dist_mtx) | any(rowSums(disim_mtx) == 0)) {
            # get required info based on the method
            seqs <- db_gp[[ifelse(cdr3, cdr3_col, junction)]]
            junc_length <- unique(stri_length(seqs))
            # find unique seqs
            df <- as.data.table(seqs)[, list(list(.I)), by = seqs]
            n_unq <- nrow(df)
            ind_unq <- df$V1
            seqs_unq <- df$seqs
            if (n_unq == 1) {
                return(list("idCluster" = rep(1, n), 
                            "n_cluster" = 1, 
                            "eigen_vals" = rep(0, n)))
            }
            # calculate unique seuences distance matrix
            disim_mtx <- pairwiseDist(seq = seqs_unq, 
                                      dist_mat = getDNAMatrix(gap = 0))
        } 
    } else if (method == "novj") {
        # get required info based on the method
        seqs <- db_gp[[ifelse(cdr3, cdr3_col, junction)]]
        junc_length <- unique(stri_length(seqs))
        # find unique seqs
        df <- as.data.table(seqs)[, list(list(.I)), by = seqs]
        n_unq <- nrow(df)
        ind_unq <- df$V1
        seqs_unq <- df$seqs
        if (n_unq == 1) {
            return(list("idCluster" = rep(1, n), 
                        "n_cluster" = 1, 
                        "eigen_vals" = rep(0, n)))
        }
        # calculate unique seuences distance matrix
        disim_mtx <- pairwiseDist(seq = seqs_unq, 
                                  dist_mat = getDNAMatrix(gap = 0))
    }
    ### pass to clustering pipeline
    result <- passToClustering_lev2(mtx = disim_mtx, 
                                    junc_length = junc_length,
                                    threshold = threshold, 
                                    base_sim = base_sim, 
                                    iter_max = iter_max, 
                                    nstart = nstart)
    idCluster_unq <- result$idCluster
    eigen_vals <- result$eigen_vals
    ### back to reality
    idCluster <- rep(NA, n)
    for (i in 1:n_unq) {
        idCluster[ind_unq[[i]]] <- idCluster_unq[i]
    }
    n_cluster <- length(unique(idCluster))
    
    ### retrun results
    return_list <- list("idCluster" = idCluster, 
                        "n_cluster" = n_cluster, 
                        "eigen_vals" = eigen_vals)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
# check special case if all distances (off diagonal elements) are the same. 
# (this also includes cases with only two sequences)
passToClustering_lev2 <- function(mtx, 
                                  junc_length = NULL, 
                                  threshold = NULL, 
                                  base_sim = 0.95, 
                                  iter_max = 1000, 
                                  nstart = 1000) {
    ### constants
    n <- nrow(mtx)
    bs <- (1 - base_sim)*junc_length
    off_diags_nuq <- unique(mtx[row(mtx) != col(mtx)])
    ### check special cases
    if (n == 1) {
        idCluster <- 1    # all in same clone
        eigen_vals <- 0
    } else if (length(off_diags_nuq) == 1) {   # seqs have equal-distances from each other
        if (off_diags_nuq > bs) { 
            idCluster <- 1:n         # all singletons
            eigen_vals <- rep(0, n)
        } else {
            idCluster <- rep(1, n)    # all in same clone
            eigen_vals <- rep(0, n)
        }
    } else if (max(mtx) <= bs) {
        idCluster <- rep(1, n)    # all in same clone
        eigen_vals <- rep(0, n)
    } else {
        results <- passToClustering_lev3(mtx = mtx, 
                                         junc_length = junc_length,
                                         threshold = threshold, 
                                         iter_max = iter_max, 
                                         nstart = nstart)
        idCluster <- results$idCluster
        eigen_vals <- results$eigen_vals
    }
    ### return list
    return_list <- list("idCluster" = idCluster, 
                        "eigen_vals" = eigen_vals)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
passToClustering_lev3 <- function(mtx, 
                                  junc_length = NULL, 
                                  threshold = NULL, 
                                  iter_max = 1000, 
                                  nstart = 1000){
    n <- nrow(mtx)
    ### set seed for reproducibility
    ### check the minimum number of data points requirement
    if (n < 3) stop("SCOPer needs at least 3 unique data points")
    ### calculate the krnl matrix
    krnl_mtx <- krnlMtxGenerator(mtx = mtx)
    ### calculate the affinity matrix
    if (!is.null(threshold)) {
        aff_mtx <- makeAffinity(mtx_o = mtx, 
                                mtx_k = krnl_mtx,
                                thd = threshold*junc_length) 
    } else {
        # aff_mtx <- round(krnl_mtx, 3)
        nearest_dist <- apply(mtx, 2,  function(x) {
            gt0 <- which(x > 0)
            if (length(gt0) != 0) { min(x[gt0]) } else { NA }
        })
        aff_mtx <- makeAffinity(mtx_o = mtx,
                                mtx_k = krnl_mtx,
                                thd = max(nearest_dist, na.rm = T))
    }
    ### affinity matrix is diagonal. Each sequence belongs to a singlton clone.
    if (all(aff_mtx[!diag(nrow(aff_mtx))] == 0)) { 
        return_list <- list("idCluster" = c(1:n),
                            "eigen_vals" = rep(0,n))
        return(return_list)
    }
    ### calculate the laplacian matrix
    L <- laplacianMtx(entry = aff_mtx)
    ### eigen-decomposition
    eigens <- eigen(L, symmetric = TRUE)
    eigen_vals <- rev(eigens$values)
    eigen_vals <- rangeAtoB(eigen_vals, 0, 1) 
    if (all(eigen_vals == 0)) {  # Each sequence belongs to a singlton clone.
        return_list <- list("idCluster" = c(1:n),
                            "eigen_vals" = rep(0,n))
        return(return_list)
    }
    ### k upper bound
    eigenValDens <- density(eigen_vals)
    k_up <- sum(eigen_vals < eigenValDens$x[which.max(eigenValDens$y)]) + 1
    if (is.na(k_up) | k_up %in% c(1,2)) k_up <- n
    # determine the number of clusters using log distance
    logEigenValsDiff <- sapply(2:k_up, function(x){ log10(eigen_vals[x]) - log10(eigen_vals[x-1]) })
    nas <- sum(is.nan(logEigenValsDiff) | is.na(logEigenValsDiff) | is.infinite(logEigenValsDiff))
    if (nas == length(logEigenValsDiff)) {
        k <- nas + 1
    } else {
        k <- nas + ifelse(nas > 0, which.max(logEigenValsDiff[-(1:nas)]), which.max(logEigenValsDiff)) + 1
    }
    
    if (k == n) { 
        idCluster <- 1:n         # all singletons
        eigen_vals <- rep(0, n)
    } else {
        ### pick k smalest eigenvectors
        # The vectors are col-normalized to unit length
        eigenVecs <- eigens$vectors
        eigenVecs <- eigenVecs[, (n-k+1):(n)]
        ### kmeans clustering
        # set.seed(12345)
        set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
        idCluster <- kmeans(x = round(eigenVecs, 6), 
                            centers = k, 
                            iter.max = iter_max, 
                            nstart = nstart)$cluster
        ### check if idclusters and affinity matrix agrees
        idCluster <- passToClustering_lev4(aff_mtx = aff_mtx, 
                                           idCluster = idCluster)
    }
    
    ### return results
    return_list <- list("idCluster" = idCluster,
                        "eigen_vals" = eigen_vals)  
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
# check affinity matrix and clusters id 
# aff_mtx_sub[which(aff_id_sub %in% gr_ls[[4]]),]
passToClustering_lev4 <- function(aff_mtx, idCluster) {
    n <- nrow(aff_mtx)
    new_idCluster <- rep(NA, n)
    new_k <- 0
    id_unq <- unique(idCluster)
    k <- length(unique(idCluster))
    for (z in 1:k) {
        aff_id_sub <- which(idCluster == id_unq[z])
        if (length(aff_id_sub) == 1) {
            new_k <- new_k + 1
            new_idCluster[aff_id_sub] <- new_k    
        } else {
            aff_mtx_sub <- aff_mtx[aff_id_sub, aff_id_sub]
            n <- nrow(aff_mtx_sub)
            rows_ls <- rep(list(NULL), n)
            for (i in 1:n) {
                rows_ls[[i]] <- aff_id_sub[aff_mtx_sub[, i] != 0]
            }
            gr_ls <- rep(list(NA), n)
            for (y in 1:n) {
                ids <- rows_ls[[y]]
                l <- length(gr_ls[!is.na(gr_ls)])
                for (x in 1:l) {
                    if (1 > l) break
                    if (length(base::intersect(ids, gr_ls[[x]])) > 0) {
                        ids <- base::union(ids, gr_ls[[x]])
                        gr_ls[[x]] <- NA
                    }
                }
                gr_ls[[y]] <- ids
                gr_ls <- c(gr_ls[!is.na(gr_ls)], gr_ls[is.na(gr_ls)])
            }
            gr_ls <- gr_ls[!is.na(gr_ls)]
            for (i in 1:length(gr_ls)) {
                new_k <- new_k + 1
                new_idCluster[gr_ls[[i]]] <- new_k    
            }
        }
    }
    return(new_idCluster)
}
# TO CHECK USE aff_mtx_sub[, which(aff_id_sub == ???)]
