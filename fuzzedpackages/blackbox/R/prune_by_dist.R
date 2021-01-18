## potential arternative to greedyMAXMINwithFixed()
# If distance threshold given, computes clusters separated by larger distances; 
# If no threshold given, tries to adjust threshold to obtain (finalSize-length(fixedRows)) clusters to sample from;
# Next, sample once from each cluster;
# Finally, sample exactly (finalSize-length(fixedRows)) rows to remove any excess one. 

prune_by_dist <- function(rownamedarray,
                          finalSize, ## including fixed rows
                          scales=NULL,
                          threshold=NULL,
                          fixedRows=NULL ## any row indices or rownames from rownamedarray
                          ) {
  if ( ! requireNamespace("igraph",quietly=TRUE)) {
    stop("Package 'igraph' is not installed, hence prune_by_dist() cannot be used.")
  }
  if ( ! is.null(fixedRows)) {
    allrownames <-  rownames(rownamedarray)
    if ( ! is.character(fixedRows)) fixedRows <- allrownames[fixedRows]
    issampleRow <- ( ! (allrownames %in% fixedRows))
  }
  sampleSize <- finalSize -length(fixedRows)
  if (sampleSize<0L) {stop("Incorrect arguments such that sampleSize < length(fixedRows)")}
  if (sampleSize==0L) {
    message.redef("Suspect arguments such that sampleSize = length(fixedRows)")
    return(fixedRows)
  }
  if (!is.null(scales)) rownamedarray <- sweep(rownamedarray,2L,sqrt(scales),FUN=`/`) # t(t(rownamedarray)/sqrt(scales))
  nr <- nrow(rownamedarray)
  distContainer <- proxy::dist(rownamedarray)
  #
  adjmat <- matrix(TRUE,nrow=nrow(rownamedarray),ncol=nrow(rownamedarray))
  distrange <- range(distContainer[])
  clu <- function(threshold) {
    adjvec <- (distContainer[] < threshold) ## adjvec is a boolean matrix
    adjmat[lower.tri(adjmat)] <- adjvec 
    adjmat[upper.tri(adjmat)] <- adjvec 
    ess <- igraph::graph_from_adjacency_matrix(adjmat,mode="undirected") ## network::network ?
    igraph::clusters(ess) ## clusters separated by distance > threshold 
    ## sna::clique.census(mids_1993, mode="graph", tabulate.by.vertex=FALSE, enumerate=FALSE) ?    
  }
  get_sample_clu_nb <- function(clustersObject) {
    if (is.null(fixedRows)) {
      return(clustersObject$no)
    } else {
      ## find clusters without fixed points:
      nofixed <- setdiff(names(table(clustersObject$membership)),names(table(clustersObject$membership[ ! issampleRow])))
      return(length(nofixed))
    }
  }
  goodrows <- NULL
  ## determine Th_clu(threshold) or else a trivial value of goodrows
  if (is.null(threshold)) {
    lowTh <- distrange[1L] + .Machine$double.eps
    upTh <- distrange[2L] - .Machine$double.eps
    lowTh_clu <- clu(lowTh)
    if (get_sample_clu_nb(lowTh_clu)<sampleSize) {
      message.redef("All distances appear similar: will use random sampling")
      goodrows <- rownames(rownamedarray)
    } else {
      upTh_clu <- clu(upTh)
      if (get_sample_clu_nb(upTh_clu)>sampleSize) {
        message.redef("All distances appear similar: will use random sampling")
        goodrows <- rownames(rownamedarray)
      }
    }
    #
    if (is.null(goodrows)) { ## non trivial threshold needs to be found
      lowQ <- 0
      upQ <- 1
      sample_clu_nb <- -1
      prev_sample_clu_nb <- -2
      while (TRUE) {
        if (sample_clu_nb==prev_sample_clu_nb) {
          thresholdQ <- (upQ+lowQ)/2
        } else {
          lowTh_clu_nb <- get_sample_clu_nb(lowTh_clu)
          thresholdQ <- lowQ+ (upQ-lowQ)*(lowTh_clu_nb-sampleSize)/(lowTh_clu_nb-get_sample_clu_nb(upTh_clu))
        }
        threshold <- quantile(distContainer[],thresholdQ)
        Th_clu <- clu(threshold)
        prev_sample_clu_nb <- sample_clu_nb
        sample_clu_nb <- get_sample_clu_nb(Th_clu)
        if (sample_clu_nb==sampleSize) {
          break
        } else if (sample_clu_nb>sampleSize) {
          lowQ <- thresholdQ
          lowTh_clu <- Th_clu
        } else if (sample_clu_nb<sampleSize) {
          upQ <- thresholdQ
          upTh_clu <- Th_clu
        } 
        if (upQ-lowQ<.Machine$double.eps) {
          Th_clu <- lowTh_clu
          break
        }
        #print(c(lowQ,upQ,sample_clu_nb,threshold))
      }
    } 
  } else Th_clu <- clu(threshold)
  #
  if (is.null(goodrows)) { # the uses Th_clu
    distContainer <- as.matrix(distContainer)
    findMostDistant <- function(id) {
      inrows <- (Th_clu$membership==id) ## $membership does not have rownames
      if ( ! is.null(fixedRows)) inrows <- (inrows & issampleRow)
      outrows <- (Th_clu$membership!=id)
      ### Will select point most distant to out-cluster fixed points: 
      # if ( ! is.null(fixedRows)) outrows <- (outrows | (! issampleRow) ) 
      subdistmat <- distContainer[inrows,outrows,drop=FALSE]
      rowmins <- apply(subdistmat,1L,min) 
      maxmin <- max(rowmins)
      mostdistantrow <- names(which(rowmins==maxmin))
      ##crude tie-breaking:
      if (length(mostdistantrow)>1L) mostdistantrow <- sample(mostdistantrow,size=1L)
      return(mostdistantrow)
    }
    if (is.null(fixedRows)) {
      goodrows <- sapply(seq(get_sample_clu_nb(Th_clu)), findMostDistant)
    } else {
      ## find clusters without fixed points:
      nofixed <- setdiff(names(table(Th_clu$membership)),names(table(Th_clu$membership[ ! issampleRow])))
      goodrows <- sapply(as.integer(nofixed), findMostDistant)
    }
  }
  if (length(goodrows)>sampleSize) goodrows <- sample(goodrows,size=sampleSize)
  goodrows <- c(goodrows,fixedRows)
  return(goodrows) ## rownames
}