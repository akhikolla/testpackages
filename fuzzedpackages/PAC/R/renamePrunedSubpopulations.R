
#' Prune away specified subpopulations in clades that are far away.
#' @param pacman_results PAC-MAN analysis result matrix that contains network annotation, clade IDs and mean (centroid) clade expression levels.
#' @param subpopulationsToPrune A vector of clade IDs; these clades will be pruned.

#' @return Returns PAC-MAN analysis result matrix with pruned clades. The pruning process creates new clades to replace the original clade ID of the specified subpopulations.
#' @export
#' 
#' 
renamePrunedSubpopulations<-function(pacman_results,subpopulationsToPrune){
  pacman_results_cladePruned<-pacman_results
  pacman_results_cladePruned[,2]<-as.character(pacman_results_cladePruned[,2])
  cladeCount<-length(grep("clade", unique(pacman_results[,2])))
  
  for(i in 1:nrow(subpopulationsToPrune)){
    pacman_results_cladePruned[which(pacman_results_cladePruned[,2]==subpopulationsToPrune[i,1] & pacman_results_cladePruned[,3]==subpopulationsToPrune[i,2]),2]<-paste0("clade", cladeCount+1)
    cladeCount<-cladeCount+1
  }
  return(pacman_results_cladePruned)
}
