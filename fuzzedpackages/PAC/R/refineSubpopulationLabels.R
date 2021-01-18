#' Refines the subpopulation labels from PAC using network alignment and small subpopulation information. Outputs a new set of files containing the representative labels.
#' 
#' 
#' @param sampleIDs sampleID vector
#' @param dim_subset a string vector of string names to subset the data columns for PAC; set to NULL to use all columns
#' @param clades_network_only the alignment results from MAN; used to translate the original sample-specific labels into clade labels
#' @param expressionGroupClamp clamps the subpopulations into desired number of expression groups for assigning small subpopulations into larger groups or their own groups.
#' @export

refineSubpopulationLabels<-function(sampleIDs, dim_subset, clades_network_only, expressionGroupClamp){
  
  subpopulationLabels<-NULL
  
  newlabels<-matrix(0,0,1)
  
  for(s in 1:length(sampleIDs)){
    #print(s)
    sampleID<-sampleIDs[s]
    load(paste0(sampleID,"_subpopulationLabels.Rdata"))
    
    new_subpopulationLabels<-subpopulationLabels
    
    #cluster bigger subpopulations by networks
    for(i in 1:length(clades_network_only)){
      newlabel<-unique(clades_network_only[[i]][grep(sampleID, names(clades_network_only[[i]]))])
      
      oldsubpopname<-names(clades_network_only[[i]][grep(sampleID, names(clades_network_only[[i]]))])
      oldlabel<-gsub("^.*\\_","", oldsubpopname)
      
      oldlabel<-as.numeric(oldlabel)
      
      for(j in 1:length(oldlabel)){
        index<-which(subpopulationLabels == oldlabel[j])
        new_subpopulationLabels[index]<-paste0("clade", newlabel)
      }
    }
    
    #get signal or expression values
    load(paste0(sampleID,"_dataMatrix.Rdata"))
    
    #subset for dimensions involved in clustering of subpopulations
    if(!is.null(dim_subset)){
      inputMatrix_withSampleName<-cbind(inputMatrix_withSampleName[,1], inputMatrix_withSampleName[,dim_subset])
    }
    
    inputMatrix_withSampleName<-as.data.frame(inputMatrix_withSampleName)
    
    data_agg<-aggregateData(inputMatrix_withSampleName,new_subpopulationLabels)
    
    expValues<-data_agg[,3:(ncol(data_agg)-1)]
    data_agg_names<-paste0(data_agg[,2],"_", data_agg[,1])
    rownames(expValues)<-data_agg_names
    
    #cluster by expression should be keeping previous labels for 'established' branches as they contain the clade information
    
    hc<-hclust(dist(expValues, method="euclidean"), "ave")
    groupB<-cutree(hc, k = expressionGroupClamp)
    clades_avg<-split(groupB, groupB) 
    
    newer_subpopulationLabels<-new_subpopulationLabels   
    e_index<-1
    
    for(e in 1:length(clades_avg)){
      newlabel<-unique(clades_avg[[e]][grep(sampleID, names(clades_avg[[e]]))])
      oldsubpopname<-names(clades_avg[[e]][grep(sampleID, names(clades_avg[[e]]))])
      oldlabel<-gsub("^.*\\_","", oldsubpopname)
      
      if(sum(grepl("clade", oldlabel))==0){
        #if there is no clade subpopulation, establish the expression-sample-specific clade using the non-clade subpopulations
        
        newlabel<-paste0("esc_", s, "_", e_index)        
        for(b in 1:length(oldlabel)){
          index<-which(newer_subpopulationLabels == oldlabel[b])
          newer_subpopulationLabels[index]<-newlabel
        }
        e_index<-e_index+1
        
      } else if (sum(grepl("clade", oldlabel))==1){
        #if there is only one clade subpopulation in the expression clade, merge the non-clade subpopulations in this expression clade into the clade subpopulation
        newlabel<-oldlabel[grepl("clade", oldlabel)]
        for(b in 1:length(oldlabel)){
          index<-which(newer_subpopulationLabels == oldlabel[b])
          newer_subpopulationLabels[index]<-newlabel
        } 
      } else if (sum(grepl("clade", oldlabel))==(length(oldlabel)-1)){
        #if one non-clade subpopulation clusters by expression with multiple clade subpopulations, merge the non-clade subpopulation into closest clade supopulation by expression level
        
        distance_expValues<-as.matrix(dist(expValues))
        oldlabel_b<-oldlabel[1-grepl("clade", oldlabel)]
        oldlabel_aug<-paste0(sampleID, "_",oldlabel[1-grepl("clade", oldlabel)])
        
        
        distance_expValues_subset<-distance_expValues[which(grepl("_clade", rownames(distance_expValues))==1),oldlabel_aug]      
        minExp_newlabel<-min(distance_expValues_subset)
        
        newlabel<-names(which(distance_expValues_subset==minExp_newlabel))
        newlabel<-gsub("^.*\\_","", newlabel)
        
        index<-which(newer_subpopulationLabels == oldlabel_b)
        newer_subpopulationLabels[index]<-newlabel
        
      } else if (sum(grepl("clade", oldlabel))<(length(oldlabel)-1)){
        #if non-clades subpopulation clusters by expression with multiple clade subpopulations, merge the non-clade subpopulations into closest clade supopulation by expression level
        
        distance_expValues<-as.matrix(dist(expValues))
        oldlabel_b<-oldlabel[!grepl("clade", oldlabel)]
        oldlabel_aug<-paste0(sampleID, "_",oldlabel_b)
        
        distance_expValues_subset<-distance_expValues[which(grepl("_clade", rownames(distance_expValues))==1),oldlabel_aug]      
        
        for(b in 1:ncol(distance_expValues_subset)){
          
          minExp_newlabel<-min(distance_expValues_subset[,b])
          
          
          newlabel<-names(which(distance_expValues_subset[,b]==minExp_newlabel))
          newlabel<-gsub("^.*\\_","", newlabel)
          
          index<-which(newer_subpopulationLabels == oldlabel_b[b])
          newer_subpopulationLabels[index]<-newlabel
        } 
      }
    }
    save(newer_subpopulationLabels, file=paste0(sampleID,"_new_subpopulations_Representative_subpops.Rdata"))
  }
}
