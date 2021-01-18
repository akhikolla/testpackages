#' @import Rtsne
#' @import ggplot2 
#' @importFrom ggrepel geom_text_repel
NULL

#' Makes constellation plot, in which the centroids are clusters are embedded in the t-SNE 2D plane and the cross-sample relationships are plotted as lines connecting related sample clusters (clades).

#' @param pacman_results PAC-MAN analysis result matrix that contains network annotation, clade IDs and mean (centroid) clade expression levels.
#' @param perplexity perplexity setting for running t-SNE
#' @param max_iter max_iter setting for running t-SNE
#' @param seed set seed to make t-SNE and consetllation plot to be reproducible
#' @param plotTitle max_iter setting for running t-SNE
#' @param nudge_x nudge on x coordinate of centroid labels
#' @param nudge_y nudge on y coordinate of centroid labels
#' 
#' @export

constellationPlot<-function(pacman_results, perplexity,max_iter, seed,plotTitle="Constellations of Clades",nudge_x=0.3, nudge_y=0.3){
  set.seed(seed)

  Rtsne_result<-Rtsne(pacman_results[,-c(1,2,3,ncol(pacman_results))], perplexity=perplexity,max_iter=max_iter)
  Rtsne_pac_result<-cbind(as.data.frame(Rtsne_result$Y), pacman_results[,2])
  colnames(Rtsne_pac_result)<-c("x", "y", "clade")
  
  recurrentClades<-names(which(table(Rtsne_pac_result[,3]) >= 2))
  
  recurrentClades_subset<-Rtsne_pac_result[(Rtsne_pac_result$clade %in% recurrentClades),]
  recurrentClades_subset$clade<-gsub("clade", "", recurrentClades_subset$clade)
  
  non_recurrentClades_subset<-Rtsne_pac_result[(!Rtsne_pac_result$clade %in% recurrentClades),]
  non_recurrentClades_subset$clade<-gsub("clade", "", non_recurrentClades_subset$clade)
  
  p<-ggplot(recurrentClades_subset[,c(1,2)], color=recurrentClades_subset[,3],
            aes(x = recurrentClades_subset$x, y = recurrentClades_subset$y,
                group=recurrentClades_subset$clade)) +
    geom_point(aes(color=recurrentClades_subset[,3]))+
    geom_line(aes(color=recurrentClades_subset[,3])) + labs(x = "bh_tsne1") + labs(y = "bh_tsne2")
  
  p2<-p+geom_point(data=non_recurrentClades_subset[,c(1,2)], 
                   aes(x = non_recurrentClades_subset$x, y = non_recurrentClades_subset$y,group=non_recurrentClades_subset$clade), color= "grey")
  
  unique_clade_indexes<-match(unique(recurrentClades_subset[,3]), recurrentClades_subset[,3])
  p3<-p2+geom_text_repel(data=recurrentClades_subset[unique_clade_indexes,c(1,2)],
                         nudge_x = nudge_x, nudge_y = nudge_y, 
                         aes(x = recurrentClades_subset[unique_clade_indexes,]$x, 
                             y=recurrentClades_subset[unique_clade_indexes,]$y,
                             group=recurrentClades_subset[unique_clade_indexes,]$clade, 
                             label=recurrentClades_subset[unique_clade_indexes,3], 
                             color=recurrentClades_subset[unique_clade_indexes,3]))
  
  p4 <-p3 + geom_text_repel(data=non_recurrentClades_subset[,c(1,2)], 
                            nudge_x = nudge_x, nudge_y = nudge_y, 
                            aes(x = non_recurrentClades_subset$x, 
                                y=non_recurrentClades_subset$y, group=non_recurrentClades_subset$clade, 
                                label=non_recurrentClades_subset[,3]), color="grey")
  
  p4 + theme_bw() +
    theme(axis.line = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) + theme(legend.position="none") +
    ggtitle(plotTitle) + theme(plot.title = element_text(hjust = 0.5))
  
}