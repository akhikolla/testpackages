#' Calculates the Jaccard similarity matrix.
#' 
#' @param network1 first network matrix input
#' @param network2 second network matrix input
#' @return the alignment/co-occurene score

JaccardSM<-function(network1, network2){
  cooccur<-sum(network1+network2==2)
  activeNodes<-sum((network1+network2)>0)
  alignmentScore<-cooccur/activeNodes
  return(alignmentScore)
}