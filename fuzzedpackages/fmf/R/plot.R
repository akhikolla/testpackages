#'PCA Plot of the Noise Score of Each Individual
#'
#' This function plots the noise score for each observation
#' @param score a vector of values indicating the optential of being a noise.
#' @param data  matrix or data frame with no label.
#' @param cl factor of true classifications of data set.
#' @param geom.ind as geom for observations, which can be set to "text", "point" and "none". The default is "text". 
#' @param labelsize size of geom_text.
#' @param geom_point_size size of geom_point and geom_none.
#' @param ... optional parameters to be passed to other methods.
#' @return an plot of PCA with the noise score of each observation

#' @author Wanwan Zheng
#' @import dplyr
#' @import FactoMineR
#' @import factoextra
#' @import ggplot2
#'
#' @examples
#' 
#' data(iris)
#' out = fmf(Species~.,iris)
#' plot(out$noise_score, iris[,-1], iris[,1])
#' 
#'@name plot
#'@export

plot = function(score,
                data,
                cl,
                geom.ind = "text",
		labelsize = 3, 
 		geom_point_size = 3,
                       ...)
{
  dev.new()
  pca.res = data %>% FactoMineR::PCA(graph = FALSE)
  p = pca.res %>% 
      factoextra::fviz_pca_ind(geom.ind = geom.ind,
                 col.ind = score,
                 gradient.cols = c("#6495ed", "#ff8c00", "#ff0000"),
                 labelsize = labelsize,
                 legend.title = "Noise score") 
  p = p  + geom_point(shape = as.numeric(cl),aes(colour = score),size = geom_point_size)
  p
}
