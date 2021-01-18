#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
error <-function(x,y){
  rf <- randomForest(x, y, prox=TRUE)
  return((rf$confusion[1,2]+rf$confusion[2,1])/sum(rf$confusion[,-3]))
}
