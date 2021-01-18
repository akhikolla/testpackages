

#' Image data of handwritten digits.
#'
#' A dataset containing 100 images of handwritten digits.
#'
#'
#'#'@format A list containing a matrix of image data and a vector of labels:
#' \describe{
#'   \item{Images}{100-by-784 matrix of image data of handwritten digits.}
#'   \item{Labels}{100-by-1 vector of labels of handwritten digits.}
#'   
#' }
#' @source \url{http://yann.lecun.com/exdb/mnist/}
#' 
#' 
#' 
#'@examples
#'data(mnist_data)
#'
#'Img_Mat = mnist_data$Images
#'Img_Label = mnist_data$Labels
#'
#'digit_data = Img_Mat[1, ]      ### image data (784-by-1 vector) of the first handwritten digit (=5) 
#'label = Img_Label[1]           ### label of the first handwritten digit (=5)
#'imgmat = matrix(digit_data, 28, 28)    ### transform the vector of image data to a matrix 
# image(imgmat, axes = FALSE, col = grey(seq(0, 1, length = 256)))   ### convert data to a real image
#'    
#' 
#'@docType data
#'@keywords datasets 
#'@name mnist_data
#'@usage data(mnist_data)
#'@export 
#' 
#' 
NULL
