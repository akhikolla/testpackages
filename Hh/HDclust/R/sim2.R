#' Synthetic dataset used in section 5.1.2 of the reference paper.
#' 
#' Dataset used for testing clustering with HMM-VB.
#' The data dimension is 5. Data points were drawn from a 10-component Gaussian Mixture Model.
#' By specific choice of the means, the data contains 10 distinct clusters. For details see the 
#' references. 
#' 
#' @format A data frame with 5000 rows and 5 variables. Last column contains ground truth cluster labels.
#' @name sim2
#' @docType data
#' @references Lin Lin and Jia Li, "Clustering with hidden Markov model on variable blocks," \strong{Journal of Machine Learning Research}, 18(110):1-49, 2017.
"sim2"