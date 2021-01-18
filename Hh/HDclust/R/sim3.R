#' Synthetic dataset used in section 5.1.3 of the reference paper
#' 
#' Dataset used for testing clustering with HMM-VB.
#' The data dimension is 40. The first 10 dimensions were generated from a 3-component Gaussian
#' Mixture Model (GMM). The remaining 30 dimensions were generated from a 5-component GMM. By specific design
#' of the means, covariance matrices and transition probabilities, the data contain 5 distinct clusters. 
#' For details see the references.
#' 
#' @format A data frame with 1000 rows and 40 variables. Last column contains ground truth cluster labels.
#' @name sim3
#' @docType data
#' @references Lin Lin and Jia Li, "Clustering with hidden Markov model on variable blocks," \strong{Journal of Machine Learning Research}, 18(110):1-49, 2017.
"sim3"