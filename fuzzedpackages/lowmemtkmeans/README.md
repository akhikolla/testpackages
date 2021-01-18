
The tkmeans package attempts to implement the trimmed k-means algorithm of García-Escudero, et. al.(2008) using as little memory as possible. Data is editted in place, the trimming is implemented using a priority queue structure in C++ trhough Rcpp and low memory use versions of utility functions are provided.

An extremely simple example:
1. Convert the iris dataset to a matrix and rescale matrix columns.

    iris_mat <- as.matrix(iris[,1:4])
    scale_params<-scale_mat_inplace(iris_mat)

1.  Cluster with 2 and 3 clusters, 10% trimming

        iris_cluster_2<- tkmeans(iris_mat, 2 , 0.1, c(1,1,1,1), 1, 10, 0.001)  
        iris_cluster_3<- tkmeans(iris_mat, 2 , 0.1, c(1,1,1,1), 1, 10, 0.001)

2.  Calculate BIC

        BIC_2 <-cluster_BIC(iris_mat, iris_cluster_2)  
        BIC_3 <-cluster_BIC(iris_mat, iris_cluster_3)

3.  Allocate using 3 clustering

        clustering <- nearest_cluster(iris_mat, iris_cluster_3)

4.  Plot results using reconstructed matrix

        library(lattice) 
        orig_matrix <- sweep(sweep(m,2,scale_params[2,],'*'),2,scale_params [1,], '+')  
        xyplot(orig_matrix[,1]~orig_matrix[,2], group=clustering) 

To install the latest version:

    install.packages("devtools")  
    library(devtools)  
    install_github("andrewthomasjones/tkmeans")
