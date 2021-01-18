# RJcluster

Maintained and written by: Rachael Shudde

Authors of the main manuscript: Shahina Rahman and Valen E. Johnson

Clustering algorithm for big data where the number of observations << the number of covariates. Implementation can be found here: https://arxiv.org/abs/1811.00956

Supports a scalable version of RJ clust.  

In the RJclust function, if the num_cut variable = 1 or is not passed, the non-scaled version of RJclust will be used. if num_cut > 1, the scaled version of RJclust will be used.
