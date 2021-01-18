## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(ClustVarLV)
#library(clv3w)

## ------------------------------------------------------------------------
data(coffee)
# 12 coffee aromas rated by 84 consumers on 15 emotion terms

## ---- results="hide"-----------------------------------------------------
resclv3w_coffee<-CLV3W(coffee,mode.scale=2,NN=TRUE,moddendoinertie=TRUE,graph=TRUE,gmax=11,cp.rand=1)
# option NN=TRUE means that consumers within a group must be positively correlated with its latent component, otherwise its loading is set to 0 

# Print of the 'clv3W' object 
print(resclv3w_coffee)

## ---- fig.width=3, fig.height=3------------------------------------------
# Dendrogram of the CLV3W hierarchical clustering algorithm :
plot(resclv3w_coffee,"dendrogram")
# Graph of the variation of the associated clustering criterion
plot(resclv3w_coffee,"delta")

## ------------------------------------------------------------------------
# Summary the CLV3W results for a partition into 2 groups
summary(resclv3w_coffee,K=2)

## ---- fig.width=4, fig.height=4------------------------------------------
# Representation of the group membership for a partition into 4 groups
plot_var.clv3w(resclv3w_coffee,K=2,labels=TRUE,cex.lab=0.8,beside=TRUE,mode3=FALSE)

## ---- fig.width=6, fig.height=6------------------------------------------
plot_var.clv3w(resclv3w_coffee,K=2,labels=TRUE,cex.lab=0.8,beside=FALSE,mode3=TRUE)

## ---- results="hide"-----------------------------------------------------
# Extract the group membership of each variable
get_partition(resclv3w_coffee,2)

# Extract the group latent variables 
get_comp(resclv3w_coffee,2)

# Extract the vector of loadings of the variables 
get_loading(resclv3w_coffee,2)

# Extract the vector of weights associated with mode3
get_weight(resclv3w_coffee,2)


## ------------------------------------------------------------------------
res.clv3wkm.rd<-CLV3W_kmeans(coffee,2,mode.scale=2,NN=TRUE,init=20,cp.rand=2)

## ------------------------------------------------------------------------
table(get_partition(resclv3w_coffee,K=2),get_partition(res.clv3wkm.rd,K=2)) 

