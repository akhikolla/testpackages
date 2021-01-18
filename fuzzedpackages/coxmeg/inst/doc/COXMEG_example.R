## ---- message=FALSE, warning=FALSE,echo=FALSE---------------------------------
# library(knitcitations)
# cleanbib()
# options("citation_format" = "pandoc")
# r<-citep("10.1534/genetics.119.302940")
# write.bibtex(file="references.bib")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("coxmeg", repos="http://R-Forge.R-project.org")

## ----echo=TRUE----------------------------------------------------------------
library(coxmeg)
library(MASS)
library(Matrix)
n_f <- 200
mat_list <- list()
size <- rep(5,n_f)
offd <- 0.5
for(i in 1:n_f)
{
  mat_list[[i]] <- matrix(offd,size[i],size[i])
  diag(mat_list[[i]]) <- 1
}
sigma <- as.matrix(bdiag(mat_list))
sigma = as(sigma,'dgCMatrix')


## ----echo=TRUE----------------------------------------------------------------
n = nrow(sigma)
tau_var <- 0.2
x <- mvrnorm(1, rep(0,n), tau_var*sigma)
pred = rnorm(n,0,1)
myrates <- exp(x+0.1*pred-1)
y <- rexp(n, rate = myrates)
cen <- rexp(n, rate = 0.02 )
ycen <- pmin(y, cen)
outcome <- cbind(ycen,as.numeric(y <= cen))
head(outcome)
sigma[1:5,1:5]


## ----echo=TRUE----------------------------------------------------------------
re = coxmeg(outcome,sigma,type='bd',X=pred,order=1,detap='diagonal')
re

## ----echo=TRUE----------------------------------------------------------------
library(coxme)
bls <- c(1)
for(i in (size[1]-1):1)
{bls <- c(bls, c(rep(offd,i),1))}
tmat <- bdsmatrix(blocksize=size, blocks=rep(bls,n_f),dimnames=list(as.character(1:n),as.character(1:n)))
re_coxme = coxme(Surv(outcome[,1],outcome[,2])~as.matrix(pred)+(1|as.character(1:n)), varlist=list(tmat),ties='breslow')
re_coxme

## ----echo=TRUE----------------------------------------------------------------
re = coxmeg(outcome,sigma,type='bd',order=1,detap='diagonal')
tau = re$tau
print(tau)
re2 = fit_ppl(pred,outcome,sigma,type='bd',tau=tau,order=1)
re2

## ----echo=TRUE----------------------------------------------------------------
library(coxmeg)
bed = system.file("extdata", "example_null.bed", package = "coxmeg")
bed = substr(bed,1,nchar(bed)-4)
pheno = system.file("extdata", "ex_pheno.txt", package = "coxmeg")
cov = system.file("extdata", "ex_cov.txt", package = "coxmeg")

## building a relatedness matrix
n_f <- 600
mat_list <- list()
size <- rep(5,n_f)
offd <- 0.5
for(i in 1:n_f)
{
  mat_list[[i]] <- matrix(offd,size[i],size[i])
  diag(mat_list[[i]]) <- 1
}
sigma <- as.matrix(bdiag(mat_list))

re = coxmeg_plink(pheno,sigma,type='bd',bed=bed,tmp_dir=tempdir(),cov_file=cov,verbose=FALSE)
re

## ----eval=FALSE---------------------------------------------------------------
#  re = coxmeg_plink(pheno,sigma,type='bd',cov_file=cov,verbose=FALSE)
#  re = coxmeg_plink(pheno,sigma,type='bd',bed=bed,tmp_dir=tempdir(),tau=re$tau,cov_file=cov,verbose=FALSE)

## ----echo=TRUE----------------------------------------------------------------
geno = matrix(rbinom(nrow(sigma)*10,2,runif(nrow(sigma)*10,0.05,0.5)),nrow(sigma),10)
pheno_m = read.table(pheno)
re = coxmeg_m(geno,pheno_m[,3:4],sigma,type='bd',verbose=FALSE)
re

## ----echo=TRUE----------------------------------------------------------------
re = coxmeg_plink(pheno,sigma,type='dense',bed=bed,tmp_dir=tempdir(),cov_file=cov,detap='slq',verbose=FALSE,solver=2)
re

## ----echo=TRUE----------------------------------------------------------------
re = coxmeg_plink(pheno,sigma,type='dense',bed=bed,tmp_dir=tempdir(),tau=re$tau,cov_file=cov,detap='slq',verbose=FALSE,solver=2,score=TRUE)
re

## ----echo=TRUE----------------------------------------------------------------
sigma[2,1] = sigma[1,2] = 1
re = coxmeg_plink(pheno,sigma,type='bd',cov_file=cov,verbose=FALSE,spd=FALSE)
re

