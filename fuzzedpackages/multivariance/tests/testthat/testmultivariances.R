set.seed(123412)
N = 50
n = 3
x = matrix(rnorm(N*n),nrow = N)

y = matrix(rnorm(N*4),nrow = N)
vec = c(1,1,2,3)

# centered distance matrices ####
context("centered distance matrices")
test_that("cdms",{
  expect_warning(cdm(matrix(1,nrow= N,ncol = n)),"constant")
  expect_equal(cdm(x),cdm(x,psi = function(x,y) sqrt(sum((x-y)^2))))
  expect_equivalent(fastdist(as.matrix(x[,1])),as.matrix(dist(x[,1])))
  expect_equivalent(fastdist(as.matrix(x)),as.matrix(dist(x)))
  x = rnorm(100)

  expect_equal(
    double.center(fastdist(as.matrix(x)),normalize = TRUE),
    cdm(x)
  )

  expect_equal(
    double.center(fastdist(as.matrix(x)),normalize = FALSE),
    cdm(x,normalize = FALSE)
  )

# not public
#  all.equal(cdms(x,normalize = FALSE),
#    cdms.mu.bcd(x,normalize = FALSE)$list.cdm)

# not implemented
#  all.equal(cdms(x,normalize = TRUE),
#    cdms.mu.bcd(x,normalize = TRUE)$list.cdm)
})

# definition of multivariances ####
context("definition of multivariances")
test_that("multivariance, total.multivariance, m.multivariance", {
  expect_warning(multivariance(matrix(1,nrow= N,ncol = n)),"constant")

  expect_equal(multivariance(x), multivariance(x[,c(2,3,1)]))
  expect_equal(multivariance(x), m.multivariance(x,m=3))
  expect_equal(multivariance(x[,c(1,2)]), m.multivariance(x[,c(1,2)],m=2))
  expect_equal((multivariance(x[,c(1,2)]) + multivariance(x[,c(1,3)]) + multivariance(x[,c(3,2)]))/3, m.multivariance(x,m=2))
  expect_equal(multivariance(x[,c(1,2)]), total.multivariance(x[,c(1,2)]))
  expect_equal(total.multivariance(x), (m.multivariance(x,m=3)+m.multivariance(x,m=2)*3)/4)

  expect_equivalent(multivariances.all(x),c(multivariance(x),total.multivariance(x),m.multivariance(x,m=2),m.multivariance(x,m=3)))
  expect_equivalent(multivariances.all(y,vec),c(multivariance(y,vec),total.multivariance(y,vec),m.multivariance(y,vec,m=2),m.multivariance(y,vec,m=3)))
})

# resampling ####
context("resampling")


set.seed(1)
quant = quantile(resample.multivariance(x)$resampled,0.95)
pval = sum(resample.multivariance(x)$resampled>=multivariance(x))/300

set.seed(1)
test_that("resampling p-values and quantiles",{
  expect_equal(resample.rejection.level(0.05,x), quant)
  expect_equal(resample.pvalue(multivariance(x),x), pval)
})

set.seed(1)

mat = matrix( c(1,2,3,1,1,2,1,2,1,2,1,1),nrow = 4,byrow = TRUE)
  for (i in 1:4) {

    for (type in c("multi","total","m.multi.2","m.multi.3")) {
     # print(paste(type,i)) for debugging
    ma = multivariances.all(x,vec = mat[i,])

    expect_equal(
      unname(ma[type]),
      as.numeric(resample.multivariance(x,vec = mat[i,],type = type,times = 2)$original) # here as.numeric is required since in the case of "NA" the return value would be of type logical
    )
    }

    expect_equal(
      ma,
      resample.multivariance(x,vec = mat[i,],type = "all",times = 2)$original
    )

  }


set.seed(1)
x = matrix(rnorm(10*10),10)
vec = c(1:5,1:5)
for (re in c(FALSE))
  for (inc in c(TRUE,FALSE)) {
    #print(paste("replace:",re,", include first:",inc))
    set.seed(1234)
    a = sample.cdms(cdms(x,vec),replace = re,incl.first = inc)
    set.seed(1234)
    b = cdms(sample.cols(x,vec ,replace = re,incl.first = inc),vec)
    #print(all.equal(a,b))
    expect_equal(a,b)
  }


# comparisons in dimensions ####
for (n in c(2,5)) {
# * function arguments ####
context(paste0("function arguments, n = ",n))
set.seed(123412)
N = 5
x = matrix(rnorm(N*n),nrow = N)

for (ty in c("total","m.multi.2","m.multi.3","multi"))
  for (pvt in c("distribution_free","resample","pearson_approx","pearson_unif"))
    multivariance.test(x,type=ty,p.value.type = pvt)

if (n == 5) { # testing argument 'vec'
  for (ty in c("total","m.multi.2","m.multi.3","multi"))
    for (pvt in c("distribution_free","resample","pearson_approx"))
      multivariance.test(x,vec= c(1,2,2,1,3),type=ty,p.value.type = pvt)

  for (ty in c("total","m.multi.2","m.multi.3","multi"))
    for (pvt in c("distribution_free","resample","pearson_approx"))
      multivariance.test(x,vec= c(1,2,2,1,2),type=ty,p.value.type = pvt)

}

for (ty in c("total","m.multi.2","m.multi.3","multi"))
  for (pvt in c("distribution_free","resample","pearson_approx"))
    multivariance.test(x,type=ty,p.value.type = pvt, psi = function(x,y) sum(abs(x-y)))

for (ty in c("total","m.multi.2","m.multi.3","multi"))
  for (pvt in c("distribution_free","resample","pearson_approx"))
    multivariance.test(x,type=ty,p.value.type = pvt, psi = function(x) abs(x), isotropic = TRUE)

for (ty in c("total","m.multi.2","m.multi.3","multi"))
  for (pvt in c("distribution_free","resample","pearson_approx"))
    expect_warning(multivariance.test(x,type=ty,p.value.type = pvt, p = 0.5),"p is not in")
# note: for "pearson_approx" there is more than one warning

for (ty in c("total","m.multi.2","m.multi.3","multi"))
  for (pvt in c("distribution_free","resample","pearson_approx"))
    multivariance.test(x,type=ty,p.value.type = pvt, p = 1.5)

for (ty in c("total"))
  for (pvt in c("distribution_free","resample","pearson_approx"))
    multivariance.test(x,type=ty,p.value.type = pvt, lambda = 2)

for (ty in c("total","m.multi.2","m.multi.3","multi"))
  for (pvt in c("resample"))
    multivariance.test(x,type=ty,p.value.type = pvt, times = 10)

for (ty in c("total","m.multi.2","m.multi.3","multi"))
  for (pvt in c("resample"))
    expect_warning(multivariance.test(x,type=ty,p.value.type = pvt, resample.type = "bootstrap"),"bootstrap")



# * equality of distances ####
context(paste0("equality of distances, n = ",n))

set.seed(123412)
N = 50
x = matrix(rnorm(N*n),nrow = N)

for (ty in c("total","m.multi.2","m.multi.3","multi"))
  for (pvt in c("distribution_free","pearson_approx")) {
    expect_equal( multivariance.test(x,type=ty,p.value.type = pvt),
      multivariance.test(x,type=ty,p.value.type = pvt,psi = function(x,y)sum(abs(x-y)) ))
    expect_equal( multivariance.test(x,type=ty,p.value.type = pvt),
      multivariance.test(x,type=ty,p.value.type = pvt,external.dm.fun = fastdist) )
    expect_equal( multivariance.test(x,type=ty,p.value.type = pvt),
      multivariance.test(x,type=ty,p.value.type = pvt,psi = function(x) abs(x),isotropic = TRUE) )
    expect_equal( multivariance.test(x,type=ty,p.value.type = pvt)[c("statistic","p.value")],
      multivariance.test(apply(x,2,function(y) (abs(rnorm(1))+1)*y+rnorm(1)),type=ty,p.value.type = pvt)[c("statistic","p.value")])

  }

}

# pearson ####
context("pearsons approximation")

set.seed(123)
x = coins(20)
cmb = multivariance:::cdms.mu.bcd(x)

# for cmb list or matrix
expect_equal(
pearson.pvalue(x),
pearson.pvalue(cmb)
)

expect_equal(
pearson.pvalue(x[,1:2]),
pearson.pvalue(cmb,1:2)
)

expect_equal(
  pearson.pvalue(x[,c(1,3)]),
  pearson.pvalue(cmb,c(1,3))
)

expect_equal(
  pearson.pvalue(x[,c(2,3)]),
  pearson.pvalue(cmb,c(2,3))
)

set.seed(123)
x = coins(20,3)
cmb = multivariance:::cdms.mu.bcd(x)

expect_equal(
  pearson.pvalue(cmb,type = "all"),
  pearson.pvalue(x,type = "all")
)

expect_equal(
  pearson.pvalue(x,type = "all"),
c(multi=pearson.pvalue(x,type = "multi"),
  total=pearson.pvalue(x,type = "total"),
  m.multi.2=pearson.pvalue(x,type = "m.multi.2"),
  m.multi.3=pearson.pvalue(x,type = "m.multi.3"))
)

# multicorrelation ####
context("multicorrelation")
set.seed(1213)

y = rnorm(10)
expect_equivalent(multicorrelation(cbind(y,2*y,1-y,y*5-pi,y+1),type="pairwise",multicorrelation.type = "unnormalized",estimator.type = "biased"),1)
expect_equivalent(multicorrelation(cbind(y,2*y,1-y,y*5-pi,y+1),type = "m.multi.3",multicorrelation.type = "unnormalized",estimator.type = "biased"),1)
expect_equivalent(multicorrelation(cbind(y,2*y,1-y,y*5-pi,y+1),type = "multi",multicorrelation.type = "unnormalized",estimator.type = "biased"),1)

x = matrix(rnorm(10*4),10)
expect_equivalent(multicorrelation(x,type = "multi",multicorrelation.type = "unnormalized",estimator.type = "biased"),
  multicorrelation(x,type = "multi",multicorrelation.type = "normalized",estimator.type = "biased"))
expect_equivalent(multicorrelation(x,type = "m.multi.2",multicorrelation.type = "unnormalized",estimator.type = "biased"),
  multicorrelation(x,type = "m.multi.2",multicorrelation.type = "normalized",estimator.type = "biased"))
# expected to be not equal
# expect_equivalent(multicorrelation(x,type = "m.multi.3",multicorrelation.type = "unnormalized"), multicorrelation(x,type = "m.multi.3",multicorrelation.type = "normalized"))

suppressWarnings(
expect_equivalent(
  multicorrelation(matrix(rep(1,10*3),ncol = 3),type = "total.lower",estimator.type = "biased"),
  0)
)

expect_warning(
  multicorrelation(matrix(rep(1,10*3),ncol = 3),type = "total.lower",estimator.type = "biased"),"Constant")


#TODO check that lower bound is a lower bound!

# * unbiased ####

suppressWarnings( expect_equivalent(
  multicorrelation(matrix(rep(1,10*3),ncol = 3)),
  c(0,0)))

# expect_warning( multicorrelation(matrix(rep(1,10*3),ncol = 3)),"zero entries")
# warning removed in 2.3.1

expect_equivalent(multicorrelation(cbind(y,2*y,1-y,y*5-pi,y+1),type="pairwise",multicorrelation.type = "normalized",estimator.type = "bias.corrected"),1)

expect_equivalent(multicorrelation(cbind(y,2*y,1-y,y*5-pi,y+1),type="all",multicorrelation.type = "normalized",estimator.type = "bias.corrected")["unnormalized"],1)


# identity with other measures ####
context("identities with other measures")

set.seed(1234)
# comparison with covariance
n = 2
N = 100
x = matrix(rnorm(N*n),N)

xx = cbind(x[,1],x[,1])

expect_equal(multivariance(xx,psi = function(x) x^2, isotropic = TRUE, Nscale = FALSE,normalize = FALSE),
(2*(N-1)/N*var(x[,1]))^2)

expect_equal(multivariance(x,psi = function(x) x^2, isotropic = TRUE, Nscale = FALSE,normalize = FALSE),
(2*(N-1)/N*cov(x[,1],x[,2]))^2)

expect_equal(
multivariance(x,psi = function(x) x^2, isotropic = TRUE, Nscale = FALSE),
(cor(x[,1],x[,2]))^2
)

expect_equivalent(
  (cor(x[,1],x[,2]))^2,
  multicorrelation(x, type = "multi",estimator.type = "biased",psi = function(x) x^2, isotropic = TRUE,squared = TRUE)
)


expect_equal(
multicorrelation(x, type = "multi",estimator.type = "biased",psi = function(x) x^2, isotropic = TRUE),
multicorrelation(x, type = "multi",estimator.type = "biased", multicorrelation.type = "normalized", psi = function(x) x^2, isotropic = TRUE)
)



## RV
set.seed(1234)
n = 2
N = 100
x = matrix(rnorm(N*n),N)
y = matrix(rnorm(N*n),N)

expect_equal(
sum((2*(N-1)/N*cov(x,y))^2), # as in formula (28)
multivariance(cbind(x,y),vec = c(rep(1,n),rep(2,n)),psi = function(x) x^2, isotropic = TRUE, Nscale = FALSE,normalize = FALSE)
)

expect_equal(
sum((2*(N-1)/N*cov(x,x))^2), # as in formula (28)
multivariance(cbind(x,x),vec = c(rep(1,n),rep(2,n)),psi = function(x) x^2, isotropic = TRUE, Nscale = FALSE,normalize = FALSE)
)

expect_equivalent(
sum((2*(N-1)/N*cov(x,y))^2) /sqrt( sum((2*(N-1)/N*cov(x,x))^2)*sum((2*(N-1)/N*cov(y,y))^2) ), # RV
multicorrelation(cbind(x,y),vec = c(rep(1,n),rep(2,n)), type = "multi",estimator.type = "biased",psi = function(x) x^2, isotropic = TRUE) # thus RV is the squared multicorrelation
)

expect_equal(
  multicorrelation(cbind(x,y),vec = c(rep(1,n),rep(2,n)), type = "multi",estimator.type = "biased",psi = function(x) x^2, isotropic = TRUE),
  multicorrelation(cbind(x,y),vec = c(rep(1,n),rep(2,n)), type = "multi",estimator.type = "biased",psi = function(x) x^2, isotropic = TRUE,multicorrelation.type = "normalized")
)

# dependence structures ####
context("dependence structures")
# just some code, to see if some errors are produced...

set.seed(1023)
x = coins(10,5)
vec = 1:ncol(x)

verbose.output = FALSE

for (sty in c("clustered","full"))
  for (ty in c("conservative","resample","pearson_approx","consistent"))
    dependence.structure(x,vec,type = ty, structure.type = sty,list.cdm = NULL, alpha = 0.05,stop.too.many = 100000,verbose = verbose.output)

for (sty in c("clustered"))
  for (ty in c("conservative","resample","pearson_approx","consistent"))
    dependence.structure(x,vec,type = ty, structure.type = sty,list.cdm = NULL, alpha = 0.05,stop.too.many = 100000,verbose = verbose.output)

vec = c(1:3,1)
for (sty in c("clustered","full"))
  for (ty in c("conservative","resample","pearson_approx","consistent"))
    dependence.structure(x,vec,type = ty, structure.type = sty,list.cdm = NULL, alpha = 0.05,stop.too.many = 100000,verbose = verbose.output)

# input problems ####

if (FALSE) {
  x = 1:10
  y = rnorm(10)

  multicorrelation(cbind(x,y))
  multicorrelation(data.frame(x,y)) # note: a data.frame is a list!
  multivariance(cbind(x,y))
  multivariance(data.frame(x,y))
}

