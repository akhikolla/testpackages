## Tests for vcpen


context("Testing the vcpen output")

data(vcexample)
nvc <- 1+length(unique(doseinfo[,2]))
id <- 1:nrow(dose)
## vcs for genetic kernel matrices
Kerns <- vector("list", length=nvc)
for(i in 1:(nvc-1)){
    Kerns[[i]] <- kernel_linear(dose[,grep(i, doseinfo[,2])])
    rownames(Kerns[[i]]) <- id
    colnames(Kerns[[i]]) <- id
}
## vc for residual variance
Kerns[[nvc]] <- diag(nrow(dose))
rownames(Kerns[[nvc]]) <- id
colnames(Kerns[[nvc]]) <- id

vcfit.save <- readRDS("vcfit.rds")
vcfit  <- vcpen(response, covmat, Kerns, frac1 = .6)
#summary(vcfit, digits=1)


###########################################################################################################
#### Basic functionality
###########################################################################################################

test_that("Basic vcpen", {
  expect_equal(vcfit$beta.grid, expected=vcfit.save$beta.grid, tolerance=1e-4)
  expect_equal(vcfit$vc_grid, expected=vcfit.save$vc_grid, tolerance=1e-2)
  expect_equal(vcfit$vc, expected=vcfit.save$vc, tolerance=1e-5)
  expect_equal(vcfit$beta, expected=vcfit.save$beta, tolerance=1e-5)
  })

