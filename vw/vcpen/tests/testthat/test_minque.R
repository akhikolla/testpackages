## Tests for minque


context("Testing the minque output")

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

prefit.save <- readRDS("prefit_minque.rds")
prefit <- minque(response, covmat, Kerns, n.iter=2)

###########################################################################################################
#### Basic functionality
###########################################################################################################

## note, this fails if you add 1e-4 to prefit.save$vc, but not if you add 1e-8
test_that("Basic minque", {
  expect_equal(prefit$vc, expected=prefit.save$vc, tolerance=1e-6)
  expect_equal(prefit$beta, expected=prefit.save$beta, tolerance=1e-6)
  expect_equal(prefit$residuals, expected=prefit.save$residuals, tolerance=1e-6)
  })

