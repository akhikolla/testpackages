context("gRain misc")

##############################

data(chest_cpt); bn <- grain(compileCPT(chest_cpt))

test_that("nodeNames()", {
    vn <- c("asia", "tub", "smoke", "lung", "bronc", "either", "xray", "dysp")
    expect_setequal(nodeNames(bn), vn)
})

## A <- matrix(1:16, nrow = 4, ncol = 4)
## a <- 1:4

## test_that("ysym()", {
##   B <- ysym(A)
##   b <- ysym(a)
  
##   expect_s3_class(B, "yac_symbol")
##   expect_s3_class(b, "yac_symbol")
  
##   expect_equal(A, as_r(B))
  
##   x <- B %*% b
##   expect_s3_class(x, "yac_symbol")
##   expect_equal(x$yacas_cmd, "{90,100,110,120}")
  
##   expect_equal(eval(yac_expr(x)), c(A %*% a))
## })
