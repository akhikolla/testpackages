

# =============================================================================
# Use the same example from <https://cran.r-project.org/web/packages/
#                            SimMultiCorrData/vignettes/workflow.html>.
# =============================================================================
rm(list = ls()); gc()
set.seed(123)
N = 10000L # Sample size.
K = 10L # 10 marginals.
# 3 PDFs, 2 nonparametric PMFs, 5 parametric PMFs:
PMFs = c(
  apply(cbind(rnorm(N), rchisq(N, 4), rbeta(N, 4, 2)), 2, function(x)
    data.frame(val = sort(x), P = 1.0 / N)),
  list(data.frame(val = 1:3 + 0.0, P = c(0.3, 0.45, 0.25))),
  list(data.frame(val = 1:4 + 0.0, P = c(0.2, 0.3, 0.4, 0.1))),
  apply(cbind(rpois(N, 1), rpois(N, 5), rpois(N, 10),
              rnbinom(N, 3, 0.2), rnbinom(N, 6, 0.8)), 2, function(x)
    data.frame(val = as.numeric(sort(x)), P = 1.0 / N))
  )


# Create the target correlation matrix `Rey`:
set.seed(11)
Rey <- diag(1, nrow = 10)
for (i in 1:nrow(Rey)) {
  for (j in 1:ncol(Rey)) {
    if (i > j) Rey[i, j] <- runif(1, 0.2, 0.7)
    Rey[j, i] <- Rey[i, j]
  }
}


system.time({result = SimJoint::SJpearsonPMF(
  PMFs = PMFs, sampleSize = N, cor = Rey, errorType = "meanSquare",
  seed = 456, maxCore = 7, convergenceTail = 8, verbose = T)})


# Check relative errors.
summary(as.numeric(abs(result$cor / Rey - 1)))




# =============================================================================
# Play with random nonparametric PMFs.
# =============================================================================
rm(list = ls()); gc()
set.seed(123)
N = 2000L
K = 20L


# Create totally random nonparametric PMFs:
PMFs = lapply(1L : K, function(x)
{
  p = runif(2, 1, 10)
  result = data.frame(
    val = sort(rnorm(200)), P = runif(200))
  result$P = result$P / sum(result$P)
  result
})


# Create a valid correlation matrix upper-bounded by `frechetUpperCor`.
while(T)
{
  targetCor = matrix(runif(K * K, -0.1, 0.3), ncol = K)
  targetCor[lower.tri(targetCor)] = t(targetCor)[lower.tri(t(targetCor))]
  diag(targetCor) = 1
  if(min(eigen(targetCor)$values) >= 0) break # Break once the correlation
  # matrix is semi-positive definite. This loop could be running for quite
  # a long time if we do not bound `runif()`.
}


result = SimJoint::SJpearsonPMF(
  PMFs = PMFs, sampleSize = N, cor = targetCor, stochasticStepDomain = c(0, 1),
  errorType = "meanSquare", seed = 456, maxCore = 7, convergenceTail = 8)


# Visualize errors and correlation matrices.
par(mfrow = c(2, 2))
hist(result$cor - targetCor, breaks = K * K, main = NULL,
     xlab = "Error", cex.lab = 1.5, cex.axis = 1.25)
hist(result$cor / targetCor - 1, breaks = K * K, main = NULL,
     xlab = "Relative error", ylab = "", cex.lab = 1.5, cex.axis = 1.25)
zlim = range(range(targetCor[targetCor < 1]), range(result$cor[result$cor < 1]))
col = colorRampPalette(c("blue", "red", "yellow"))(K * K)
tmp = targetCor[, K : 1L]
image(tmp, xaxt = "n", yaxt = "n", zlim = zlim, bty = "n",
      main = "Target cor", col = col)
tmp = result$cor[, K : 1L]
image(tmp, xaxt = "n", yaxt = "n", zlim = zlim, bty = "n",
      main = "Cor reached", col = col)
par(mfrow = c(1, 1))














