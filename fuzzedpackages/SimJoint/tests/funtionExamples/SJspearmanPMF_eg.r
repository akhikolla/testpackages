

# =============================================================================
# Play with completely random nonparametric PMFs.
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


result = SimJoint::SJspearmanPMF(
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











