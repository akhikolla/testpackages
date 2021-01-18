# =============================================================================
# Benchmark against R package `SimMultiCorrData`. Use the same example
# from <https://cran.r-project.org/web/packages/SimMultiCorrData/
#       vignettes/workflow.html>.
# =============================================================================

set.seed(123)
N = 10000L # Sample size.
K = 10L    # 10 marginals.
# Sample from 3 PDFs, 2 nonparametric PMFs, 5 parametric PMFs:
marginals = cbind(
  rnorm(N), rchisq(N, 4), rbeta(N, 4, 2),
  SimJoint::LHSpmf(data.frame(val = 1:3, P = c(0.3, 0.45, 0.25)), N,
               seed = sample(1e6L, 1)),
  SimJoint::LHSpmf(data.frame(val = 1:4, P = c(0.2, 0.3, 0.4, 0.1)), N,
               seed = sample(1e6L, 1)),
  rpois(N, 1), rpois(N, 5), rpois(N, 10),
  rnbinom(N, 3, 0.2), rnbinom(N, 6, 0.8))
# The seeding for `LHSpmf()` is unhealthy, but OK for small examples.


marginals = apply(marginals, 2, function(x) sort(x))


# Create the example target correlation matrix `Rey`:
set.seed(11)
Rey <- diag(1, nrow = K)
for (i in 1:nrow(Rey)) {
  for (j in 1:ncol(Rey)) {
    if (i > j) Rey[i, j] <- runif(1, 0.2, 0.7)
    Rey[j, i] <- Rey[i, j]
  }
}


system.time({result = SimJoint::SJpearson(
  X = marginals, cor = Rey, errorType = "meanSquare", seed = 456,
  maxCore = 1, convergenceTail = 8, verbose = FALSE)})
# user  system elapsed
# 0.30    0.00    0.29
# One the same platform, single-threaded speed (Intel i7-4770 CPU
# @ 3.40GHz, 32GB RAM, Windows 10, g++ 4.9.3 -Ofast, R 3.5.2) is more
# than 50 times faster than `SimMultiCorrData::rcorrvar()`:
# user  system elapsed
# 16.05   0.34   16.42


# Check error statistics.
summary(as.numeric(round(cor(result$X) - Rey, 6)))
# Min.          1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.000365   -0.000133  -0.000028  -0.000047  0.000067    0.000301


# Post simulation optimization further reduce the errors:
resultOpt = SimJoint::postSimOpt(
  X = result$X, cor = Rey, convergenceTail = 10000)
summary(as.numeric(round(cor(resultOpt$X) - Rey, 6)))
# Min.        1st Qu.    Median      Mean   3rd Qu.      Max.
# -7.10e-05 -3.10e-05 -1.15e-05 -6.48e-06  9.00e-06  7.10e-05


# Max error magnitude is less than 1% of that from
# `SimMultiCorrData::rcorrvar()`:
# Min.          1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.008336   -0.001321          0  -0.000329   0.001212    0.00339
# This table is reported in Step 4, correlation methods 1 or 2.




# =============================================================================
# Use the above example and benchmark against John Ruscio & Walter
# Kaczetow (2008) iteration. The R code released with their paper was
# erroneous. A corrected version is given by Github user "nicebread":
# <https://gist.github.com/nicebread/4045717>, but his correction was
# incomprehensive and can only handle 2-dimensional instances. Please change
# Line 32 to `Target.Corr <- rho` and source the file.
# =============================================================================
# \donttest{# Test Ruscio-Kaczetow's code.
#   set.seed(123)
#   RuscioKaczetow = GenData(Pop = as.data.frame(marginals),
#                            Rey, N = 1000) # By default, the function takes 1000
#   # samples from each marginal population of size 10000.
#   summary(round(as.numeric(cor(RuscioKaczetow) - Rey), 6))
#   # Min.         1st Qu.    Median      Mean   3rd Qu.      Max.
#   # -0.183274 -0.047461  -0.015737 -0.008008  0.027475  0.236662
# }


result = SimJoint::SJpearson(
  X = apply(marginals, 2, function(x) sort(sample(x, 1000, replace = TRUE))),
  cor = Rey, errorType = "maxRela", maxCore = 2) # CRAN does not allow more
# than 2 threads for running examples.
summary(round(as.numeric(cor(result$X) - Rey), 6))
# Min.           1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.0055640  -0.0014850 -0.0004810 -0.0007872  0.0000000  0.0025920
resultOpt = SimJoint::postSimOpt(
  X = result$X, cor = Rey, convergenceTail = 10000)
summary(as.numeric(round(cor(resultOpt$X) - Rey, 6)))
# Min.          1st Qu.     Median       Mean    3rd Qu.       Max.
# -6.240e-04 -2.930e-04 -2.550e-05 -6.532e-05  1.300e-04  5.490e-04




# =============================================================================
# Benchmark against R package `GenOrd`
# <https://cran.r-project.org/web/packages/GenOrd/index.html> using the
# example above Statistics cannot be collected because it has been running
# for more than 10 hours.
# =============================================================================
# \donttest{# Library `GenOrd` should have been installed and attached.
  # system.time({resultGenOrd = ordsample(
  #   N, marginal = lapply(1L : K, function(x) (1 : (N - 1)) / N), Rey,
  #   support = as.data.frame(marginals))})
# }




# =============================================================================
# Benchmark against R package `EnvStats` using its manual example on Page 1156
# of <https://cran.r-project.org/web/packages/EnvStats/EnvStats.pdf>. The
# function `simulateVector()` imposes rank correlations.
# =============================================================================
# \donttest{# Library `EnvStats` should have been installed and attached.
  cor.mat = matrix(c(1, 0.8, 0, 0.5, 0.8, 1, 0, 0.7, 0, 0, 1, 0.2, 0.5,
                     0.7, 0.2, 1), 4, 4)
  pareto.rns <- simulateVector(100, "pareto", list(location = 10, shape = 2),
                               sample.method = "LHS", seed = 56)
  mat <- simulateMvMatrix(
    1000, distributions = c(Normal = "norm", Lognormal = "lnormAlt",
                            Beta = "beta", Empirical = "emp"),
    param.list = list(Normal = list(mean=10, sd=2),
                      Lognormal = list(mean=10, cv=1),
                      Beta = list(shape1 = 2, shape2 = 3),
                      Empirical = list(obs = pareto.rns)),
    cor.mat = cor.mat, seed = 47, sample.method = "LHS")


  round(cor(mat, method = "spearman"), 2)
  #          Normal Lognormal    Beta Empirical
  #Normal      1.00    0.78     -0.01      0.47
  #Lognormal   0.78    1.00     -0.01      0.67
  #Beta       -0.01   -0.01      1.00      0.19
  #Empirical   0.47    0.67      0.19      1.00


  # Imposing rank correlations is equivalent to imposing Pearson correlations
  # on ranks.
  set.seed(123)
  marginals = cbind(sort(rnorm(1000, 10, 2)),
                    sort(rlnormAlt(1000, 10, 1)),
                    sort(rbeta(1000, 2, 3)),
                    sort(sample(pareto.rns, 1000, replace = TRUE)))
  marginalsRanks = cbind(1:1000, 1:1000, 1:1000, 1:1000)
  # Simulate the joint for ranks:
  tmpResult = SimJoint::SJpearson(
    X = marginalsRanks, cor = cor.mat, errorType = "meanSquare", seed = 456,
    maxCore = 2, convergenceTail = 8, verbose = TRUE)$X
  # Reorder `marginals` by ranks.
  result = matrix(mapply(function(x, y) y[as.integer(x)],
                         as.data.frame(tmpResult),
                         as.data.frame(marginals), SIMPLIFY = TRUE), ncol = 4)
  round(cor(result, method = "spearman"), 2)
  # 1.0  0.8  0.0  0.5
  # 0.8  1.0  0.0  0.7
  # 0.0  0.0  1.0  0.2
  # 0.5  0.7  0.2  1.0
# }




# ============================================================================
# Play random numbers.
# ============================================================================
set.seed(123)
N = 2000L
K = 20L
# The following essentially creates a mixture distribution.
marginals = c(runif(10000L, -2, 2), rgamma(10000L, 2, 2), rnorm(20000L))
marginals = matrix(sample(marginals, length(marginals)), ncol = K)
# This operation made the columns comprise samples from the same
# mixture distribution.
marginals = apply(marginals, 2, function(x) sort(x))


# May take a while to generate valid correlation matrix.
while(TRUE)
{
  targetCor = matrix(runif(K * K, -0.1, 0.4), ncol = K)
  targetCor[lower.tri(targetCor)] = t(targetCor)[lower.tri(t(targetCor))]
  diag(targetCor) = 1
  if(all(eigen(targetCor)$values >= 0)) break
}


result = SimJoint::SJpearson(
  X = marginals, cor = targetCor, errorType = "meanSquare", seed = 456,
  maxCore = 2, convergenceTail = 8, verbose = TRUE)
resultOpt = SimJoint::postSimOpt(
  X = result$X, cor = targetCor, convergenceTail = 10000)


# Visualize errors and correlation matrices.
par(mfrow = c(2, 2))
hist(resultOpt$cor - targetCor, breaks = K * K, main = NULL,
     xlab = "Error")
hist(resultOpt$cor / targetCor - 1, breaks = K * K, main = NULL,
     xlab = "Relative error")
zlim = range(range(targetCor[targetCor < 1]),
             range(resultOpt$cor[resultOpt$cor < 1]))
col = colorRampPalette(c("blue", "red", "yellow"))(K * K)
tmp = targetCor[, K : 1L]
image(tmp, xaxt = "n", yaxt = "n", zlim = zlim, bty = "n",
      main = "Target cor", col = col)
tmp = resultOpt$cor[, K : 1L]
image(tmp, xaxt = "n", yaxt = "n", zlim = zlim, bty = "n",
      main = "Cor reached", col = col)
par(mfrow = c(1, 1))




# =============================================================================
# An example where the functional relationships between marginals are highly
# nonlinear and the target correlations are hard to impose. Other packages
# would fail or report theoretical infeasibility.
# =============================================================================
set.seed(123)
N = 10000L
K = 10L


# Several 2-parameter PDFs in R:
marginals = list(rbeta, rcauchy, rf, rgamma, rnorm, runif, rweibull)
Npdf = length(marginals)


if(Npdf >= K) chosenMarginals =
  marginals[sample(Npdf, K, replace = TRUE)] else chosenMarginals =
  marginals[c(1L : Npdf, sample(Npdf, K - Npdf, replace = TRUE))]


# Sample from the marginal PDFs.
marginals = as.matrix(as.data.frame(lapply(chosenMarginals, function(f)
{
  para = sort(runif(2, 0.1, 10))
  rst = f(N, para[1], para[2])
  sort(rst)
})))
dimnames(marginals) = NULL


frechetUpperCor = cor(marginals) # The correlation matrix should be
# upper-bounded by that of the perfectly rank-correlated
# joint (Frechet upper bound). The lower bound is characterized by
# d-countercomonotonicity and depends not only on marginals.
cat("Range of maximal correlations between marginals:",
    range(frechetUpperCor[frechetUpperCor < 1]))
# Two perfectly rank-correlated marginals can have a Pearson
# correlation below 0.07. This is due to highly nonlinear functional
# relationships between marginal PDFs.


# Create a valid correlation matrix upper-bounded by `frechetUpperCor`.
while(TRUE)
{
  targetCor = sapply(frechetUpperCor, function(x)
    runif(1, -0.1, min(0.5, x * 0.8)))
  targetCor = matrix(targetCor, ncol = K)
  targetCor[lower.tri(targetCor)] = t(targetCor)[lower.tri(t(targetCor))]
  diag(targetCor) = 1
  if(min(eigen(targetCor)$values) >= 0) break # Stop once the correlation
  # matrix is semi-positive definite. This loop could run for
  # a long time if we do not bound the uniform by 0.3.
}


system.time({result = SimJoint::SJpearson(
  X = marginals, cor = targetCor, stochasticStepDomain = c(0, 1),
  errorType = "meanSquare", seed = 456, maxCore = 7, convergenceTail = 8, verbose = F)})
# \donttest{
resultOpt = SimJoint::postSimOpt( # Could take many seconds.
  X = result$X, cor = targetCor, convergenceTail = 10000)


# Visualize errors and correlation matrices.
par(mfrow = c(2, 2))
hist(resultOpt$cor - targetCor, breaks = K * K, main = NULL,
     xlab = "Error")
hist(resultOpt$cor / targetCor - 1, breaks = K * K, main = NULL,
     xlab = "Relative error")
zlim = range(range(targetCor[targetCor < 1]),
             range(resultOpt$cor[resultOpt$cor < 1]))
col = colorRampPalette(c("blue", "red", "yellow"))(K * K)
tmp = targetCor[, K : 1L]
image(tmp, xaxt = "n", yaxt = "n", zlim = zlim, bty = "n",
      main = "Target cor", col = col)
tmp = resultOpt$cor[, K : 1L]
image(tmp, xaxt = "n", yaxt = "n", zlim = zlim, bty = "n",
      main = "Cor reached", col = col)
par(mfrow = c(1, 1))
# }


# Different `errorType` could make a difference.
result = SimJoint::SJpearson(
  X = marginals, cor = targetCor, stochasticStepDomain = c(0, 1),
  errorType = "maxRela", seed = 456, maxCore = 2, convergenceTail = 8)
# \donttest{
resultOpt = SimJoint::postSimOpt(
X = result$X, cor = targetCor, convergenceTail = 10000)


# Visualize errors and correlation matrices.
par(mfrow = c(2, 2))
hist(resultOpt$cor - targetCor, breaks = K * K, main = NULL,
     xlab = "Error")
hist(resultOpt$cor / targetCor - 1, breaks = K * K, main = NULL,
     xlab = "Relative error")
zlim = range(range(targetCor[targetCor < 1]),
             range(resultOpt$cor[resultOpt$cor < 1]))
col = colorRampPalette(c("blue", "red", "yellow"))(K * K)
tmp = targetCor[, K : 1L]
image(tmp, xaxt = "n", yaxt = "n", zlim = zlim, bty = "n",
      main = "Target cor", col = col)
tmp = resultOpt$cor[, K : 1L]
image(tmp, xaxt = "n", yaxt = "n", zlim = zlim, bty = "n",
      main = "Cor reached", col = col)
par(mfrow = c(1, 1)) # }






















