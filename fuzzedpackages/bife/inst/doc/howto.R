## ---- eval = FALSE------------------------------------------------------------
#  data(psid, package = "bife")
#  head(psid)

## ---- eval = FALSE------------------------------------------------------------
#  ##    ID LFP KID1 KID2 KID3     INCH AGE TIME
#  ## 1:  1   1    1    1    1 58807.81  26    1
#  ## 2:  1   1    1    0    2 41741.87  27    2
#  ## 3:  1   1    0    1    2 51320.73  28    3
#  ## 4:  1   1    0    1    2 48958.58  29    4
#  ## 5:  1   1    0    1    2 53634.62  30    5
#  ## 6:  1   1    0    0    3 50983.13  31    6

## ---- eval = FALSE------------------------------------------------------------
#  library(bife)
#  stat <- bife(LFP ~ KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID, psid, "probit")
#  summary(stat)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID
#  ##
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.7144667  0.0562414 -12.704   < 2e-16 ***
#  ## KID2      -0.4114554  0.0515524  -7.981  1.45e-15 ***
#  ## KID3      -0.1298776  0.0415477  -3.126   0.00177 **
#  ## log(INCH) -0.2417657  0.0541720  -4.463  8.08e-06 ***
#  ## AGE        0.2319724  0.0375351   6.180  6.40e-10 ***
#  ## I(AGE^2)  -0.0028846  0.0004989  -5.781  7.41e-09 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 6058.88,
#  ## null deviance= 8152.05,
#  ## nT= 5976, N= 664
#  ##
#  ## ( 7173 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6
#  ##
#  ## Average individual fixed effect= -1.121

## ---- eval = FALSE------------------------------------------------------------
#  apes_stat <- get_APEs(stat)
#  summary(apes_stat)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -9.278e-02  8.034e-03 -11.549   < 2e-16 ***
#  ## KID2      -5.343e-02  7.228e-03  -7.393  1.44e-13 ***
#  ## KID3      -1.687e-02  6.009e-03  -2.807     0.005 **
#  ## log(INCH) -3.140e-02  7.515e-03  -4.178  2.95e-05 ***
#  ## AGE        3.012e-02  5.306e-03   5.677  1.37e-08 ***
#  ## I(AGE^2)  -3.746e-04  7.071e-05  -5.298  1.17e-07 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE------------------------------------------------------------
#  stat_bc <- bias_corr(stat)
#  summary(stat_bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID
#  ##
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.6308839  0.0555073 -11.366   < 2e-16 ***
#  ## KID2      -0.3635269  0.0511325  -7.110  1.16e-12 ***
#  ## KID3      -0.1149869  0.0413488  -2.781   0.00542 **
#  ## log(INCH) -0.2139549  0.0536613  -3.987  6.69e-05 ***
#  ## AGE        0.2052708  0.0373054   5.502  3.75e-08 ***
#  ## I(AGE^2)  -0.0025520  0.0004962  -5.143  2.70e-07 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 6062.8,
#  ## null deviance= 8152.05,
#  ## nT= 5976, N= 664
#  ##
#  ## ( 7173 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6

## ---- eval = FALSE------------------------------------------------------------
#  apes_stat_bc <- get_APEs(stat_bc)
#  summary(apes_stat_bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -9.127e-02  7.830e-03 -11.657   < 2e-16 ***
#  ## KID2      -5.259e-02  7.146e-03  -7.359  1.85e-13 ***
#  ## KID3      -1.664e-02  5.962e-03  -2.790   0.00526 **
#  ## log(INCH) -3.095e-02  7.406e-03  -4.180  2.92e-05 ***
#  ## AGE        2.970e-02  5.274e-03   5.632  1.79e-08 ***
#  ## I(AGE^2)  -3.692e-04  7.031e-05  -5.251  1.51e-07 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE------------------------------------------------------------
#  library(data.table)
#  setDT(psid)
#  setkey(psid, ID, TIME)
#  psid[, LLFP := shift(LFP), by = ID]

## ---- eval = FALSE------------------------------------------------------------
#  dyn <- bife(LFP ~ LLFP + KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID, psid, "probit")
#  dyn_bc <- bias_corr(dyn, L = 1L)
#  summary(dyn_bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ LLFP + KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) |
#  ##     ID
#  ##
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## LLFP       1.0025625  0.0473066  21.193   < 2e-16 ***
#  ## KID1      -0.4741275  0.0679073  -6.982  2.91e-12 ***
#  ## KID2      -0.1958365  0.0625921  -3.129  0.001755 **
#  ## KID3      -0.0754042  0.0505110  -1.493  0.135482
#  ## log(INCH) -0.1946970  0.0621143  -3.134  0.001722 **
#  ## AGE        0.2009569  0.0477728   4.207  2.59e-05 ***
#  ## I(AGE^2)  -0.0024142  0.0006293  -3.836  0.000125 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 4822.99,
#  ## null deviance= 6549.14,
#  ## nT= 4792, N= 599
#  ##
#  ## ( 1461 observation(s) deleted due to missingness )
#  ## ( 6896 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6
#  ##
#  ## Average individual fixed effect= -1.939

## ---- eval = FALSE------------------------------------------------------------
#  apes_dyn_bc <- get_APEs(dyn_bc)
#  summary(apes_dyn_bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## LLFP       0.1537523  0.0072038  21.343   < 2e-16 ***
#  ## KID1      -0.0617373  0.0078862  -7.828  4.94e-15 ***
#  ## KID2      -0.0255003  0.0072685  -3.508  0.000451 ***
#  ## KID3      -0.0098186  0.0058864  -1.668  0.095314 .
#  ## log(INCH) -0.0253520  0.0070041  -3.620  0.000295 ***
#  ## AGE        0.0261671  0.0054330   4.816  1.46e-06 ***
#  ## I(AGE^2)  -0.0003144  0.0000714  -4.403  1.07e-05 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

