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
#  library(alpaca)
#  stat <- feglm(LFP ~ KID1 + KID2 + KID3 + log(INCH) | ID + TIME, psid, binomial("probit"))
#  summary(stat)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ KID1 + KID2 + KID3 + log(INCH) | ID + TIME
#  ##
#  ## Estimates:
#  ##            Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.676905   0.056301 -12.023   < 2e-16 ***
#  ## KID2      -0.344383   0.049896  -6.902  5.13e-12 ***
#  ## KID3      -0.007039   0.035344  -0.199     0.842
#  ## log(INCH) -0.234136   0.054403  -4.304  1.68e-05 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 6069.65,
#  ## null deviance= 8152.05,
#  ## n= 5976, l= [664, 9]
#  ##
#  ## ( 7173 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6

## ---- eval = FALSE------------------------------------------------------------
#  apes.stat <- getAPEs(stat)
#  summary(apes.stat)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.0880151  0.0179164  -4.913  8.99e-07 ***
#  ## KID2      -0.0447764  0.0106707  -4.196  2.71e-05 ***
#  ## KID3      -0.0009169  0.0050091  -0.183   0.85476
#  ## log(INCH) -0.0304425  0.0095164  -3.199   0.00138 **
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE------------------------------------------------------------
#  stat.bc <- biasCorr(stat)
#  summary(stat.bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ KID1 + KID2 + KID3 + log(INCH) | ID + TIME
#  ##
#  ## Estimates:
#  ##            Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.596290   0.056301 -10.591   < 2e-16 ***
#  ## KID2      -0.303340   0.049896  -6.079  1.21e-09 ***
#  ## KID3      -0.006124   0.035344  -0.173  0.862431
#  ## log(INCH) -0.207057   0.054403  -3.806  0.000141 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 6069.65,
#  ## null deviance= 8152.05,
#  ## n= 5976, l= [664, 9]
#  ##
#  ## ( 7173 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6

## ---- eval = FALSE------------------------------------------------------------
#  apes.stat.bc <- getAPEs(stat.bc)
#  summary(apes.stat.bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##            Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.086459   0.016344  -5.290  1.22e-07 ***
#  ## KID2      -0.043983   0.010017  -4.391  1.13e-05 ***
#  ## KID3      -0.000888   0.005044  -0.176   0.86026
#  ## log(INCH) -0.030022   0.009222  -3.255   0.00113 **
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE------------------------------------------------------------
#  library(data.table)
#  setDT(psid)
#  psid[, LLFP := shift(LFP), by = ID]

## ---- eval = FALSE------------------------------------------------------------
#  dyn <- feglm(LFP ~ LLFP + KID1 + KID2 + KID3 + log(INCH) | ID + TIME, psid, binomial("probit"))
#  dyn.bc <- biasCorr(dyn, L = 1L)
#  summary(dyn.bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ LLFP + KID1 + KID2 + KID3 + log(INCH) | ID + TIME
#  ##
#  ## Estimates:
#  ##           Estimate Std. error z value Pr(> |z|)
#  ## LLFP       1.01608    0.04695  21.643   < 2e-16 ***
#  ## KID1      -0.45389    0.06787  -6.687  2.27e-11 ***
#  ## KID2      -0.15737    0.06031  -2.610   0.00907 **
#  ## KID3       0.01562    0.04331   0.361   0.71839
#  ## log(INCH) -0.18834    0.06165  -3.055   0.00225 **
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 4777.58,
#  ## null deviance= 6549.14,
#  ## n= 4792, l= [599, 8]
#  ##
#  ## ( 1461 observation(s) deleted due to missingness )
#  ## ( 6896 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6

## ---- eval = FALSE------------------------------------------------------------
#  apes.dyn.bc <- getAPEs(dyn.bc)
#  summary(apes.dyn.bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##            Estimate Std. error z value Pr(> |z|)
#  ## LLFP       0.156044   0.029129   5.357  8.46e-08 ***
#  ## KID1      -0.059108   0.013078  -4.520  6.19e-06 ***
#  ## KID2      -0.020493   0.007684  -2.667   0.00765 **
#  ## KID3       0.002033   0.004875   0.417   0.67663
#  ## log(INCH) -0.024526   0.008165  -3.004   0.00267 **
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

