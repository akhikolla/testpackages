## ---- eval = FALSE------------------------------------------------------------
#  # Import the data set
#  library(haven)
#  library(data.table)
#  cudata <- read_dta("dataaxj1.dta")
#  setDT(cudata)
#  
#  # Subsetting relevant variables
#  var.nms <- c("exp1to2", "custrict11", "ldist", "comlang", "border", "regional",
#               "comcol", "curcol", "colony", "comctry", "cuwoemu", "emu", "cuc",
#               "cty1", "cty2", "year", "pairid")
#  cudata <- cudata[, var.nms, with = FALSE]
#  
#  # Generate identifiers required for structural gravity
#  cudata[, pairid := factor(pairid)]
#  cudata[, exp.time := interaction(cty1, year)]
#  cudata[, imp.time := interaction(cty2, year)]
#  
#  # Generate dummies for disaggregated currency unions
#  cudata[, cuau := as.numeric(cuc == "au")]
#  cudata[, cube := as.numeric(cuc == "be")]
#  cudata[, cuca := as.numeric(cuc == "ca")]
#  cudata[, cucf := as.numeric(cuc == "cf")]
#  cudata[, cucp := as.numeric(cuc == "cp")]
#  cudata[, cudk := as.numeric(cuc == "dk")]
#  cudata[, cuea := as.numeric(cuc == "ea")]
#  cudata[, cuec := as.numeric(cuc == "ec")]
#  cudata[, cuem := as.numeric(cuc == "em")]
#  cudata[, cufr := as.numeric(cuc == "fr")]
#  cudata[, cugb := as.numeric(cuc == "gb")]
#  cudata[, cuin := as.numeric(cuc == "in")]
#  cudata[, cuma := as.numeric(cuc == "ma")]
#  cudata[, cuml := as.numeric(cuc == "ml")]
#  cudata[, cunc := as.numeric(cuc == "nc")]
#  cudata[, cunz := as.numeric(cuc == "nz")]
#  cudata[, cupk := as.numeric(cuc == "pk")]
#  cudata[, cupt := as.numeric(cuc == "pt")]
#  cudata[, cusa := as.numeric(cuc == "sa")]
#  cudata[, cusp := as.numeric(cuc == "sp")]
#  cudata[, cuua := as.numeric(cuc == "ua")]
#  cudata[, cuus := as.numeric(cuc == "us")]
#  cudata[, cuwa := as.numeric(cuc == "wa")]
#  cudata[, cuwoo := custrict11]
#  cudata[cuc %in% c("em", "au", "cf", "ec", "fr", "gb", "in", "us"), cuwoo := 0L]
#  
#  # Set missing trade flows to zero
#  cudata[is.na(exp1to2), exp1to2 := 0.0]
#  
#  # Construct binary and lagged dependent variable for the extensive margin
#  cudata[, trade := as.numeric(exp1to2 > 0.0)]
#  cudata[, ltrade := shift(trade), by = pairid]

## ---- eval = FALSE------------------------------------------------------------
#  mod <- feglm(exp1to2 ~ emu + cuwoo + cuau + cucf + cuec + cufr + cugb + cuin + cuus +
#                 regional + curcol | exp.time + imp.time + pairid | cty1 + cty2 + year, cudata,
#               family = poisson())
#  summary(mod, "sandwich")

## ---- eval = FALSE------------------------------------------------------------
#  ## poisson - log link
#  ##
#  ## exp1to2 ~ emu + cuwoo + cuau + cucf + cuec + cufr + cugb + cuin +
#  ##     cuus + regional + curcol | exp.time + imp.time + pairid |
#  ##     cty1 + cty2 + year
#  ##
#  ## Estimates:
#  ##           Estimate Std. error z value Pr(> |z|)
#  ## emu       0.048895   0.010277   4.758  1.96e-06 ***
#  ## cuwoo     0.765988   0.053272  14.379   < 2e-16 ***
#  ## cuau      0.384469   0.118832   3.235   0.00121 **
#  ## cucf     -0.125608   0.099674  -1.260   0.20760
#  ## cuec     -0.877318   0.083451 -10.513   < 2e-16 ***
#  ## cufr      2.095726   0.062952  33.291   < 2e-16 ***
#  ## cugb      1.059957   0.034680  30.564   < 2e-16 ***
#  ## cuin      0.169745   0.147029   1.154   0.24830
#  ## cuus      0.018323   0.021530   0.851   0.39473
#  ## regional  0.159181   0.008714  18.267   < 2e-16 ***
#  ## curcol    0.386882   0.046827   8.262   < 2e-16 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 35830778.98,
#  ## null deviance= 2245707302.17,
#  ## n= 1610165, l= [11227, 11277, 34104]
#  ##
#  ## ( 1363003 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 13

## ---- eval = FALSE------------------------------------------------------------
#  summary(mod, "clustered", cluster = ~ cty1 + cty2 + year)

## ---- eval = FALSE------------------------------------------------------------
#  ## poisson - log link
#  ##
#  ## exp1to2 ~ emu + cuwoo + cuau + cucf + cuec + cufr + cugb + cuin +
#  ##     cuus + regional + curcol | exp.time + imp.time + pairid |
#  ##     cty1 + cty2 + year
#  ##
#  ## Estimates:
#  ##          Estimate Std. error z value Pr(> |z|)
#  ## emu       0.04890    0.09455   0.517   0.60507
#  ## cuwoo     0.76599    0.24933   3.072   0.00213 **
#  ## cuau      0.38447    0.22356   1.720   0.08547 .
#  ## cucf     -0.12561    0.35221  -0.357   0.72137
#  ## cuec     -0.87732    0.29493  -2.975   0.00293 **
#  ## cufr      2.09573    0.30625   6.843  7.75e-12 ***
#  ## cugb      1.05996    0.23766   4.460  8.19e-06 ***
#  ## cuin      0.16974    0.30090   0.564   0.57267
#  ## cuus      0.01832    0.05092   0.360   0.71898
#  ## regional  0.15918    0.07588   2.098   0.03593 *
#  ## curcol    0.38688    0.15509   2.495   0.01261 *
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 35830778.98,
#  ## null deviance= 2245707302.17,
#  ## n= 1610165, l= [11227, 11277, 34104]
#  ##
#  ## ( 1363003 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 13

## ---- eval = FALSE------------------------------------------------------------
#  library(car)
#  cus <- c("cuwoo", "cuau", "cucf", "cuec", "cufr", "cugb", "cuin", "cuus")
#  linearHypothesis(mod, cus, vcov. = vcov(mod, "clustered", cluster = ~ cty1 + cty2 + year))

## ---- eval = FALSE------------------------------------------------------------
#  ## Linear hypothesis test
#  ##
#  ## Hypothesis:
#  ## cuwoo = 0
#  ## cuau = 0
#  ## cucf = 0
#  ## cuec = 0
#  ## cufr = 0
#  ## cugb = 0
#  ## cuin = 0
#  ## cuus = 0
#  ##
#  ## Model 1: restricted model
#  ## Model 2: exp1to2 ~ emu + cuwoo + cuau + cucf + cuec + cufr + cugb + cuin +
#  ##     cuus + regional + curcol | exp.time + imp.time + pairid |
#  ##     cty1 + cty2 + year
#  ##
#  ## Note: Coefficient covariance matrix supplied.
#  ##
#  ##   Df  Chisq Pr(>Chisq)
#  ## 1
#  ## 2  8 96.771  < 2.2e-16 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE------------------------------------------------------------
#  mods <- feglm(trade ~ emu + cuwoo + cuau + cucf + cuec + cufr + cugb + cuin + cuus +
#                  regional + curcol | exp.time + imp.time + pairid, cudata,
#                family = binomial("probit"))
#  summary(mods)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## trade ~ emu + cuwoo + cuau + cucf + cuec + cufr + cugb + cuin +
#  ##     cuus + regional + curcol | exp.time + imp.time + pairid
#  ##
#  ## Estimates:
#  ##          Estimate Std. error z value Pr(> |z|)
#  ## emu       0.28677    0.12842   2.233  0.025546 *
#  ## cuwoo     0.57036    0.07675   7.432  1.07e-13 ***
#  ## cuau      1.09776    0.28863   3.803  0.000143 ***
#  ## cucf      0.14035    0.05250   2.673  0.007507 **
#  ## cuec      0.02973    0.63120   0.047  0.962432
#  ## cufr      2.31013    0.24334   9.493   < 2e-16 ***
#  ## cugb      0.18260    0.04252   4.294  1.75e-05 ***
#  ## cuin      0.25333    0.15696   1.614  0.106519
#  ## cuus     -0.36223    0.07908  -4.580  4.64e-06 ***
#  ## regional  0.06259    0.01782   3.513  0.000444 ***
#  ## curcol   -0.93748    0.23913  -3.920  8.84e-05 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 681026.48,
#  ## null deviance= 1947737.97,
#  ## n= 1406505, l= [11227, 11273, 29630]
#  ##
#  ## ( 1566663 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 8

## ---- eval = FALSE------------------------------------------------------------
#  modsbc <- biasCorr(mods, panel.structure = "network")
#  summary(modsbc)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## trade ~ emu + cuwoo + cuau + cucf + cuec + cufr + cugb + cuin +
#  ##     cuus + regional + curcol | exp.time + imp.time + pairid
#  ##
#  ## Estimates:
#  ##          Estimate Std. error z value Pr(> |z|)
#  ## emu       0.26936    0.12842   2.097  0.035952 *
#  ## cuwoo     0.50427    0.07675   6.571  5.01e-11 ***
#  ## cuau      0.99328    0.28863   3.441  0.000579 ***
#  ## cucf      0.12976    0.05250   2.472  0.013448 *
#  ## cuec     -0.02330    0.63120  -0.037  0.970558
#  ## cufr      2.07974    0.24334   8.547   < 2e-16 ***
#  ## cugb      0.16612    0.04252   3.907  9.35e-05 ***
#  ## cuin      0.20516    0.15696   1.307  0.191182
#  ## cuus     -0.33674    0.07908  -4.258  2.06e-05 ***
#  ## regional  0.06273    0.01782   3.520  0.000431 ***
#  ## curcol   -0.72433    0.23913  -3.029  0.002453 **
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 681026.48,
#  ## null deviance= 1947737.97,
#  ## n= 1406505, l= [11227, 11273, 29630]
#  ##
#  ## ( 1566663 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 8

## ---- eval = FALSE------------------------------------------------------------
#  apesbc <- getAPEs(modsbc, n.pop = nrow(cudata))
#  summary(apesbc)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##           Estimate Std. error z value Pr(> |z|)
#  ## emu       0.018890   0.008526   2.216  0.026713 *
#  ## cuwoo     0.035655   0.005908   6.035  1.59e-09 ***
#  ## cuau      0.070750   0.014709   4.810  1.51e-06 ***
#  ## cucf      0.009045   0.003783   2.391  0.016796 *
#  ## cuec     -0.001612   0.038730  -0.042  0.966810
#  ## cufr      0.154181   0.023722   6.500  8.05e-11 ***
#  ## cugb      0.011596   0.002997   3.869  0.000109 ***
#  ## cuin      0.014349   0.011038   1.300  0.193634
#  ## cuus     -0.022868   0.005267  -4.342  1.41e-05 ***
#  ## regional  0.004360   0.001238   3.521  0.000430 ***
#  ## curcol   -0.047850   0.018521  -2.584  0.009778 **
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE------------------------------------------------------------
#  modd <- feglm(trade ~ ltrade + emu + cuwoo + cuau + cucf + cuec + cufr + cugb + cuin + cuus +
#                  regional + curcol | exp.time + imp.time + pairid, cudata,
#                family = binomial("probit"))
#  moddbc <- biasCorr(modd, L = 2L, panel.structure = "network")
#  summary(moddbc)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## trade ~ ltrade + emu + cuwoo + cuau + cucf + cuec + cufr + cugb +
#  ##     cuin + cuus + regional + curcol | exp.time + imp.time + pairid
#  ##
#  ## Estimates:
#  ##          Estimate Std. error z value Pr(> |z|)
#  ## ltrade    1.24347    0.00443 280.716   < 2e-16 ***
#  ## emu       0.18254    0.13543   1.348   0.17772
#  ## cuwoo     0.41716    0.08064   5.173  2.30e-07 ***
#  ## cuau      0.79027    0.30405   2.599   0.00935 **
#  ## cucf      0.10150    0.05446   1.864   0.06238 .
#  ## cuec      0.17233    0.68493   0.252   0.80135
#  ## cufr      1.46437    0.24241   6.041  1.53e-09 ***
#  ## cugb      0.10085    0.04452   2.265   0.02348 *
#  ## cuin      0.17074    0.16391   1.042   0.29757
#  ## cuus     -0.25676    0.08216  -3.125   0.00178 **
#  ## regional  0.04644    0.01867   2.488   0.01284 *
#  ## curcol   -0.43458    0.25141  -1.729   0.08388 .
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 603824.81,
#  ## null deviance= 1930471.73,
#  ## n= 1394201, l= [11133, 11178, 29538]
#  ##
#  ## ( 45048 observation(s) deleted due to missingness )
#  ## ( 1533919 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 22

## ---- eval = FALSE------------------------------------------------------------
#  apedbc <- getAPEs(moddbc, n.pop = nrow(cudata[!is.na(ltrade)]))
#  summary(apedbc)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##            Estimate Std. error z value Pr(> |z|)
#  ## ltrade    0.1148321  0.0005004 229.459   < 2e-16 ***
#  ## emu       0.0114102  0.0077264   1.477  0.139737
#  ## cuwoo     0.0265827  0.0047666   5.577  2.45e-08 ***
#  ## cuau      0.0520368  0.0151715   3.430  0.000604 ***
#  ## cucf      0.0063053  0.0031961   1.973  0.048516 *
#  ## cuec      0.0107638  0.0369895   0.291  0.771055
#  ## cufr      0.1019601  0.0246482   4.137  3.52e-05 ***
#  ## cugb      0.0062640  0.0024262   2.582  0.009827 **
#  ## cuin      0.0106628  0.0094163   1.132  0.257478
#  ## cuus     -0.0155832  0.0043727  -3.564  0.000366 ***
#  ## regional  0.0028744  0.0010406   2.762  0.005739 **
#  ## curcol   -0.0261449  0.0154963  -1.687  0.091571 .
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

