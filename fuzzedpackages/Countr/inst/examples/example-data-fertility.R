## make sure the countr pkg is loaded first, then load the fertility data
data(fertility)

## specify the model as in McShane(2008)
form <- children ~ german + years_school + voc_train + university + religion +
    year_birth + rural + age_marriage


## fit the weibull model
wei <- renewalCount(formula = form, data = fertility, dist = "weibull",
                    computeHessian = TRUE, weiMethod = "conv_dePril",
                    control = renewal.control(trace = 0,
                        method = "nlminb"),
                    convPars = list(convMethod = "dePril")
                    )

pois <- glm(formula = form, data = fertility, family = poisson())

## compare residuals: nothing much you can deduce
par(mfrow = c(1, 2))
res_wei <- residuals(wei, type = "pearson")
qqnorm(res_wei, ylim = range(res_wei), main = "Weibull Renewal Model")
qqline(res_wei, ylim = range(res_wei))
grid()
qqnorm(residuals(pois), ylim = range(res_wei), main = "GLM Poisson")
qqline(residuals(pois), ylim = range(res_wei))
grid()

## comparing expected and predicted frequencies
## inspired from Cameron (2013) Chapter 5.3.4
breaks_ <- c(0:5, 7, 9)
pears <- compareToGLM(poisson_model = pois,
                      breaks = breaks_, weibull = wei)


frequency_plot(pears$Counts, pears$Actual,
               dplyr::select(pears, contains("_predicted"))
               )

## run the formal chi-sq test gof
test <- chiSq_gof(wei, breaks = breaks_)

