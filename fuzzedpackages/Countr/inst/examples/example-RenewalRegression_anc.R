library(testthat)
form_main <- children ~ german + years_school + voc_train + university + religion +
    year_birth + rural + age_marriage

form_anc <- children ~ german + years_school + voc_train

IV_scale <- glm(form_main, family = poisson(), data = fertility)
IV_shape <- glm(form_anc, family = poisson(), data = fertility)
startR <- renewalCoef(IV_scale, target = "gamma")
startS <- renewalCoef(IV_shape, target = "gamma")
start <- c(startR[!grepl("shape",  names(startR))], 
           startS[grepl("shape",  names(startS))])

## add regresiion on the shape parameter
anc <- list(shape = form_anc)
## =========================== flexsurv API =====================================
print("............ flexsurv API ............")
res_old <- renewalCount(formula = form_main, data = fertility, dist = "gamma",
                        computeHessian = FALSE, anc = anc,
                        standardise = FALSE,
                        control = renewal.control(trace = 0, start = start)
                        )

## =========================== FORMULA API ======================================
print("............ FORMULA API  ............")
form_new <- children ~ german + years_school + voc_train + university + religion +
    year_birth + rural + age_marriage | german + years_school + voc_train

res_new <- renewalCount(formula = form_new, data = fertility, dist = "gamma",
                        computeHessian = FALSE,  standardise = FALSE,
                        control = renewal.control(trace = 0, start = start)
                        )
     
expect_equal(coef(res_new), coef(res_old))
                 
                 
