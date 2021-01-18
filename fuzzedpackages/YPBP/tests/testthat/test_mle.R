
library(YPBP)

data("gastric")

# fiting the models:
fit1 <- ypbp(Surv(time, status)~trt, data=gastric, baseline = "hazard")
fit2 <- ypbp(Surv(time, status)~trt, data=gastric, baseline = "odds")

# summarizing:
summary(fit1)
summary(fit2)

# extracting regression coefficients:
coef(fit1)
coef(fit2)

# extracting covariance matrix associated with the regression coefficients:
vcov(fit1)
vcov(fit2)

# CI for regression coefficients:
confint(fit1)
confint(fit2)


# extracting the design matrix:
model.matrix(fit1)
model.matrix(fit2)

# obtaining the survival probatilities:
ekm <- survival::survfit(Surv(time, status)~trt, data=gastric)
newdata <- data.frame(trt=0:1)
St <- survfit(fit1, newdata)
plot(ekm, col=1:2)
with(St, lines(time, surv[[1]]))
with(St, lines(time, surv[[2]], col=2))


# computing the crossing survival time:
newdata1 <- data.frame(trt=0)
newdata2 <- data.frame(trt=1)
tcross <- crossTime(fit1, newdata1, newdata2, nboot = 100)
tcross
